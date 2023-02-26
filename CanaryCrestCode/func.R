## shedding_parameters
##' @export
shedding_parameters <- function(M=NA, frac_ctdna=NA, tumor_diameter_cm=NA, sample_blood_volume_mL=15, normal_ng_per_mL=2.9, total_blood_volume_mL=5000, e=33, r=0.004, q=1.4e-4, death_rate=0.136) {

    ## LUAD/LUSC values:
    ## r=0.004, q=1.4e-4, death_rate=0.136

    ## expected number of normal-derived fragments
    hge_normal_dna <- normal_ng_per_mL * sample_blood_volume_mL * 1000 / 6.6
    normal_fragments_per_bp <- hge_normal_dna*2

    if(!is.na(M)) {

        ## expected number of tumor-derived fragments
        tumor_volume_cm3 <- M/(10^9)
        tumor_diameter_cm <- 2*((3*tumor_volume_cm3)/(4*pi))^(1/3)
        C <- M*death_rate*q/(e+r)
        hge_tumor_dna <- C*sample_blood_volume_mL/total_blood_volume_mL ## scale C by 15mL out of 5000mL
        tumor_fragments_per_bp <- hge_tumor_dna*2
        total_fragments_per_bp <- normal_fragments_per_bp + tumor_fragments_per_bp

        ## overall fraction of tumor-derived cfDNA
        frac_ctdna <- tumor_fragments_per_bp / total_fragments_per_bp

    } else if(!is.na(frac_ctdna)) {

        ## we know f, calculate tumor fragments per bp from it
        tumor_fragments_per_bp <- (frac_ctdna / (1-frac_ctdna)) * normal_fragments_per_bp
        total_fragments_per_bp <- normal_fragments_per_bp + tumor_fragments_per_bp

        ## calculate tumor attributes that give the specified f
        hge_tumor_dna <- tumor_fragments_per_bp / 2
        C <- hge_tumor_dna * total_blood_volume_mL / sample_blood_volume_mL
        M <- C*(e+r)/(death_rate*q)
        tumor_volume_cm3 <- M/(10^9)
        tumor_diameter_cm <- 2*((3*tumor_volume_cm3)/(4*pi))^(1/3)

    } else if(!is.na(tumor_diameter_cm)) {
        tumor_volume_cm3 <- (4/3)*pi*(tumor_diameter_cm/2)^3
        M <- tumor_volume_cm3 * 10^9
        C <- M * death_rate * q / (e+r)

        hge_tumor_dna <- C*sample_blood_volume_mL/total_blood_volume_mL ## scale C by 15mL out of 5000mL
        tumor_fragments_per_bp <- hge_tumor_dna*2
        total_fragments_per_bp <- normal_fragments_per_bp + tumor_fragments_per_bp

        ## overall fraction of tumor-derived cfDNA
        frac_ctdna <- tumor_fragments_per_bp / total_fragments_per_bp
    }

    list(normal_fragments_per_bp = normal_fragments_per_bp, tumor_fragments_per_bp = tumor_fragments_per_bp, total_fragments_per_bp = total_fragments_per_bp, tumor_diameter_cm = tumor_diameter_cm, tumor_volume_cm3 = tumor_volume_cm3, frac_ctdna = frac_ctdna, M=M)
}



##' simulate_fragments
##' @export
simulate_fragments <- function(M, ext_dir, covariate_binsize=1000) {

    ## get parameters from shedding model
    params <- shedding_parameters(M=M)

    ## load covariates (1000kb bins with annotated mappability/GC%, excludes blacklisted/telomere/centromere)
    d <- fread(file.path(ext_dir,'covariates_hg19.txt.gz'))

    ## get the total number of shed DNA fragments for tumor and normal
    hGE_normal <- params$normal_fragments_per_bp
    hGE_tumor <- params$tumor_fragments_per_bp

    ## get allele1/allele2 DNA fragments from Normal and Tumor
    ## the number of fragments should be scaled by the mean number overall
    d$normal_fragments <- rpois(n=nrow(d),lambda=hGE_normal*covariate_binsize)
    d$normal_fragments_allele1 <- rbinom(size=d$normal_fragments,prob=0.5,n=nrow(d))
    d$normal_fragments_allele2 <- d$normal_fragments - d$normal_fragments_allele1
    d$tumor_fragments <- rpois(n=nrow(d),lambda=hGE_tumor*covariate_binsize)
    d$tumor_fragments_allele1 <- rbinom(size=d$tumor_fragments,prob=0.5,n=nrow(d))
    d$tumor_fragments_allele2 <- d$tumor_fragments - d$tumor_fragments_allele1

    d
}


##' sample_ascn_profiles
##' @export
sample_ascn_profiles <- function(tcga_project, ext_dir, n_tumors=NA) {
    project_file <- file.path(ext_dir,paste0('ASCN/',toupper(tcga_project),'_segments.txt.gz'))
    message(project_file)

    if(!file.exists(project_file)) stop('Allele-specific copy number segments for [',tcga_project,'] not found!')

    ## load 
    seg <- fread(project_file)

    ## randomly select a number of profiles with replacement
    if(!is.na(n_tumors)) {
        tumors <- sort(unique(seg$PATIENT_ID))
        if(n_tumors > length(tumors)) {
            message('Requested ',n_tumors,' (greater than ',length(tumors),' available), sampling with replacement.')
            tumors <- sample(x=tumors, replace=T, size=n_tumors) 
        } else {
            message('Requested ',n_tumors,' (less than ',length(tumors),' available), sampling without replacement.')
            tumors <- sample(x=tumors, replace=F, size=n_tumors)        
        }
        sample_tumors <- function(i, tumors, seg) {
            sample_segments <- seg[PATIENT_ID %in% tumors[i],]
            sample_segments <- cbind(sample=i, sample_segments)
            sample_segments
        }
        seg <- rbindlist(lapply(1:length(tumors), sample_tumors, tumors, seg))
    
    } else {
        message('Number of tumors not specified, returning all segments for this project.')
        sample_num <- as.integer(factor(seg$PATIENT_ID, levels=unique(seg$PATIENT_ID)))
        seg <- cbind(sample=sample_num, seg)
    }

    seg
}


##' sample_hetsnps
##' @export
sample_hetsnps <- function(ext_dir, n_genomes=NA) {
    genomes <- sort(dir(file.path(ext_dir,'1kg'),full.names=T))

    ## randomly select a number of profiles with replacement
    if(!is.na(n_genomes)) {
        if(n_genomes > length(genomes)) {
            message('Requested ',n_genomes,' (greater than ',length(genomes),' available), sampling with replacement.')
            genomes <- sample(x=genomes, replace=T, size=n_genomes) 
        } else {
            message('Requested ',n_genomes,' (less than ',length(genomes),' available), sampling without replacement.')
            genomes <- sample(x=genomes, replace=F, size=n_genomes)        
        }
    } else {
        message('Number of tumors not specified, returning all segments for this project.')
    }
    
    ## load the genotypes for the (sampled) genomes
    get_hetsnps <- function(snp_file) {
        snp <- fread(snp_file)
        names(snp) <- c('chr','pos','rsid','ref','alt','gt')
        snp <- snp[ref %in% c('A','C','G','T') & alt %in% c('A','C','G','T'),]
        snp$rsid <- paste0(snp$chr,':',snp$ref,snp$pos,snp$alt)
        snp[gt==paste(ref,alt,sep='|'),allele_order:='FWD']
        snp[gt==paste(alt,ref,sep='|'),allele_order:='REV']
        snp$file <- snp_file
        snp
    }

    snps <- rbindlist(lapply(genomes, get_hetsnps))
    snps
}


##' simulate_reads
##' @export
simulate_reads <- function(gt, seg, fragments, C, adjusted_coverage, poisson_coverage, seq_error_rate) {
    require(mc2d)
    ## seq_error_rate=0.0024 (PMID: 30026539)

    ## combine fragment counts per bins with the het-snp genotypes
    gt <- copy(gt)
    gt$pos2 <- gt$pos + 1
    setkey(gt,'chr','pos','pos2')
    fragments <- fragments[,c('covariate_bin','chr','arm','start','end','mappability','gc_percent','median_scaling_factor','normal_fragments_allele1','normal_fragments_allele2','tumor_fragments_allele1','tumor_fragments_allele2'),with=F]
    setkey(fragments,'chr','start','end')
    gt <- foverlaps(gt, fragments, type='within')
    gt <- gt[!is.na(covariate_bin),] ## NB: this bin refers to the bin in the covariates file
    gt[,c('ref','alt','gt'):=NULL]

    ## merge the CNV segments onto the snps
    seg <- copy(seg)
    #seg[,c('num_probes','seg_length','arm'):=NULL]
    setkey(gt,'chr','pos','pos2')
    setkey(seg,'chr','seg_start','seg_end')
    if('arm' %in% names(seg)) seg[,arm:=NULL]
    snp <- foverlaps(gt, seg, type='any')
    tumor <- unique(snp$sample[!is.na(snp$sample)])
    PROJECT_CODE <- unique(snp$PROJECT_CODE[!is.na(snp$PROJECT_CODE)])
    PATIENT_ID <- unique(snp$PATIENT_ID[!is.na(snp$PATIENT_ID)])
    snp$sample <- tumor
    snp$PATIENT_ID <- PATIENT_ID
    snp$PROJECT_CODE <- PROJECT_CODE
    snp[,c('seg_start','seg_end','pos2'):=NULL]
    snp[is.na(cn1),cn1:=1]
    snp[is.na(cn2),cn2:=1]
    snp[is.na(ccf1),ccf1:=1]
    snp[is.na(ccf2),ccf2:=1]

    ## adjust allele1/2 fragments due to allele-specific copy number, but retain the SNP phase
    snp_FWD <- snp[allele_order=='FWD']
    snp_FWD$tumor_copies_allele1 <- snp_FWD$ccf1 * snp_FWD$cn1 + (1-snp_FWD$ccf1) * 1
    snp_FWD$tumor_copies_allele2 <- snp_FWD$ccf2 * snp_FWD$cn2 + (1-snp_FWD$ccf2) * 1
    snp_REV <- snp[allele_order=='REV']
    snp_REV$tumor_copies_allele1 <- snp_REV$ccf2 * snp_REV$cn2 + (1-snp_REV$ccf2) * 1
    snp_REV$tumor_copies_allele2 <- snp_REV$ccf1 * snp_REV$cn1 + (1-snp_REV$ccf1) * 1
    snp <- rbind(snp_FWD, snp_REV)
    snp$tumor_fragments_allele1 <- round(snp$tumor_fragments_allele1 * snp$tumor_copies_allele1)
    snp$tumor_fragments_allele2 <- round(snp$tumor_fragments_allele2 * snp$tumor_copies_allele2)
    snp$tumor_fragments <- snp$tumor_fragments_allele1 + snp$tumor_fragments_allele2
    snp$normal_fragments <- snp$normal_fragments_allele1 + snp$normal_fragments_allele2
    snp$total_fragments <- snp$normal_fragments + snp$tumor_fragments
    snp$depth <- C
    if(adjusted_coverage) snp$depth <- snp$total_fragments / mean(snp$total_fragments) * snp$depth * snp$median_scaling_factor

    ## get the probability of a read mapping to each class
    snp$p_allele1_tumor_reads <- snp$tumor_fragments_allele1 / snp$total_fragments
    snp$p_allele2_tumor_reads <- snp$tumor_fragments_allele2 / snp$total_fragments
    snp$p_allele1_normal_reads <- snp$normal_fragments_allele1 / snp$total_fragments
    snp$p_allele2_normal_reads <- snp$normal_fragments_allele2 / snp$total_fragments

    ## calculate the depth at each SNP = based on depth C x the median-scaling-fractor
    if(poisson_coverage) snp$depth <- rpois(lambda=snp$depth, n=nrow(snp)) 

    ## simulate the 4 classes of reads based on a multinomial distribution with total trials = depth
    prob_mat = cbind(snp$p_allele1_tumor_reads,snp$p_allele2_tumor_reads,snp$p_allele1_normal_reads,snp$p_allele2_normal_reads)
    tbl <- rmultinomial(nrow(snp), size = snp$depth, prob=prob_mat)
    tbl <- as.data.table(tbl)
    names(tbl) <- c('allele1_tumor_reads','allele2_tumor_reads','allele1_normal_reads','allele2_normal_reads')
    snp <- cbind(snp, tbl)
    snp <- snp[order(chr,pos),]

    ## add in effects of sequencing error (err=0.001?)
    allele1_tumor_reads_dropped <- rbinom(n=nrow(snp),size=snp$allele1_tumor_reads,prob=0.6667*seq_error_rate)
    allele1_tumor_reads_added_to_allele2 <- rbinom(n=nrow(snp),size=snp$allele1_tumor_reads,prob=0.3333*seq_error_rate)
    allele2_tumor_reads_dropped <- rbinom(n=nrow(snp),size=snp$allele2_tumor_reads,prob=0.6667*seq_error_rate)
    allele2_tumor_reads_added_to_allele1 <- rbinom(n=nrow(snp),size=snp$allele2_tumor_reads,prob=0.3333*seq_error_rate)
    allele1_normal_reads_dropped <- rbinom(n=nrow(snp),size=snp$allele1_normal_reads,prob=0.6667*seq_error_rate)
    allele1_normal_reads_added_to_allele2 <- rbinom(n=nrow(snp),size=snp$allele1_normal_reads,prob=0.3333*seq_error_rate)
    allele2_normal_reads_dropped <- rbinom(n=nrow(snp),size=snp$allele2_normal_reads,prob=0.6667*seq_error_rate)
    allele2_normal_reads_added_to_allele1 <- rbinom(n=nrow(snp),size=snp$allele2_normal_reads,prob=0.3333*seq_error_rate)

    ## update the read counts including the reads dropped due to error and switched from between allele1/2 due to error
    snp$allele1_normal_reads <- snp$allele1_normal_reads - allele1_normal_reads_dropped - allele1_normal_reads_added_to_allele2 + allele2_normal_reads_added_to_allele1
    snp$allele2_normal_reads <- snp$allele2_normal_reads - allele2_normal_reads_dropped - allele2_normal_reads_added_to_allele1 + allele1_normal_reads_added_to_allele2
    snp$allele1_tumor_reads <- snp$allele1_tumor_reads - allele1_tumor_reads_dropped - allele1_tumor_reads_added_to_allele2 + allele2_tumor_reads_added_to_allele1
    snp$allele2_tumor_reads <- snp$allele2_tumor_reads - allele2_tumor_reads_dropped - allele2_tumor_reads_added_to_allele1 + allele1_tumor_reads_added_to_allele2

    ## get the total allele1-reads and allele2-reads in each snp. NB: these will actually be the average number of reads for a1/a2 in the snp among any het-SNPs in the snp
    snp$allele1_reads <- snp$allele1_tumor_reads + snp$allele1_normal_reads
    snp$allele2_reads <- snp$allele2_tumor_reads + snp$allele2_normal_reads
    snp$adjusted_coverage <- adjusted_coverage
    snp$poisson_coverage <- poisson_coverage
    snp$seq_error_rate <- seq_error_rate
    snp
}


##' bin_genome
##' @export
bin_genome <- function(binsize=1e6) {
    expand_to_regions <- function(chr,binsize) {
        l <- chr$length
        bins <- ceiling(l/binsize)
        region_start <- seq(1,l,by=binsize)
        region_end <- region_start + binsize - 1
        out <- data.table(bin_start=region_start,bin_end=region_end)
        out <- out[bin_end <= l]
        out
    }
    chr <- fread(system.file("extdata", "chr_lengths_hg19.txt", package = "AIdetector"))
    regions <- chr[,expand_to_regions(.SD,binsize),by=chr]
    regions <- cbind(bin=1:nrow(regions),regions)
    regions
}


##' bin_snps
##' @export
bin_snps <- function(d, binsize) { 
    ## use bin_genome() to split the genome into desired binsize, then split SNP data into these bins

    d$pos2 <- d$pos + 1
    bins <- bin_genome(binsize)
    old_fields <- intersect(c('bin','start','end'),names(d))
    if(length(old_fields) > 0) d[,(old_fields):=NULL]
    names(bins) <- c('bin','chr','start','end')
    # 279, 280, 281 to get overlaps
    setkey(bins,'chr','start','end') # Exons table I generated
    setkey(d,'chr','pos','pos2') # snps
    d <- foverlaps(d, bins, type='within')
    d <- d[!is.na(bin),]
    d[,pos2:=NULL]
    d
}





##' prepare_hetsnp_counts
##' @export
prepare_hetsnp_counts <- function(seg, gt, C, M=NA, frac_ctdna=NA, tumor_diameter_cm=NA, remove_acrocentric=T, adjusted_coverage=T, poisson_coverage=T, seq_error_rate=0.0024) { 
    if(!is.na(M)) {
        param <- shedding_parameters(M=M); param   
    } else if(!is.na(frac_ctdna)) {
        param <- shedding_parameters(frac_ctdna=frac_ctdna); param
    } else if(!is.na(tumor_diameter_cm)) {
        param <- shedding_parameters(tumor_diameter_cm=tumor_diameter_cm); param
    }
    
    M <- param$M
    frac_ctdna <- param$frac_ctdna
    tumor_diameter_cm <- param$tumor_diameter_cm

    ## simulate allele-specific fragment counts (the numbers in each bin correspond to the average expected number per basepair)
    message('Simulating cfDNA fragments with M=',M,' ...')
    fragments <- simulate_fragments(M=M, ext_dir='ext')

    ## simulate read counts bassd on the avg fragments shed in each bin
    message('Simulating cfDNA reads with C=',C,' ...')
    # Modify gt
    d <- simulate_reads(gt, seg, fragments, C=C, adjusted_coverage=adjusted_coverage, poisson_coverage=poisson_coverage, seq_error_rate=seq_error_rate)
    d$allele1_reads_unphased <- d$allele1_reads
    d$allele2_reads_unphased <- d$allele2_reads
    d$M <- M
    d$frac_ctdna <- frac_ctdna
    d$tumor_diameter_cm <- tumor_diameter_cm


    ## remove bad chromosome arms
    if(remove_acrocentric==T) {
        message('Removing short arms from acrocentric chromosomes: 13p, 14p, 15p, 21p, 22p.')
        d$charm <- paste0(d$chr,d$arm)
        bad_arms <- c('13p','14p','15p','21p','22p')
        d <- d[!charm %in% bad_arms,]
    }

    ## apply haplotype phasing
    message('Adjusting counts for haplotype phasing ...')
    get_phased <- function(d) {
        ## nb: i confirmed phasing is working (c=100,000, m from 0 to 10^10)
        d_fwd <- d[allele_order=='FWD']
        d_rev <- d[allele_order=='REV']
        setnames(d_rev,
                 c('allele1_tumor_reads','allele2_tumor_reads','allele1_normal_reads','allele2_normal_reads'),
                 c('allele2_tumor_reads','allele1_tumor_reads','allele2_normal_reads','allele1_normal_reads'))
        d_rev <- d_rev[,(names(d_fwd)),with=F]
        d <- rbind(d_fwd, d_rev)
        d <- d[order(chr,pos),]
        d$allele1_reads <- d$allele1_tumor_reads + d$allele1_normal_reads
        d$allele2_reads <- d$allele2_tumor_reads + d$allele2_normal_reads
        d$total_reads <- d$allele1_reads + d$allele2_reads
        d
    }
    d <- get_phased(d)

    ## correct gc/mappability bias. NB: total_reads_adj maintains the origianl grand mean, mu
    message('Adjusting counts to remove %GC and mappability biases ...')
    mu <- mean(d$total_reads)
    b <- mgcv::gam(total_reads ~ 0 + s(gc_percent, bs='cr', sp=0.6) + s(mappability, bs='cr', sp=0.3), data=d)
    reads_expected_from_bias <- predict(b)
    adjustment_factor <- (d$total_reads - reads_expected_from_bias) / d$total_reads
    d$total_reads_adj <- d$total_reads * adjustment_factor
    d$allele1_reads_adj <- d$allele1_reads * adjustment_factor
    d$allele2_reads_adj <- d$allele2_reads * adjustment_factor

    d
}


##' betareg_shapes
##' @export
betareg_shapes <- function(y) {
    require(betareg)
    tryCatch({
        fit <- betareg(y ~ 1)
        coef <- summary(fit)$coef
        mu <- plogis(coef$mean[1])
        phi <- coef$precision[1]
        shape1 <- mu*phi
        shape2 <- phi*(1-mu)
        list(shape1=shape1,shape2=shape2)
    }, error=function(e) {
        list(shape1=as.numeric(NA),shape2=as.numeric(NA))
    },warnings=function(w) {
        list(shape1=as.numeric(NA),shape2=as.numeric(NA))   
    })
}


##' ctda_fraction_mle
##' @export
ctdna_fraction_mle <- function(VAF,c1,c2) {
    require(boot)
    require(stats)

    loglikelihood <- function(f, c1, c2, VAF) {
        G <- shedding_parameters(frac_ctdna=f)$total_fragments_per_bp
        shape1 <- (1-f+c1*f)*G
        shape2 <- (1-f+c2*f)*G
        ll <- -1*sum(dbeta(VAF, shape1, shape2, log=T))
        ifelse(f < 0 | f > 1, as.numeric(NA), ll)
    }
    
    opt <- optim(0, loglikelihood, c1=c1, c2=c2, VAF=VAF, method="Brent", lower=0, upper=1)
    opt <- data.table(f=opt$par[1], likelihood=-1*opt$value, c1=c1, c2=c2)
    opt

}


##' write.tsv
##' @export
write.tsv <- function(d, file, sep = "\t", quote = F, row.names = F, ...) {
    write.table(d, file = file, sep = sep, quote = quote, row.names = row.names, ...)
}


##' test_segments
##' @export
test_segments <- function(d,c1=NA,c2=NA,min_f=1e-4,R, n_cpus=1) {
    charm <- d$charm[1]
    message(charm)


    if(is.na(c1) | is.na(c2)) {
        ## decide which c-vals we should use. Take c-vals from the most likely run. If there's a tie, choose cvals that give f closest to the median value.
        #pre_run0 <- ctdna_fraction_mle(d$VAF,1,1) ## balanced
        pre_run1 <- ctdna_fraction_mle(d$VAF,0,1) ## Del1
        pre_run2 <- ctdna_fraction_mle(d$VAF,1,0) ## Del2
        #pre_run3 <- ctdna_fraction_mle(d$VAF,1,2) ## Dup1
        #pre_run4 <- ctdna_fraction_mle(d$VAF,2,1) ## Dup2
        #pre_run5 <- ctdna_fraction_mle(d$VAF,0,2) ## Del1+Dup2
        #pre_run6 <- ctdna_fraction_mle(d$VAF,2,0) ## Del2+Dup1
        pre_run <- rbind(pre_run1, pre_run2) 
        #pre_run <- rbind(pre_run0, pre_run1, pre_run2, pre_run3, pre_run4, pre_run5, pre_run6)
        pre_run <- pre_run[order(likelihood,f,decreasing=T),]
        pre_run <- pre_run[likelihood==max(likelihood),]
        mid_val <- median(pre_run$f)
        pre_run$diff_from_median <- abs(pre_run$f - mid_val)
        pre_run <- pre_run[order(diff_from_median,decreasing=F),]
        c1 <- pre_run$c1[1]
        c2 <- pre_run$c2[1]
    }

    ## now run for the best cvals with bootstrapped CIs
    fn <- function(data, c1, c2, indices) {
        d <- data[indices,]
        res <- ctdna_fraction_mle(d$VAF, c1, c2)
        return(res$f)
    }

    if(n_cpus > 1) {
        bt <- boot(data=d, statistic=fn, R=R, c1=c1, c2=c2, parallel='multicore', ncpus=n_cpus)
    } else {
        bt <- boot(data=d, statistic=fn, R=R, c1=c1, c2=c2)   
    }

    L <- boot::empinf(bt, data=d, statistic=fn, c1=c1, c2=c2, type='jack')
    var <- var.linear(unlist(L))
    mu <- mean(bt$t)

    estimate <- bt$t0
    boot_vals <- bt$t 
    boot_vals <- boot_vals[boot_vals > 0 & boot_vals < 1 & !is.na(boot_vals)]
    if(length(boot_vals) >= 30 & estimate > min_f) {
        shapes <- suppressMessages(betareg_shapes(boot_vals))
    } else {
        shapes <- list(shape1=as.numeric(NA),shape2=as.numeric(NA))
    }


    bt$t[bt$t < min_f] <- 0
    bt$t0[bt$t0 < min_f] <- 0
    estimate <- bt$t0
    ci <- as.numeric(quantile(bt$t,c(0.025,0.975)))
    prop_gtr0 <- mean(bt$t >= min_f)

    list(charm=charm, estimate=estimate, lwr=ci[1], upr=ci[2], c1=c1, c2=c2, mu=mu, var=var, R=R, prop_gtr0=prop_gtr0, shape1=shapes[[1]], shape2=shapes[[2]])
}


##' get_fraction_estimate
##' @export
get_fraction_estimate <- function(result,k=NA,min_f) {

    ## estimating overall tumor DNA fraction
    ## 1. If no segments had lower-95%CI > 0, then use the estimate and CI from whichever segment has the median estimate.
    ## 2. Otherwise, use 1D clustering of the estimates to find a cluster of at least two segments with the highest center. Use the estimate/CI from the segment in this cluster closest to its center.
    ## 3. If no clusters had > 1 segments, then take the median center among the clusters, then use the estimate/CI from the segment closest to this center.

    tryCatch({
        if(is.na(k)) {
            k <- round(sqrt(sum(result$lwr >= min_f)))
        }

        if(sum(result$lwr >= min_f) > k) {
            message('Attempting clustering with k=',k,' ...')
            ## there are at least 4 segments with evidence for f > 0. Attempt clustering them to choose the best representative subset.
            tmp <- result[result$lwr >= min_f]
            clus <- Ckmeans.1d.dp(tmp$estimate, k=k) ## require 3 clusters (low,mid,high)
            tmp$cluster <- clus$cluster
            km <- data.table(center=clus$center, size=clus$size, cluster=unique(clus$cluster))
            if(any(km$size >= k)) {
                ## if any clusters have enough segments, use segments from the cluster of this min-size with the greatest mean
                km <- km[size >= k]
                km <- km[order(center,decreasing=T),]
                km <- km[1,]
                tmp <- tmp[cluster %in% km$cluster,]
            }
            ## NB: if no clusters with 2+ segments, then clustering fails. Use all segments with lwr > min_f

        } else if(sum(result$lwr > min_f) <= k) {
            ## too few segments with f > min_f for clustering. Just use all of these.
            tmp <- result[result$lwr > min_f]    

        } else {
            ## no clusters had significant evidence of f > 0, so we will consider all segments
            tmp <- copy(result)
        }

        ## metaanalysis?
        ## - https://stats.stackexchange.com/questions/3652/combining-two-confidence-intervals-point-estimates
        ## - https://en.wikipedia.org/wiki/Pooled_variance
        pooled_mean <- mean(tmp$estimate)

        ## pooled variance
        pooled_var <- sum((tmp$R - 1)*tmp$var) / sum(tmp$R-1)

        ## pooled 95% CIs
        pooled_sd <- sqrt(pooled_var)
        pooled_N <- nrow(tmp)
        pooled_lwr <- pooled_mean - 1.96*pooled_sd / sqrt(pooled_N)
        pooled_upr <- pooled_mean + 1.96*pooled_sd / sqrt(pooled_N)

        list(estimate=pooled_mean, lwr=pooled_lwr, upr=pooled_upr) ## round to 4 dec.

    },error=function(e) {
        na <- as.numeric(NA)
        list(estimate=na, lwr=na, upr=na)
    })

}



