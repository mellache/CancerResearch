import pandas as pd
import os
import numpy as np
import random

import matplotlib.pyplot as plt

from sksurv.datasets import load_breast_cancer
from sksurv.linear_model import CoxPHSurvivalAnalysis, CoxnetSurvivalAnalysis
from sksurv.preprocessing import OneHotEncoder

from sklearn.model_selection import GridSearchCV, KFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

import warnings
from sklearn.exceptions import FitFailedWarning

from sksurv.ensemble import RandomSurvivalForest
from sksurv.metrics import cumulative_dynamic_auc
from sksurv.metrics import concordance_index_censored
from sksurv.metrics import concordance_index_ipcw
from sksurv.metrics import as_concordance_index_ipcw_scorer
from sksurv.metrics import as_cumulative_dynamic_auc_scorer

data_dir = "/u/project/xjzhou/huran/BIG_Summer/"

sample_sheet = pd.read_csv(data_dir+"LIHC_sample_sheet.csv", sep=",", index_col=0)

m_mat = pd.read_csv(data_dir+"processed_m_values.csv", sep=",", index_col=0)
primary_tumor_samples = sample_sheet[sample_sheet["Sample Type"]=="Primary Tumor"]
m_mat_primary_tumor = m_mat[primary_tumor_samples.index]
m_mat_primary_tumor.columns = sample_sheet.loc[m_mat_primary_tumor.columns, "Case ID"]

univariate_m_values = pd.read_csv(data_dir+"univariate_Cox_selected_features_m.txt", sep=",", index_col=0)
select_probes = univariate_m_values[univariate_m_values["pvalue"]<0.01].index #narrow down the number of features by using lower p-value

annot = pd.read_csv(os.path.join(data_dir, "train_set.csv"), index_col=0)

X = m_mat_primary_tumor.loc[select_probes, annot.index] #5562*218
X = X.T #218*5562

def get_y_label(annot):
    label = annot.iloc[:,:2]
    y = []
    for i in range(len(label)):
        data = label.iloc[i]
        if data["status"]:
            y.append((True, data["time"]))
        else:
            y.append((False, data["time"]))
    y = np.array(y, dtype=[('cens', '?'), ('time', '<f8')])
    return y

y = get_y_label(annot)



coxnet_pipe = make_pipeline(
    StandardScaler(),
    CoxnetSurvivalAnalysis(l1_ratio=0.9, alpha_min_ratio=0.01, max_iter=1000)
)
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", FitFailedWarning)
coxnet_pipe.fit(X, y)


estimated_alphas = coxnet_pipe.named_steps["coxnetsurvivalanalysis"].alphas_
cv = KFold(n_splits=5, shuffle=True, random_state=0)
gcv = GridSearchCV(
    make_pipeline(StandardScaler(), CoxnetSurvivalAnalysis(l1_ratio=0.9)),
    param_grid={"coxnetsurvivalanalysis__alphas": [[v] for v in estimated_alphas]},
    cv=cv,
    error_score=0.5,
    n_jobs=4).fit(X, y)

cv_results = pd.DataFrame(gcv.cv_results_)


alphas = cv_results.param_coxnetsurvivalanalysis__alphas.map(lambda x: x[0])
mean = cv_results.mean_test_score
std = cv_results.std_test_score

fig, ax = plt.subplots(figsize=(9, 6))
ax.plot(alphas, mean)
ax.fill_between(alphas, mean - std, mean + std, alpha=.15)
ax.set_xscale("log")
ax.set_ylabel("concordance index")
ax.set_xlabel("alpha")
ax.axvline(gcv.best_params_["coxnetsurvivalanalysis__alphas"][0], c="C1")
ax.axhline(0.5, color="grey", linestyle="--")
ax.grid(True)
plt.savefig('TCGA_Coxnet_1_plot.pdf')

best_model = gcv.best_estimator_.named_steps["coxnetsurvivalanalysis"]
best_coefs = pd.DataFrame(
    best_model.coef_,
    index=X.columns,
    columns=["coefficient"]
)

non_zero = np.sum(best_coefs.iloc[:, 0] != 0)
print("Number of non-zero coefficients: {}".format(non_zero))

non_zero_coefs = best_coefs.query("coefficient != 0")
coef_order = non_zero_coefs.abs().sort_values("coefficient").index

_, ax = plt.subplots(figsize=(6, 8))
non_zero_coefs.loc[coef_order].plot.barh(ax=ax, legend=False)
ax.set_xlabel("coefficient")
ax.grid(True)
plt.savefig('TCGA_Coxnet_2_plot.pdf')

out_dir = "/u/project/xjzhou/huran/BIG_Summer/results/"
non_zero_coefs.to_csv(out_dir+"TCGA_Coxnet_nonzero_features.csv")

##Validate the model on validation dataset

annot_validate = pd.read_csv(os.path.join(data_dir, "validate_set.csv"), index_col=0)

X_validate = m_mat_primary_tumor.loc[select_probes, annot_validate.index] #5562*108
X_validate = X_validate.T #108*5562

y_validate = get_y_label(annot_validate) 

coxnet_pred = make_pipeline(
    StandardScaler(),
    CoxnetSurvivalAnalysis(l1_ratio=0.9, fit_baseline_model=False)
)
coxnet_pred.set_params(**gcv.best_params_)
coxnet_pred.fit(X, y)

risk_score = coxnet_pred.predict(X_validate)

time_list = [int(i) for i in np.quantile([i[1] for i in y_validate],[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])]
#lower, upper = np.percentile(y["time"], [10, 90])
#gbsg_times = np.arange(lower, upper + 1)

#tAUROC = cumulative_dynamic_auc(survival_train=y, survival_test=y_validate, estimate=risk_score, times=gbsg_times)
tAUROC = cumulative_dynamic_auc(survival_train=y, survival_test=y_validate, estimate=risk_score, times=time_list)

AUC = tAUROC[0]
metrics = []
for i in range(len(time_list)):
    metrics.append(round(AUC[i],4))

#cindex = concordance_index_ipcw(survival_train=y, survival_test=y_validate, estimate=risk_score, tau=gbsg_times[-1])[0]
cindex = concordance_index_censored(event_indicator=[i[0] for i in y_validate], 
                               event_time=[i[1] for i in y_validate], 
                               estimate=risk_score)[0]

metrics.append(round(cindex,4))
metrics = pd.DataFrame(metrics).T

metrics.to_csv(out_dir+"Coxnet_validate_results2.csv")

risk_score = pd.DataFrame(risk_score, index=annot_validate.index)
risk_score.to_csv(out_dir+"Coxnet_validate_risk_score.txt")


