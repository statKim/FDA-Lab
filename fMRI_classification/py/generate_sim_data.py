import pandas as pd
import numpy as np
import random
import glob
import time
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings(action='ignore')

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score

from sklearn.pipeline import Pipeline
from sklearn.decomposition import SparsePCA, PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.naive_bayes import GaussianNB
from sklearn.linear_model import LogisticRegression, Lasso, Ridge
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier

from scipy.linalg import cholesky, solve
from scipy.stats import multivariate_normal

import os
print(os.getcwd())

def data_load(path):

    file_paths = glob.glob(path + '/*.csv')

    data_list = []

    for file_path in file_paths:
        data = pd.read_csv(file_path, index_col=0)
        data_list.append(data)

    return data_list

data = data_load("./fMRI_Classification/Clean_fMRI_PekingUniversity")
data

# num of region
region_data = data[0]["Region_ID"].value_counts().reset_index()
region_data.columns = ['Region_ID', 'Count']

region_data = region_data.sort_values(by='Region_ID').reset_index(drop=True)
region_data

# Generate input X, label y
X_ = [d.iloc[:,16:] for d in data]
y_ = [0]*len(data)

for i in range(len(data)):
    if data[i]['Phenotype_Diagnosis_L1'][0] == 'ADHD Diagnosed':
        y_[i] = 1

# Exclude outliers
ADHD_list = [i for i in range(len(y_)) if y_[i] == 1]
control_list = [i for i in range(len(y_)) if (y_[i] == 0)&(i!=16)&(i!=69)]

alignedsub = {}

for i in range(191):
    file_path = f'./fMRI_Classification/AlignedSubject/AlignedSub{i}.csv'
    alignedsub[i] = pd.read_csv(file_path, index_col=0)

# Mean curves
ADHD_meancurve = pd.DataFrame(alignedsub[ADHD_list[0]])
control_meancurve = pd.DataFrame(alignedsub[control_list[0]])

for num in ADHD_list[1:]:
    ADHD_meancurve=ADHD_meancurve+alignedsub[num]
    
for num in control_list[1:]:
    control_meancurve=control_meancurve+alignedsub[num]
    
ADHD_meancurve = ADHD_meancurve/len(ADHD_list)
control_meancurve = control_meancurve/len(control_list)

# mean function 50*200
group1_meanfunc = ADHD_meancurve.iloc[:50,:200]
group2_meanfunc = pd.concat([ADHD_meancurve.iloc[:10,:200], control_meancurve.iloc[10:50,:200]])

# Covariance function
def calc_sigma(X1, X2, l=1):
    Sigma = np.zeros((len(X1), len(X2)))
    for i in range(len(X1)):
        for j in range(len(X2)):
            Sigma[i, j] = np.exp(-0.5 * (np.abs(X1[i] - X2[j]) / l) ** 2)
    return Sigma

# Parameters
n1=40
n2=60
sigma_n1=2
sigma_n2=0.5

time_points = np.linspace(0, 199, 200)
x_star = np.linspace(0, 199, 200)
y_ = [1]*n1 + [0]*n2

k_xx = calc_sigma(time_points, time_points)
k_xxs = calc_sigma(time_points, x_star)
k_xsx = calc_sigma(x_star, time_points)
k_xsxs = calc_sigma(x_star, x_star)
k_xx_inv1 = solve(k_xx + sigma_n1**2 * np.eye(len(k_xx)), np.eye(len(k_xx)))
k_xx_inv2 = solve(k_xx + sigma_n2**2 * np.eye(len(k_xx)), np.eye(len(k_xx)))


# Covariance function for each group
f_bar_star_g1 = k_xsx @ k_xx_inv1 @ group1_meanfunc.T
f_bar_star_g2 = k_xsx @ k_xx_inv2 @ group2_meanfunc.T
cov_f_star_g1 = k_xsxs - k_xsx @ k_xx_inv1 @ k_xxs
cov_f_star_g2 = k_xsxs - k_xsx @ k_xx_inv2 @ k_xxs

# for i in range(100):
    
#     # Set seed number for reproducible result
#     np.random.seed(i)

#     # Generate per each region
#     func_data = np.zeros((100, 200, 50))  
#     for k in range(50):
#         # Group 1
#         func_data[:n1, :, k] = multivariate_normal.rvs(mean=f_bar_star_g1.iloc[:, k], cov=cov_f_star_g1, size = n1)

#         # Group 2
#         func_data[n1:, :, k] = multivariate_normal.rvs(mean=f_bar_star_g2.iloc[:, k], cov=cov_f_star_g2, size = n2)


# Generate simulation data: n-m-p array (n = 100, m = 200, p = 50)
def gen_sim_data(f_bar_star_g1, f_bar_star_g2, cov_f_star_g1, cov_f_star_g2, n1, n2, y, seed = 1):
    # Set seed number for reproducible result
    np.random.seed(seed)

    # Generate per each region
    func_data = np.zeros((100, 200, 50))  
    for k in range(50):
        # Group 1
        func_data[:n1, :, k] = multivariate_normal.rvs(mean=f_bar_star_g1.iloc[:, k], cov=cov_f_star_g1, size = n1)

        # Group 2
        func_data[n1:, :, k] = multivariate_normal.rvs(mean=f_bar_star_g2.iloc[:, k], cov=cov_f_star_g2, size = n2)
    

    X_train, X_test, y_train, y_test = train_test_split(func_data, y, test_size=0.3)

    return X_train, X_test, y_train, y_test
