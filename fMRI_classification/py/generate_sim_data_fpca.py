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

from skfda.representation.grid import FDataGrid
from skfda.preprocessing.dim_reduction.projection import FPCA


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


# Crop dataframe to 50*200
alignedsub = {key: df.iloc[:50, :200] for key, df in alignedsub.items()}


n1 = 40
n2 = 60
y_ = [1]*n1 + [0]*n2

def generate_samples(data, n_samples, n_components=50):
    
    fd = FDataGrid(data_matrix=data.values, grid_points=np.linspace(0, 1, 200))
    fpca = FPCA(n_components=n_components)
    fpca.fit(fd)

    meanfunc = pd.DataFrame(fpca.mean_.data_matrix[0])
    eigenfuncs = np.array([ef.data_matrix[0] for ef in fpca.components_])
    eigenfuncs = pd.DataFrame(eigenfuncs[:,:,0])
    eigenvals = fpca.explained_variance_
    
    ksi = np.random.normal(0, np.sqrt(eigenvals), (n_samples, len(eigenvals)))
    samples = meanfunc.values.flatten() + np.dot(ksi, eigenfuncs)

    return samples

ADHDsub = [alignedsub[key] for key in alignedsub if key in ADHD_list]
controlsub = [alignedsub[key] for key in alignedsub if key in control_list]    



# Generate simulation data: n-m-p array (n = 100, m = 200, p = 50)
def gen_sim_data(n1, n2, y, seed = 1):
    # Set seed number for reproducible result
    np.random.seed(seed)

    # Generate per each region
    func_data = np.zeros((100, 200, 50))  
    for j in range(50):
        ADHD_data = pd.DataFrame([df.iloc[j, :] for df in ADHDsub])
        control_data = pd.DataFrame([df.iloc[j, :] for df in controlsub]) 

        func_data[:n1, :, j] = generate_samples(ADHD_data, n1)
        func_data[n1:, :, j] = generate_samples(control_data, n2)

    X_train, X_test, y_train, y_test = train_test_split(func_data, y, test_size=0.3)

    return X_train, X_test, y_train, y_test
