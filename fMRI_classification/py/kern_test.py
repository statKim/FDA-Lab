# All of these (except TQDM) are standard imports and are included for example in Anaconda
# TQDM is a package for progress bars and can be easily installed, see https://pypi.org/project/tqdm/
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics.pairwise import rbf_kernel
from sklearn.metrics import pairwise_distances
from sklearn.metrics import pairwise_kernels
import scipy.linalg as lin
from scipy.special import zeta
from sklearn.decomposition import PCA
import pandas as pd
import sklearn.gaussian_process as gp
from scipy.stats import t
from scipy.stats import percentileofscore 
from tqdm import tqdm_notebook as tqdm
import os
import pandas as pd
import random
from matplotlib import image


# Norm Functions and T Maps
# The next two cells contain the code for the five different settings used in the numerical tests: ID, CEXP, FPCA, POLY and SQR.
# To see how to combine these to recreate the numerical results see comments in the cells which perform the tests.

def K_ID(X,Y,gamma=1):
    """
    Forms the kernel matrix K for the two sample test using the SE-T kernel with bandwidth gamma
    where T is the identity operator
    
    Parameters:
    X - (n_samples,n_obs) array of samples from the first distribution 
    Y - (n_samples,n_obs) array of samples from the second distribution 
    gamma - bandwidth for the kernel, if -1 then median heuristic is used to pick gamma
    
    Returns:
    K - matrix formed from the kernel values of all pairs of samples from the two distributions
    """
    n_obs = X.shape[1]
    XY = np.vstack((X,Y))
    dist_mat = (1/np.sqrt(n_obs))*pairwise_distances(XY, metric='euclidean')
    if gamma == -1:
        gamma = np.median(dist_mat[dist_mat > 0])
   
    K = np.exp(-0.5*(1/gamma**2)*(dist_mat**2))
    return K


def FPCA(X,n_comp = 0.95):
    """
    Computes principal components of given data up to a specified explained variance level
    
    Parameters:
    X - (n_samples,n_obs) array of function values
    n_comp - number of principal components to compute. If in (0,1) then it is the explained variance level
    
    Returns:
    Normalised eigenvalues and eigenfunctions
    """
    n_points = np.shape(X)[1]
    pca = PCA(n_components = n_comp)
    pca.fit(X)
    return (1/n_points)*pca.explained_variance_,pca.components_


def K_COV(X,Y,gamma=1):
    """
    Forms the kernel matrix K for the two sample test using the COV kernel
    
    Parameters:
    X - (n_samples,n_obs) array of samples from the first distribution 
    Y - (n_samples,n_obs) array of samples from the second distribution 
    gamma - dummy variable noot used in function, is an input for ease of compatability with other kernels
    
    Returns:
    K - matrix formed from the kernel values of all pairs of samples from the two distributions
    """    
    n_obs = X.shape[1]
    XY = np.vstack((X,Y))
    return ((1/n_obs)*np.dot(XY,XY.T))**2


def K_FPCA(X,Y,gamma = 1,n_comp = 0.95):
    """
    Forms the kernel matrix K for the two sample test using the SE-T kernel with bandwidth gamma
    where T is the FPCA decomposition operator
    
    Parameters:
    X - (n_samples,n_obs) array of samples from the first distribution 
    Y - (n_samples,n_obs) array of samples from the second distribution 
    gamma - bandwidth for the kernel, if -1 then median heuristic is used to pick gamma
    n_comp - number of principal components to compute. If in (0,1) then it is the explained variance level
    
    Returns:
    K - matrix formed from the kernel values of all pairs of samples from the two distributions
    """
    n_obs = X.shape[1]
    XY = np.vstack((X,Y))
    e_vals,e_funcs = FPCA(XY,n_comp = n_comp)
    scaled_e_funcs = e_funcs*np.sqrt(e_vals[:,np.newaxis])
    XY_e = (1/n_obs)*np.dot(XY,scaled_e_funcs.T)
    dist_mat = pairwise_distances(XY_e,metric='euclidean')
    if gamma == -1:
        gamma = np.median(dist_mat[dist_mat > 0])
    K = np.exp(-0.5*(1/gamma**2)*(dist_mat**2))
    return K


def K_SQR(X,Y,gamma = 1):
    """
    Forms the kernel matrix K for the two sample test using the SE-T kernel with bandwidth gamma
    where T is the map which sends x -> (x,x^{2}) in the Cartesian product of L^{2} with itself.
    
    Parameters:
    X - (n_samples,n_obs) array of samples from the first distribution 
    Y - (n_samples,n_obs) array of samples from the second distribution 
    gamma - bandwidth for the kernel to be used on the two norms, if -1 then median heuristic 
            is used to pick a different gamma for each norm, if gamma = 0 then median heuristic
            is used to pick a single gamma for each norm.
            
    Returns:
    K - matrix formed from the kernel values of all pairs of samples from the two distributions
    """
    n_obs = X.shape[1]
    XY = np.vstack((X,Y))
    dist_mat_1 = (1/np.sqrt(n_obs))*pairwise_distances(XY, metric='euclidean')
    dist_mat_2 = (1/np.sqrt(n_obs))*pairwise_distances(XY**2, metric='euclidean')
    dist_mat = dist_mat_1 + dist_mat_2
    if gamma == 0:
        gamma = np.median(dist_mat[dist_mat > 0])
        K = np.exp(-0.5*(1/gamma**2)*dist_mat**2)
        return K
    if gamma == -1:
        gamma_1 = np.median(dist_mat_1[dist_mat_1 > 0])
        gamma_2 = np.median(dist_mat_2[dist_mat_2 > 0])
        K = np.exp(-0.5*((1/gamma_1**2)*dist_mat_1**2 + (1/gamma_2**2)*dist_mat_2**2))
        return K
    K = np.exp(-0.5*((1/gamma**2)*(dist_mat**2)))
    return K


def cos_exp_kernel(x,y,n_freqs = 5,l=1):
    """
    The c-exp kernel
    
    Parameters:
    x,y - inputs 
    n_freqs - number of frequencies to include in the sum
    l- bandwidth of the kernel
    
    Returns:
    Kernel values given x,y
    """
    
    cos_term = np.sum([np.cos(2*np.pi*n*(x-y)) for n in range(n_freqs)])
    return cos_term*np.exp(-(0.5/(l**2))*(x-y)**2)


def CEXP(X,n_freqs = 20,l=np.sqrt(10)):
    """
    Transforms an array of function values using the integral operator induced by the cos-exp kernel. 
    The function values are assumed to be on [0,1]
    
    Parameters:
    X - (n_samples,n_obs) array of function values
    n_freqs - number of frequencies to include in the sum
    l- bandwidth of the kernel
    
    Returns:
    cos_exp_X - (n_samples,n_obs) array of function values where each function has been passed
                through the integral operator induced by the cos-exp kernel
    """
    n_obs = X.shape[1]
    obs_grid = np.linspace(0,1,n_obs)
    T_mat = pairwise_kernels(obs_grid.reshape(-1,1), metric = cos_exp_kernel, n_freqs = n_freqs,l=l)
    cos_exp_X = (1./n_obs)*np.dot(X,T_mat)
    return cos_exp_X


# Testing Functions
# The next three cells contain the functions used to conduct the two sample test.
# MMD_K calculates an empirical MMD quantity
# two_sample_test performs a two-sample test, using MMD_K multiple times
# power_test performs two_sample_test numerous times to calculate an estimate of test power

def MMD_K(K,M,N):
    """
    Calculates the empirical MMD^{2} given a kernel matrix computed from the samples and the sample sizes of each distribution.
    
    Parameters:
    K - kernel matrix of all pairwise kernel values of the two distributions
    M - number of samples from first distribution
    N - number of samples from first distribution
    
    Returns:
    MMDsquared - empirical estimate of MMD^{2}
    """
    
    Kxx = K[:N,:N]
    Kyy = K[N:,N:]
    Kxy = K[:N,N:]
    
    t1 = (1./(M*(M-1)))*np.sum(Kxx - np.diag(np.diagonal(Kxx)))
    t2 = (2./(M*N)) * np.sum(Kxy)
    t3 = (1./(N*(N-1)))* np.sum(Kyy - np.diag(np.diagonal(Kyy)))
    
    MMDsquared = (t1-t2+t3)
    
    return MMDsquared


def two_sample_test(X,Y,gamma,n_perms,z_alpha = 0.05,make_K = K_ID,return_p = False, seed = 1):
    """
    Performs the two sample test and returns an accept or reject statement
    
    Parameters:
    X - (n_samples,n_obs) array of samples from the first distribution 
    Y - (n_samples,n_obs) array of samples from the second distribution 
    gamma - bandwidth for the kernel
    n_perms - number of permutations performed when bootstrapping the null
    z_alpha - rejection threshold of the test
    return_p - option to return the p-value of the test
    make_K - function called to construct the kernel matrix used to compute the empirical MMD
    
    Returns:
    rej - 1 if null rejected, 0 if null accepted
    p-value - p_value of test
    
    """
    
    np.random.seed(seed)
    
    # Number of samples of each distribution is identified and kernel matrix formed
    M = X.shape[0]
    N = Y.shape[0]
    K = make_K(X,Y,gamma = gamma)
    
    # Empirical MMD^{2} calculated
    MMD_test = MMD_K(K,M,N)
    
    # For n_perms repeats the kernel matrix is shuffled and empirical MMD^{2} recomputed
    # to simulate the null
    shuffled_tests = np.zeros(n_perms)
    for i in range(n_perms):
            idx = np.random.permutation(M+N)
            K = K[idx, idx[:, None]]
            shuffled_tests[i] = MMD_K(K,M,N)
    
    # Threshold of the null calculated and test is rejected if empirical MMD^{2} of the data
    # is larger than the threshold
    q = np.quantile(shuffled_tests, 1.0-z_alpha)
    rej = int(MMD_test > q)
    
    if return_p:
        p_value = 1-(percentileofscore(shuffled_tests,MMD_test)/100)
        return rej, p_value
    else:
        return rej    


def power_test(X_samples,Y_samples,gamma,n_tests,n_perms,z_alpha = 0.05,make_K = K_ID,return_p = False):
    """
    Computes multiple two-sample tests and returns the rejection rate
    
    Parameters:
    X_samples - (n_samples*n_tests,n_obs) array of samples from the first distribution 
    Y_samples - (n_samples*n_tests,n_obs) array of samples from the second distribution 
    gamma - bandwidth for the kernel
    n_tests - number of tests to perform
    n_perms - number of permutations performed when bootstrapping the null
    z_alpha - rejection threshold of the test
    make_K - function called to construct the kernel matrix used to compute the empirical MMD
    return_p - option to return the p-value of the test
    
    Returns:
    power - the rate of rejection of the null
    """
    
    # Number of samples of each distribution is identified
    M = int(X_samples.shape[0]/n_tests)
    N = int(Y_samples.shape[0]/n_tests)
    rej = np.zeros(n_tests)
    
    # For each test, extract the data to use and then perform the two-sample test
    for t in tqdm(range(n_tests)):
        X_t = X_samples[t*M:(t+1)*M,:]
        Y_t = Y_samples[t*N:(t+1)*N,:]
        rej[t] = two_sample_test(X_t,Y_t,gamma,n_perms,z_alpha = z_alpha,make_K = make_K,return_p = return_p)
    
    # Compute average rate of rejection
    power = np.mean(rej)
    return power
