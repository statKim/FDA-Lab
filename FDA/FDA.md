# Functional Data Analysis



## 1. What is "Functional Data"?

- **Definition**

  > "Observations on subjects that you can imagine as $Xi_(s_i)​$,
  > where $s_i​$ is continuous"

  $$
  y_i(t) = X_i(t) + \epsilon_i(t) \\
  \epsilon_i : error
  $$

  

  예를 들어, **time-series**에서는 데이터를 seasonality, cyclical component 등을 차분과 같은 방법으로 상쇄하여 stationary한 데이터로 변환하여 분석한다.

  반면, **FDA**에서는 각 time마다 discrete하게 관찰된 데이터를 실수 집합에서 연속이 되게 smoothing 과정을 통해 사이의 값들을 interpolate한다. smoothing된 데이터의 형태를 함수라 가정하여 이를 추정한다. 추정된 함수를 통해 이후에 분석에 사용하게 된다.

- **Functional data 특징**

  - High dimensional
  - Temporal and/or spatial structure (시간 또는 공간의 구조)
  - Interpretability across subject domains (특정 주제를 통한 해석)

- **Functional data의 종류**

  - Dense functional data : 일정한 간격(a fine regular grid)을 가지는 데이터, i.e., $ x_i = ( X_i(\frac{1}{N}), X_i(\frac{2}{N}), ..., X_i(1) ) ​$

    ex) spectral data, imaging data, accelerometry, ...

  - Sparse functional data : 각 observation의 간격이 불규칙(irregular)하고 random인 경우

    ex) CD4(에이즈 바이러스 항원) count, blood pressure, etc.

  ##### => 각 case에서 특정한 interpolation이 필요!!!

- **왜 multivariate techniques (MANOVA, clustering, multiple regression, etc.)을 바로 적용하지 않을까?**

  - functional data techniques은 **데이터의 구조**를 고려해야 함

    > MDA는 순서에 불변(permutation-invariant)하지만,
    >
    > FDA는 그렇지 않다!!

  - FDA 방법론의 발전은 MDA의 연장(extension)이다.



## 2. Functional data to smoothing function

- ##### A basis function system is a set of $K$ <u>known</u> functions $\phi_k(t)$

  - ###### linearly independent of each other

  - ###### can be extended to include any number $K$ in the system

- ##### 함수 $X(t)$는 basis function의  linear combination

$$
X(t) = \sum_{k=1}^{K} c_k \phi_k(t)
\\
X(t) = \bold C^T \Phi(t) \ \ \ (matrix \ form)
$$

- ##### basis function에 따라 smoothing 방법이 달라짐

- ##### 간단히 미분값을 계산할 수 있지만 모든 basis가 미분가능하지는 않다. 



### 2-1. Basis functions

- #### Monomial(Polynomial) Basis : $ X(t) = \sum_{k=1}^{K} c_k t^{k-1} $

  - ###### basis functions are the monomials: $ 1, t, t^2, ... $

  - ###### $K=5$까지는 간단히 계산 가능하지만, 정밀한 지역화된 특징(sharp localized features)을 잡는데 심각한 문제가 있고 unequally spaced data에 대해서 computation 문제가 발생할 수 있다.

  - ###### For a polynomial of degree $m$, the derivative of order $m + 1$ is zero.

- #### Fourier Basis

  - ###### basis functions are sine and cosine functions of increasing frequency: 

    $$
    1, sin(\omega t), cos(\omega t), sin(2 \omega t), cos(2 \omega t), ...,
    \\
    sin(m \omega t), cos(m \omega t), ...
    $$

  - 

- #### Spline Basis

  - ##### B-spline

  - ##### P-spline

- #### Power Basis : $ t^{\lambda_1}, t^{\lambda_2}, t^{\lambda_3}, ... where \ \lambda_i s \ are \ distinct$

- #### Exponential Basis : $ e^{\lambda_1 t}, e^{\lambda_2 t}, e^{\lambda_3 t}, ...  where \ \lambda_i s \ are \ distinct $