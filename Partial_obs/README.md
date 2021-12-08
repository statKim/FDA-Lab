### Functional snippets

- IRLS algorithm -> 문제 없음!!!
  - `MASS::rlm()` 과 `IRLScpp()` 을 비교했을 떄, 결과가 거의 동일
    - 즉, 알고리즘 자체는 맞는 방법이나, 전체적으로 under-estimate되는 경향이 보임...
  - 기존 IRLScpp에서 weight의 분모 계산에서 0이 나옴(residual = 0인 경우가 있었음)
    - 0이 되는 경우는 임의의 작은 수인 0.0001로 대체
    - C++ Eigen으로 해보려고 며칠을 헤매다가, 걍 for 문으로 해결함
  - 기존의 IRLScpp에서 scale updating은 잘못된 방법임
    - initial로 fixed해서 하면 extreme noise에서 잘 작동 안함
    - M-scale (Bisquare)로 수정 (iterative method)
      - 일단 너무 느리다... 중첩 for문이 너무 길어짐
        - pmax?? pmin?? 그거로 책에 있는 exact form의 element-wise min으로 계산하자
- Outlier setting
  - Extreme value distribution에서 generate??
  - 여기서는 snippet이라 extreme outlier여도 될듯 (현재 out_type=2 처럼)
- Real data
  - CD4 데이터 (Boente(2021) 참고)
    - snippet이라고 해도 될지 애매... sparse 케이스이긴 함
      - 전체 domain으로 보면, long-snippet 형태라 적용할 수 있긴 할듯
    - matern corr structure를 따를지를 모르겠음
    - 문제점...
      - noise variance의 bw가 너무 크게 나와서 안돌아감...
      - log변환한 다음에 하면, noise가 NaN이 나옴...
  - Spinal bone mineral density
    - Lin & Wang(2020)에서는 raw cov를 그려보니까 diagonal과 인접한 부분부터 급격하게 0으로 떨어지기 때문에 matern covariance structure를 사용함
    - outlier가 없음...
- 결과 요약
  - Eigenfunction, var, col + intra, extra??
    - 근데 cov 요약은 좀 애매해보임.... 아예 true PC로부터 generate해서 PC를 얼마나 잘 estimate하는지 해볼까?
  - Eigen, PC score
    - 30번 한 결과를 보니까, var, cov 요약하지 말고 PCA 결과 위주로 요약하는 것이 좋아보임
  - ROBPCA 논문에서는 eigenvalue로도 요약함 (maxsub, PVE, MSE)
- Variance smoothing이 계속 over-smoothing됨 (bw가 너무 큰 것이 선택됨)
- cov structure가 matern corr를 만족하는지 확인하는 test 같은 것이 있나??
- 시뮬레이션 세팅을 score로부터 generate하되, cov형태가 matern corr 형태를 따르는 경우로 찾아보기
  - Delaigle 세팅 model 2를 사용해서 score를 Normal에서 generate함
    - outlier는 t-dist에서 generate?? 아니면 score 자체를 t-dist로??
  - Boente(2021) 논문 세팅 model 2 확인!!!
  - 근데 cov structure도 matern 말고 fourier(Lin 논문)나 다른 corr structure도 적용해서 해봐야할듯
- Boente에서 bw 작을 때 에러나는 이유가 epan kernel의 support 때문이었음(내가 짠 거나 mcfda는 따로 안해줌)
  - gauss 추가해보자!!
    - 근데 이거 하려면, CV 함수랑 모두 바꿔줘야됨....
    - 그냥 옵션에서 k.epan 함수를 다르게 지정하도록 함 
      - kern=="gauss"일 떄는 k.epan 함수를 guess kernel로 지정 (전역변수로 다시 지정)
- CV 해결하기 => 현재는 너무 over-smoothing 되는 경향이 있음...
  - curve-out CV 말고 그냥 obs-wise CV 하면 잘 working하는 거 같음 (mcfda에서도 이렇게 함, 근데 논문에서는 아님)
- IRLS tolerence 계산식 수정?? 지금은 absolute sum으로 되어있음 (근데 크게 중요하지는 않은듯)
- bisquare일 떄 bw_loss 확인 => 지금 Huber로 돌아간거 같음
- Biweight loss 추가해보기 (CV 포함)
  - Huber와 크게 다르지 않음
- Lin paper 처럼 diagonal part, off-diagonal part 차이 계산해보기
- completion은 의미가 좀 없어보이고, clustering이나 classification으로??
- irregular grid로 바꾸기
- 우선 fixed bw로 결과 좀 확인해보기 (Yao랑 Boente가 계속 말썽임...)
  - epan kernel로는 에러 많이 떴었음!! Gauss kernel로 해보자!!
  - 근데 Boente는 epan만 되게 해놨음
- IRLS.cpp에 locpolysmooth에서 maxit, tol 옵션 추가하고 maxit = 50으로 수정
- `robfpca`에서 delta = 1.345로 fixed하고, bw = range/5로 fixed 하기
- simulation setting을 irregular(`synfd`)와 regular(Delaigle)로 하면 될듯

### Robust covariance estimation for partially observed functional data

- Theory

    - Boente et al.(2014)는 curve 자체를 elliptical process로부터 나왔다고 정의
    - 우리 페이퍼는 각 curve의 element가 elliptical process로부터 나왔다고 정의
    - 그래서 나중에 covariance를 element-wise estimate하고 얘의 consistency 보일 수 있기 때문

- tdist의 경우는 pre-smoothiing하고 해도 되지 않을까??

- t distribution with df 3인 경우의 M-estimator(MLE) 구하는 코드 짜기

    - `MASS::fitdistr()`

- OGK에 hard rejection 추가했는데 Boente에서만 이상하고, Delaigle에서는 막 이상하게 나오지는 않음

    - 근데 여전히 eigenfunction은 별로임

- Boente 세팅에서 yao 제외한 나머지 방법들이 2nd PC var > 1st PC var

    - 또한, PC function 순서는 또 나머지들이 맞는거 같고 yao가 틀린거 같음

- OGK 문제

    - cov_ogk에서 trim 옵션은 현재 사용하면 안됨 (OGK 과정에서 M으로만 돌도록 되어있음)

    - PCA 과정에서 missing이 있다보니 score의 variance가 매우 이상 (1st PC var이 가장 커야되는데 아님...)

        ```
         [1]  41.603410 198.423465 106.069488  11.904357  23.407688  22.137286 ...
        ```

        - 왜냐하면 모든 obs를 이용해야하는데, NA가 있어서 그게 안됨 => 1st score의 값과 variation이 너무 작아짐

- GK에서 missing 없는 데이터로 해도 corr > 1인 경우가 생김

    - 다른 식을 사용해서 해결함

- OGK 다시 짜보기
    - Missing 제외해서 계산하는 방법으로 적용 (따라서, 중간 과정에서의 정확한 PC score 계산이 불가능...)
    - 근데 현재, gk에서 cor > 1 인 문제가 발생하고 이를 사용한 결과로 완전히 확신하기 어려움...

- GK에서의 문제
    - NA 제외해서 하니까 cor이 1이 넘는 경우 발생... (full data로 하면 max가 정확히 1이 나옴)
        - `is.na(z1)`을 하면서 obs 몇개가 빠짐으로 인해, corr=1의 가정이 무너짐
        - trimmed_sd로는 고쳤으나, M-estimator로는 안됨.... (알고보니 trim도 안됨... )
            - extreme한 값 때문에 dispersion estimator가 크게 추정되서 1이 넘는거 같음...
            - non-robust sd로 하면 max = 1
    - 이건 psd 만족 안하기 때문에, 옵션이 필요할듯... (pm은 psd 만족해서 함수에서 psd 옵션 제거함)
    - Boente model 4에서 MAD=0이 되어버림...

- PM에서 NA 제외하고 cov 계산하더라도 psd랑 affine equivariance 만족하는거 보일 수 있나??
    - imputation된 상태에서 cov 계산하면 psd, affine equivaraince 모두 확보됨
    - 근데 NA 존재하면 정의 자체가 불가능한듯...

- PM에서 location, dispersion estimator 방법 결정하기
    - 이건 너무 튀는 방법만 고르지 않으면 크게 상관은 없을듯

- Score 또는 score distance로 adjbox 해보기 / robMah 구체적으로 찾아보기
    - `foutliers`에서 `robMah` 방법은 first score 2개를 가지고 함 (즉, completion이 필요없음)
        - 현재 `fun_outlier()` 함수 새로 짜서 PCA에서 select된 개수만큼 사용(중간 단계는 foutlier의 코드 가져옴)
        - 이건 한 번 비교해보자!!
    - 또한, chisquare cutoff 쓰는 것이 적절해보이지 않음 (df = p 인데 functional data에서는 infinite임...)
        - 근데 Hyndman(2010) bagplot paper에서는 dense grid 개수로 cut-off하는 것이 언급되어 있음

- Raymaekers(2021) 방법론에서 missing 처리 부분을 바꿔보기?
    - imputation하지 말고 NA 그대로 두고, missing 제외하고 cov 계산
        - 평균 imputation보다 결과 좋음 (Boente에서는 압도적으로 좋음!!)
        - 근데 평균으로 대체하면 psd 되지만, 다른 경우에는 만족이 안될듯...(psd 정의 생각해보면 cov는 항상 psd)
    - missing 제외 distance 가까운 것들 몇개로 평균?
        - 결과가 좋긴 하지만, 기준을 결정하기가 애매....(몇 개를 사용할지? 비율?)
        - 이걸 잘 정리하면 contribution이 될 수 있을까??

- PM10 eigenfunction, outlier detection 해보기
    - smoothing과정에서 knots selection이 중요!! => Xiao(2013), Remark 3
        - 논문에서 $p/2$ 를 recommend함

- Yao paper처럼 noise 계산해보기(diagonal은 가까운 값으로 채워서 평균 하든가 해서 하기)

- bandwidth=0.2랑 0.3 차이가 나는 것은 `.Machine$double.eps` 때문이었음

- Trimmed noise variance가 bandwidth에 따라 값이 매우 달라짐...
    - L-estimate
    - 75% quantile?? Boxplot 계열의 outlier detection??
    - adjbox나 boxplot은 bimodal 형태처럼 나오는 것을 잘 잡지 못함...

- Boente(2015) 처럼 oultier detect하고 completion 해보기 (score로 density 확인??)
    - 근데 문제는 Yao가 어중간한 것이 아니라 그냥 전체적으로 나쁨...(즉, robust 방법들도 안 좋은 것은 얘도 안좋음)
    - adjbox는 bi-modal인 경우를 잘 detect 못함...
    - ROBPCA 논문의 score, orthogonal distance에 대한 cut-off로 outlier detection(PC score based 방법)
      - 근데 이 경우는 score가 normally distributed 가정을 했기 때문에 chi-squared cut-off를 사용함
      - 만약 heavy tail dist(t or Elliptical dist) 사용하게 되면 cut-off 달라져야 할듯
    - completion한 후에 functional data outlier detection도 적용
      - `rainbow` 패키지 함수를 사용했고, Mah와 HU를 제외한 나머지 방법은 outlier가 없는 것으로 결과가 나옴...

- Boente(2015) 세팅 시뮬레이션
    - M-est(sm)가 eigenvalue, vector가 complex가 나옴
        - 우선 imagenary part가 0에 가까운 경우에 한해 real part만 뽑아서 사용함
    - Model 1 : M-est, M-est(smooth)가 좋게 나옴 (noise var이 너무 커서 오히려 안좋음 => Huber에서 var 모두 0 나옴)
    - Model 2 : Kraus-M이 압도적으로 좋음 (noise var 매우매우 큼)
    - Model 3 : noise var 작음, noise 빼준 것이 매우 좋음 (안빼준 것은 매우 나쁨)
    - 뭔가 diagonal part를 adjust 해주어야 더 좋게 나오는듯... 문제는 그 이유를 모르겠음
    - Model 1, 2는 measurement error 더해줬지만 딱히 얘를 추정하지는 않음

- outlier 아닌 것만 completion하지 말고, spike 추가되기 전을 잘 맞추는지 확인??

- mean smoothing spline으로 하는거 penalized spline으로 수정하기(하나로 통일)
    - 둘의 차이는 Knots의 개수인데, smoothing spline에서는 모든 데이터 point가 knots이라 가정함

- 시뮬레이션 t dist로 heavy tail로 생성한 다음에 해보기?? 별로임....

- Under normal, PACE be BLUP, Under Elliptical, Boente(2020) derive the equation

- time series처럼 뾰족한 데이터에 적용해보기

- 지금 방법은 오히려 형태만 달라지는 경우는 매우 못맞힘... 이때는 Huber가 매우 좋음(이건 아마 structure 가정 때문인듯)

- 현재 sandwich smoothing한 후에 noise var 빼면 negative 나옴... (smoothing 전에 빼주고 psd 까지 해주기?)

- Kraus 세팅에서 seed=23일 때 아래 에러 발생 (complex eigenvalue 발생) 

    - `invalid comparison with complex values`

- local을 kernel weight 줘서 한 번 해보기

    - 식으로 적어보니 Nadaraya-watson M-estimator와 exactly same

- noise var trimmed 말고 winsorization으로 바꿔서 해보기

    - 근데 이건 문제가 $A_0$에서 winsorize한 index를 $A_1$에 어떻게 적용하지??

- `MASS::huber`를 `robustbase::huberM`으로 바꾸기

- noise var에서의 trimmed하는 기준을 outlying measure로 일정 curve를 제외하는 것으로??

- Covariance smoothing
    1. M-est + 2D smoothing
        - Eigenvalue가 왜 천천히 작아질까??
            - smooth해질수록 PVE가 점점 커짐
            - 아마 smooth할수록 low rank approximation이 가능해져서 그런가??
            - smooth할수록 eigenfunction이 True와 더 비슷해짐
            - 메인을 하나 정해야 하나?? (PCA를 잘할지? 아니면 cov estimate을 잘할지??)
            - 내가 생각하기에는 surface에 spike가 있으면 그 축에서의 variance가 커지고, 1st eigen direction의 분산과 2nd direction의 분산이 큰 차이가 안나게 됨... 그렇기  때문에 eigenvalue가 빠르게 감소하지 않는듯
              - 그럼 smoothing을 해주긴 해야하는데 지금처럼 kernel smoothing이나 smoothing spline은 적용 못할듯...(CV에서 점점 값이 축소되는데 그럼 당연히 CV error는 커질 수 밖에 없음)
        - 밑의 Yao처럼 diagonal 제외하고 해봤더니 안좋네...ㅠㅠ 다른 방법 찾아보기
            - diagonal 제외하고 해보되, bw 좀 늘려서 괜찮은지 테스트만 해보기(raw cov를 M-est라고 가정하고 해보자!!)
            - 거기에 cardinality 곱해서도 해보기(보니까 Yao 페이퍼에는 average가 아닌 곱한 term 사용)
        - 이 방법의 문제점
            - 각 grid point마다 값이 1개씩만 존재하기 때문에 이 상태에서 smoothing을 하게되면 전체적으로 값이 down scale되게 됨 (이건 다른 smoothing에서도 그렇긴 하지만...)
            - 또한 값이 1개씩 밖에 없기 때문에 CV를 하기가 어려움
        
    2. M-est + Eigenfunction smoothing
        - Eigenfunction smoothing spline 적용하고 Gram-schmit로 orthogonalization하는 것 고려
        - `stats::smooth.spline()` 함수가 CV or GCV 기준으로 smoothing para 결정해줌 
        - 근데 이거 그램슈미츠가 안돌아감... (Not full-rank)
            - 아마 first eigenvector들은 잘 smoothing이 될텐데 뒤쪽으로 갈 수록 noise형태로 되다보니 smoothing하면 대부분 한 직선으로 되는 경우가 생기는듯
        
    3. Raw cov + Robust 2D smoothing
        - 형태가 엉망인듯 하고 eigenfunction도 1st가 매우 크게 벗어남
        
    4. Raw cov + Local M-est

        - `get_raw_cov_cpp`에서 같은 timepoint 반복 안하게 설정하기

        - Cov는 그냥 그런 것 같긴한데 completion은 괜찮은편
        - 근데 이 경우에는 지금 diagonal 제외하고 local 방식으로 계산하기 떄문에 어쩌면 variance가 음수가 될 수도 있을듯... 이건 한 번 확인해보자!!
        - seed=31(Delaigle 세팅), 35(Kraus 세팅) 일 떄, Kraus(sm)이 엄청 커짐....
          - 근데 psd로 바꿔준 후에는 안그럼
        - 근데 이건 결국 local regrssion과 비교해봐야할 듯... (대부분 regression type을 쓰는 데에는 이유가 있을 것임...)
          - 단순 local이 의미가 있나?? 이걸 improvement한 것이 Nadaraya-Watson (Local constant)나 Local polynomial regression인데.....
        - CV로 했을 때가 bw=0.05로 fix했을 때에 비해 훨씬 안좋음...

- Local smoothing이랑 M+smooth.spline은 거의 비슷함

- `.Machine$double.eps`로 아주 작은 값은 0으로 처리할 것 (`get_eigen()`)

    - 어쩌면 이게 completion에서의 극단적인 spike를 보여주는 것일수도??

- Yao et al.(2005) covariance estimation 과정
    - cov estimation에서 raw cov 계산하고, diagonal 제외한 상태로 2D smoothing해서 cov estimation
    - 그리고 diagonal term만 따로 1D smoothing해서 $G(t,t) + \sigma^2$ obtain
    - 여기서 둘을 빼어서 $\hat\sigma^2$ 계산

- Assumptions
    - observed on regular girds but partially observed
    - few curves are contaminated by extreme spikes

- Fisher consistency vs Consistency
    - Fisher consistency : design(or method) consistency라고도 불리며 sample version의 estimator 대신 true value(population으로 부터 얻어지는 실제 값)을 넣었을 때, 실제 estimator가 obtain되면 이를 fisher consistency라고 함 (따라서 asymptotic consistency와는 개념이 다름)
    - https://stats.stackexchange.com/questions/67879/sample-variance-fisher-consistency
    - 즉, variance의 MLE는 Fisher consistent but unbiased estimator는 Not Fisher consistent
    - 그런데 둘 다 (method) consistent

- Theory  전개
    - PACE 식에서 Gaussian 가정하지 말고, Elliptical 가정?? (Heavy tail dist이기 때문)
    - Consistency
        - mean, cov : Kraus(2015)
        - PC score : Yao(2005) under Gaussian assumption, Boente(2020) under Elliptical dist
    - Identifiable?? 즉, unique한 솔루션이 존재하냐는 의미!! => canonical form이 존재 (예를 들어 X_bar처럼~~)

- noise 먼저 빼고 psd로 만드는 거랑 psd 먼저 하고 noise 빼는 거랑 별 차이 없음

- funPCA vs FPCA 비교해보기
    - off-diagonal은 똑같고 diagonal에서 아주 조금 차이 나는데 왜 inverse 구하면 값이 엄청나게 커지지???
        - ginv로 해볼까??
        - 지금 robfpca 패키지에 ginv로 수정해놓음... (혹시 이상하면 나중에 수정할 것)
        - Singular 문제를 ginv로 하고 noise var = 0으로 두고 해도 결과는 비슷
    - fitted cov 계산해서 하니까 문제 없는듯??
        - 근데 이렇게 psd로 만들어주는 방식이 일반적인 방법인가??
    - 그럼 M-est의 경우는 2번하는 셈인데, 1번만 하는게 맞는건가??

- Delaigle setting에서 True eigenfunction도 알고있음!!! (즉, True score 개수도 알고 있음!!)

- elliptical distribution이 normal, t-dist 등을 포함하고 있는 더 포괄적인 개념

- 근데 애초에 적은 비율만 extreme noise에 영향을 받기 때문에, noise var을 가정하지 않아도 될듯??

- `fields::smooth.2d`에서 bw 옵션인 theta가 최신버전에서는 aRange로 변경됨
    - 이거 CV를 어떻게 해야하지??

- 시뮬레이션 해보기
    - seed=6일 때 스파이크 너무 큼...
    - 3가지 경우 (partially obs, snippet, sparse)
        - partially obs는 결과 있으니까, 나머지 2가지 해보기
        - Partially obs 일때는 epanichnikov kernel
        - snippets 일때는 gauss kernel (근데 M-est는 계산이 불가능함...)
    - 비교 방법론
        - Yao: sparse에 적합 but non-robust
        - Huber: Lin(snippet에 적합) + robust
        - Kraus
        - Kraus-M
        - Boente
        - M-est
        - M-est-noise
        - M-est(sm)
        - M-est(sm)-noise
    - 전체 reconstruction error랑 eigenfunction error 계산
        - 이거는 함수로 만들기(애초에 matrix로 받아서 나중에 label 추가해서 보기 좋게 만드는 거로)

- 2단계로 robust하게 바꾸기 (Cov estimation => PCA)
    - PACE를 robust하게 바꾸기

- dopper처럼 너무 wiggly한 것 보다는 덜 fluctuate한 경우로 시뮬레이션 하는 것이 좋아보임

- Lounici 방법이 robust한 경향이 있어서(cauchy나 t dist contamination), PACE를 robust하게 하는 방법을 좀 찾아보기
    - Projection-pursuit : variance가 아닌 robust scale estimator를 maximize하는 PC score를 obtain

- Lounici (2014) Bernoulli paper 참고
    - Oracle property 3가지?? (Lasso paper인가 sparse matrix paper인가에 나와있음)
    - sparse covariance estimation (lasso penalty 넣어서)
        => 이것도 reference 찾아보기
        나중에 시뮬레이션에서도 sparse한 경우, 아닌 경우 2가지로 할 수도??
    - eqation (1.2)에서 \Sigma_n에 지금 M-estimator로 구한 것 넣어서 구하는 방법 (또는 absolute loss??)
        => 이거 estimate하는 방법이나 알고리즘 더 찾아보기
    - eq (1.4)에서 \Sgima_n^\delta에 M-est로 바꾼거 넣고 eq (1.5)로 최종 estimate
    - eq (1.5) proof 1040쪽에 exact 솔루션 있음 (convex여서 그런듯??)

- Oracle
    - Oracle estimator : variable의 영향이 없는, 즉, 어떤 coef=0인지를 알고 있는 estimator (True sparsity를 알고 있는 estimator)
    - Oracle property [참고](https://www.math.arizona.edu/~hzhang/math574m/2020F_Lect11_shrink_lasso.pdf) 
        - selection consistency : true parameter space가 아닌 경우의 estimation = 0, else != 0
        - estimation consistency : asymptotic normality (convergence in distribution)
        - Weak : Oracle estimator의 coef=0인 set이 estimator의 regualarization para에 따라 달라지는 set에 포함될 확률이 1로 수렴
        - Strong : Oracle estimator의 coef=0인 것과 estimator가 같을 확률이 1로 수렴
    - Oracle inequality : 예를 들어, LASSO는 Oracle property를 가지는데, LASSO estimator와 Oracle estimator의 차이가 bounded

- 애초에 PACE에서도 noise var 보정(?)을 해주기 때문에, 이론적으로 M-est-noise를 사용하는 것이 맞음
   단, Kraus일 때가 문제긴한데....

- robust FPCA 논문 찾아보기 (pre-smoothing을 하는지 ref 찾아보기!!!)

- Trimming은 5%만큼의 데이터를 제외하지만, Winsorization은 5%의 outlier들을 5, 95th quantile 값으로 대체 

- elliptical distribution : MVN에 비해 heavy tail을 가지거나 매우 좁은 tail을 가진 symmetric distribution (MVN은 elliptical dist의 special case!!!)

- 결과요약
    - dense한 경우에만 돌아가는 것들이 있으니까, 이거랑 partially observed인 경우랑 둘 다 요약
    - Cov estimation

- Delaigle 세팅의 경우, noise를 가정하지 않았기 때문에 pre-smooth를 안해도 될듯??
    - 근데 outlier는 noise 있는데 smoothing해야됨??

- 시행착오!!
    - raw cov 계산하고, 이후에 robust smoothing함으로써 noise 제거
    - score에 noise가 있다고 가정??

- Kraus의 경우, pre-smooth를 해서 noise 제거(penalized spline)

- noise var 뺀 값이 0이 되는 바깥쪽은 우선 그냥 제일 근접한 값으로 대체함

- huber loss delta가 큰 영향을 주는지 시뮬레이션 ISE로 확인하기
        => 왜 그런지는 모르겠으나, 1.345로 fix했을 때, 더 좋음...
        => bandwidth는 CV하는 것이 좋은 경우도 있음 (이건 하는 것이 좋아보임)

- snippets보다 partially observed로 방법 찾아보기(cov est, PCA)
    - Lin은 데이터 더 많아져도 비슷함... 즉, 데이터가 많아져도 그 정보를 활용 잘 못하는듯...
    - PCA도 PACE 말고 정보를 더 활용할 수 있는 방법 찾아보기

- Boente(2020) sparseFPCA 시뮬레이션 세팅으로 해보기

- data adaptive method(지금 방법론은 너무 밋밋함...)

- 근데 variance를 robust하게 estimate하면, 사실상 noise를 무시하게 되는거니까 굳이 estimate 안해도 되지 않나??
        => 애초에 outlier가 껴있어서 크게 estimate되는 것이 당연한데, robust하게 구하는 것에 이것이 필요할까??

- Clustering simulation
    - normal noise 대신 다른 noise 줘서 outlier로 정의하고 그럼에도 불구하고 얘네들이 true로 잘 clustering 되는지?
    - outlier 없을 때, Kraus가 굉장히 좋음.....(다 맞히는 경우도 많음)
    - outlier 있을 때, 모두 안좋음
        => robust하게 score를 얻었기 때문에, outlier detection problem이 되어버림. 따라서, clustering이 잘 안됨
        => 그러면, 1st stage에서 outlier detect를 하고, 2nd stage에서 clustering?? (trimmed k-means??)
    - completion 후에는 대부분 경우가 CCR = 1이 되어버림....

- Application
    - Completion
        - 표 요약은 1st curve만 하던지? 아니면 outlier 아닌 경우만 값으로 요약
        - 앞뒤만 missing인 경우로 바꿔보기
        - 단순한 curve 위주로 해보기(변곡점 2개인 경우 최대한 줄여서...)
    - AMI clustering
        - missing이 있더라도 비슷한 형태의 curve들과 같이 묶임
        - 3년치 다 쓰지 말고 1년치만 가져다가 쓰고, 우선 3일치 평균 => 결과 괜찮으면 1일치로 쓰기
    - AMI clustering after completion
        - imputation을 mean(x_i)랑 completion이랑 비교
        - kCFC도 적용가능

- outlier 만들 때, 기존 generate된 값에 추가로 더해주는 것으로 하기

- Robust covariance for functional snippets
    - title은 snippet으로 가고 simulation에서 missing 정도 다양하게 해서 결과 요약
        - missing 정도를 snippets, fragments, partially observed 정도로 나누어서 모두 simulation하기
        - outlyingness도 다양하게 해서 표 요약하기
	- 확실히 variance가 under-estimate 되다보니 양 끝 부분에서의 prediction이 좋지 않음...
            - noise variance
                - outlier 들어가게되면 A_1 term이 커지게 되니까, 이 값이 0.75 quantile보다 커지는 것들 제외한 값들로만 계산 
                - 25% trimmed mean(Winsorized mean)으로 바꾸는거 어떻게 생각함?? 값 자체는 크게 차이나지 않음
                    => 이전 Lin setting 다시 해보기 (원래는 너무 크게 estimate되어서 var 모두 0로 되었었음...)
	- robust smoothing
		- 다른 robust smoothing이랑 비슷비슷함 (LOWESS, WRM)
		- Huber랑 biweight 등 다른 M-type도 같이 시뮬레이션하기
		- 이거를 바꾸기보다는 PCA 부분을 더 손보기???
		- natural cubic spline처럼 boundary 부분에서 linear하게
	- robfpca 패키지
		- error = F or noise = F 추가하기
			- varfunc에서는 sig2 = 0
			- covfunc에서는 sig2x = 0
		- Lt, LY <=> matrix 이렇게 변환해주는 함수 만들기
		- simulation data 코드 추가
		- CV 코드 이상... parallel 제대로 안돌아가고 있음...
			- 특히, wrapper function의 경우 에러가 발생
			- @import Rcpp 하니까 제대로 돌아감~~
		- package에서 print되는 부분 cat으로 다 바꾸기
		- match.arg 써서 조건 처리하기
		- https://adv-r.hadley.nz/rcpp.html
		- C++ 코드 패키지에서 build될 떄, inline 부분때문에 계속 에러 떴었음...
	- IRLS에서 scale parameter도 interative하게 update해줌
	- 다른 robust covariance method (비교방법론)
		- LOWESS
			- Huber 쓴거랑 거의 비슷
			- Cleveland, W. S. (1979). Robust locally weighted regression and smoothing scatterplots. JASA
		- WRM
			- WRM.cpp에서 Line 88 부분에 xdat != xdat[i] 추가함
			- 분모가 0이 되어서 Inf 또는 NaN이 나오게 되어 이 부분을 없애기 위함
			- 하지만 robfilter 안의 함수와 약간 값이 달라지게 되는 문제 발생
			- 그림 그려보니까 거의 비슷함
		- Boente(2020)
			- robust covariance estimation + reconstruction
			- robust loss 적용한 bivariate smoother로 cov estimate
			- 뭐 비슷해보이긴 한데 Yao 논문이랑은 다르다고 함		
	- Completion
		- curve imputation; 
		- Kraus setting에서 covaraince가 positive하게 나오도록 수정해서 해야할 듯
		- cov가 잘 estimate되면, reconstruction에서 interpolation은 매우 좋음(PVE 0.99) But, completion은 그닥인 것처럼 보임...
		- curve alignment
			- 기울기를 맞추고 평행이동
			- A_i에서 B_i까지의 변환하는 비율을 linear한 변화량을 가지도록 함
	- AMI data clustering : 일단 보류...
		- 데이터 잘 확인해봐야됨...(NA, 0 모두 포함되어있음)
		- pre-smoothing을 할까? 너무 spike가 많기도 하고 obs 너무 많은듯..
			- 일별 평균 내던가 뭔가 처리를 먼저 하고 해야할 듯...
		- Yao, Lin, Huber, Kraus, 기존 kCFC 비교
		- 단, kCFC는 missing 포함된 데이터는 제외함
		- 보통 real-data에서는 pre-smoothing을 하기도 하나?? (물론, method에서는 당연히 안하지만...)
		- 현재 dimension 매우 큼... 우선 3일치로 변환해서 했다고 함
			- 일별 변환 (na.rm = TRUE)
			- 3일치로 변환 (na.rm = FALSE)
		- gr = seq(0, 1, length.out = 개수)로 변환하고(어차피 equal grid니까), 이 grid에서 PCA 되도록 설정하기
	- visualization function 만들기
		- sample trajectories
		- var, cov, etc
		- cov surface ggplot으로 해보기
		https://homepage.divms.uiowa.edu/~luke/classes/STAT4580/threenum.html