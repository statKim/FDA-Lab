## k-Centres Riemannian Functional Clustering (kCRFC)

### 2022.01.20(목)


- Dai(2018)의 RFPCA + kCFC 해서 clustering
- sphere의 경우에는 intrinsic과 extrinsic이 동일해서 RFPCA (Dai, 2018)를 써도 괜찮을듯
- 지금 fdapace에 있는 kCFC 함수는 ith obs 빼서 ISE 계산 안하도록 되어있음
  - 이대로 하면 매우 느리고, 결과도 거의 비슷함 (그래서 function에서 제외하면서 하지 않고 한 번에 한듯...)
- RFPCA가 더 좋게 나오는 세팅을 찾긴 했는데, 자꾸 몇 가지 seed에서 매우 안좋음
  - FVEthreshold를 늘릴지 줄일지???
  - 뭔가 PC가 많이 뽑히면 CCR이 이상해지는 경우가 생기는듯 (근데 딱히 이것도 아님...)
- clustering Riemannian manifold에 대한 논문에서 시뮬레이션 세팅 찾아보기
- RFPCA의 경우, sparse하지 않으면 bw 조절이 필요 없음 (smoothing 안함)
- 시물레이션 데이터 세팅
  1. sign changed + lambda 다르게
  2. cos flucuate + lambda 다르게
  3. sin fluctuate + lambda 다르게

### 2022.01.27(목)

- Real 데이터 찾아보기
  
  - Airline trajectory data
    - paper에 나와있는 링크는 유료인 것 같아서 아마 못쓸듯... 
      - AOS published version paper 뒤쪽에 Acknowledgements에 데이터 사용허가 받았다고 함
    - TrajAir
      - https://theairlab.org/trajair/#additional-info
      - https://kilthub.cmu.edu/articles/dataset/TrajAir_A_General_Aviation_Trajectory_Dataset/14866251
      - 근데 관측 범위가 너무 좁은듯...
    - Kaggle 데이터
      - https://www.kaggle.com/open-flights/flight-route-database
      - 시간별 데이터가 아닌듯...
    - OpenSky Network
      - 데이터 받을 수 있긴 한데, 변수 설명이 잘 안되어있음
        - https://opensky-network.org/data/impala
      - https://opensky-network.org/datasets/
      - https://opensky-network.org/datasets/flights_data5_develop_sample/
    - traffic library in Python
      - https://traffic-viz.github.io/index.html
      - OpenSky Network 에서 받은 데이터를 사용하기 편하게 정리(?)해둔 라이브러리
      - 근데 라이브러리 설치가 안됨...
    - Aircraft Localization Competition
      - 데이터 다운 : https://www.aicrowd.com/challenges/cyd-campus-aircraft-localization-competition
      - 데이터가 OpenSky Network 기반이라고 함
      - 근데 여기서 timeAtServer 변수가 비행중 시간이 아닌 것 같음...
  - Bird migration data
    - `MigConnectivity` R package
    - https://methodsblog.com/2017/11/09/migratory-connectivity/
    - 근데 이건 시간별 관측이 아니라, 시작점이랑 끝점만 있는듯
  - 태풍 경로 데이터
    - application to asian typhoon trajectories 로 하면 괜찮을듯??
    - 데이터 다운로드 : https://www.ncei.noaa.gov/products/international-best-track-archive?name=gisSLD
    - 태풍 경로로 clustering한 paper들이 꽤 있어서 읽어보면 좋을듯
- kCFC property를 generalize할 수 있을지 찾아보기
  - 우선 paper의 theory 파트 잘 읽어보기
  - metric을 Riemannian metric으로 바꾸면 될까??
  - (C1)이 Assumption 1, (C2)가 Thm 1-(a), (C3)가 Thm 1-(b)
  - Thm 1이 의미하는 것은 non-idntifiability condition이 hold하면, 각 cluster component로 predicte된 curve들을 구분할 수 없다.
  - 즉, 다시 말해, identifiability condition이 hold하면, 두 L2 distance 차이가 0으로 bounded below. (즉, 항상 구분이 가능)
- CCR에서 튀는 것들 잡을 수 있는 방법 찾아보기
  - generate하는 PC 개수 줄이기??
  - 대체로 spike 튀는 경우는 보통 4개 방법 모두 안좋은 경우임...
  
    - 애초에 k-means는 항상 나쁘기 때문에 se가 작은 것임
  - funclust의 경우, C++에서 seed를 결정하다보니 결과가 계속 바뀜...

    - 이건 아예 함수 내부에서 세팅하지 않으면 계속 이럴듯...
    - 또 R에서의 seed랑 C++에서의 seed가 다를 수 있음
- 시뮬레이션 세팅 관련 reference 찾아보기
  - Dai(2022) paper의 Simulation studies에 2개 그룹으로 나눠서 generate한듯??
    - 이 페이퍼의 real data도 한 번 읽어보기
      - taxi demands 데이터인데, 아래 링크에 있음
      - https://github.com/toddwschneider/nyc-taxi-data
    - 근데 이건 내가 하고 있는 케이스와는 다른 경우인듯... 여기서는 Hilbert sphere
      - 즉, function $f \in L_2$에서 $\lVert f \rVert \le 1$를 만족하는 function들을 모아놓은 것을 의미
      - Hilbert space가 infinite space이기 때문에, 이러한 function의 family는 infinite dimensinal sphere이고, $S^\infty \subset \mathcal{H}$ 의 형태로 표현가능
- `RFPCA` 패키지에서 Euclidean 옵션이 MFPCA가 맞는지 확인이 필요
  
  - Multivariate FPCA와 동일!!
- 비교 방법론 더 찾아보기 (FPCA 기반이나 multivariate functional data에 적용 가능한 방법)
  - 아래 방법들은 교수님 paper에서도 비교한 방법론들임!!

  - funHDDC 패키지
    - Schmutz, A., Jacques, J., Bouveyron, C., Cheze, L., & Martin, P. (2020). Clustering multivariate functional data in group-specific functional subspaces. *Computational Statistics*, 1-31.

  - Funclustering 패키지
    - Jacques, J., & Preda, C. (2014). Model-based clustering for multivariate functional data. *Computational Statistics & Data Analysis*, *71*, 92-106.


### 2022.02.07(월)

- OpenSky Network 리얼 데이터 찾아보기
  - 가입하고 https://opensky-network.org/data/apply 에서 access application 제출해야됨
  - 근데 abstract 같은 것 300글자 이상 쓰라고 해서 일단 안함
- 태풍 데이터로 clustering 결과 뽑아보기
  - kCFC(M)은 제외해도 될듯 (어차피 proposed도 아니고 기존에 있던 방법도 아님)
- Multivariate functional clustering 최근 방법론 찾아보기
- 방법론 한 번 쭉 정리하기 (kCFC 책 참고)
  - notation이랑 metric 쭉 정리하고, Thm 1을 manifold 버전으로 수정해보기!!
- 시뮬레이션 세팅은 나중에 리얼데이터의 meanfunction 형태를 보고 이와 비슷하게 generate해도 좋을듯
  - spike 튀는 경우 해결할 수 있으면 해볼것!!
    - m = 51, K = 5로 조절해봐도 크게 다르지 않음...