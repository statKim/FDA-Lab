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
  - 근데 abstract 같은 것 300 단어 이상 쓰라고 해서 일단 안함
- 태풍 데이터로 clustering 결과 뽑아보기
  - 2000~2017년 442개 태풍 가지고 테스트해봄
    - 참고논문에서는 같은 기간에 432개라고 하는데 일단은 무시
    - Misumi(2019) Multivariate functional clustering and its application to typhoon data
  - kCFC(M)은 제외해도 될듯 (어차피 proposed도 아니고 기존에 있던 방법도 아님)
  - 사실상 위치 데이터만을 가지고 clustering하는 거라서, 전체적인 trajectory의 코스가 비슷한 것들끼리 묶여야할 것 같음
    - 이 관점에서 kCFC(R)의 결과가 별로임... 오히려 initial clustering이 더 좋아보임
    - iteration 과정이 오히려 결과를 나쁘게 하는 것 같기도...
      - geodesic distance 사용하는 방법이 별로인가??
- Multivariate functional clustering 최근 방법론 찾아보기
  - gmfd 결과 추가
  - Martino, A., Ghiglietti, A., Ieva, F., & Paganoni, A. M. (2019). A k-means procedure based on a Mahalanobis type distance for clustering multivariate functional data. *Statistical Methods & Applications*, 28(2), 301-322.
- 방법론 한 번 쭉 정리하기 (kCFC 책 참고)
  - notation이랑 metric 쭉 정리하고, Thm 1을 manifold 버전으로 수정해보기!!
- 시뮬레이션 세팅은 나중에 리얼데이터의 meanfunction 형태를 보고 이와 비슷하게 generate해도 좋을듯
  - spike 튀는 경우 해결할 수 있으면 해볼것!!
    - m = 51, K = 5로 조절해봐도 크게 다르지 않음...

### 2022.02.25(금)

- Opensky Network 데이터 사용해보기

  - https://opensky-network.org/data/impala
  - 데이터 파악
  
    - callsign 앞의 3개 문자로 항공사 알 수 있음
      - https://en.wikipedia.org/wiki/List_of_airline_codes
    - 데이터가 너무 크기 때문에 hour 변수로 조건을 줘야하는듯 (1시간 단위라고 함)
      - time이 초단위 관측이고, hour가 시간별임
      - 또한 반드시 hour를 where문에 넣어주어야함!!! (time만 조건 걸면 오래 걸린다고 함)
  - 데이터 가져오기

    - ICAO-24
      - 인천공항 : RKSI
      - LA : KLAX
      - 런던 : EGLL
    - `인천 -> LA` 는 너무 위도가 비슷해서 curve 형태가 잘 안보임
      - 위도가 크게 바뀌는 직항노선을 고려하자 (런던, 시드니)
    - 필요한 DB : flights_data4, state_vectors_data4
    - `flights_data4`
      - 우선 flights_data4 에서 도착지에 대한 icao24 리스트 쭉 뽑기 (예를 들어, 인천 -> LA)
      - hour 변수로 기간 조절 (근데 이 table에는 hour 변수가 없음...)
        - 대신 day 변수가 있긴 함
      - callsign 변수로 항공사 지정
      - estdepartureairport, estarrivalairport 로 출발지랑 도착지 설정
      - 근데 2019년1월~3월말까지 비행기록이 없음... (3월 말에 5일정도만 데이터 존재)
        - 2019-03-01 ~ 2019-06-30
      - firstseen이 출발, lastseen이 도착시간인듯!!!
    - `state_vectors_data4`
      - 해당 icao24 리스트에 대한 비행기록을 state_vectors_data4 에서 가져오기
      - 지금보다 조건을 더 걸어야 할듯.... 매우 느림
      - onground 변수는 좀 불확실한듯... 사용 안하는게 좋아보임
    - 출발, 도착 시간을 어떻게 정하지??
      - firstseen이 출발, lastseen이 도착시간인듯!!! 
    - 도착지를 LA 말고 다른 곳으로 하자
      - 대서양 넘어가는 루트가 중간에 관측이 안되는 것 같음
      - 근데 원래 태평양 지나가는 루트인데 왜 관측이 안되지??
  - Flightaware 데이터 사용방법
  - OpenSky Network 데이터에서 적당히 길면서 그룹 잘 나뉜 루트 찾아보기

    - 중간중간 끊긴 데이터는 interpolation해서 사용해도 괜찮을듯!!
    - 런던 -> LGAV(아테네)
    
      - 여기도 여전히 중간에 missing인 경우가 있음
    - LA(KLAX) -> 뉴욕(KJFK)
    - 런던 -> WSSS(싱가포르), OBBI(바레인), VIDP(델리), OMDB(두바이), SBGR(상파울루), YSSY(시드니), OTHH(도하, 카타르), LLBG(이스라엘), VHHH(홍콩)
    
      - 싱가포르는 중간에 끊기는 구간이 너무 길다...
      - 델리도 별로.....
      - 시드니는 애초에 출발점에서의 관측치도 잘 없음
      - 상파울루도 너무 missing 구간이 길다
      - 홍콩도 아시아대륙 부근에서 관측이 안됨...
    - Dai and Muller 페이퍼에서는 총 969개 항공경로 사용함
    - 싱가포르 -> 런던 667개 돌고 있음!!
      - 근데 이 결과도 그닥 좋지 않음...
      - clustering이다보니 가까운 것들끼리 묶이는데, 노선들이 골고루 섞여있어서 구분이 쉽지 않음

### 다른 리얼데이터

- 초파리 데이터

  - 데이터
    - Chiou and Muller (2014), Linear manifold modelling of multivariate functional data
  - compositional data를 square root 취해서 sphere-valued data로 바꿔서 함
  - 62개 trajectories
  - 각 cuve별 관측 끝나는 지점이 달라서, 적당히 day <= 37 로 잘라서 사용
  - RFPCA, MFPCA가 거의 비슷
- Bird migration 데이터
    - Movebank
      - https://www.movebank.org/cms/webapp?gwt_fragment=page=studies,path=study9651291
    - Reference paper
      - https://doi.org/10.1007/s10531-018-1538-6
    - 봄, 가을로 구분되는 것 같은데, 정확한 기준을 찾아보자
    - 이동할 때 연도 바뀌는 경우가 있을 수도 있음
      - GPS 데이터이다 보니, 매달 관측됨 (위치가 크게 달라지지 않더라도 센서를 통해 관측됨)
      - 예를 들어, 11월부터 이동 시작해서 2월까지 이동하는 경우
      - 그냥, 매뉴얼하게 하는게 빠를듯...ㅠㅠ
        - 매뉴얼하게 하니까 56개 trajectories
- DTI(Diffusion Tensor Imaging) data

    - refund 패키지에 DTI라는 데이터가 있긴 함
    - 근데 raw DTI data를 preprocessing해서 correlation matrix로 만들어주고, 얘를 SPD matrix로 간주하여 모델링
    - HCP 900 subjects release
      - 
    - ADNI
      - 근데 이건 dense는 아니고 sparse인듯
      - 만약 이 데이터 사용할 경우, 반드시 라이센스 다시 읽어보기!!!

        - draft도 먼저 보내줘야하고, 사사문구도 적어주어야 함

### 2022.03.18(금)

- kCFC 페이퍼의 property를 Riemannian mfd로 확장해보기

### 2022.07

- kCFC 논문의 Theorem 확장은 component의 rate이 필요한데, iRFPCA에서의 rate이 알려진 바가 없음
  - 이론 전개는 어려울듯...
- iRFPCA 코드 짜는 것도 시간이 꽤 걸릴듯... (내가 geometry를 이해를 잘 못하고 있다...)
  - sphere-valued data를 찾아서 하는 것이 좋을 듯
- 지구 위의 데이터 또는 compositional 데이터 찾아보기
- kCRFC가 별 차이가 없는 이유?
  - 잘 생각해보니, true와 estimated의 geodesic distance를 계산해서 비교하는데, sphere 위에서는 두 점 사이의 거리가 매우 짧을 것이기 때문에, geodisic distance나 L2 distance나 거의 차이가 없을 것임
    - 그래서 iteration step에서 MFPCA나 RFPCA나 차이가 없었던 것 같음
  - reconstruction 후에 distance 비교 말고 다른 방법 있을까?
    - 처음에 2개로 initial clustering, 다시 각각 2개씩 총 4개로 clustering, 거기서 prediction error가 reconstruction error보다 작으면 같은 cluster로 묶고, 다시 같은 과정 반복? (similarity measure)
    - Riemannian CLT 써서 normality 적용한 model based algoritm?? EM algorithm으로 optimize?
    - Fuzzy c-means clustering하고, probability 애매한 애들만 weighted reclassification?
      - assign probability < 0.8 인 경우에만 reclassification step 적용 (iteration 1번)
- 리얼 데이터
  - CO, NO2, PM10 등 대기오염 종류들 여러개로 index 만들어서 compositional 데이터로 이용?
    - Air Pollution Index (API)
    - https://link.springer.com/content/pdf/10.1007/s00477-018-1542-0.pdf
    - 구성비로 계산하다보니 값의 변화가 거의 없음...(심지어 RFPCA가 안돌아감)

### 2023.02.18(토) - 1st revision

- bird migration 전처리에서 오타 있었음... (총 57개임!! 이전에는 오타 떄문에 56개였음)

  - 이것도 draft에 수정함

- mode of variation(분산) 똑같이 하고 Euclidean manifold에서 생성하고 mean function만 다르게 한 시뮬결과 확인할 것!

  1. mode of variation 똑같고 Sphere

     - S2

       ```
             Method           CCR          aRand           CCR         aRand           CCR         aRand
       1    kCFC(R) 0.931 (0.083)  0.767 (0.247) 0.858 (0.169) 0.622 (0.387) 0.754 (0.163) 0.357 (0.358)
       2    kCFC(M) 0.938 (0.043)  0.771 (0.149) 0.823 (0.157) 0.510 (0.361) 0.770 (0.166) 0.394 (0.360)
       3 K-means(R) 0.863 (0.045)  0.531 (0.129) 0.708 (0.098) 0.204 (0.189) 0.731 (0.082) 0.233 (0.155)
       4 K-means(M) 0.863 (0.050)  0.534 (0.145) 0.712 (0.108) 0.217 (0.218) 0.725 (0.093) 0.228 (0.179)
       5   funclust 0.841 (0.070)  0.479 (0.173) 0.706 (0.112) 0.212 (0.184) 0.858 (0.093) 0.542 (0.253)
       6    funHDDC 0.827 (0.055)  0.435 (0.146) 0.595 (0.063) 0.043 (0.068) 0.665 (0.068) 0.119 (0.094)
       7       gmfd 0.539 (0.029) -0.001 (0.014) 0.709 (0.111) 0.216 (0.200) 0.582 (0.094) 0.053 (0.168)
       ```

     - S3

       ```
             Method           CCR         aRand           CCR         aRand           CCR         aRand
       1    kCFC(R) 0.935 (0.092) 0.788 (0.236) 0.812 (0.173) 0.503 (0.379) 0.904 (0.123) 0.711 (0.302)
       2    kCFC(M) 0.878 (0.129) 0.633 (0.319) 0.665 (0.133) 0.172 (0.261) 0.876 (0.122) 0.620 (0.306)
       3 K-means(R) 0.836 (0.060) 0.459 (0.158) 0.639 (0.084) 0.096 (0.121) 0.794 (0.062) 0.353 (0.141)
       4 K-means(M) 0.833 (0.071) 0.457 (0.185) 0.635 (0.082) 0.091 (0.119) 0.794 (0.063) 0.356 (0.145)
       5   funclust 0.782 (0.124) 0.373 (0.287) 0.743 (0.147) 0.315 (0.296) 0.803 (0.090) 0.395 (0.174)
       6    funHDDC 0.862 (0.052) 0.530 (0.145) 0.648 (0.075) 0.101 (0.106) 0.806 (0.046) 0.378 (0.112)
       7       gmfd 0.540 (0.031) 0.001 (0.013) 0.656 (0.135) 0.162 (0.233) 0.754 (0.142) 0.331 (0.294)
       ```

  2. mode of variation 똑같고 Euclidean manifold

     - 근데 RFPCA에서 sphere로 계산시 fretchet mean 계산 과정에서 오류

       - Sphere structure를 만족하지 않는 경우에 fretchet mean 계산이 불가능

     - E2

       ```
             Method           CCR         aRand           CCR         aRand           CCR         aRand
       1    kCFC(R) 0.995 (0.029) 0.983 (0.094) 0.967 (0.114) 0.924 (0.252) 0.788 (0.204) 0.493 (0.471)
       2    kCFC(M) 0.995 (0.029) 0.983 (0.094) 0.967 (0.114) 0.924 (0.252) 0.788 (0.204) 0.493 (0.471)
       3 K-means(R) 0.984 (0.019) 0.939 (0.073) 0.974 (0.061) 0.913 (0.191) 0.851 (0.122) 0.547 (0.321)
       4 K-means(M) 0.984 (0.019) 0.939 (0.073) 0.974 (0.061) 0.913 (0.191) 0.851 (0.122) 0.547 (0.321)
       5   funclust 0.820 (0.160) 0.506 (0.355) 0.757 (0.152) 0.350 (0.293) 0.772 (0.160) 0.391 (0.336)
       6    funHDDC 0.968 (0.030) 0.880 (0.110) 0.658 (0.077) 0.115 (0.128) 0.738 (0.094) 0.255 (0.203)
       7       gmfd 0.542 (0.032) 0.001 (0.016) 0.561 (0.075) 0.028 (0.141) 0.541 (0.033) 0.001 (0.016)
       ```

     - E3

       ```
             Method           CCR         aRand           CCR         aRand           CCR         aRand
       1    kCFC(R) 1.000 (0.000) 1.000 (0.000) 0.977 (0.095) 0.945 (0.219) 0.962 (0.111) 0.902 (0.270)
       2    kCFC(M) 1.000 (0.000) 1.000 (0.000) 0.977 (0.095) 0.945 (0.219) 0.962 (0.111) 0.902 (0.270)
       3 K-means(R) 0.994 (0.011) 0.976 (0.043) 0.982 (0.071) 0.947 (0.184) 0.957 (0.048) 0.845 (0.159)
       4 K-means(M) 0.994 (0.011) 0.976 (0.043) 0.982 (0.071) 0.947 (0.184) 0.957 (0.048) 0.845 (0.159)
       5   funclust 0.826 (0.188) 0.563 (0.425) 0.791 (0.183) 0.466 (0.384) 0.816 (0.189) 0.537 (0.413)
       6    funHDDC 0.989 (0.014) 0.958 (0.054) 0.686 (0.076) 0.153 (0.136) 0.945 (0.037) 0.794 (0.129)
       7       gmfd 0.542 (0.035) 0.002 (0.017) 0.552 (0.042) 0.008 (0.027) 0.757 (0.175) 0.381 (0.410)
       ```

       




## Appendix. ADNI data preprocessing

### 0. Preperation

- Install cmake

- dcm2niix : https://github.com/rordenlab/dcm2niix
- FSL : https://artiiicy.tistory.com/70

```
# install cmake
brew install cmake

# install dcm2niix
cd ~
git clone https://github.com/rordenlab/dcm2niix.git
cd dcm2niix
mkdir build && cd build
cmake -DZLIB_IMPLEMENTATION=Cloudflare -DUSE_JPEGLS=ON -DUSE_OPENJPEG=ON ..
make

# make python2
conda create --name py2 python=2
source activate py2

# FSL install (fslinstaller.py 먼저 다운)
python2 ~/tmp/fslinstaller.py
# Install MRtrix3 (mrcalc, dwinormalise)
conda install -c omnia eigen3
conda install -c dsdale24 qt5
conda install -c mrtrix3 mrtrix3

# intall ANTs
# https://github.com/ANTsX/ANTs/wiki/Compiling-ANTs-on-Linux-and-Mac-OS
mkdir ANTs
cd ANTs
git clone https://github.com/ANTsX/ANTs.git
# make a build and install dir
mkdir build install
cd build
ccmake ../ANTs
위의 링크 순서대로 실행
# use multicore
ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=4

# make python 3.7
conda deactivate
conda create --name py37 python=3.7
source activate py37

# install Scilpy
git clone https://github.com/scilus/scilpy.git
pip install -e .
```

### Preprocessing

```
cd ~/tmp
mkdir diffusion/
mkdir diffusion/dti

~/dcm2niix/build/bin/dcm2niix -o ~/tmp/diffusion ~/tmp/data

# 폴더 내에 nii, bval, bvec 을 하나하나씩 변환해야함!!
cd ~/tmp/diffusion
mv *.bval bval
mv *.bvec bvec
mv *.nii temp.nii
rm *.json


mrconvert -fslgrad bvec bval temp.nii temp.mif -force

dwidenoise temp.mif temp_denoised.mif -noise temp_noise.mif -force

dwifslpreproc temp_denoised.mif temp_preproc.mif -rpe_none -pe_dir j -nocleanup -force

dwiextract temp_preproc.mif - -bzero | mrmath - mean b0.mif -axis 3

mrconvert temp_preproc.mif data.nii.gz -export_grad_mrtrix temp_grad_table -export_grad_fsl bvecs bvals -force

mrconvert b0.mif b0.nii.gz

bet2 b0.nii.gz b0_bet.nii.gz -f 0.2

dwibiascorrect ants data.nii.gz data_bias_corrected.nii.gz -fslgrad bvecs bvals -mask b0_bet.nii.gz -ants.c [300x150x75x50,1e-6] -force
# 작은 따옴표로 묶어줘야함!! (zsh: no matches found: [300x150x75x50,1e-6])
dwibiascorrect ants data.nii.gz data_bias_corrected.nii.gz -fslgrad bvecs bvals -mask b0_bet.nii.gz -ants.c '[300x150x75x50,1e-6]' -force

dwinormalise individual data_bias_corrected.nii.gz b0_bet.nii.gz data_normalized.nii.gz -fslgrad bvecs bvals

dtifit --data=data_normalized.nii.gz --out=dti --mask=b0_bet.nii.gz --bvecs=bvecs --bvals=bvals --save_tensor

```

