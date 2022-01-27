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