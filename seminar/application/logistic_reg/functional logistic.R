setwd("C:\\Users\\user\\Desktop\\KHS\\Thesis\\seminar\\application\\logistic_reg")

library(fda.usc)
library(dplyr)
library(ggplot2)
library(data.table)

set.seed(100)
cell <- read.delim("mbc_9_12_3273__CDCDATA.txt")
gene <- cell[, c(1, 6:23)] # 0~120min, 18 points
gene <- na.omit(gene)

train <- gene[, -1] # remove the column of gene's names
train.fdata <- fdata(train)
train.fd <- fdata2fd(train.fdata, nbasis=5)
train.pca <- create.pc.basis(train.fdata, l=1:5, lambda=1)
train.pca$basis$data <- -train.pca$basis$data

# epsilon 1~5, group label을 성분으로 하는 100개의 datasets(trian 100, test 100) 생성
means <- seq(0.6, 0.2, -0.1)
variances <- c(2.6950, 0.8850, 0.1957, 0.1266, 0.1079)
random.sets <- list()
for (i in 1:100) {
  e1 <- append( rnorm(100, means[1], sqrt(variances[1])), rnorm(100, -means[1], sqrt(variances[1])) )
  e2 <- append( rnorm(100, means[2], sqrt(variances[2])), rnorm(100, -means[2], sqrt(variances[2])) )
  e3 <- append( rnorm(100, means[3], sqrt(variances[3])), rnorm(100, -means[3], sqrt(variances[3])) )
  e4 <- append( rnorm(100, means[4], sqrt(variances[4])), rnorm(100, -means[4], sqrt(variances[4])) )
  e5 <- append( rnorm(100, means[5], sqrt(variances[5])), rnorm(100, -means[5], sqrt(variances[5])) )
  group <- as.factor(c( rep("g1", 100), rep("g0", 100) ))
  random.set <- list(e1=e1, e2=e2, e3=e3, e4=e4, e5=e5, group=group)
  set.i <- paste("set", as.character(i), sep="")
  random.sets[[set.i]] <- random.set
}
random.sets$set1$e1[1:10]
append( rnorm(10, means[1], sqrt(variances[1])), rnorm(10, -means[1], sqrt(variances[1])) )
e1.new <- as.data.frame((random.sets$set1$e1[c(1:5)] %o% train.pca$basis$data[1,]))
e2.new <- as.data.frame(t(random.sets$set1$e2[c(1:5)] %o% train.pca$basis$data[2,]))

# mean curve에 mth eigenfunction_m과 epsilon_m을 곱한 값을 더한 X samples 생성
X.curves <- list()
mean.curve <- as.vector(train.pca$mean$data)
pca.coefs <- train.pca$basis$data
i <- 1
for (set.i in random.sets) {
  sum.e <- mean.curve # 1 datase합t당 PC X e들의 합
  for (j in 1:5) {
    ej <- paste("e", as.character(j), sep="")
    ej.new <- as.data.frame( (set.i[[ej]] %o% pca.coefs[j,]) )
    sum.e <- ej.new + sum.e
  }
  X.curves[[names(random.sets)[i]]] <- sum.e
  i <- i + 1
}

length(X.curves)
dim(X.curves[[1]])
random.set$group


#######################################
### generate된 데이터로 classification
#######################################
# time points 변수 만들기
time <- colnames(gene)[-1]
time <- as.integer(gsub(pattern="alpha|min",
                        replacement="",
                        x=time))

# reduced rank model에 맞도록 데이터 변환
library(reshape2)
sub <- cbind(gene=1:200, X.curves[[1]])
data <- melt(sub, "gene", variable.name="time", value.name="gene_exp")
# data$time <- as.numeric(data$time)
data$gene <- as.factor(data$gene)
data <- data %>% 
  arrange(gene, time)
dim(data)
head(data)
data$time <- time   # time points 생성


# PC score 계산
source("../EM_reduced_rank.R")
par(mfrow=c(2,3))
init_value <- c(.1, .1, .1, .1, .1)
for (kn in c(4, 9, 14)) {
  fit <- fpca.fit(data, iter=50, init_value=init_value, num_knots=kn, num_pc=5, grid_sep=1)
  
  # curves and mean function
  # par(mfrow=c(1,2))
  ind <- unique(data[,1])
  sub <- data[which(data[,1]==ind[1]), ]
  plot(sub[,2], sub[,3],
       type="o",
       xlim=c(range(data[,2])[1], range(data[,2])[2]),
       ylim=c(range(data[,3])[1], range(data[,3])[2]),
       xlab="minutes",
       ylab="gene expression level")
  for (i in ind) {
    sub <- data[which(data[,1]==i), ]
    lines(sub[,2], sub[,3], type="o")
  }
  lines(fit$Timegrid, fit$MeanFunction, col="red", lwd=3, type="l")
  
  # PC function
  for (i in 1:5){
    plot(fit$Timegrid, fit$FPCscore[,i], type="l", lwd=3, ylim=c(-.4, .2),
         xlab="minutes", ylab="Princ. comp.", main=paste("PC", i))
  }
  
  # loglikelihood 비교
  print( paste("Reduced rank(", kn, " knots) : ", round(fit$Loglik, 2), sep="") )
  
  # object 저장
  # save(fit, file=paste("rrm_", kn, "knots.RData", sep=""))
}

load("rrm_9knots.RData")

# functional logistic regression
fit$FPCscore[which(fit$Timegrid %in% time),]

y <- factor( c(rep(1, 50), rep(0, 50)) )
data <- dcast(data, gene~time, value.val=gene_exp)[, -1]
dim(data)
train <- cbind(y=y, data[1:100, ])
test <- cbind(y=y, data[101:200, ])
fit.logit <- glm(y ~ ., data=train, family = "binomial")
summary(fit.logit)
predict(fit.logit, test)
