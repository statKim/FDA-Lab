##################################
### bagging classifier with FPCA
##################################

### majority voting functon
majority_vote <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


### train, test split
##  - input data form : cbind(id, time, val, y)
train_test_split <- function(data, train.prop=2/3) {
  # # train, test split => range(test set) > range(train set) 인 경우의 id 제거
  # min.max.grid <- data %>% 
  #   filter(time %in% range(time)) %>% 
  #   dplyr::select(id) %>% 
  #   unique
  # id <- setdiff(unique(data$id), min.max.grid$id)
  # 
  # N <- length(unique(data$id))
  # N.train <- ceiling(N * train.prop)
  # id.train <- c(sample(id, N.train-length(unique(min.max.grid$id))),
  #               min.max.grid$id)
  # id.test <- id[-which(id %in% id.train)]
  
  ### train, test split => train set에 없는 time point 포함시 에러
  # range(test set) > range(train set) => 넘는 부분에 해당하는 id 제거
  id <- unique(data$id)
  N <- length(unique(data$id))
  N.train <- ceiling(N * train.prop)
  id.train <- sample(id, N.train)
  id.test <- NULL
  range.train <- range(data$time[which(data$id %in% id.train)])
  range.test <- range(data$time[-which(data$id %in% id.train)])
  if (range.test[1] < range.train[1]) {
    over.ind <- which(data$time[-which(data$id %in% id.train)] < range.train[1] )
    over.ind <- data$id[-which(data$id %in% id.train)][over.ind]
    id.test <- id[-which(id %in% c(id.train, over.ind))]
  }
  if (range.test[2] > range.train[2]) {
    over.ind <- which(data$time[-which(data$id %in% id.train)] > range.train[2] )
    over.ind <- data$id[-which(data$id %in% id.train)][over.ind]
    if (is.numeric(id.test)) {
      id.test <- intersect(id.test,
                           id[-which(unique(data$id) %in% c(id.train, over.ind))])
    } else {
      id.test <- id[-which(id %in% c(id.train, over.ind))]
    }
  }
  if (!is.numeric(id.test)) {
    id.test <- id[-which(id %in% id.train)]
  }
  
  # transform to FPCA input
  train <- data[which(data$id %in% id.train), ]
  test <- data[which(data$id %in% id.test), ]
  
  X.train <- MakeFPCAInputs(IDs = train$id,
                            tVec = train$time,
                            yVec = train$val)
  X.test <- MakeFPCAInputs(IDs = test$id,
                           tVec = test$time,
                           yVec = test$val)
  y.train <- unique(train[, c("id", "y")])[, "y"]
  y.test <- unique(test[, c("id", "y")])[, "y"]
  
  return(list(X.train = X.train,
              X.test = X.test,
              y.train = y.train,
              y.test = y.test))
}




### Single classifier
get_single_err <- function(X.train, X.test, y.train, y.test) {
  fpca.fit <- FPCA(X.train$Ly, 
                   X.train$Lt, 
                   optns=list(dataType="Sparse",
                              methodXi="CE",
                              # methodSelectK="BIC",
                              FVEthreshold=0.99,
                              verbose=F))
  k <- fpca.fit$selectK   # optimal number of PCs
  fpc.score <- fpca.fit$xiEst
  
  # train FPC scores
  train.fpc <- data.frame(y = y.train, 
                          x = fpc.score)
  colnames(train.fpc) <- c("y", paste("FPC", 1:k, sep=""))
  
  # tune the hyperparmeters
  tune.linear <- tune.svm(y ~ .,
                          data = train.fpc,
                          kernel = "linear",
                          cost = c(10^(-3:1), 2^(5:10)))
  tune.radial <- tune.svm(y ~ .,
                          data = train.fpc,
                          kernel = "radial",
                          cost = c(10^(-3:1), 2^(5:10)),
                          gamma = c(10^(-3:1), 2^(5:10)))
  
  # fit classifiers
  fit.logit <- glm(y~., train.fpc, family=binomial)
  fit.svm.linear <- svm(y ~ .,
                        data = train.fpc, 
                        kernel = "linear", 
                        cost = tune.linear$best.model$cost)
  fit.svm.radial <- svm(y ~ .,
                        data = train.fpc, 
                        kernel = "radial", 
                        cost = tune.radial$best.model$cost,
                        gamma = tune.radial$best.model$gamma)
  fit.lda <- lda(y~., train.fpc)
  fit.qda <- qda(y~., train.fpc)
  fit.nb <- naiveBayes(y~., train.fpc)
  
  # test FPC scores
  fpc.score.test <- predict(fpca.fit, X.test$Ly, X.test$Lt, K=k, xiMethod="CE")
  test.fpc <- data.frame(y = y.test, 
                         x = fpc.score.test)
  colnames(test.fpc) <- c("y", paste("FPC", 1:k, sep=""))
  
  # predict
  pred <- list()
  pred.logit <- predict(fit.logit, test.fpc, type="response")
  pred[[1]] <- factor(ifelse(pred.logit > 0.5, 1, 0), levels=c(0, 1))
  pred[[2]] <- predict(fit.svm.linear, test.fpc)
  pred[[3]] <- predict(fit.svm.radial, test.fpc)
  pred[[4]] <- predict(fit.lda, test.fpc)$class
  pred[[5]] <- predict(fit.qda, test.fpc)$class
  pred[[6]] <- predict(fit.nb, test.fpc, type="class")
  
  # save the error rate
  err.single <- sapply(pred, function(x){ mean(x != y.test) })
  
  return(list(err.single = err.single,
              hyper.para = list(tune.linear = tune.linear,
                                tune.radial = tune.radial)))
}



## Bagging
# Bootstrap aggregating
get_bag_err <- function(X.train, X.test, y.train, y.test, B = 100, packages = c("fdapace","e1071","MASS"), hyper.para) {
  tune.linear <- hyper.para$tune.linear
  tune.radial <- hyper.para$tune.radial
  
  start.time <- Sys.time()
  N <- length(X.train$Ly)
  boot.mat <- matrix(sample(1:N, N*B, replace=T), N, B)
  y.pred <- foreach(b=1:B, .packages=packages) %dopar% {
    boot.ind <- boot.mat[, b]
    
    # bootstrap sample의 time range를 벗어나는 test set sample 제거
    id.test <- NULL
    max.boot <- max(unlist(X.train$Lt[boot.ind]))
    min.boot <- min(unlist(X.train$Lt[boot.ind]))
    if (min(unlist(X.test$Lt)) < min.boot) {
      # train set에서의 index(Not id)
      over.ind <- which(sapply(X.test$Lt, 
                               function(x){ sum(x < min.boot) }) != 0)
      id.test <- unlist(X.test$Lid)[-over.ind]
    }
    if (max(unlist(X.test$Lt)) > max.boot) {
      # train set에서의 index(Not id)
      over.ind <- which(sapply(X.test$Lt, 
                               function(x){ sum(x > max.boot) }) != 0)
      if (is.numeric(id.test)) {
        id.test <- intersect(id.test,
                             unlist(X.test$Lid)[-over.ind])
      } else {
        id.test <- unlist(X.test$Lid)[-over.ind]
      }
    }
    if (!is.numeric(id.test)) {
      id.test <- unlist(X.test$Lid)
    }
    
    # bootstrap sample의 time range를 벗어나는 OOB sample 제거
    id.oob <- NULL
    if (min(unlist(X.train$Lt[-unique(boot.ind)])) < min(unlist(X.train$Lt[boot.ind]))) {
      # train set에서의 index(Not id)
      over.ind <- which(sapply(X.train$Lt[-unique(boot.ind)], 
                               function(x){ sum( x < min(unlist(X.train$Lt[boot.ind])) ) }) != 0)
      id.oob <- unlist(X.train$Lid)[-unique(boot.ind)][-over.ind]
    }
    if (max(unlist(X.train$Lt[-unique(boot.ind)])) > max(unlist(X.train$Lt[boot.ind]))) {
      # train set에서의 index(Not id)
      over.ind <- which(sapply(X.train$Lt[-unique(boot.ind)], 
                               function(x){ sum( x > max(unlist(X.train$Lt[boot.ind])) ) }) != 0)
      if (is.numeric(id.oob)) {
        id.oob <- intersect(id.oob,
                            unlist(X.train$Lid)[-unique(boot.ind)][-over.ind])
      } else {
        id.oob <- unlist(X.train$Lid)[-unique(boot.ind)][-over.ind]
      }
    }
    if (!is.numeric(id.oob)) {
      id.oob <- unlist(X.train$Lid)[-unique(boot.ind)]
    }
    
    # fit FPCA for bootstrapped data
    fpca.fit <- FPCA(X.train$Ly[boot.ind], 
                     X.train$Lt[boot.ind], 
                     optns=list(dataType="Sparse",
                                methodXi="CE",
                                # methodSelectK="BIC",
                                FVEthreshold=0.99,
                                verbose=F))
    k <- fpca.fit$selectK   # optimal number of PCs
    fpc.score <- fpca.fit$xiEst
    
    # train FPC scores
    train.fpc <- data.frame(y = y.train[boot.ind], 
                            x = fpc.score)
    colnames(train.fpc) <- c("y", paste("FPC", 1:k, sep=""))
    
    # classifiers with bootstrapped FPC scores
    fit.logit <- glm(y~., train.fpc, family=binomial)
    fit.svm.linear <- svm(y ~ ., 
                          data = train.fpc,
                          kernel = "linear",
                          cost = tune.linear$best.model$cost)
    fit.svm.radial <- svm(y ~ ., 
                          data = train.fpc,
                          kernel = "radial",
                          cost = tune.radial$best.model$cost,
                          gamma = tune.radial$best.model$gamma)
    fit.lda <- lda(y~., train.fpc)
    fit.qda <- qda(y~., train.fpc)
    fit.nb <- naiveBayes(y~., train.fpc)
    
    # test FPC scores
    ind <- which(unlist(X.test$Lid) %in% id.test)
    fpc.score.test <- predict(fpca.fit, X.test$Ly[ind], X.test$Lt[ind], K=k, xiMethod="CE")
    test.fpc <- data.frame(y=y.test[ind], 
                           x=fpc.score.test)
    colnames(test.fpc) <- c("y", paste("FPC", 1:k, sep=""))
    
    # predict
    pred.logit <- predict(fit.logit, test.fpc, type="response")
    pred.logit <- factor(ifelse(pred.logit > 0.5, 1, 0), levels=c(0, 1))
    pred.svm.linear <- predict(fit.svm.linear, test.fpc)
    pred.svm.radial <- predict(fit.svm.radial, test.fpc)
    pred.lda <- predict(fit.lda, test.fpc)$class
    pred.qda <- predict(fit.qda, test.fpc)$class
    pred.nb <- predict(fit.nb, test.fpc, type="class")
    
    # OOB FPC scores
    oob.ind <- which(X.train$Lid %in% id.oob)
    fpc.score.oob <- predict(fpca.fit,
                             X.train$Ly[oob.ind],
                             X.train$Lt[oob.ind],
                             K=k,
                             xiMethod="CE")
    oob.fpc <- data.frame(y=y.train[oob.ind],
                          x=fpc.score.oob)
    colnames(oob.fpc) <- c("y", paste("FPC", 1:k, sep=""))
    
    # OOB error rate
    oob.error <- c(mean(factor(ifelse(predict(fit.logit, oob.fpc, type="response") > 0.5, 1, 0), levels=c(0, 1)) != oob.fpc$y),
                   mean(predict(fit.svm.linear, oob.fpc) != oob.fpc$y),
                   mean(predict(fit.svm.radial, oob.fpc) != oob.fpc$y),
                   mean(predict(fit.lda, oob.fpc)$class != oob.fpc$y),
                   mean(predict(fit.qda, oob.fpc)$class != oob.fpc$y),
                   mean(predict(fit.nb, oob.fpc, type='class') != oob.fpc$y) )
    
    # train error rate
    train.error <- c(mean(factor(ifelse(predict(fit.logit, train.fpc, type="response") > 0.5, 1, 0), levels=c(0, 1)) != train.fpc$y),
                     mean(predict(fit.svm.linear, train.fpc) != train.fpc$y),
                     mean(predict(fit.svm.radial, train.fpc) != train.fpc$y),
                     mean(predict(fit.lda, train.fpc)$class != train.fpc$y),
                     mean(predict(fit.qda, train.fpc)$class != train.fpc$y),
                     mean(predict(fit.nb, train.fpc, type='class') != train.fpc$y) )
    
    return( list(oob.error = oob.error,
                 train.error = train.error,
                 boot = data.frame(id = id.test,
                                   y = y.test[ind],
                                   # predicted value
                                   logit = pred.logit,
                                   svm.linear = pred.svm.linear,
                                   svm.radial = pred.svm.radial,
                                   lda = pred.lda,
                                   qda = pred.qda,
                                   nb = pred.nb,
                                   # predicted value * OOB accuracy
                                   logit.oob = as.numeric(pred.logit) * (1 - oob.error[1]),
                                   svm.linear.oob = as.numeric(pred.svm.linear) * (1 - oob.error[2]),
                                   svm.radial.oob = as.numeric(pred.svm.radial) * (1 - oob.error[3]),
                                   lda.oob = as.numeric(pred.lda) * (1 - oob.error[4]),
                                   qda.oob = as.numeric(pred.qda) * (1 - oob.error[5]),
                                   nb.oob = as.numeric(pred.nb) * (1 - oob.error[6]))) )
  }
  end.time <- Sys.time()
  print(end.time - start.time)
  
  ## save the accuracy
  res <- as.data.frame(rbindlist(lapply(y.pred, function(x){ x$boot })))
  # majority voting
  pred <- res %>% 
    group_by(id, y) %>% 
    summarise(logit = majority_vote(logit),
              svm.linear = majority_vote(svm.linear),
              svm.radial = majority_vote(svm.radial),
              lda = majority_vote(lda),
              qda = majority_vote(qda),
              nb = majority_vote(nb))
  err.majority <- 1 - apply(pred[, 3:8], 2, function(x){ mean(x == pred$y) })
  
  # oob error
  oob.acc <- colSums( 1 - t(sapply(y.pred, function(x){ x$oob.error })) )
  pred <- res %>% 
    group_by(id, y) %>% 
    summarise(logit = factor(ifelse(sum(logit.oob)/oob.acc[1] > 1.5, 1, 0), 
                             levels=c(0, 1)),
              svm.linear = factor(ifelse(sum(svm.linear.oob)/oob.acc[2] > 1.5, 1, 0), 
                                  levels=c(0, 1)),
              svm.radial = factor(ifelse(sum(svm.radial.oob)/oob.acc[3] > 1.5, 1, 0), 
                                  levels=c(0, 1)),
              lda = factor(ifelse(sum(lda.oob)/oob.acc[4] > 1.5, 1, 0), 
                           levels=c(0, 1)),
              qda = factor(ifelse(sum(qda.oob)/oob.acc[5] > 1.5, 1, 0), 
                           levels=c(0, 1)),
              nb = factor(ifelse(sum(nb.oob)/oob.acc[6] > 1.5, 1, 0), 
                          levels=c(0, 1)))
  err.oob <- 1 - apply(pred[, 3:8], 2, function(x){ mean(x == pred$y) })
  
  return(list(err.majority = err.majority,
              err.oob = err.oob,
              y.pred = y.pred))
}


### generate sparse functional data
sim.curve <- function(n, sparsity=5:10, model="A") {
  ## generate curves
  y <- factor(c(rep(0, n/2), rep(1, n/2)), levels=c(0, 1))
  data <- list()
  # random sparsify => generate curve
  for (i in 1:n) {
    # random sparsify
    num.obs <- sample(sparsity, 1)
    t <- sort(runif(num.obs, 0, 10))
    
    # mean function and variance
    if (model == "A") {
      if (i < n/2) {
        mu <- t + sin(t)
        lambda <- c(4, 2, 1)
      } else {
        mu <- t + cos(t)
        lambda <- c(16, 8, 4)
      }
    } else if (model == "B") {
      if (i < n/2) {
        mu <- t + sin(t)
        lambda <- c(4, 2, 1)
      } else {
        mu <- t + cos(t)
        lambda <- c(4, 2, 1)
      }
    } else if (model == "C") {
      if (i < n/2) {
        mu <- t + sin(t)
        lambda <- c(4, 2, 1)
      } else {
        mu <- t + sin(t)
        lambda <- c(16, 8, 4)
      }
    } else {
      stop("Model is not existed. 3 models(A, B, C) can be applied.")
    }
    
    # eigenfunctions
    phi <- cbind(cos(pi*t/5) / sqrt(5),
                 sin(2*pi*t/5) / sqrt(5),
                 cos(3*pi*t/5) / sqrt(5))
    
    # FPC scores
    xi <- cbind(rnorm(1, 0, sqrt(lambda[1])),
                rnorm(1, 0, sqrt(lambda[2])),
                rnorm(1, 0, sqrt(lambda[3])))
    
    # measurement error
    eps <- rnorm(num.obs, 0, 0.5)
    
    # generate the curve with measurement error
    X <- mu + xi %*% t(phi) + eps
    
    data$id[[i]] <- rep(i, num.obs)
    data$y[[i]] <- rep(y[i], num.obs)
    data$Lt[[i]] <- t
    data$Ly[[i]] <- X
  }
  
  # return type of data
  id <- unlist( mapply(function(x, y){ rep(x, y) },
                       1:n,
                       sapply(data$Lt, length)) )
  time <- unlist(data$Lt)
  val <- unlist(data$Ly)
  y <- unlist( mapply(function(x, z){ rep(x, z) },
                      y,
                      sapply(data$Lt, length)) )
  data <- data.frame(id, time, val, y)

  return(data)
}
