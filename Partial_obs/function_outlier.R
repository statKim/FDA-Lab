### Outlier detection by robust Mahalanobis distance based on first 2 FPC scores
# input: "funPCA" object which is the output of function robfpca::funPCA()
# return: length n vector (if outlier, return 1, else 0)
fun_outlier <- function(funPCA.obj, method = c("robMah","outmap")) {
  if (method == "robMah") {
    # # first 2 FPC scores
    # pc.score <- funPCA.obj$pc.score[, 1:2]
    
    # selected K FPC scores
    k <- funPCA.obj$K
    pc.score <- funPCA.obj$pc.score[, 1:k]
    
    n <- nrow(pc.score)
    rownames(pc.score) <- 1:n
    s <- cbind(pc.score, 
               rep(1, n))
    
    # function from "foutliers" function in "cran/rainbow"
    robout <- function (data, nclass, meth = c("mve", "mcd"), rep = 10) {
      ncol = dim(data)[2]
      tempo = data[data[, ncol] == nclass, 1:(ncol - 1)]
      # rownames(tempo) <- 1:100
      namestempo = rownames(tempo)
      nrow = dim(tempo)[1]
      roboutl = NULL
      roboutall = matrix(0, nrow, rep)
      rownames(roboutall) = namestempo
      for (i in 1:rep) {
        mcdc = cov.rob(tempo, method = meth)
        mbc = sqrt(mahalanobis(tempo, mcdc$center, mcdc$cov, 
                               to = 1e-14))
        roboutl = c(roboutl, boxplot(mbc, plot = FALSE)$out)
        roboutall[, i] = mbc
      }
      a = as.matrix(roboutl)
      b = apply(roboutall, 1, mean)
      outme = rev(sort(b))
      topo = rev(sort(b))[1:10]
      top = table(as.numeric(rownames(a)))
      top1 = top[top > ceiling(rep / 2)]
      topout = as.numeric(names(top1))
      ntops = length(topout)
      outly = rep(0, ntops)
      for (i in 1:ntops) {
        outly[i] = mean(a[as.numeric(rownames(a)) == topout[i]])
      }
      return(list(outliers = topout, depth.total = outme))
    }
    
    out <- robout(s, 1, "mcd")$outliers
  } else if (method == "outmap") {
    # Outlier map (Score dist VS Orthogonal dist), See Hubert(2005)
    SD <- score_dist(funPCA.obj)   # score distance
    OD <- orthogonal_dist(funPCA.obj)   # orthogonal distance
    # mcd_fit <- covMcd(OD^(2/3))
    mcd_fit <- MASS::cov.mcd(matrix(OD^(2/3)))
    cut_y <- (mcd_fit$center + sqrt(as.numeric(mcd_fit$cov))*qnorm(0.975))^(3/2)
    cut_x <- sqrt(qchisq(0.975, k))
    
    out <- which(SD > cut_x & OD > cut_y)
  }
  
  out_y <- rep(0, n)
  out_y[out] <- 1
  
  return(out_y)
}




# score distance
score_dist <- function(funPCA.obj) {
  K <- funPCA.obj$K
  score <- funPCA.obj$pc.score
  lambda <- funPCA.obj$lambda
  
  SD <- apply(score, 1, function(row){
    sqrt(sum(row^2 / lambda))
  })
  
  return(SD)
}

# orthogonal distance
orthogonal_dist <- function(funPCA.obj) {
  X <- list2matrix(funPCA.obj$data)
  K <- funPCA.obj$K
  X_hat <- predict(funPCA.obj, K = K)
  
  OD <- (X - X_hat)^2
  OD <- apply(OD, 1, function(row){
    sqrt(sum(row, na.rm = T))
  })
  
  return(OD)
}