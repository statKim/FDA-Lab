# devtools::install_github("rizbicki/FlexCoDE")
# devtools::install_github("rizbicki/predictionBands")
# library(predictionBands)
# library(FlexCoDE)

# finds k nearest neighbors in xTrain of each xTest
which_neighbors <- function(xTrain,xTest,k){
  return(FNN::get.knnx(data=xTrain,query = xTest,k=k)$nn.index)
}

profile_density <- function(t_grid,y_grid,cde_estimate){
  v2 <- cde_estimate[order(cde_estimate)]
  v2s <- rev(cumsum(rev(v2)))*(y_grid[2]-y_grid[1])
  v2s <- v2s[findInterval(t_grid, v2) + 1]
  v2s[which(is.na(v2s))] <- 0
  return(v2s)
}

cum_dist <- function(y_grid,cde_estimates,y_values){
  which_closest <- FNN::get.knnx(data = y_grid,
                                 query=y_values,k=1)$nn.index
  apply(as.matrix(1:nrow(cde_estimates)),1,function(xx){
    return(sum(cde_estimates[xx,1:which_closest[xx]])*diff(y_grid)[1])
  })
}


#' Fit conformal prediction bands based on density estimation for regression
#'
#' @param x Matrix with covariates of training set
#' @param y Vector with the (continuous) responses of training set
#' @param per_train # percentage of samples used for traning density estimator (defaults to 40\%)
#' @param per_val # percentage of samples used for tuning density estimator (defaults to 10\%)
#' @param per_ths  # percentage of samples used for computeing thresholds for the conformal method (defaults to 50\%)
#' @param k # Number of clusters for cd-split. Default to round(per_ths*nrow(as.matrix(x))/100) so that each cluster has on average 100 samples
#' @param regressionFunction # regression function to be used for FlexCode. Defaults to Random Forests. See FlexCode documentation for additional regression methods.
#' @param ... Additional arguments to FlexCoDE::fitFlexCoDE
#'
#' @return Returns an object of the class predictionBands with the following components:
#' \item{density_fit}{Object of the class FlexCoDE with the estimated density}
#' \item{cum_dist_evaluated_train}{Cumulative conditional distribution functions on the training set (for dist-split)}
#' \item{conformity_score_train}{Conformal scores on the training set (for cd-split)}
#' \item{conformity_score_train_hpd}{Conformal scores on the training set (for hpd-split)}
#' \item{t_grid}{Level sets of the densities}
#' \item{g_train}{Profiles of the training sample}
#' \item{center_kmeans}{The center of the clusters found by kmeans (in the profile space)}
#' @export
#'
#' @examples
#'
#' # generate data
#' n <- 1000
#' n_new <- 50
#' d <- 10
#' data <- matrix(NA,n,d+1)
#' x <- matrix(rnorm(n*d),n,d)
#' y <- x[,1]+rnorm(n,0,0.1)
#' fit <- fit_predictionBands(x,y,0.5,0.4,0.1)
#'
#' xnew <- matrix(rnorm(n_new*d),n_new,d)
#' ynew <- xnew[,1]+rnorm(n_new,0,0.1)
#'
#'  # Dist-split
#'  bands <- predict(fit,xnew,type="dist")
#'  bands[[1]]
#'  bands[[2]]
#'  plot(bands)
#'  plot(bands,ynew)
#'
#'  # CD-split
#'  bands <- predict(fit,xnew,type="cd")
#'  bands[[1]]
#'  bands[[2]]
#'  plot(bands)
#'  plot(bands,ynew)
fit_predictionBands <- function(x,y,
                                per_train=0.4,
                                per_val=0.1,
                                per_ths=1-per_train-per_val,
                                splits = NULL,
                                k=max(round(per_ths*nrow(as.matrix(x))/100),1),
                                regressionFunction=FlexCoDE::regressionFunction.Forest,
                                ...) {
  x <- as.matrix(x)
  n_levels <- 1000
  if (is.null(splits)) {
    splits <- sample(c("Train","Validation","Threshold"),
                     size = nrow(x),
                     prob = c(per_train,per_val,per_ths),
                     replace = TRUE)
  }
  
  
  fit <- FlexCoDE::fitFlexCoDE(xTrain=x[splits=="Train",],
                               zTrain=y[splits=="Train"],
                               xValidation=x[splits=="Validation",],
                               zValidation=y[splits=="Validation"],
                               regressionFunction = regressionFunction,
                               ...)
  
  pred_train_cde <- FlexCoDE::predict.FlexCoDE(fit,x[splits!="Threshold",])
  t_grid <- seq(0,max(pred_train_cde$CDE),length.out = n_levels)
  g_train_cde <- matrix(NA,nrow(pred_train_cde$CDE),
                        length(t_grid))
  for(ii in 1:nrow(pred_train_cde$CDE))
  {
    g_train_cde[ii,] <- profile_density(t_grid,pred_train_cde$z,
                                        pred_train_cde$CDE[ii,])
  }
  
  kmeans_result <- try(kmeanspp(g_train_cde,k=k),
                       silent = TRUE)
  if(class(kmeans_result)=="try-error")
  {
    kmeans_result <- kmeans(g_train_cde,centers = k)
  }
  centers_kmeans <- kmeans_result$centers
  rm(g_train_cde)
  rm(pred_train_cde)
  
  
  pred_train <- FlexCoDE::predict.FlexCoDE(fit, x[splits=="Threshold",])  ## done as FlexCoDE didn't use .S3method
  
  # hpd-split
  which_select <- cbind(1:length(y[splits=="Threshold"]),
                        which_neighbors(as.matrix(pred_train$z),
                                        as.matrix(y[splits=="Threshold"]),
                                        k=1))
  which_smaller <- apply(pred_train$CDE<=pred_train$CDE[which_select],1,which)
  conformity_score_train_hpd <- rep(NA,nrow(pred_train$CDE))
  for(ii in 1:nrow(pred_train$CDE))
  {
    conformity_score_train_hpd[ii] <- sum(pred_train$CDE[ii,which_smaller[[ii]]])
  }
  band <- diff(pred_train$z)[1]
  conformity_score_train_hpd <- conformity_score_train_hpd*band
  
  
  t_grid <-seq(0,max(pred_train$CDE),length.out = n_levels)
  # observed densities:
  which_select <- cbind(1:length(y[splits=="Threshold"]),
                        which_neighbors(as.matrix(pred_train$z),
                                        as.matrix(y[splits=="Threshold"]),
                                        k=1))
  conformity_score_train <- pred_train$CDE[which_select]
  
  # Profiles
  g_train <- matrix(NA,length(conformity_score_train),
                    length(t_grid))
  for(ii in 1:length(conformity_score_train))
  {
    g_train[ii,] <- profile_density(t_grid,pred_train$z,
                                    pred_train$CDE[ii,])
  }
  
  
  cum_dist_evaluated_train <- cum_dist(pred_train$z,pred_train$CDE,
                                       y[splits=="Threshold"])
  
  return_value <- NULL
  return_value$density_fit <- fit
  return_value$cum_dist_evaluated_train <- cum_dist_evaluated_train
  return_value$conformity_score_train <- conformity_score_train
  return_value$conformity_score_train_hpd <- conformity_score_train_hpd
  return_value$t_grid <- t_grid
  return_value$band <- band
  return_value$g_train <- g_train
  return_value$centers_kmeans <- centers_kmeans
  return_value$splits <- splits
  return_value$pred_calib <- pred_train
  class(return_value) <-"predictionBands"
  return(return_value)
}




# Fit Mixture Density Network
mdn_pred_band <- function(x_train, y_train, x_calib, y_calib,
                          k = max(round(nrow(as.matrix(x_calib))/100), 1),
                          n_components = 4, 
                          hidden_dim = 30,
                          seed = 1000) {
  n_levels <- 1000

  # Mixture density network
  fit <- py$fit_mixture_density_network(x = py$np$array(x_train),
                                        y = py$np$array(matrix(y_train, ncol = 1)),
                                        n_components = py$np$int16(n_components), 
                                        hidden_dim = py$np$int16(hidden_dim),
                                        seed = py$np$int16(seed))
  pred_train_cde <- py$predict_mixture_density_network(fit,
                                                       py$np$array(x_train), 
                                                       py$np$array(matrix(y_train, ncol = 1)))
  t_grid <- seq(0, max(pred_train_cde), length.out = n_levels)
  
  y_grid <- seq(min(y_train), max(y_train), length.out = 1000)  # candidate of y
  g_train_cde <- matrix(NA, nrow(x_train), length(t_grid))
  pred_train_cde <- list(
    z = y_grid,
    CDE = matrix(NA, length(y_train), 1000)
  )
  for (i in 1:nrow(x_train)) {
    x_temp <- matrix(rep(x_train[i, ], 1000), ncol = ncol(x_train), byrow = T)
    pred_train_cde$CDE[i, ] <- py$predict_mixture_density_network(fit,
                                                                  x_temp, 
                                                                  matrix(pred_train_cde$z, ncol = 1))
    g_train_cde[i, ] <- profile_density(t_grid, pred_train_cde$z, pred_train_cde$CDE[i, ])
  }
  
  kmeans_result <- try(kmeanspp(g_train_cde, k=k),
                       silent = TRUE)
  if(class(kmeans_result)=="try-error")  {
    kmeans_result <- kmeans(g_train_cde,centers = k)
  }
  centers_kmeans <- kmeans_result$centers
  rm(g_train_cde)
  rm(pred_train_cde)
  
  
  # Conditional densities in calibration set
  pred_train <- list(
    z = y_grid,
    CDE = matrix(NA, length(y_calib), 1000)
  )
  for (i in 1:length(y_calib)) {
    # Estimate of f(y|x_i) in calibration set
    x_temp <- matrix(rep(x_calib[i, ], 1000), ncol = ncol(x_calib), byrow = T)
    pred_train$CDE[i, ] <- py$predict_mixture_density_network(fit,
                                                              x_temp, 
                                                              matrix(pred_train$z, ncol = 1))
  }
  
  
  # pred_train <- FlexCoDE::predict.FlexCoDE(fit,x[splits=="Threshold",])  ## done as FlexCoDE didn't use .S3method
  # 
  # # hpd-split
  # which_select <- cbind(1:length(y[splits=="Threshold"]),
  #                       which_neighbors(as.matrix(pred_train$z),
  #                                       as.matrix(y[splits=="Threshold"]),
  #                                       k=1))
  # which_smaller <- apply(pred_train$CDE<=pred_train$CDE[which_select],1,which)
  # conformity_score_train_hpd <- rep(NA,nrow(pred_train$CDE))
  # for(ii in 1:nrow(pred_train$CDE)) {
  #   conformity_score_train_hpd[ii] <- sum(pred_train$CDE[ii,which_smaller[[ii]]])
  # }
  # band <- diff(pred_train$z)[1]
  # conformity_score_train_hpd <- conformity_score_train_hpd*band
  # 
  # 
  # t_grid <-seq(0,max(pred_train$CDE),length.out = n_levels)
  # # observed densities:
  # which_select <- cbind(1:length(y[splits=="Threshold"]),
  #                       which_neighbors(as.matrix(pred_train$z),
  #                                       as.matrix(y[splits=="Threshold"]),
  #                                       k=1))
  # conformity_score_train <- pred_train$CDE[which_select]
  
  # Observed conditional densities f(y_i|x_i) in calibration set
  conformity_score_train <- py$predict_mixture_density_network(fit,
                                                               x_calib, 
                                                               matrix(y_calib, ncol = 1))
  
  # Profiles
  t_grid <- seq(0, max(pred_train$CDE), length.out = n_levels)
  # t_grid <- seq(0, max(sapply(pred_train, function(x){ max(x$CDE) })), length.out = n_levels)
  g_train <- matrix(NA, length(conformity_score_train), length(t_grid))
  # cum_dist_evaluated_train <- rep(NA, length(conformity_score_train))
  for (i in 1:length(conformity_score_train)) {
    g_train[i, ] <- profile_density(t_grid, pred_train$z, pred_train$CDE[i, ])
  }
  cum_dist_evaluated_train <- cum_dist(pred_train$z, pred_train$CDE,
                                       y_calib)
  
  out <- list(
    density_fit = fit,
    cum_dist_evaluated_train = cum_dist_evaluated_train,
    conformity_score_train = conformity_score_train,
    t_grid = t_grid,
    g_train = g_train,
    centers_kmeans = centers_kmeans,
    pred_calib = pred_train
  )
  class(out) <- "mdn_pred_band"
  
  return(out)
}



#' Risk Controlling Prediction Set (RCPS) + CD-split
#'
rcps_partition <- function(x_train, y_train, x_calib, y_calib, x_test, 
                           cde = "flexcode", alpha = 0.1, delta = 0.05, seed, ...) {
  if (cde == "flexcode") {
    # Flexcode
    
    set.seed(seed)
    x <- rbind(x_train, x_calib) 
    y <- c(y_train, y_calib)
    
    # Conditional density estimates for the each neighborhood A_i in training data
    # It also contains the calibration estimate.
    splits <- sample(c("Train","Validation"),
                     size = nrow(x_train),
                     prob = c(0.8, 0.2),
                     replace = TRUE)
    splits <- c(splits, rep("Threshold", nrow(x_calib)))
    # Too long time...
    cd_split_fit <- fit_predictionBands(x, y, 
                                        # per_train = 0.4,
                                        # per_val = 0.1,
                                        # per_ths = 0.5, 
                                        splits = splits,
                                        # regressionFunction.extra = list(nCores = 10),
                                        # regressionFunction = FlexCoDE::regressionFunction.NW,
                                        ...)
    # Warning message in FlexCoDE::fitFlexCoDE(xTrain = x[splits == "Train", ], zTrain = y[splits == :
    # “the optimal I found was exactly nIMax; try increasing nIMax if you want to improve performance”
    # Warning message:
    #   “did not converge in 10 iterations”
    # user   system  elapsed 
    # 7434.482 3856.160 9108.065 (서버 기준)
    
    # Predict CDE for test data
    pred_test <- FlexCoDE::predict.FlexCoDE(cd_split_fit$density_fit, x_test)
  } else if (cde == "mdn") {
    # Mixture density network with 3 hidden layers
    cd_split_fit <- mdn_pred_band(x_train, y_train, x_calib, y_calib, seed = seed, ...)
                                  # n_components = 4, hidden_dim = 30)
    # user   system  elapsed 
    # 1004.735   34.873  148.326 
    # Warning message:
    # Quick-TRANSfer stage steps exceeded maximum (= 3247500) 
    
    # Predict CDE for test data
    pred_test <- list(
      z = cd_split_fit$pred_calib$z,
      CDE = matrix(NA, length(y_test), 1000)
    )
    # pred_test <- list()
    for (i in 1:length(y_test)) {
      x_temp <- matrix(rep(x_test[i, ], 1000), ncol = ncol(x_test), byrow = T)
      pred_test$CDE[i, ] <- py$predict_mixture_density_network(cd_split_fit$density_fit,
                                                               x_temp, 
                                                               matrix(pred_test$z, ncol = 1))
    }
  }
  
  
  # Profile density estimate for test data
  g_test <- matrix(NA, nrow(x_test), length(cd_split_fit$t_grid))
  for(i in 1:nrow(x_test)){
    g_test[i, ] <- profile_density(cd_split_fit$t_grid,
                                   pred_test$z,
                                   pred_test$CDE[i, ])
  }
  # Find A_i which belongs to new data X_{n+1}
  # It uses k nearest neighbors in data of each query
  which_partition_test <- FNN::get.knnx(data = cd_split_fit$centers_kmeans,
                                        query = g_test,
                                        k = 1)$nn.index
  which_partition_train <- FNN::get.knnx(data = cd_split_fit$centers_kmeans,
                                         query = cd_split_fit$g_train,
                                         k = 1)$nn.index
  # length(which_partition_test)
  # length(which_partition_train)
  
  # # par(mfrow = c(1, 2))
  # # plot(x_test, y_test, col = which_partition_test)
  # plot(x[cd_split_fit$splits == "Threshold"], 
  #      y[cd_split_fit$splits == "Threshold"], col = which_partition_train)
  
  
  
  # Predict CDE for calibration data
  pred_calib <- cd_split_fit$pred_calib
  
  # Find optimal lambda for nested set-valued predcitor from calibration set
  if (cde == "flexcode") {
    lambda_list <- -seq(1e-8, min(apply(pred_calib$CDE, 1, max))-(1e-8), length.out = 20)  # candidates of lambda
  } else if (cde == "mdn") {
    lambda_max <- apply(pred_calib$CDE, 1, max)
    # print(lambda_max %>% sort() %>% head())
    lambda_max <- min(lambda_max[lambda_max > 1e-4])
    lambda_list <- -seq(1e-8, lambda_max-(1e-8), length.out = 20)  # candidates of lambda
  }
  
  R_hat_lambda <- rep(0, length(lambda_list))
  for (j in 1:length(which_partition_train)) {
    # Calibration set belongs to A_i
    idx_calib_in_A_i <- which(which_partition_train == which_partition_train[j])
    y_calib_A_i <- y_calib[idx_calib_in_A_i]
    
    # CDE for i-th test data
    f_hat <- data.frame(x = pred_calib$z,
                        y = pred_calib$CDE[j, ])
    
    # Nested set-valued predictors
    # Do not need to use Greedy algorithm for tolerance interval
    T_lambda <- lapply(lambda_list, function(lambda){
      sort(f_hat$x[f_hat$y > -lambda])
    })
    
    # # Nested set-valued predictors using Greedy algorithm
    # T_lambda <- list()
    # d <- 0.1   # step size rate
    # zeta <- 3  # bound of random variables in Hoeffding's inequality
    # for (i in 1:length(lambda_list)) {
    #   lambda <- lambda_list[i]
    #   T_lambda_i <- c()
    #   
    #   # iter <- 0
    #   zeta_update <- zeta
    #   while (zeta_update > -lambda) {
    #     # iter <- iter + 1
    #     zeta_update <- zeta_update - d*zeta_update
    #     T_lambda_i <- c(T_lambda_i, 
    #                     f_hat$x[(f_hat$y > zeta_update) & !(f_hat$x %in% T_lambda_i)])
    #   }
    #   T_lambda[[i]] <- sort(T_lambda_i)
    #   # print(iter)
    # }
    # # sapply(T_lambda, length)
    
    
    # Split the dis-connected region
    grid_size_density <- round(diff(f_hat$x)[1], 5)   # equal grid size of density estimate
    T_lambda_split <- lapply(T_lambda, function(tol_band){
      if (length(tol_band) == 0) {
        # No covered interval
        return( list(c(0, 0)) )
      } else if (length(tol_band) == 1) {
        # Only 1 data contained => make region
        return( list(c(tol_band - grid_size_density, 
                       tol_band + grid_size_density)) )
      }
      
      # Find connected components of bands
      dif <- round(diff(tol_band), 5)
      split_idx <- which(dif > grid_size_density)
      if (length(split_idx) > 0) {
        tol_band_split <- list()
        for (i in 1:length(split_idx)) {
          if (i == 1) {
            tol_band_split[[1]] <- tol_band[1:split_idx[1]]
          } else {
            tol_band_split[[i]] <- tol_band[(split_idx[i-1]+1):split_idx[i]]
          }
        }
        tol_band_split[[length(split_idx)+1]] <- tol_band[(split_idx[length(split_idx)]+1):length(tol_band)]
      } else {
        tol_band_split <- list(tol_band)
      }
      
      # Obtain the range of each band
      tol_band_split <- lapply(tol_band_split, range)
      
      return(tol_band_split)
    })
    
    # # Visualization
    # par(mfrow = c(1, 2))
    # plot(f_hat, type = "l")
    # for (i in 1:length(lambda_list)) {
    #   tol_band_split <- T_lambda_split[[i]]
    #   for (j in 1:length(tol_band_split)) {
    #     lines(tol_band_split[[j]], rep(-lambda_list[[i]], 2), col = i, lwd = 3)
    #   }
    # }
    # plot(f_hat[which(f_hat$y > 0), ], type = "l")
    
    
    # Empirical risk
    R_hat_lambda <- R_hat_lambda + sapply(1:length(T_lambda_split), function(i){
      number <- sapply(T_lambda_split[[i]], function(interval){
        (y_calib[j] >= interval[1] & y_calib[j] <= interval[2])
      })
      # (1 - sum(number)) / length(y_calib_A_i)
      (1 - sum(number)) / length(y_calib)
    })
  }
  
  ## Find UCB using the concentration inequality
  # Bound of concentration inequality
  bound <- sqrt(log(1/delta)/(2*length(y_calib)))  # Hoeffding's inequality
  
  # Find lambda_hat (We use max since lambda has decreasing order!)
  ucb_rule <- which(R_hat_lambda < alpha - bound)
  if (length(ucb_rule) == 0) {
    lambda_hat_order <- NA
  } else {
    lambda_hat_order <- max(ucb_rule)
  }
  lambda_hat <- lambda_list[lambda_hat_order]
  
  
  # Risk-controliling Prediction Set (RCPS)
  rcps <- list()
  for (j in 1:length(which_partition_test)) {
    # Calibration set belongs to A_i
    idx_calib_in_A_i <- which(which_partition_train == which_partition_test[j])
    y_calib_A_i <- y_calib[idx_calib_in_A_i]
    
    # CDE for i-th test data
    f_hat <- data.frame(x = pred_test$z,
                        y = pred_test$CDE[j, ])
    
    # Nested set-valued predictors with lambda_hat
    # Do not need to use Greedy algorithm for tolerance interval
    T_lambda <- sort(f_hat$x[f_hat$y > -lambda_hat])
    if (length(T_lambda) == 0) {
      # No covered interval
      T_lambda <- c(0, 0)
    } else if (length(T_lambda) == 1) {
      # Only 1 data contained => make region
      T_lambda <- c(T_lambda - grid_size_density, 
                    T_lambda + grid_size_density)
    }
    
    
    # # Nested set-valued predictors using Greedy algorithm with lambda_hat
    # d <- 0.1   # step size rate
    # zeta <- 3  # bound of random variables in Hoeffding's inequality
    # T_lambda <- c()
    # # iter <- 0
    # zeta_update <- zeta
    # while (zeta_update > -lambda_hat) {
    #   # iter <- iter + 1
    #   zeta_update <- zeta_update - d*zeta_update
    #   T_lambda <- c(T_lambda, 
    #                 f_hat$x[(f_hat$y > zeta_update) & !(f_hat$x %in% T_lambda)])
    # }
    # T_lambda <- sort(T_lambda)
    # # print(iter)
    # # length(T_lambda)
    
    # Find connected components of bands
    dif <- round(diff(T_lambda), 5)
    split_idx <- which(dif > min(dif))
    if (length(split_idx) > 0) {
      T_lambda_split <- list()
      for (i in 1:length(split_idx)) {
        if (i == 1) {
          T_lambda_split[[1]] <- T_lambda[1:split_idx[1]]
        } else {
          T_lambda_split[[i]] <- T_lambda[(split_idx[i-1]+1):split_idx[i]]
        }
      }
      T_lambda_split[[length(split_idx)+1]] <- T_lambda[(split_idx[length(split_idx)]+1):length(T_lambda)]
    } else {
      T_lambda_split <- list(T_lambda)
    }
    
    # Obtain the range of each band
    T_lambda_split <- lapply(T_lambda_split, range)  
    
    # RCPS with lambda_hat
    rcps[[j]] <- T_lambda_split
  }

  out <- list(
    cd_split_fit = cd_split_fit,
    rcps = rcps,
    pred_calib = pred_calib,
    pred_test = pred_test,
    partition = list(train = which_partition_train,
                     test = which_partition_test)
  )
  
  return(out)  
}



summary_rcps <- function(rcps, y_test, idx_bright, idx_faint) {
  df <- data.frame(
    Total = c(
      Coverage = mapply(function(rcps_i, y_i){
        sapply(rcps_i, function(interval){
          y_i >= interval[1] & y_i <= interval[2]
        }) %>% 
          sum()
      }, rcps, y_test) %>% 
        mean(),
      Size = sapply(rcps, function(rcps_i){
        sapply(rcps_i, function(interval){
          diff(interval)
        }) %>% 
          sum()
      }) %>%
        mean()
    ),
    Bright = c(
      Coverage = mapply(function(rcps_i, y_i){
        sapply(rcps_i, function(interval){
          y_i >= interval[1] & y_i <= interval[2]
        }) %>% 
          sum()
      }, rcps[idx_bright], y_test[idx_bright]) %>% 
        mean(),
      Size = sapply(rcps[idx_bright], function(rcps_i){
        sapply(rcps_i, function(interval){
          diff(interval)
        }) %>% 
          sum()
      }) %>%
        mean()
    ),
    Faint = c(
      Coverage = mapply(function(rcps_i, y_i){
        sapply(rcps_i, function(interval){
          y_i >= interval[1] & y_i <= interval[2]
        }) %>% 
          sum()
      }, rcps[idx_faint], y_test[idx_faint]) %>% 
        mean(),
      Size = sapply(rcps[idx_faint], function(rcps_i){
        sapply(rcps_i, function(interval){
          diff(interval)
        }) %>% 
          sum()
      }) %>%
        mean()
    )
  )
  return(df)
}