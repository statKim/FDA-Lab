library(MASS)

### penalized spline for each curve
pspline_curve <- function(t, x) {
  grid.length <- length(x)
  x_sm <- rep(NA, grid.length)
  NA_ind <- which(is.na(x))   # check NA
  
  if (length(NA_ind) > 0) {
    # Partially observed curve
    ind_1 <- min(NA_ind)
    ind_2 <- max(NA_ind)
    if (ind_1 == 1 || ind_2 == grid.length) {
      # missing for start or end
      t_obs <- t[-NA_ind]
      x_obs <- x[-NA_ind]
      x_sm[-NA_ind] <- pspline::sm.spline(t_obs, x_obs)$ysmth
    } else {
      # missing for middle
      # concatenate for 2 parts of observed periods, smoothing, and seperate
      t_obs <- t[-NA_ind]
      x_obs <- x[-NA_ind]
      sm_fit <- pspline::sm.spline(t_obs, x_obs)
      x_sm[1:(ind_1-1)] <- sm_fit$ysmth[1:(ind_1-1)]
      x_sm[(ind_2+1):grid.length] <- sm_fit$ysmth[-(1:(ind_1-1))]
      
      # t_obs_1 <- t[1:(ind_1-1)]
      # x_obs_1 <- x[1:(ind_1-1)]
      # sm_fit <- sm.spline(t_obs_1, x_obs_1)
      # x_sm[1:(ind_1-1)] <- as.numeric(sm_fit$ysmth)
      # 
      # t_obs_2 <- t[(ind_2+1):grid.length]
      # x_obs_2 <- x[(ind_2+1):grid.length]
      # sm_fit <- sm.spline(t_obs_2, x_obs_2)
      # x_sm[(ind_2+1):grid.length] <- as.numeric(sm_fit$ysmth)
    }
  } else {
    # Fully observed curve
    sm_fit <- pspline::sm.spline(t, x)
    x_sm <- as.numeric(sm_fit$ysmth)
  }
  
  return(x_sm)
}


## main programs for principal component analysis and reconstruction
## of incomplete functional data using regularised linear prediction

# functional data are evaluated on an equidistant grid of points
# in an interval
# replications are stored in rows of a matrix
# d.t is the distance between the grid points (i.e., length of the
# domain divided by the number of columns)
# if d.t is unspecified, the the length of the domain is assumed
# to be 1 and d.t is set to 1/ncol(x)

#####################################################################
### (1) mean function, covariance operator and its eigendecomposition
### for possibly incomplete functional data
#####################################################################
### M-estimator for covaraince function
mean.rob.missfd <- function(x, smooth = F, work.grid = NULL) {
  if (smooth == T) {
    x.2 <- matrix2list(x)
    if (is.null(work.grid)) {
      work.grid <- seq(0, 1, length.out = ncol(x))
    }
    mu.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, method = "huber", 
                           kernel = "epanechnikov",
                           bw = 0.2, delta = 1.345)
    mu <- predict(mu.obj, work.grid)
  } else {
    if (is.list(x)) {
      x <- list2matrix(x)
    }
    mu <- apply(x, 2, function(t){ MASS::huber(t)$mu })
  }

  return(mu)
}

var.rob.missfd <- function(x, smooth = F, 
                           bw = 0.1,
                           make.pos.semidef = TRUE, noise.var = 0) {
  n <- nrow(x)
  p <- ncol(x)
  rob.var <- matrix(0, p, p)
  
  if (sum(is.na(x)) == 0) {
    mu <- mean.rob.missfd(x, smooth = smooth)
    # complete curves
    for (s in 1:p) {
      for (t in 1:p) {
        if (s <= t) {
          A <- (x[, s] - mu[s])*(x[, t] - mu[t])
          rob.var[s, t] <- huber(A)$mu
        } else {
          rob.var[s, t] <- rob.var[t, s]
        }
      }
    }
  } else {
    # include partially observed curves
    for (s in 1:p) {
      for (t in 1:p) {
        if (s <= t) {
          NA_ind <- apply(x, 1, function(xi){ sum(is.na(xi[c(s, t)])) > 0 })
          if (sum(NA_ind) > 0) {
            ind <- which(NA_ind == FALSE)
          } else {
            ind <- 1:n
          }
          
          mu <- mean.rob.missfd(x[ind, ], smooth = smooth, work.grid = c(s, t))
          A <- (x[ind, s] - mu[1])*(x[ind, t] - mu[2])
          rob.var[s, t] <- huber(A)$mu
        } else {
          rob.var[s, t] <- rob.var[t, s]
        }
      }
    }
  }
  
  # rob.var <- matrix(0, p, p)
  # for (s in 1:p) {
  #   for (t in 1:p) {
  #     if (s <= t) {
  #       NA_ind <- apply(x, 1, function(xi){ sum(is.na(xi[c(s, t)])) > 0 })
  #       if (sum(NA_ind) > 0) {
  #         ind <- which(NA_ind == FALSE)
  #       } else {
  #         ind <- 1:n
  #       }
  #       
  #       mu <- mean.rob.missfd(x[ind, ], smooth = smooth, work.grid = c(s, t))
  #       A <- (x[ind, s] - mu[1])*(x[ind, t] - mu[2])
  #       rob.var[s, t] <- huber(A)$mu
  #     } else {
  #       rob.var[s, t] <- rob.var[t, s]
  #     }
  #   }
  # }
  
  # subtract noise variance
  diag(rob.var) <- diag(rob.var) - noise.var
  
  # 2-dimensional smoothing
  if (smooth == T) {
    gr <- seq(0, 1, length.out = p)
    rob.var <- fields::smooth.2d(as.numeric(rob.var),
                                 x = expand.grid(gr, gr), surface = F,
                                 theta = bw, nrow = p, ncol = p)
  }
  
  # make positive-semi-definite
  if (isTRUE(make.pos.semidef)) {
    eig <- eigen(rob.var)
    k <- which(eig$values > 0)
    rob.var <- eig$vectors[, k] %*% diag(eig$values[k]) %*% t(eig$vectors[, k])
  }
  
  return(rob.var)
}


eigen.missfd = function(R,d.t,...)
{
  # d.t is the distance between the points of the equidistant grid
  out = eigen(R,symmetric=TRUE,...)
  if (missing(d.t)) d.t = 1/ncol(R)
  out$values = out$values*d.t
  out$vectors = out$vectors/sqrt(d.t)
  out
}



#############################################################
### (2) prediction of scores of a partially observed function
#############################################################
pred.score.rob.missfd = function(x1,phi,x,R,mu,n,C,ncomp,alpha,d.t,
                             gcv.df.factor=1,gcv.plot=FALSE,
                             gcv.print=FALSE,
                             smooth = FALSE)
{
  # x1 is one partially observed curve
  # this function predicts the Fourier scores of x1
  # with respect to phi (i.e., inner products of x1 and phi)
  # phi is a matrix containing functions in columns
  # if missing, phi defaults to the first four principal components
  # x is a functional data set (functions in rows)
  # R is the covar oper, mu the mean fun, n the sample size of x
  # C, ncomp are the covar oper and sample size for the set of
  # complete functions
  # R, mu, n, C, ncomp are optional, may be precomputed
  # beforehand, e.g., to save time when this function is to be
  # called repeatedly (for different incomplete x1)
  # alpha is the ridge regularisation parameter, determined by gcv
  # from complete curves if missing
  # d.t is the distance between the points of the equidistant grid
  miss = is.na(x1)
  obs = !miss
  k.O = sum(obs)
  k.M = sum(miss)
  k = k.O + k.M
  if (missing(d.t)) d.t = 1/k
  if (missing(phi)) {
    R = var.rob.missfd(x, smooth = smooth)
    phi = eigen.missfd(R)$vectors[,1:4]
  }
  phi = as.matrix(phi)
  nscores = ncol(phi)
  if (missing(mu)) {
    mu = mean.rob.missfd(x, smooth = smooth)
    n = nrow(x)
  }
  if (missing(R)) {
    R = var.rob.missfd(x, smooth = smooth)
    n = nrow(x)
  }
  if (missing(C)) {
    comp = apply(!is.na(x),1,all)
    ncomp = sum(comp)
    C = var.rob.missfd(x[comp,], smooth = smooth) # covariance operator of complete curves
  }
  
  if (k.M==0) { # complete observation, no prediction needed
    b = rep(NA,nscores)
    for (j in 1:nscores) {
      b[j] = sum((x1[obs]-mu[obs])*phi[obs,j])*d.t
    }
  } else {
    b = rep(NA,nscores)
    b.details = matrix(NA,7,nscores)
    rownames(b.details) = c("score.obs","score.miss","se","relerr","alpha","df.alpha","varprop.alpha")
    
    # precompute quantities that don't depend on alpha
    # to speed up the optimisation in gcv
    cMM = colSums(phi[miss,,drop=FALSE]*(C[miss,miss,drop=FALSE]%*%phi[miss,,drop=FALSE]))*(d.t^2)
    eig.ROO = eigen.missfd(R[obs,obs,drop=FALSE],d.t=d.t)
    lambda.ROO = eig.ROO$values
    psi.ROO = eig.ROO$vectors
    # rO = crossprod(R[obs,miss,drop=FALSE],phi[miss,,drop=FALSE])*d.t # columns = righ-hand sides of the inverse problem
    psirO = crossprod(psi.ROO,(R[obs,miss,drop=FALSE]%*%phi[miss,,drop=FALSE])*d.t)*d.t # projection of the right-hand sides on psi.ROO
    psirO.psicO = psirO*crossprod(psi.ROO,(C[obs,miss,drop=FALSE]%*%phi[miss,,drop=FALSE])*d.t)*d.t
    psirOrOpsi.psiCOOpsi = array(0,c(k.O,k.O,nscores))
    for (j in 1:nscores) {
      psirOrOpsi.psiCOOpsi[,,j] = tcrossprod(psirO[,j])*crossprod(psi.ROO,C[obs,obs,drop=FALSE]%*%psi.ROO)*d.t^2
    }
    
    if (missing(alpha)) { # gcv selection of alpha for each score
      alpha = rep(NA,nscores)
      if (gcv.plot) {
        a = ceiling(sqrt(nscores)) # number of columns in the plot
        op = par(mfrow=c(ceiling(nscores/a),a))
      }
      for (j in 1:nscores) {
        # gcv function for scores and for functions is the same
        # just use it with the right input parameters
        alpha[j] = alpha.gcv.pred.rob.missfd(trace.CMM=cMM[j],
                                         diag.psiROMCMOpsi=psirO.psicO[,j],
                                         psiROMMOpsi.psiCOOpsi=psirOrOpsi.psiCOOpsi[,,j],
                                         lambda.ROO=lambda.ROO,n=n,
                                         ncomp=ncomp,df.factor=gcv.df.factor,
                                         gcv.plot=gcv.plot,gcv.print=gcv.print)
        if (gcv.plot) title(main=j)
      }
      if (gcv.plot) par(op)
    } else {
      if (length(alpha)==1) alpha = rep(alpha,nscores)
    }
    
    rMM = colSums(phi[miss,,drop=FALSE]*(R[miss,miss,drop=FALSE]%*%phi[miss,,drop=FALSE]))*(d.t^2)
    for (j in 1:nscores) {
      b.details["score.obs",j] = sum((x1[obs]-mu[obs])*phi[obs,j])*d.t
      b.details["score.miss",j] = sum(psirO[,j]*crossprod(psi.ROO,x1[obs]-mu[obs])/(lambda.ROO+alpha[j]))*d.t
      b[j] = b.details["score.obs",j] + b.details["score.miss",j]
      b.details["se",j] = sqrt( rMM[j] - sum(psirO[,j]^2*lambda.ROO/(lambda.ROO+alpha[j])^2) )
      b.details["relerr",j] = b.details["se",j]/sqrt(sum(phi[,j]*(R%*%phi[,j]))*(d.t)^2)
      b.details["alpha",j] = alpha[j]
      b.details["df.alpha",j] = sum(lambda.ROO/(lambda.ROO+alpha[j]))
      b.details["varprop.alpha",j] = sum(lambda.ROO^3/(lambda.ROO+alpha[j])^2)/sum(lambda.ROO)
    }
    attr(b,"details") = b.details
  }
  
  b  
}





################################################################
### (3) prediction of the missing part of an incomplete function
################################################################
################################################################
pred.rob.missfd = function(x1,x,R,mu,n,C,ncomp,alpha,d.t,
                       gcv.df.factor=1,gcv.plot=FALSE,gcv.print=FALSE,
                       smooth = FALSE)
{
  # x1 is one partially observed curve
  # this function predicts the missing part of x1
  # x is a functional data set
  # R is the covar oper, mu the mean fun, n the sample size of x
  # C, ncomp covar oper and sample size for the set of
  # complete functions
  # R, mu, n, C, ncomp are optional, may be precomputed
  # beforehand, e.g., to save time when this function is to be
  # called repeatedly (for different incomplete x1)
  # alpha is the ridge regularisation parameter, determined by gcv
  # from complete curves if missing
  # d.t is the distance between the points of the equidistant grid
  miss = is.na(x1)
  obs = !miss
  k.O = sum(obs)
  k.M = sum(miss)
  k = k.O + k.M
  if (missing(d.t)) d.t = 1/k
  if (missing(mu)) {
    mu = mean.rob.missfd(x, smooth = smooth)
    n = nrow(x)
  }
  if (missing(R)) {
    R = var.rob.missfd(x, smooth = smooth)
    n = nrow(x)
  }
  if (missing(C)) {
    comp = apply(!is.na(x),1,all)
    ncomp = sum(comp)
    C = var.rob.missfd(x[comp,], smooth = smooth) # covariance operator of complete curves
  }
  
  x1.pred = rep(NA,k) # this will contain predicted x1 on miss and NAs on obs
  x1.pred.covar = matrix(NA,k,k) # covar oper of the predictive distribution
  
  if (k.M>0) {
    # precompute quantities that don't depend on alpha
    # to speed up the optimisation in gcv
    trace.CMM = sum(diag(C)[miss])*d.t
    eig.ROO = eigen.missfd(R[obs,obs,drop=FALSE],d.t=d.t)
    lambda.ROO = eig.ROO$values
    psi.ROO = eig.ROO$vectors
    psiROM = crossprod(psi.ROO,R[obs,miss,drop=FALSE])*d.t
    diag.psiROMCMOpsi = rowSums(psiROM*crossprod(psi.ROO,C[obs,miss,drop=FALSE])*d.t)*d.t
    psiROMMOpsi.psiCOOpsi = tcrossprod(psiROM)*d.t*crossprod(psi.ROO,C[obs,obs,drop=FALSE]%*%psi.ROO)*d.t^2
    
    # if alpha not provided on input, use gcv
    if (missing(alpha)) {
      alpha = alpha.gcv.pred.rob.missfd(trace.CMM=trace.CMM,
                                    diag.psiROMCMOpsi=diag.psiROMCMOpsi,
                                    psiROMMOpsi.psiCOOpsi=psiROMMOpsi.psiCOOpsi,
                                    lambda.ROO=lambda.ROO,n=n,
                                    ncomp=ncomp,df.factor=gcv.df.factor,
                                    gcv.plot=gcv.plot,gcv.print=gcv.print)
    }
    
    # compute the regularised prediction and predictive covariance operator
    x1.pred[miss] = crossprod(psiROM,crossprod(psi.ROO,x1[obs]-mu[obs])/(lambda.ROO+alpha))*d.t + mu[miss]
    x1.pred.covar[miss,miss] = R[miss,miss] - crossprod(psiROM*(lambda.ROO/(lambda.ROO+alpha)^2),psiROM)
    
    attributes(x1.pred) = list(covar=x1.pred.covar,se=sqrt(diag(x1.pred.covar)),
                               relerr=sqrt(sum(diag(x1.pred.covar)[miss])/sum(diag(R))),
                               alpha=alpha,df.alpha=sum(lambda.ROO/(lambda.ROO+alpha)),
                               varprop.alpha=sum(lambda.ROO^3/(lambda.ROO+alpha)^2)/sum(lambda.ROO))
  } else { # complete curve, no prediction
    x1.pred = x1
  }
  
  x1.pred
}

gcv.pred.rob.missfd = function(log.alpha, trace.CMM, diag.psiROMCMOpsi,
                           psiROMMOpsi.psiCOOpsi, lambda.ROO, ncomp,
                           df.factor=1, prn=FALSE)
{
  # gcv on complete curves
  # this version is fully based on eigendecomposition
  alpha = exp(log.alpha)
  gof.term = trace.CMM - 2*sum(diag.psiROMCMOpsi/(lambda.ROO+alpha)) + sum(psiROMMOpsi.psiCOOpsi*tcrossprod(1/(lambda.ROO+alpha)))
  df.alpha = sum(lambda.ROO/(lambda.ROO+alpha))
  out = log(gof.term) - 2*log(1-df.factor*df.alpha/ncomp)
  if (prn) print(c(alpha=alpha,gof.term=gof.term,df.alpha=df.alpha,log.gcv=out))
  out
}

alpha.gcv.pred.rob.missfd = function(trace.CMM, diag.psiROMCMOpsi,
                                 psiROMMOpsi.psiCOOpsi, lambda.ROO,
                                 n, ncomp, df.factor=1, gcv.plot=FALSE,
                                 gcv.print=FALSE)
{
  # scale (standardise) everything
  std = lambda.ROO[1]
  trace.CMM = trace.CMM/std
  lambda.ROO = lambda.ROO/std
  diag.psiROMCMOpsi = diag.psiROMCMOpsi/(std^2)
  psiROMMOpsi.psiCOOpsi = psiROMMOpsi.psiCOOpsi/(std^3)
  # find minimum alpha so that df <= df.max
  nn = sum(lambda.ROO>.Machine$double.eps^.5) # rank of ROO
  if (nn>=length(lambda.ROO)) { # ROO full rank, alpha small OK
    log.alpha.min = log(.Machine$double.eps^.5)
  } else {
    df.max = min(nn/4) # maximum df that we will allow in optim
    # now need to find alpha corresponding to df.max
    log.alpha.min = log(.Machine$double.eps^.5)
    m1 = log.alpha.min
    m2 = log(n/df.max)
    if (df.alpha.eq(log.alpha=log.alpha.min,deg.fr=df.max,lambda=lambda.ROO)>0) {
      # find log alpha between m1, m2 such that df equals df.max
      # df at m2 certainly < df.max
      # if df at m1 > df.max, there exists alpha such that df equals df.max
      # (opposite signs at end points), we find it using uniroot
      # (if same signs, we keep m1)
      log.alpha.min = ( uniroot(df.alpha.eq,c(m1,m2),deg.fr=df.max,lambda=lambda.ROO)$root )
    }
  }
  # minimise gcv
  log.alpha = optim(max(log(mean(lambda.ROO)),log.alpha.min),
                    gcv.pred.rob.missfd,NULL,
                    trace.CMM=trace.CMM,
                    diag.psiROMCMOpsi=diag.psiROMCMOpsi,
                    psiROMMOpsi.psiCOOpsi=psiROMMOpsi.psiCOOpsi,
                    lambda.ROO=lambda.ROO,ncomp=ncomp,
                    df.factor=df.factor,method="Brent", upper=99, #"L-BFGS-B",
                    # control=list(trace=0),
                    lower=log.alpha.min)$par
  # plot gcv on a grid of alpha values
  if (gcv.plot) {
    log.alpha.max = log.alpha+3
    log.alpha.grid = seq(log.alpha.min,log.alpha.max,len=40)
    log.gcv.grid = double(length(log.alpha.grid))
    for (i in 1:length(log.gcv.grid)) {
      log.gcv.grid[i] = gcv.pred.rob.missfd(log.alpha=log.alpha.grid[i],
                                        trace.CMM=trace.CMM,
                                        diag.psiROMCMOpsi=diag.psiROMCMOpsi,
                                        psiROMMOpsi.psiCOOpsi=psiROMMOpsi.psiCOOpsi,
                                        lambda.ROO=lambda.ROO,ncomp=ncomp,
                                        df.factor=df.factor,
                                        prn=gcv.print) #+log(std)
    }
    plot(log.alpha.grid+log(std),log.gcv.grid+log(std),type="l",xlab="log(alpha)",ylab="log(gcv(alpha))")
    abline(v=log.alpha.min+log(std),lty=3)
    abline(v=log.alpha+log(std),col=2)
  }
  exp(log.alpha)*std
}

df.alpha.eq = function(log.alpha,deg.fr,lambda)
{
  sum(lambda/(lambda+exp(log.alpha))) - deg.fr
}

#######################################################################
### (4) prediction bands for the missing part of an incomplete function
#######################################################################

qsupgp2 = function(p,R,h=pmax(sqrt(diag(R)),.2*sqrt(max(diag(R)))),nsim=2500)
{
  # simulate the p-th quantile of sup(|X|/h) where X is a mean zero
  # Gaussian process with covariance function R
  eig.R = eigen(R,symmetric=TRUE)
  npos = sum(eig.R$values>0)
  y = (eig.R$vectors[,1:npos]*rep(sqrt(eig.R$values[1:npos]),each=nrow(R))) %*% matrix(rnorm(npos*nsim),npos,nsim) # simulated curves are in columns
  quantile(apply(abs(y)/h,2,max), probs=p)
}

gpband = function(mu,R,h=pmax(sqrt(diag(R)),.2*sqrt(max(diag(R),na.rm=TRUE))),coverage=.95)
{
  # for a Gaussian process with mean mu and covariance R
  # this function computes a band of the form mu +- u*h that
  # contains a trajectory with probability coverage,
  # h is the boundary function for the band (h=1 for const width,
  # default h reflects the sd but is bounded away from 0)
  # the band is constructed only in regions where mu is not NA
  w = !is.na(mu)
  if (length(h)==1) h = rep(h,length(mu))
  u = qsupgp2(p=coverage,R=R[w,w,drop=FALSE],h=h[w])
  out = cbind(lower=mu-u*h,upper=mu+u*h)
}
