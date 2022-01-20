# This file includes the following functions (see the Documentation file):

# mvquantile
# qpc
# outlqpc

# =========================================================
# mvquantile: function to compute directional quantiles
# =========================================================

# fd ............ (n * d) matrix of sample points
# alpha ......... vector of values of alpha
# u ............. (m * d) matrix of projection directions (m>1)


mvquantile <- function(fd, alpha, u) {
  d <- dim(fd)[2]
  dire <- dim(u)[1]
  lalp <- length(alpha)
  q <- array(dim = c(dire, d, lalp))
  for (i in 1:dire) {
    univ <-
      apply((fd - matrix(
        apply(fd, 2, mean),
        nr = dim(fd)[1],
        nc = d,
        byr = T
      )) * matrix(
        u[i, ],
        nr = dim(fd)[1],
        nc = d,
        byr = T
      ), 1, sum)
    qq <- quantile(univ, alpha)
    q[i, , ] <- outer(u[i, ], qq) + apply(fd, 2, mean)
  }
  return(q)
}


# =========================================================
# qpc: Function to compute the first principal quantile direction
# =========================================================

# fd ...... data (n x m)
# alpha ... alpha-principal component
# np ...... number of proyections
# ct ...... center "mean" or "median"
# depth ... deepest curve "FM", "RP", "RPD" (see library fda.usc)

require(fda.usc)

qpc <- function(fd,
                alpha = 0.5,
                np = 20000,
                ct = "median",
                depth = "FM") {
  if (ct == "mean") {
    fd_ct <- func.mean(fd)$data
  }
  else{
    func <- paste("depth.", depth, sep = "")
    fd_ct <- do.call(func, list(fd))$median$data
  }
  n <- nrow(fd)
  m <- ncol(fd)
  data_ct <- fd - matrix(as.numeric(fd_ct),
                         nr = n,
                         nc = m,
                         byr = T)
  
  # Random projections
  qq <- numeric()
  zr <- matrix(nrow = np, ncol = m)
  for (i in 1:np) {
    z <- rnorm(m)
    modulo <- sum(z ^ 2)
    z <- z / sqrt(modulo)
    zr[i, ] <- z
    dfproj <- data_ct %*% z
    qq[i] <- quantile(dfproj, alpha)
  }
  
  # Principal quantile direction
  umax <- which.max(abs(qq))
  pqc <- zr[umax, ]
  dfproj <- data_ct %*% pqc
  return(list(
    "pqc" = pqc,
    "proj" = dfproj,
    "data" = fd,
    "fd_ct" = fd_ct,
    "alpha" = alpha
  ))
}


# =========================================================
# outlqpc: Function to compute the ouliers from the quantile principal direction
# =========================================================

# fd ...... data (n x m)
# alpha ... alpha-principal component
# np ...... number of proyections
# ct ...... center "mean" or "median"
# depth ... deepest curve "FM", "RP", "RPD"


outlqpc <- function(fd,
                    alpha = 0.5,
                    np = 20000,
                    ct = "median",
                    depth = "FM") {
  qpc_out <- qpc(fd, alpha, np, ct, depth)
  
  qA <- quantile(qpc_out$proj, 0.25)
  qB <- quantile(qpc_out$proj, 0.75)
  IR <- qB - qA
  k <- 1.5
  outL <- (qpc_out$proj < max(min(qpc_out$proj), qA - k * IR))
  outU <- (qpc_out$proj > min(max(qpc_out$proj), qB + k * IR))
  outl <- outL | outU
  outl <- which(outl)
  outl <- row.names(fd)[unique(outl)]
  return(
    list(
      "outliers" = outl,
      "pqc" = qpc_out$pqc,
      "proj" = qpc_out$proj,
      "data" = fd,
      "fd_ct" = qpc_out$fd_ct,
      "alpha" = alpha
    )
  )
}
