## functions for generating random functional data and missing periods

simul.fd = function(n = 200, grid = seq(0,1,len=200), lambda.cos = 3^(-(2*(1:300)-1)), lambda.sin = 3^(-(2*(1:300))), randcoef = norm.randcoef)
{
  x = matrix(0,n,length(grid))
  R = matrix(0,length(grid),length(grid))
  for (j in 1:length(lambda.cos)) {
    f = sqrt(lambda.cos[j])*sqrt(2)*cos(2*pi*j*grid)
    x = x + randcoef(n)%*%t(f)
    R = R + f%*%t(f)
  }
  for (j in 1:length(lambda.sin)) {
    f = sqrt(lambda.sin[j])*sqrt(2)*sin(2*pi*j*grid)
    x = x + randcoef(n)%*%t(f)
    R = R + f%*%t(f)
  }
  attr(x,"R") = R
  x
}

norm.randcoef = function(n) rnorm(n,0,1)
unif.randcoef = function(n) runif(n,-1,1)*sqrt(3)
t5.randcoef = function(n) rt(n,5)/sqrt(5/3)

simul.obs = function(n = 100, grid = seq(0,1,len=200),d=1.4,f=.2)
{
  out = matrix(TRUE,n,length(grid))
  for (i in 1:n) {
    cen = d*sqrt(runif(1))
    e = f*runif(1)
    out[i,(cen-e<grid)&(grid<cen+e)] = FALSE # missing interval = (cen-u,cen+u)
  }
  out
}
