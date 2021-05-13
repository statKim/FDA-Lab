##########################
#### CD4 data example
##########################
require(refund)
library(face)
data(cd4)
n <- nrow(cd4)
T <- ncol(cd4)
id <- rep(1:n,each=T)
t <- rep(-18:42,times=n)
y <- as.vector(t(cd4))
sel <- which(is.na(y))

## organize data and apply FACEs
data <- data.frame(y=log(y[-sel]),
                   argvals = t[-sel],
                   subj = id[-sel])
data <- data[data$y>4.5,]
fit_face <- face.sparse(data,argvals.new=(-20:40))
data.h <- data
tnew <- fit_face$argvals.new


## scatter plots
Xlab <- "Months since seroconversion"
Ylab <- "log (CD4 count)"
par(mfrow=c(1,1),mar = c(4.5,4.5,3,2))
id <- data.h$subj
uid <- unique(id)
plot(data.h$argvals,data.h$y,
     type = "n", ylim = c(4.5,8),
     xlab = Xlab, ylab = Ylab,
     cex.lab = 1.25,cex.axis=1.25,cex.main = 1.25)
for(i in 1:10){
  seq <- which(id==uid[i])
  lines(data.h$argvals[seq],data.h$y[seq],lty=1,col="gray",lwd=1,type="l")
  #points(data.h$argvals[seq],data.h$y[seq],col=1,lty=1,pch=1)
}
Sample <- seq(10,50,by=10)
for(i in Sample){
  seq <- which(id==uid[i])
  lines(data.h$argvals[seq],data.h$y[seq],lty=1,col="black",lwd=1,type="l")
}
lines(tnew,fit_face$mu.new,lwd=2,lty=2,col="red")


## plots of variance/correlation functions
Cov <- fit_face$Chat.new
Cov_diag <- diag(Cov)
Cor <- fit_face$Cor.new
par(mfrow=c(1,2),mar=c(4.5,4.1,3,4.5))
plot(tnew,Cov_diag,type="l",
     xlab = Xlab, ylab="",main= "CD4: variance function",
     #ylim = c(0.8,1.5),
     cex.axis=1.25,cex.lab=1.25,cex.main=1.25,lwd=2)
require(fields)
image.plot(tnew,tnew,Cor,
           xlab=Xlab, ylab = Xlab,
           main = "CD4: correlation function",
           cex.axis=1.25,cex.lab=1.25,cex.main=1.25,
           axis.args = list(at = c(0,0.2,0.4,0.6,0.8,1.0)),
           legend.shrink=0.75,legend.line=-1.5)


## prediction of several subjects
par(mfrow=c(2,2),mar=c(4.5,4.5,3,2))
Sample <- c(30,40,50,60)
for(i in 1:4){
  sel <- which(id==uid[Sample[i]])
  dati <- data.h[sel,]
  seq <- -20:40
  k <- length(seq)
  dati_pred <- data.frame(y = rep(NA,nrow(dati) + k ),
                          argvals = c(rep(NA,nrow(dati)),seq),
                          subj=rep(dati$subj[1],nrow(dati) + k )
  )
  dati_pred[1:nrow(dati),] <- dati
  yhat2 <- predict(fit_face,dati_pred)
  data3 <- dati
  Ylim <- range(c(data3$y,yhat2$y.pred))
  plot(data3$argvals,data3$y,xlab=Xlab,ylab=Ylab, main = paste("Male ",i,sep=""),
       ylim = c(4,8.5),
       cex.lab=1.25,cex.axis = 1.25,cex.main = 1.25,pch=1,xlim=c(-20,40))
  Ord <- nrow(dati) + 1:k
  lines(dati_pred$argvals[Ord],yhat2$y.pred[Ord],col="red",lwd=2)
  lines(dati_pred$argvals[Ord],
        yhat2$y.pred[Ord] - 1.96*yhat2$se.pred[Ord], col="red",lwd=1,lty=2)
  lines(dati_pred$argvals[Ord],
        yhat2$y.pred[Ord] + 1.96*yhat2$se.pred[Ord], col="red",lwd=1,lty=2)
  lines(tnew,fit_face$mu.new,lty=3,col="black",lwd=2)
  legend("bottomleft",c("mean","prediction"),lty=c(3,1),col=1:2,lwd=2,bty="n")
}




require(refund)
data(cd4)
n <- nrow(cd4)
T <- ncol(cd4)
id <- rep(1:n,each=T)
t <- rep(-18:42,times=n)
y <- as.vector(t(cd4))
sel <- which(is.na(y))
## organize data
data <- data.frame(y=log(y[-sel]),
                   argvals = t[-sel],
                   subj = id[-sel])
data <- data[data$y>4.5,]
## smooth
fit <- pspline(data)
## plot
plot(data$argvals,fit$mu.new,type="p")
## prediction
pred <- predict(fit,quantile(data$argvals,c(0.2,0.6)))
pred
