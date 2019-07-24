# ------------------------------
#    Canadian weather data
# ------------------------------
setwd("C:\\Users\\hyunsungKim\\Desktop\\Thesis\\FDA\\fda_packages\\data")

library(fda)
names(CanadianWeather)
names(daily)

# sse <- vector()
# for (i in 1:15) {
#   n <- 2*i+1   # 홀수
#   daybasis <- create.fourier.basis(rangeval=c(0, 365), nbasis=n)
#   tt <- smooth.basis(day.5,
#                      CanadianWeather$dailyAv[,,"Temperature.C"][,c("Pr. Rupert", "Montreal", "Edmonton", "Resolute")],
#                      daybasis, fdnames=list("Day", "Station", "Deg C"))
#   sse[i] <- tt$SSE
# }

basis <- create.fourier.basis(rangeval=c(0, 365), nbasis=11)

tt <- smooth.basis(monthMid,
                   CanadianWeather$monthlyTemp[,c("Pr. Rupert", "Montreal", "Edmonton", "Resolute")],
                   basis, fdnames=list("Month", "Station", "Temperature"))
plot(tt, xaxt="n", main="Temperature of 4 stations in canada")
axis(side=1, at=c(0, tt$argvals[c(2,4,6,8,10,12)]), labels=c(0,2,4,6,8,10,12))
text(20, 5.7, adj=0.5, 'PR.rupert')
text(20, -8, adj=0.5, 'Montreal', col=2)
text(20, -16, adj=0.5, 'Edmonton', col=3)
text(20, -30, adj=0.5, 'Resolute', col=4)

tt <- smooth.basis(monthMid,
                   CanadianWeather$monthlyTemp,
                   basis, fdnames=list("Month", "Station", "Temperature"))

# set up the harmonic acceleration operator
harmaccelLfd365 <- vec2Lfd(c(0,(2*pi/365)^2,0), c(0, 365))

harmfdPar <- fdPar(basis, harmaccelLfd365, 1e5)

daytemppcaobj <- pca.fd(tt$fd, nharm=4, harmfdPar)

op <- par(mfrow=c(2,2), pty="m")
plot.pca.fd(daytemppcaobj)

names(daytemppcaobj)
daytemppcaobj$harmonics
daytemppcaobj$values
daytemppcaobj$scores
daytemppcaobj$meanfd

par(mfrow=c(2,1))
plot(daytemppcaobj$meanfd, main="Mean function", xaxt="n")
axis(side=1, at=c(0, tt$argvals[c(2,4,6,8,10,12)]), labels=c(0,2,4,6,8,10,12))
plot(daytemppcaobj$harmonics, main="4 PC weight functions", ylab="weight", xaxt="n")
axis(side=1, at=c(0, tt$argvals[c(2,4,6,8,10,12)]), labels=c(0,2,4,6,8,10,12))
legend('bottomleft', c("PC1","PC2","PC3","PC4"), lty=1:4, col=1:4, cex=0.5)

## varimax rotation
daytemppcaobjVM <- varmx.pca.fd(daytemppcaobj, nharm=4)
str(daytemppcaobjVM)
dimnames(daytemppcaobjVM$scores)[[2]] <- paste("PC", 1:4, sep="")
round(daytemppcaobjVM$scores)

#  plot harmonics
op <- par(mfrow=c(2,2), pty="m")
plot.pca.fd(daytemppcaobjVM)
par(op)

# plotscores(daytemppcaobj, daytemppcaobj$harmonics)
plot.basisfd(daytemppcaobj$harmonics$coefs[,1])


# #  plot log eigenvalues
# daytempeigvals <- daytemppcaobjVM[[2]]
# 
# plot(1:10, log10(daytempeigvals[1:10]), type="b",
#      xlab="Eigenvalue Number", ylab="Log 10 Eigenvalue")
# abline(lsfit(5:10, log10(daytempeigvals[5:10])), lty=2)

#  plot factor scores
harmscr <- daytemppcaobj[[3]]

plot(harmscr[,1], harmscr[,2], xlab="PC 1", ylab="PC 2", pch=".")
text(harmscr[,1], harmscr[,2], CanadianWeather$place)





#### 1~12로 했을 때
basis <- create.fourier.basis(rangeval=c(0, 12), nbasis=11)

tt <- smooth.basis(1:12,
                   CanadianWeather$monthlyTemp[,c("Pr. Rupert", "Montreal", "Edmonton", "Resolute")],
                   basis, fdnames=list("Month", "Station", "Temperature"))
plot(tt)
# axis(side=1, at=tt$argvals, labels=1:12)
# axis(side=1, at=c(0, tt$argvals[c(2,4,6,8,10,12)]), labels=c(0, 2,4,6,8,10,12))

tt <- smooth.basis(1:12,
                   CanadianWeather$monthlyTemp,
                   basis, fdnames=list("Month", "Station", "Temperature"))

# set up the harmonic acceleration operator
harmaccelLfd365 <- vec2Lfd(c(0,(2*pi/12)^2,0), c(0, 12))

harmfdPar <- fdPar(basis, harmaccelLfd365, 1e5)

daytemppcaobj <- pca.fd(tt$fd, nharm=4, harmfdPar)

op <- par(mfrow=c(2,2), pty="m")
plot.pca.fd(daytemppcaobj)

names(daytemppcaobj)
daytemppcaobj$harmonics
daytemppcaobj$values
daytemppcaobj$scores
daytemppcaobj$meanfd
plot(daytemppcaobj$meanfd)
plot(daytemppcaobj$harmonics)


## varimax rotation
daytemppcaobjVM <- varmx.pca.fd(daytemppcaobj, nharm=4)
str(daytemppcaobjVM)
dimnames(daytemppcaobjVM$scores)[[2]] <- paste("PC", 1:4, sep="")
round(daytemppcaobjVM$scores)

#  plot harmonics
op <- par(mfrow=c(2,2), pty="m")
plot.pca.fd(daytemppcaobjVM)
par(op)

#  plot log eigenvalues
daytempeigvals <- daytemppcaobjVM[[2]]

plot(1:20, log10(daytempeigvals[1:20]), type="b",
     xlab="Eigenvalue Number", ylab="Log 10 Eigenvalue")
abline(lsfit(5:20, log10(daytempeigvals[5:20])), lty=2)

#  plot factor scores
harmscr <- daytemppcaobj[[3]]

plot(harmscr[,1], harmscr[,2], xlab="PC 1", ylab="PC 2", pch=".")
text(harmscr[,1], harmscr[,2], CanadianWeather$place)

# plot(harmscr[,3], harmscr[,4], xlab="Harmonic 3", ylab="Harmonic 4")
# text(harmscr[,3], harmscr[,4], CanadianWeather$place)








# install.packages("fda")
# install.packages("fpca")
# install.packages("fdapace")
library(fda)
library(fpca)
library(fdapace)
weather<-CanadianWeather
daily<-daily
monthtemp<-as.data.frame(weather[["monthlyTemp"]])
montemp<-as.matrix(monthtemp)
ID.m<-rep(c(1:35),each=12)
IDa<-rep(c(1:20),each=12)
IDc<-rep(c(1:12),each=12)
IDAr<-rep(c(1:3),each=12)
mon<-rep(c(1:12),times=35)
mon1<-rep(c(1:12),times=20)   ##해양성
mon2<-rep(c(1:12),times=12)   ##대륙성
mon3<-rep(c(1:12),times=3)    ##극지방
Atl<-montemp[,1:15]
conti1<-montemp[,16:24]
conti2<-montemp[,30:32]
conti<-cbind(conti1,conti2)
pac<-montemp[,25:29]
arc<-montemp[,33:35]
marit<-cbind(Atl,pac)
#######################
for (i in 1:5){
  plot(montemp[,i],type='l',ylim=c(-20,15),ylab="Temperature")
  par(new=T)
}
#######################
input<- MakeFPCAInputs(IDs=ID.m,tVec=mon,yVec=montemp)
input1<-MakeFPCAInputs(IDs=IDa,tVec=mon1,yVec=marit)
input2<-MakeFPCAInputs(IDs=IDc,tVec=mon2,yVec=conti)
input3<-MakeFPCAInputs(IDs=IDAr,tVec=mon3,yVec=arc)
fpca<-FPCA(input$Ly,input$Lt,optns=list(FVEthreshold=0.999))
fpca1<-FPCA(input1$Ly,input1$Lt,optns=list(FVEthreshold=0.9999))
fpca2<-FPCA(input2$Ly,input2$Lt,optns=list(FVEthreshold=0.9999))
fpca3<-FPCA(input3$Ly,input3$Lt,optns=list(FVEthreshold=0.9999))
FPCs<-fpca[["xiEst"]]
FPC1<-FPCs[,1]
FPC2<-FPCs[,2]
FPC3<-FPCs[,3]
FPC4<-FPCs[,4]
plot(FPC1,FPC2,xlim=c(-60,60),ylim=c(-10,20))
text(FPC1,FPC2, CanadianWeather$place, col=4)
# identity(FPC1)
##################################
e<-fpca[["phi"]]
e1<-e[,1]
e2<-e[,2]
e3<-e[,3]
e4<-e[,4]
par(mfrow=c(2,2))
plot(e1,type='l',ylim=c(-0.6,0.6))
plot(e2,type='l',ylim=c(-0.6,0.6))
plot(e3,type='l',ylim=c(-0.6,0.6))
plot(e4,type='l',ylim=c(-0.6,0.6))
mtemp<-as.matrix(weather[["monthlyTemp"]])
