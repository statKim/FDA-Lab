### SPD case test

gr <- seq(0, 1, length.out = 11)
mu <- lapply(gr, function(t){
    matrix(c(t^(0.4), 0.5*t, 0.1*t^(1.5), 
             0.5*t, t^(0.5), 0.5*t, 
             0.1*t^(1.5), 0.5*t, t^(0.6)), 
           nrow = 3, byrow = T)
})

i <- 2
plot3d(ellipse3d(mu[[i]]), col = "red", alpha = 0.5, aspect = TRUE)




### Ellipsoid mean function
clear3d()   # remove graph
mfrow3d(1, 10)   # par(mfrow = c(2, 1))
for (i in 2:length(gr)) {
  plot3d(ellipse3d(mu[[i]]), 
         col = "red",
         xlab = '', ylab = '', zlab = '', axes = FALSE,
         xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1))
}
rgl.close()



### DTI dataset
library(refund)
data(DTI)
class(DTI)
head(DTI)
dim(DTI)
str(DTI)

length(unique(DTI$ID))

dim(DTI$cca)  # 382  93
dim(DTI$rcst) # 382  55

DTI$cca %>% is.na() %>% sum()
DTI$rcst %>% is.na() %>% sum()


head(DTI[, 1:4])

DTI$Nscans


str(DTI2)

plot3d(ellipse3d(matrix(DTI$cca[1, ], ncol = 1)), 
       col = "red",
       xlab = '', ylab = '', zlab = '', axes = FALSE,
       xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1))
