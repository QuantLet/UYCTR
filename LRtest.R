#------------------------------------------------------------------------
# Testing the significant function of yield curve factors in forecasting the exchange rates
#------------------------------------------------------------------------

rm(list = ls(all = TRUE))
graphics.off()
setwd("")

# Read the data
library(readxl)
e      = read_excel("../TBR.xlsx")
e.1    = as.matrix(e[1:dim(e)[1],  2:dim(e)[2]])
ee     = apply(e.1, c(1,2), as.numeric)
b0     = 1
b1     = length(ee[, 1])
y      = matrix(NA, nrow = b1, ncol = 10)
x      = matrix(NA, nrow = 10, ncol = 3)
beta1  = rep(NA, b1)
beta2  = rep(NA, b1)
beta3  = rep(NA, b1)
re     = matrix(NA, nrow = b1, ncol = 10)
tao0   = 1
lamda  = 0.0609

# NS fitting
R1 = c(1, (1-exp(-lamda*tao0))/(lamda*tao0), ((1-exp(-lamda*tao0))/(lamda*tao0)-exp(-lamda*tao0)))
for(tao in c(3, 6, 12, 24, 36, 60, 84, 120, 240)){
R  = c(1, (1-exp(-lamda*tao))/(lamda*tao), ((1-exp(-lamda*tao))/(lamda*tao)-exp(-lamda*tao)))
R1 = rbind(R1,R)}

for(t in b0:b1){  
y[t, ]   = ee[t, 1:10]
x        = R1
fit      = lm(y[t, ] ~ 0+ x)
beta1[t] = fit$coefficients[1]
beta2[t] = fit$coefficients[2]
beta3[t] = fit$coefficients[3]
re[t, ]  = y[t, ]-x%*%c(beta1[t], beta2[t], beta3[t])
}
Rsqure1  = rep(NA, 10)
for(j in c(1:10)){
Rsqure1[j] = 1 - (t(re[, j])%*%re[, j]/((y[, j]-mean(y[, j]))%*%(y[, j]-mean(y[, j]))))
}
Rsqure1

## The unrestricted model
xM     = ee[, 11:13]
beltaM = solve(t(xM)%*%xM)%*%(t(xM)%*%re)
ro0    = xM%*%beltaM
Ro0    = (re[1, ] - ro0[1, ])%*%t(re[1, ] - ro0[1, ])
for(j in c(2:b1)){
Ro0    = Ro0 + (re[j, ] - ro0[j, ])%*%t(re[j, ] - ro0[j, ])
}
urv    = abs(det(Ro0))

# The in sample R^2
Rsqure2 = 1 - ((t(re[, 6]-ro0[, 6])%*%(re[, 6]-ro0[, 6]))/((re[, 6]-mean(re[, 6]))%*%(re[, 6]-mean(re[, 6]))))
Rsqure2

# The restricted model
Ro1   = (re[1, ])%*%t(re[1, ])
for(j in c(2:b1)){
  Ro1 = Ro1 + (re[j, ])%*%t(re[j, ])
}

rv    = abs(det(Ro1))


## The likelihood ratio test
q     = (b1-b0)*(log(rv)-log(urv))
pchisq(q, df=3, lower.tail = F)
