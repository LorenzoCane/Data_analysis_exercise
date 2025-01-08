#Data Analysis course A.A. 2024/2025 - Part 1 - Exam Exercise
#Lorenzo Cane


#overall consideration: in all the model considered the Ljung-Box test p-value decreases (but still p-value > 0.05) if lags = 15 or 20 is selected
#also the acf plot of residuals shows peak betweeen 10 and 20 lags. 
#A correlation in medium-lags range can be hard to be understand but it seems to be present
#*************************************************************************
#Install and Load libraries
if (!requireNamespace("forecast", quietly = TRUE)) {
  install.packages("forecast", dependencies = TRUE)
}
if (!requireNamespace("urca", quietly = TRUE)) {
  install.packages("urca", dependencies = TRUE)
}
library(forecast)
library(urca)

#*************************************************************************
#Import data
input_file <- "Sun_data.txt" #data input file

data <- read.table(input_file, header = TRUE,  sep = '\t',  fill = TRUE)
flux <- data$flux_1GHz
flux <-flux[flux > 0] #Remove non-positive values 
flux <- log(flux)
min(flux)
max(flux)
#*************************************************************************
#Preliminary study

#plot
plot.ts(flux, main ="Sun Flux at 1GHz", ylab= ("Flux [SFU]"))
pdf("img/Flux_plot.pdf")
plot.ts(flux, main ="Sun Flux at 1GHz", ylab= ("Flux [SFU]"))
dev.off()

#acf and pacf of the flux
par(mfrow=c(2,1))
pdf("img/acf_pacf_flux.pdf", width = 15)
acf(flux, lag.max = 150, main="ACF of original Flux")
pacf(flux, lag.max = 150, main = "PACF of original Flux")
dev.off()
par(mfrow=c(2,1))

#At this stage (no differenciations) series seems not stationary

#KPSS test (confermed)
kpss_test_or <- ur.kpss(flux, type="tau", lags = "short")
print("Before differencing:")
summary(kpss_test_or)

#***************************************************************************
#Try 1 diff
flux_diff <- diff(flux, differences = 1)

#ACF/PACF plots for diff. flux
par(mfrow=c(2,1))
pdf("img/acf_pacf_diff_flux.pdf", width = 15)
acf(flux_diff, lag.max = 150, main = "ACF Flux diff")
pacf(flux_diff, lag.max = 150, main ="PACF Flux diff")
dev.off()
par(mfrow=c(2,1))

#Now ts seems stationary ->test (conferm stat)
kpss_test_diff <- ur.kpss(flux_diff, type = "tau", lags = "short")
print("After differencing:")
summary(kpss_test_diff)

#***************************************************************************
#Fit ARIMA(2,1,3) model 
print("****************************************************")
print("MODEL 1: ARIMA(2,1,3)")
#previous plots suggest for an ARMA(2,4), ARMA(2,3) or ARMA(2,5) or ARMA(2,2) after differencing


#model <- auto.arima(flux)
#summary(model)

model1 = arima(flux, order = c(2,1,3))
model1$coef

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#residual analysis for model 1
res1 <- residuals(model1)
res1 <- res1[2:length(res1)] #avoid strange (possible) behaviour, differencing the 1st res may be artificial

plot.ts(res1, main = "Residuals plot for ARIMA(2,1,3)" , ylab = "Residuals")
pdf("img/res_arima213.pdf")
plot.ts(res1, main = "Residuals plot for ARIMA(2,1,3)" , ylab = "Residuals")
dev.off()

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#acf and pacf of residuals
par(mfrow=c(2,1))
pdf("img/acf_pacf_arima213.pdf")
acf(res1, main = "ACF of Residuals for ARIMA(2,1,3)")
pacf(res1, main = "PACF of Residuals for ARIMA(2,1,3)")
dev.off()
par(mfrow=c(2,1))

#Graphically uncorrelated (now check)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#Test

#Indipendence test (Box test)
Box.test(res1, lag = 10, type = "Ljung-Box", fitdf = 5) # seem uncorrelated (p-value = 0.72) at short lags level

#Normality plot + test (Shapiro-Wilk and Kolmogorov-Smirnov)
pdf("img/qqplot_res_arima213.pdf")
qqnorm(res1, main = "Normal Q-Q Plot (ARIMA(2,1,3))", xlab = "Th. Quantiles", ylab = "Sample Quantiles")
qqline(res1, distribution = qnorm, col = "red")
dev.off()

ks.test(res1, y = pnorm)
shapiro.test(res1)
#not normally distr. for both test

#***************************************************************************
#Fit ARIMA(2,1,5) model 
print("****************************************************")
print("MODEL 2: ARIMA(2,1,5)")

model2 = arima(flux, order = c(2,1,5))
model2$coef

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#residual analysis for model 2
res2 <- residuals(model2)
res2 <- res2[2:length(res2)] #avoid strange (possible) behaviour, differencing the 1st res may be artificial

plot.ts(res2, main = "Residuals plot for ARIMA(2,1,5)" , ylab = "Residuals")
pdf("img/res_arima215.pdf")
plot.ts(res2, main = "Residuals plot for ARIMA(2,1,5)" , ylab = "Residuals")
dev.off()

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#acf and pacf of residuals
par(mfrow=c(2,1))
pdf("img/acf_pacf_arima215.pdf")
acf(res2, main = "ACF of Residuals for ARIMA(2,1,5)")
pacf(res2, main = "PACF of Residuals for ARIMA(2,1,5)")
dev.off()
par(mfrow=c(2,1))

#Graphically uncorrelated (now check)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#Test

#Indipendence test (Box test)
Box.test(res2, lag = 10, type = "Ljung-Box", fitdf = 7) # seem uncorrelated (p-value = 0.59) at short lags level

#Normality plot + test (Shapiro-Wilk and Kolmogorov-Smirnov)
pdf("img/qqplot_res_arima215.pdf")
qqnorm(res2, main = "Normal Q-Q Plot (ARIMA(2,1,5))", xlab = "Th. Quantiles", ylab = "Sample Quantiles")
qqline(res2, distribution = qnorm, col = "red")
dev.off()

ks.test(res2, y = pnorm)
shapiro.test(res2)
#not normally distr. for both test

#***************************************************************************
#Fit ARIMA(2,1,4) model 
print("****************************************************")
print("MODEL 3: ARIMA(2,1,4)")

model3 = arima(flux, order = c(2,1,4))
model3$coef

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#residual analysis for model 3
res3 <- residuals(model3)
res3 <- res3[2:length(res3)] #avoid strange (possible) behaviour, differencing the 1st res may be artificial

plot.ts(res3, main = "Residuals plot for ARIMA(2,1,4)" , ylab = "Residuals")
pdf("img/res_arima214.pdf")
plot.ts(res3, main = "Residuals plot for ARIMA(2,1,4)" , ylab = "Residuals")
dev.off()

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#acf and pacf of residuals
par(mfrow=c(2,1))
pdf("img/acf_pacf_arima214.pdf")
acf(res3, main = "ACF of Residuals for ARIMA(2,1,4)")
pacf(res3, main = "PACF of Residuals for ARIMA(2,1,4)")
dev.off()
par(mfrow=c(2,1))

#Graphically uncorrelated (now check)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#Test

#Indipendence test (Box test)
Box.test(res3, lag = 10, type = "Ljung-Box", fitdf = 6) # seem uncorrelated (p-value = 0.59) at short lags level

#Normality plot + test (Shapiro-Wilk and Kolmogorov-Smirnov)
pdf("img/qqplot_res_arima214.pdf")
qqnorm(res3, main = "Normal Q-Q Plot (ARIMA(2,1,4))", xlab = "Th. Quantiles", ylab = "Sample Quantiles")
qqline(res3, distribution = qnorm, col = "red")
dev.off()

ks.test(res3, y = pnorm)
shapiro.test(res3)
#not normally distr. for both test

#***************************************************************************
#Fit ARIMA(2,1,2) model 
print("****************************************************")
print("MODEL 4: ARIMA(2,1,2)")

model4 = arima(flux, order = c(2,1,2))
model4$coef

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#residual analysis for model 4
res4 <- residuals(model4)
res4 <- res4[2:length(res4)] #avoid strange (possible) behaviour, differencing the 1st res may be artificial

plot.ts(res4, main = "Residuals plot for ARIMA(2,1,2)" , ylab = "Residuals")
pdf("img/res_arima212.pdf")
plot.ts(res4, main = "Residuals plot for ARIMA(2,1,2)" , ylab = "Residuals")
dev.off()

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#acf and pacf of residuals
par(mfrow=c(2,1))
pdf("img/acf_pacf_arima212.pdf")
acf(res4, main = "ACF of Residuals for ARIMA(2,1,2)")
pacf(res4, main = "PACF of Residuals for ARIMA(2,1,2)")
dev.off()
par(mfrow=c(2,1))

#Graphically uncorrelated (now check)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#Test

#Indipendence test (Box test)
Box.test(res4, lag = 10, type = "Ljung-Box", fitdf = 4) # seem uncorrelated (p-value = 0.87) at short lags level

#Normality plot + test (Shapiro-Wilk and Kolmogorov-Smirnov)
pdf("img/qqplot_res_arima212.pdf")
qqnorm(res4, main = "Normal Q-Q Plot (ARIMA(2,1,2))", xlab = "Th. Quantiles", ylab = "Sample Quantiles")
qqline(res4, distribution = qnorm, col = "red")
dev.off()

ks.test(res4, y = pnorm)
shapiro.test(res4)
#not normally distr. for both test#

#***************************************************************************
#Model Validation through comparison
#Based on res analysis ARIMA(2,1,2) seems to be the best one (based on L-B test results), residuals plots are all very similar

model5 = arima(flux, order=c(1,1,1))

AIC(model1,model2,model3,model4,model5)