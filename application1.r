rm(list=ls(all=TRUE))

#Packages
library(readxl)
library(gamlss.cens)
library(xtable)
library(knitr)
library(dplyr)
library(ggplot2)
library(survival)
library(survminer)
library(xtable)

#Gamlss model

source("~GOLLGR_gamlss.R")

#Dataset

dados <- read_excel("~dados.xlsx")

#Descriptive analysis

#Percentage of censorship
(100*summary(as.factor(dados$censura))/sum(summary(as.factor(dados$censura))))

km1 <- survfit(Surv(tempos, censura) ~diabete, data = dados)
ggsurvplot(km1, data = dados, 
           pval = T,   
           legend.title="",
           palette =c(1, 'darkblue'),
           ylab='Survival probability',
           xlab='Time (years)',
           risk.table = T,
           ggtheme =theme_light())

hist(dados$idade, xlab='Age (years)',30,main = '')

#Getting initial values

fit.gr <- gamlss(tempos~1,family = GOLLGR(mu.link = "own",sigma.link = "log"),n.cyc=200,c.crit=0.01,
                 nu.start = 1,nu.fix = T,sigma.start = sd(dados$tempos)/1000,
                 tau.start = 1,tau.fix = T, data=dados)

fit1 <- gamlss(Surv(tempos,as.numeric(censura))~1,sigma.formula =~ 1,family=cens(GOLLGR)(sigma.link="log", mu.link="own"),data=dados,nu.start = 1,nu.fix = T, tau.start = 1,tau.fix = T,
               n.cyc=200, sigma.start = fit.gr$sigma.fv, mu.start = fit.gr$mu.fv)


# Fit regression

# GR model
fit1r <- gamlss(Surv(tempos,(censura))~  idade+ as.factor(dados$diabete),sigma.formula =~ idade+ as.factor(dados$diabete),family=cens(GOLLGR)(sigma.link="log", mu.link="own"),data=dados,nu.start = 1,nu.fix = T, tau.start = 1,tau.fix = T,c.crit=0.01,
                n.cyc=200, sigma.start = fit1$sigma.fv, mu.start = fit1$mu.fv)


# EGR model
fit2r <- gamlss(Surv(tempos,(censura))~ idade+ as.factor(dados$diabete),sigma.formula =~  idade+ as.factor(dados$diabete),family=cens(GOLLGR)(sigma.link="log", mu.link="own"),data=dados,
                tau.start = 1,tau.fix = T,c.crit=0.01,
                n.cyc=200, sigma.start = fit1r$sigma.fv, mu.start = fit1r$mu.fv)


#OLLGR model
fit3r <- gamlss(Surv(tempos,censura)~  idade+ as.factor(dados$diabete),sigma.formula =~  idade+ as.factor(dados$diabete),family=cens(GOLLGR)(sigma.link="log", mu.link="own"),data=dados,
                c.crit=0.01,nu.start = 1,nu.fix = T,
                n.cyc=200, sigma.start = fit1r$sigma.fv, mu.start = fit1r$mu.fv)

#GOLLGR model
fit4r <- gamlss(Surv(tempos,censura)~ idade+ as.factor(dados$diabete),sigma.formula =~  idade+ as.factor(dados$diabete),family=cens(GOLLGR)(sigma.link="log", mu.link="own"),data=dados,
                c.crit=0.01,
                n.cyc=2000, sigma.start = fit1r$sigma.fv, mu.start = fit1r$mu.fv)


##Tables

aic <- c(AIC(fit4r),AIC(fit3r),AIC(fit2r),AIC(fit1r))
bic <- c(BIC(fit4r),BIC(fit3r),BIC(fit2r),BIC(fit1r))
n <- length(dados$tempos)
CAIC <- c((-2*logLik(fit4r))+(4*(log(n)+1)),
          (-2*logLik(fit3r))+(5*(log(n)+1)),
          (-2*logLik(fit2r))+(5*(log(n)+1)),
          (-2*logLik(fit1r))+(6*(log(n)+1)))
t2 <- data.frame(aic,bic,CAIC)
t2


# LR test
LR.test(fit1r,fit4r)
LR.test(fit2r,fit4r)
LR.test(fit3r,fit4r)

# MLEs
summary(fit4r,type='qr')

# Residuals

par(mfrow=c(1,2)) 

plot(fit4r$residuals , ylim =c(-4 ,4) , ylab =" Quantile residuals ",pch=20)
abline(h=c(-3,3),xlab='Index',lty=2)

qqnorm(fit4r$residuals,xlim=c(-3.5,3.5),ylim=c(-4,4),main="(b)",pch=20)
abline(a=0,b=1,lwd=2,col='red')

