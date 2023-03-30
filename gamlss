require(gamlss)

#Link function (log shifted to -1)

own.linkfun = function(mu) log(mu + 1 - 1e-05)
own.linkinv = function(eta) -1 + pmax(.Machine$double.eps, exp(eta))
own.mu.eta = function(eta) pmax(.Machine$double.eps, exp(eta))
own.valideta = function(eta) TRUE   
make.link.gamlss("own")

#GOLLGR


GOLLGR <- function (mu.link = "own", sigma.link="log", tau.link = "log",nu.link = "log")
{
  mstats <- checklink("mu.link", "generalized odd log logistic GR", substitute(mu.link), 
                      c("logit", "log", "identity", "own"))
  dstats <- checklink("sigma.link", "generalized odd log logistic GR", substitute(sigma.link), 
                      c("inverse", "log", "identity", "own"))
  vstats <- checklink("nu.link", "generalized odd log logistic GR", substitute(nu.link),    
                      c('inverse', "log", "identity", "own"))
  tstats <- checklink("tau.link", "generalized odd log logistic GR", substitute(tau.link),   
                      c("inverse", "log", "identity", "own")) 
  structure(
    list(family = c("GOLLGR", "generalized odd log logistic GR"),
         parameters = list(mu=TRUE, sigma=TRUE, nu=TRUE, tau=TRUE), 
         nopar = 4, 
         type = "Continuous",
         
         mu.link = as.character(substitute(mu.link)),  
         sigma.link = as.character(substitute(sigma.link)), 
         nu.link = as.character(substitute(nu.link)), 
         tau.link = as.character(substitute(tau.link)), 
         
         mu.linkfun = mstats$linkfun, 
         sigma.linkfun = dstats$linkfun, 
         nu.linkfun = vstats$linkfun,
         tau.linkfun = tstats$linkfun,  
         
         mu.linkinv = mstats$linkinv, 
         sigma.linkinv = dstats$linkinv,
         nu.linkinv = vstats$linkinv,
         tau.linkinv = tstats$linkinv, 
         
         mu.dr = mstats$mu.eta, 
         sigma.dr = dstats$mu.eta, 
         nu.dr = vstats$mu.eta,
         tau.dr = tstats$mu.eta, 
         
         dldm = function(y,mu,sigma,nu,tau){ #----------------------------------------------------- ok
           nd1 = gamlss:::numeric.deriv(dGOLLGR(y, mu, sigma, nu,tau,log = T),"mu", delta = 1e-04)
           dldm = as.vector(attr(nd1, "gradient")) 
           dldm
         },
         d2ldm2 = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok
           nd1 = gamlss:::numeric.deriv(dGOLLGR(y, mu, sigma, nu,tau,log = T),"mu", delta = 1e-04)
           dldm = as.vector(attr(nd1, "gradient")) 
           d2ldm2 = -dldm * dldm
         },     
         dldd = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok  
           nd1 = gamlss:::numeric.deriv(dGOLLGR(y, mu, sigma, nu,tau,log =T),"sigma", delta = 1e-04)
           dldd = as.vector(attr(nd1, "gradient"))
           dldd
         } ,
         d2ldd2 = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok
           nd1 = gamlss:::numeric.deriv(dGOLLGR(y, mu, sigma, nu,tau,log =T),"sigma", delta = 1e-04)
           dldd = as.vector(attr(nd1, "gradient"))
           d2ldd2 = -dldd*dldd
           d2ldd2 
         },   
         dldv = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok
           nd1 = gamlss:::numeric.deriv(dGOLLGR(y, mu, sigma, nu,tau,log =T),"nu", delta = 1e-04)
           dldv = as.vector(attr(nd1, "gradient"))
           dldv 
         },
         d2ldv2 = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok 
           nd1 = gamlss:::numeric.deriv(dGOLLGR(y, mu, sigma, nu,tau,log =T),"nu", delta = 1e-04)
           dldv = as.vector(attr(nd1, "gradient"))
           d2ldv2 = -dldv * dldv
         },
         dldt = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok
           nd1 = gamlss:::numeric.deriv(dGOLLGR(y, mu, sigma, nu,tau,log =T),"tau", delta = 1e-04)
           dldt = as.vector(attr(nd1, "gradient"))           
           dldt
         } ,
         d2ldt2 = function(y,mu,sigma,nu,tau){ #----------------------------------------------------- ok
           nd1 = gamlss:::numeric.deriv(dGOLLGR(y, mu, sigma, nu,tau,log =T),"tau", delta = 1e-04)
           dldt = as.vector(attr(nd1, "gradient"))  
           d2ldt2 = -dldt * dldt
           d2ldt2
         },
         d2ldmdd = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok
           nd1 = gamlss:::numeric.deriv(dGOLLGR(y, mu, sigma, nu, tau,log= TRUE), "mu", delta = 1e-04)
           dldm = as.vector(attr(nd1, "gradient"))
           nd1 = gamlss:::numeric.deriv(dGOLLGR(y, mu, sigma, nu, tau,log=TRUE), "sigma", delta = 1e-04)
           dldd = as.vector(attr(nd1, "gradient"))           
           d2ldmdd = -dldm * dldd
           d2ldmdd               
         },
         d2ldmdv = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok
           nd1 = gamlss:::numeric.deriv(dGOLLGR(y, mu, sigma, nu,tau,log = TRUE), "mu", delta = 1e-04)
           dldm = as.vector(attr(nd1, "gradient"))
           nd1 = gamlss:::numeric.deriv(dGOLLGR(y, mu, sigma, nu,tau,log = TRUE), "nu", delta = 1e-04)
           
           dldv = as.vector(attr(nd1, "gradient"))
           d2ldmdv = -dldm * dldv
           d2ldmdv			
         },
         
         d2ldmdt = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok
           nd1 = gamlss:::numeric.deriv(dGOLLGR(y, mu, sigma, nu,tau,log = TRUE), "mu", delta = 1e-04)
           dldm = as.vector(attr(nd1, "gradient"))
           nd1 = gamlss:::numeric.deriv(dGOLLGR(y, mu, sigma, nu,tau,log = TRUE), "tau", delta = 1e-04)
           
           dldv = as.vector(attr(nd1, "gradient"))
           d2ldmdv = -dldm * dldv
           d2ldmdv         },
         
         d2ldddv = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok
           nd1 = gamlss:::numeric.deriv(dGOLLGR(y, mu, sigma, nu,tau,log = TRUE), "sigma", delta = 1e-04)
           dldm = as.vector(attr(nd1, "gradient"))
           nd1 = gamlss:::numeric.deriv(dGOLLGR(y, mu, sigma, nu,tau,log = TRUE), "nu", delta = 1e-04)
           
           dldv = as.vector(attr(nd1, "gradient"))
           d2ldmdv = -dldm * dldv
           d2ldmdv	
         },
         d2ldddt = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok
           nd1 = gamlss:::numeric.deriv(dGOLLGR(y, mu, sigma, nu,tau,log = TRUE), "sigma", delta = 1e-04)
           dldm = as.vector(attr(nd1, "gradient"))
           nd1 = gamlss:::numeric.deriv(dGOLLGR(y, mu, sigma, nu,tau,log = TRUE), "tau", delta = 1e-04)
           
           dldv = as.vector(attr(nd1, "gradient"))
           d2ldmdv = -dldm * dldv
           d2ldmdv
         },
         d2ldvdt = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok
           nd1 = gamlss:::numeric.deriv(dGOLLGR(y, mu, sigma, nu,tau,log = TRUE), "nu", delta = 1e-04)
           dldm = as.vector(attr(nd1, "gradient"))
           nd1 = gamlss:::numeric.deriv(dGOLLGR(y, mu, sigma, nu,tau,log = TRUE), "tau", delta = 1e-04)
           dldv = as.vector(attr(nd1, "gradient"))
           d2ldmdv = -dldm * dldv
           d2ldmdv 
         },
         #----------------------------------------------------- ok
         G.dev.incr  = function(y,mu,sigma,nu,tau,...) 
         { 
           -2*dGOLLGR(y,mu,sigma,nu,tau,log=TRUE)
         } ,                     
         rqres = expression(   
           rqres(pfun="pGOLLGR", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu, tau=tau)) ,
         
         mu.initial = expression(mu <- rep(median(y), length(y))),
         sigma.initial = expression(sigma <- rep(sd(y),length(y))), 
         nu.initial = expression( nu <- rep(2, length(y))), 
         tau.initial = expression(tau <-rep(2, length(y))), 
         
         mu.valid = function(mu) all(mu > -1)  , 
         sigma.valid = function(sigma)  all(sigma > 0),
         nu.valid = function(nu) all(nu > 0), 
         tau.valid = function(tau) all(tau > 0), 
         y.valid = function(y)  all(y > 0) 
    ),
    class = c("gamlss.family","family"))
}
#-----------------------------------------------------------------
# Probability Density Function
dGOLLGR <- function(x,mu=-0.5,sigma=0.0001,nu=0.5,tau=0.5,log = FALSE){
  if (any(mu <= -1))  stop(paste("mu must > -1", "\n", "")) 
  if (any(nu <= 0))  stop(paste("nu must be positive", "\n", "")) 
  if (any(tau <= 0))  stop(paste("tau must be positive", "\n", "")) 
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
  if (any(x < 0))  stop(paste("x must be positive", "\n", "")) 
  g <- ((2*(sigma^(mu+1)))*(x^((2*mu)+1))*exp(-sigma*(x^2)))/(gamma(mu+1))
  a <- mu+1
  b <- sigma*(x^2)
  G <- pgamma(q=b,shape=a,rate=1)
  f <- (nu*tau*g*(G^((nu*tau)-1))*((1-G^tau)^(nu-1)))/((G^(nu*tau)+((1-G^tau)^(nu)))^2)
  if(log==FALSE) fx  <- f else fx <- log(f) 
  fx
}  

#-----------------------------------------------------------------  
# Cumulative Density Function
pGOLLGR <- function(q,mu=-0.5,sigma=0.0001,nu=0.5,tau=0.5,lower.tail = TRUE, log.p = FALSE){
  if (any(nu <= 0))  stop(paste("nu must be positive", "\n", "")) 
  if (any(tau <= 0))  stop(paste("tau must be positive", "\n", "")) 
  if (any(mu <= -1))  stop(paste("mu must > -1", "\n", "")) 
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
  if (any(q < 0))  stop(paste("q must be positive", "\n", ""))
  a <- mu+1
  b <- sigma*(q^2)
  G <-pgamma(q=b,shape=a,rate=1)
  cdf <- G^(nu*tau)/(G^(nu*tau)+(1-G^tau)^nu)
  if(lower.tail==TRUE) cdf  <- cdf else  cdf <- 1-cdf 
  if(log.p==FALSE) cdf  <- cdf else  cdf <- log(cdf) 
  cdf
}
#-----------------------------------------------------------------  
# Quantile Function
qGOLLGR <-  function(p,mu=-0.5,sigma=0.0001,nu=0.5,tau=0.5, lower.tail = TRUE, log.p = FALSE){   
  if (any(mu <= -1))  stop(paste("mu must be > -1", "\n", "")) 
  if (any(sigma <= 0))  stop(paste("sigma must  be positive", "\n", "")) 
  if (any(nu < 0))  stop(paste("nu must be positive", "\n", ""))
  if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))
  if (any(p < 0) | any(p > 1))  stop(paste("p must be between 0 and 1", "\n", ""))    
  if (log.p==TRUE) p <- exp(p) else p <- p
  if (lower.tail==TRUE) p <- p else p <- 1-p
  e1 <- 1/nu
  e2 <- 1/tau
  p1 <- (p^(e1*e2))/((((1-p)^e1)+p^e1)^e2)
  u1 <- NULL
  for(i in 1:length(p1)){
    if(p1[i]==0){u1[i]=0.00001}
    if(p1[i]==1){u1[i]=0.99999}
    else {u1[i]=p1[i]}
  }
  q <- sqrt((qgamma(p=u1,shape=mu+1,rate=1))/sigma)
  q
}
#----------------------------------------------------------------- 
# Random generating function
rGOLLGR <- function(n, mu=-0.5,sigma=0.0001,nu=0.5,tau=0.5){
  if (any(mu <= -1))  stop(paste("mu must be > -1", "\n", "")) 
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", ""))   
  if (any(nu < 0))  stop(paste("nu must be positive", "\n", ""))  
  if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))  
  if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))  
  n <- ceiling(n)
  u <- runif(n,0,1)
  r <- qGOLLGR(u, mu =mu, sigma =sigma, nu=nu, tau=tau)
  r
}
#----------------------------------------------------------------- 
#hazard function

hGOLLGR <- function(x, mu=0.5, sigma=0.01,nu=2,tau=2,lower.tail = TRUE, log.p = FALSE)
{  
  h <- dGOLLGR(x,mu,sigma,nu,tau)/(1-pGOLLGR(x,mu,sigma,nu,tau))
  h
}
#----------------------------------------------------------------- 
#survival

sGOLLGR <- function(x, mu=0.5, sigma=0.01, nu=2,tau=2,lower.tail = TRUE, log.p = FALSE)
{  
  S <- (1-pGOLLGR(x,mu,sigma,nu,tau))
  S
}
