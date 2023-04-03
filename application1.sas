data lista3;
  input t censur ;
 cards;
2	1
10	1
13	1
23	1
23	1
28	1
30	1
65	1
80	1
88	1
106	1
143	1
147	1
173	1
181	1
212	1
245	1
247	1
261	1
266	1
275	1
293	1
300	1
300	1
300	1
300	1
300	1
300	1
300	1
300	1
;
  proc print data=lista3;
 run;

//** GOLLGR distribution **//

proc nlmixed cov data=lista3;

parms alpha=0.001  beta=0.001  delta=29.001 theta=0.0001;

bounds alpha>0;
bounds beta>0;
bounds theta>0;
bounds delta>-1;

GG1=probgam(theta*(t**2),delta+1);
gg=2*(theta**(delta+1))*(t**(2*delta+1))*exp(-theta*(t**(2)))*((GAMMA(delta+1))**(-1));
ff=alpha*beta*(gg)*((GG1)**(alpha*beta-1))*((1-(GG1)**(beta))**(alpha-1))*(((GG1)**(alpha*beta))+(1-(GG1)**(beta))**(alpha)  )**(-2);

logp=log(ff);

model t ~ general(logp);
run;


  //** OLLGR distribution **//

proc nlmixed cov data=lista3;

parms  alpha=0.002  delta=37.01  theta=0.001;
bounds theta>0;
bounds delta>-1;
bounda alpha>0;
beta=1;

GG1=probgam(theta*(t**2),delta+1);
gg=2*(theta**(delta+1))*(t**(2*delta+1))*exp(-theta*(t**(2)))*((GAMMA(delta+1))**(-1));
ff=alpha*beta*(gg)*((GG1)**(alpha*beta-1))*((1-(GG1)**(beta))**(alpha-1))*(((GG1)**(alpha*beta))+(1-(GG1)**(beta))**(alpha)  )**(-2);

logp=log(ff);

model t ~ general(logp);
run;



  //** EGR distribution  **//

proc nlmixed cov data=lista3;

parms  beta=0.001  delta=3.01  theta=0.0001;
bounds theta>0;
bounds delta>-1;
bounda beta>0;
alpha=1;

GG1=probgam(theta*(t**2),delta+1);
gg=2*(theta**(delta+1))*(t**(2*delta+1))*exp(-theta*(t**(2)))*((GAMMA(delta+1))**(-1));
ff=alpha*beta*(gg)*((GG1)**(alpha*beta-1))*((1-(GG1)**(beta))**(alpha-1))*(((GG1)**(alpha*beta))+(1-(GG1)**(beta))**(alpha)  )**(-2);

logp=log(ff);

model t ~ general(logp);
run;



//** RG distribution **//

proc nlmixed cov data=lista3;

parms   delta=-0.5 theta=0.00001;
bounds theta>0;
bounds delta>-1;
alpha=1;
beta=1;
GG1=probgam(theta*(t**2),delta+1);
gg=2*(theta**(delta+1))*(t**(2*delta+1))*exp(-theta*(t**(2)))*((GAMMA(delta+1))**(-1));
ff=alpha*beta*(gg)*((GG1)**(alpha*beta-1))*((1-(GG1)**(beta))**(alpha-1))*(((GG1)**(alpha*beta))+(1-(GG1)**(beta))**(alpha)  )**(-2);


logp=log(gg);


model t ~ general(logp);
run;



  //**  Classical Rayleigh distribution **//

proc nlmixed cov data=lista3;

parms    lambda=50;
bounds lambda>0;
delta=0;

gg=2*((1/(lambda**2))**(delta+1))*(t**(2*delta+1))*exp(-(1/(lambda**2))*(t**(2)))*((GAMMA(delta+1))**(-1));


logp=log(gg);


model t ~ general(logp);
run;


  //** Maxwell  distribution **//

proc nlmixed cov data=lista3;

parms    lambda=100;
bounds lambda>0;
delta=1/2;

gg=2*((1/(2*(lambda)**2))**(delta+1))*(t**(2*delta+1))*exp(-(1/(2*(lambda)**2))*(t**(2)))*((GAMMA(delta+1))**(-1));




logp=log(gg);


model t ~ general(logp);
run;


  //** Chi-square distribution**//

proc nlmixed cov data=lista3;

parms   tau=200 n=1.00001;
bounds tau>0;
bounds n>0;

gg=2*((1/(2*(tau)**2))**(((n/2)-1)+1))*(t**(2*((n/2)-1)+1))*exp(-(1/(2*(tau)**2))*(t**(2)))*((GAMMA(((n/2)-1)+1))**(-1));

logp=log(gg);


model t ~ general(logp);
run;


  //** Half normal distribution **//

proc nlmixed cov data=lista3;

parms   sigma=100;
bounds sigma>0;
delta=-1/2;

gg=2*((2/(sigma**2))**(delta+1))*(t**(2*delta+1))*exp(-(2/(sigma**2))*(t**(2)))*((GAMMA(delta+1))**(-1));

logp=log(gg);


model t ~ general(logp);
run;

