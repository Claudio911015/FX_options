###########################################
### M\'etodos de interpolaci\'on de vol ###
###########################################


# Vanna-Volga method -------------------------------------------------------

VV<-function(k,S0,t,d_d,d_f,v_atm,v_rr,v_vwb){
  
  F<-S0*d_f/d_d
  
  k_atm<-F*exp(v_atm^2*t/2)
  
  v_25c<-v_atm+v_vwb+v_rr/2
  
  v_25p<-v_atm+v_vwb-v_rr/2
  
  k_25c<-getstrikefromdelta(F,D=.25,v_25c,d_f,t)
  
  k_25p<-getstrikefromdelta(F,D=.25,v_25p,d_f,t)
  
  x1<-log(k_atm/k)*log(k_25c/k)/(log(k_atm/k_25p)*log(k_25c/k_25p))
  
  x2<-log(k/k_25p)*log(k_25c/k)/(log(k_atm/k_25p)*log(k_25c/k_atm))
  
  x3<-log(k/k_25p)*log(k/k_atm)/(log(k_25c/k_25p)*log(k_25c/k_atm))
  
  return(v_25p*x1+v_atm*x2+v_25c*x3)
}

getstrikefromdelta<-function(F,D=.25,v,d_f,t){
  
  return(F*exp(-qnorm(D/d_f)*v*sqrt(t)+v^2*t/2))
  
}

S0<-1.205
tenor<-94/365
ddf<-.9902752
fdf<-.9945049
ATM<-.0905
RR<--.0050
VWB<-.0013


k<-seq(.2,2.5,length.out = 50)

vol<-sapply(k,VV,S0=S0,t=tenor,d_d=ddf,d_f=fdf,v_atm=ATM, v_rr=RR,v_vwb=VWB)

plot(k,vol)
