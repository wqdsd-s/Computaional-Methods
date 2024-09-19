#Lecture 1
#1.(a)
p<-0.64
m<-2**31-1
a<-7**5
b<-0

#generate random integer sequence
random_seq<-function(len){
  m<-2**31-1
  a<-7**5
  b<-0
  seq<-vector(l=len)
  seq[1]<-1
  for(t in 2:len) {
    seq[t]<-(a*seq[t-1])%%m
  }
  return(seq)
}
ranseq<-random_seq(44000)

#generate uniform sequence
uni<-ranseq/m

#generate bernoulli sequence
bern<-vector(l=44000)
for(i in 1:44000) {
  if(uni[i]<p) {bern[i]=1}
  else {bern[i]=0}
}
#generate binomial sequence
n<-44
bino<-vector(l=1000)
for(i in 1:1000) {
  bino[i]<-sum(bern[(i*n-43):(i*n)])
}

#calculate for theoretical value of p(X>=40)
res1<-0
for(i in 40:44) {
  res1<-res1+dbinom(i,n,p)
}
print(res1)

#calculate for sample value of p(X>=40)
print(length(bino[bino>=40])/1000)

#1.(b)
rm(list = ls())

n<-10000
#generate uniform sequence
uni_seq<-function(len){
  m<-2**31-1
  a<-7**5
  b<-0
  seq<-vector(l=len)
  seq[1]<-1
  for(t in 2:len) {
    seq[t]<-(a*seq[t-1])%%m
  }
  return(seq/m)
}
uni<-uni_seq(n)
#generate exponential sequence
lambda<-1.5
exp_seq<--lambda*log(1-uni)

#estimate P(X>=1) and P(X>=4)
print(length(exp_seq[exp_seq>=1])/n)
print(length(exp_seq[exp_seq>=4])/n)

#sample mean and std
print(mean(exp_seq))
print(sd(exp_seq))

#histogram
hist(exp_seq,100)

#1.(c)
rm(list = ls())
n<-5000
#generate uniform sequence
uni_seq<-function(len){
  m<-2**31-1
  a<-7**5
  b<-0
  seq<-vector(l=len)
  seq[1]<-1
  for(t in 2:len) {
    seq[t]<-(a*seq[t-1])%%m
  }
  return(seq/m)
}
uni<-uni_seq(n)

#generate normal sequence
norm_seq<-vector(l=n)
for(i in 1:2500) {
  norm_seq[2*i-1]<-sqrt(-2*log(uni[2*i-1]))*cos(2*pi*uni[2*i])
  norm_seq[2*i]<-sqrt(-2*log(uni[2*i-1]))*sin(2*pi*uni[2*i])
}

#1.(d)
rm(list = ls())
n<-7000
#generate uniform sequence
uni_seq<-function(len){
  m<-2**31-1
  a<-7**5
  b<-0
  seq<-vector(l=len)
  seq[1]<-1
  for(t in 2:len) {
    seq[t]<-(a*seq[t-1])%%m
  }
  return(seq/m)
}
uni<-uni_seq(n)

#generate normal sequence
V<-2*uni-1
norm_seq<-vector(l=n)
for(i in 1:3000) {
  if ((V[2*i-1]**2+V[2*i]**2)<=1){
    norm_seq[2*i-1]=V[2*i-1]*sqrt(-2*log(V[2*i-1]**2+V[2*i]**2)/(V[2*i-1]**2+V[2*i]**2))
    norm_seq[2*i]=V[2*i]*sqrt(-2*log(V[2*i-1]**2+V[2*i]**2)/(V[2*i-1]**2+V[2*i]**2))
    }
  else{
    norm_seq[2*i-1]=NaN
    norm_seq[2*i]=NaN
    }
}
norm_seq<-norm_seq[!is.nan(norm_seq)][1:5000]
print(mean(norm_seq))
print(sd(norm_seq))

#1.(e)
rm(list = ls())

#time for B-M Method
timestart1<-Sys.time()
n<-5000
#generate uniform sequence
uni_seq<-function(len){
  m<-2**31-1
  a<-7**5
  b<-0
  seq<-vector(l=len)
  seq[1]<-1
  for(t in 2:len) {
    seq[t]<-(a*seq[t-1])%%m
  }
  return(seq/m)
}
uni<-uni_seq(n)
#generate normal sequence
norm_seq<-vector(l=n)
for(i in 1:2500) {
  norm_seq[2*i-1]<-sqrt(-2*log(uni[2*i-1]))*cos(2*pi*uni[2*i])
  norm_seq[2*i]<-sqrt(-2*log(uni[2*i-1]))*sin(2*pi*uni[2*i])
}
timeend1<-Sys.time()

rm(list = ls()[ls()!='timeend1' & ls()!='timestart1'])

#time for P-M Method
timestart2<-Sys.time()
n<-7000
#generate uniform sequence
uni_seq<-function(len){
  m<-2**31-1
  a<-7**5
  b<-0
  seq<-vector(l=len)
  seq[1]<-1
  for(t in 2:len) {
    seq[t]<-(a*seq[t-1])%%m
  }
  return(seq/m)
}
uni<-uni_seq(n)
#generate normal sequence
V<-2*uni-1
norm_seq<-vector(l=n)
for(i in 1:3000) {
  if ((V[2*i-1]**2+V[2*i]**2)<=1){
    norm_seq[2*i-1]=V[2*i-1]*sqrt(-2*log(V[2*i-1]**2+V[2*i]**2)/(V[2*i-1]**2+V[2*i]**2))
    norm_seq[2*i]=V[2*i]*sqrt(-2*log(V[2*i-1]**2+V[2*i]**2)/(V[2*i-1]**2+V[2*i]**2))
  }
  else{
    norm_seq[2*i-1]=NaN
    norm_seq[2*i]=NaN
  }
}
norm_seq<-norm_seq[!is.nan(norm_seq)][1:5000]
timeend2<-Sys.time()

print(timeend1-timestart1)
print(timeend2-timestart2)

#Lecture 2
#2.(a)
rm(list = ls())

#generate normal sequence with given length
norm_seq<-function(len){
  #generate uniform sequence
  uni_seq<-function(len){
    m<-2**31-1
    a<-7**5
    b<-0
    seq<-vector(l=len)
    seq[1]<-1
    for(t in 2:len) {
      seq[t]<-(a*seq[t-1])%%m
    }
    return(seq/m)
  }
  uni<-uni_seq(len)
  
  #generate normal sequence
  norm_seq<-vector(l=len)
  for(i in 1:(len/2)) {
    norm_seq[2*i-1]<-sqrt(-2*log(uni[2*i-1]))*cos(2*pi*uni[2*i])
    norm_seq[2*i]<-sqrt(-2*log(uni[2*i-1]))*sin(2*pi*uni[2*i])
  }
  return (norm_seq)
}

#generate wiener sequence at time=t and length of len
wie_seq<-function(len,t){
  return (norm_seq(len)*sqrt(t))
}

#estimate A(t)and B(t)
t<-1
len<-10000
Wt<-wie_seq(len,t)
At<-Wt**2+sin(Wt)
Bt<-exp(t/2)*cos(Wt)
print('t=1')
print(mean(At))
print(mean(Bt))

t<-3
Wt<-wie_seq(len,t)
At<-Wt**2+sin(Wt)
Bt<-exp(t/2)*cos(Wt)
print('t=3')
print(mean(At))
print(mean(Bt))

t<-5
Wt<-wie_seq(len,t)
At<-Wt**2+sin(Wt)
Bt<-exp(t/2)*cos(Wt)
print('t=5')
print(mean(At))
print(mean(Bt))


#3.(a)
rm(list = ls()[ls()!='wie_seq' & ls()!='norm_seq'])

r<-0.055
sigma<-0.2
S0<-100
T<-5
X<-100

len<-10000
WT<-wie_seq(len,T)
#simulate for ST
ST<-S0*exp(sigma*WT+(r-sigma**2/2)*T)
#estimate for call option price
temp<-vector(l=len)
for(i in 1:len) {
  temp[i]<-max(ST[i]-X,0)
  }
c<-mean(temp)*exp(-r*T)
print(c)

#3.(b)
d1<-(log(S0/X)+(r+sigma**2/2)*T)/sigma/sqrt(T)
d2<-d1-sigma*sqrt(T)
c1<-S0*pnorm(d1)-X*exp(-r*T)*pnorm(d2)
print(c1)

#4.(a)
rm(list = ls()[ls()!='wie_seq' & ls()!='norm_seq'])
r<-0.055
sigma<-0.2
S0<-88
T<-5
X<-100

len<-1000
#do each simulation from 1 to 10
E<-vector(l=10)
for(n in 1:10) {
  WT<-wie_seq(len,n)
  ST<-S0*exp(sigma*WT+(r-sigma**2/2)*n)
  E[n]<-mean(ST)
}
#plot
plot(0,0,main='E(Sn)',xlim=c(0,10),ylim=c(50,160),xlab='n',ylab='E(Sn)')
lines(E)

#4.(b) (c)
#plot 3 paths
t<-seq(0,10,length.out=1000)
for(i in 1:3) {
  path<-vector(l=1000)
  for(j in 1:1000) {
    WT<-wie_seq(1100,t[j])[3*i+j]
    path[j]<-S0*exp(sigma*WT+(r-sigma**2/2)*t[j])
  }
  lines(path)
}

#4.(d)

sigma<-0.3

#plot E(Sn) for sigma=0.3
E1<-vector(l=10)
for(n in 1:10) {
  WT<-wie_seq(len,n)
  ST<-S0*exp(sigma*WT+(r-sigma**2/2)*n)
  E1[n]<-mean(ST)
}
lines(E1,col='red')
#plot 3 paths for sigma=0.3
t<-seq(0,10,length.out=1000)
for(i in 1:3) {
  path<-vector(l=1000)
  for(j in 1:1000) {
    WT<-wie_seq(1100,t[j])[3*i+j]
    path[j]<-S0*exp(sigma*WT+(r-sigma**2/2)*t[j])
  }
  lines(path,col='red')
}
legend('topleft',legend=c('sigma=0.2','sigma=0.3'),col=c('black','red'),lwd=1)

#Lecture 3
#5.(a)
rm(list = ls()[ls()!='wie_seq' & ls()!='norm_seq'])
call_Euler<- function( S0,T,X,r,sigma){
  N<-1000
  len<-1000
  W<-wie_seq(len*N,T)
  price_seq<-vector(l=N)
  
  for (i in 1:N){
    WT<-W[(len*i-len+1):(len*i)]
    #simulate for ST
    ST<-S0*exp(sigma*WT+(r-sigma**2/2)*T)
    #estimate for call option price
    temp<-vector(l=len)
    for(j in 1:len) {
      temp[j]<-max(ST[j]-X,0)
    }
    c<-mean(temp)*exp(-r*T)
    price_seq[i]<-c
  }
  
  return(c(mean(price_seq),sd(price_seq)))
}
  
print(call_Euler(100,5,90,0.055,0.2))

#5.(b)
rm(list = ls()[ls()!='wie_seq' & ls()!='norm_seq'])
call_Milsh<- function( S0,T,X,r,sigma){
  N<-1000
  len<-1000
  W<-wie_seq(len*N,T)
  price_seq<-vector(l=N)
  
  for (i in 1:N){
    WT<-W[(len*i-len+1):(len*i)]
    #simulate for ST
    ST<-S0*exp(sigma*WT+(r-sigma**2/2)*T)
    #estimate for call option price
    temp<-vector(l=len)
    for(j in 1:len) {
      temp[j]<-max(ST[j]-X,0)
    }
    c<-mean(temp)*exp(-r*T)
    price_seq[i]<-c
  }
  
  return(c(mean(price_seq),sd(price_seq)))
}

print(call_Milsh(100,5,100,0.055,0.2))

#5.(c)
rm(list = ls()[ls()!='wie_seq' & ls()!='norm_seq'])
BS<- function( S0,T,X,r,sigma){
  d1<-(log(S0/X)+(r+sigma**2/2)*T)/sigma/sqrt(T)
  d2<-d1-sigma*sqrt(T)
  c1<-S0*pnorm(d1)-X*exp(-r*T)*pnorm(d2)
  return (c1)
}

print(BS(100,5,100,0.055,0.2))

#5.(e)

delta<-function(dt=0.05, S0,T,X,r,sigma){
  return ((BS( S0+dt,T,X,r,sigma)-BS( S0,T,X,r,sigma))/dt)
}

gamma<-function(dt=0.05, S0,T,X,r,sigma){
  return ((delta(dt, S0+dt,T,X,r,sigma)-delta(dt, S0,T,X,r,sigma))/dt)
}

theta<-function(dt=0.05, S0,T,X,r,sigma){
  return ((BS( S0,T+dt,X,r,sigma)-BS( S0,T,X,r,sigma))/dt)
}
vega<-function(dt=0.05, S0,T,X,r,sigma){
  return ((BS( S0,T,X,r,sigma+dt)-BS( S0,T,X,r,sigma))/dt)
}

#plot
X<-100
sigma<-0.25
r<-0.055
T<-0.5

delta_seq<-vector(l=11)
gamma_seq<-vector(l=11)
theta_seq<-vector(l=11)
vega_seq<-vector(l=11)
for (S in 95:105){
  delta_seq[S-94]<-delta(dt=0.05, S,T,X,r,sigma)
  gamma_seq[S-94]<-gamma(dt=0.05, S,T,X,r,sigma)
  theta_seq[S-94]<-theta(dt=0.05, S,T,X,r,sigma)
  vega_seq[S-94]<-vega(dt=0.05, S,T,X,r,sigma)
}

plot(0,0,main='delta',xlim=c(95,105),ylim=c(0.3,0.8),xlab='S0',ylab="value")
lines(95:105,delta_seq)

plot(0,0,main='gamma',xlim=c(95,105),ylim=c(0.015,0.025),xlab='S0',ylab="value")
lines(95:105,gamma_seq)

plot(0,0,main='theta',xlim=c(95,105),ylim=c(8,10),xlab='S0',ylab="value")
lines(95:105,theta_seq)

plot(0,0,main='vega',xlim=c(95,105),ylim=c(26,28),xlab='S0',ylab="value")
lines(95:105,vega_seq)

#6
rm(list = ls()[ls()!='wie_seq' & ls()!='norm_seq'])
rho<--0.6
r<-0.055
S0<-100
X<-100
V0<-0.05
sigma<-0.42
alpha<-5.8
beta<-0.0625
dt<-0.05
N<-10000

#simulate for St and Vt
St<-vector(l=N)
Vt<-vector(l=N)
St[1]<-S0
Vt[1]<-V0
for (i in 2:N){
  St[i]<-St[i-1]+r*St[i-1]*dt+sqrt(Vt[i-1])*St[i-1]*rnorm(1,0,dt)
  Vt[i]<-alpha*(beta-Vt[i-1])*dt+sigma*sqrt(Vt[i-1])*(rnorm(1,0,dt)**2)
}

#calculate option price
c<-St[N]*exp(-r*dt*N)
print(c)

#7
rm(list = ls()[ls()!='wie_seq' & ls()!='norm_seq'])

#multi-dimensional Halton sequence function
halton_seq<- function( n, dim=1, start=0, bases=NULL ){
    
    #   Get the first 100 primes.
    first.primes <- c(  
      2,   3,   5,   7,  11,  13,  17,  19,  23,  29,  31,  37,  41,  43,  47,  53,  59,  61,
      67,  71,  73,  79,  83,  89,  97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 
      157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 
      257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 
      367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 
      467, 479, 487, 491, 499, 503, 509, 521, 523, 541)
    if( dim > length(first.primes) ){
      first.primes <- primes(dim)  # will get up to first 1e8 primes
    }
    
    if( length(start) == 1 ){
      start <- rep(start, dim)
    }
    
    if( length(start) != dim ){
      stop( "The start vector must either have length 1 or length equal to the number of dimensions") 
    }
    
    
    #    Bases of the first dim dimensions of the Halton sequence
    if( is.null(bases)){
      bases <- first.primes[ 1:dim ]  
    }
    
    
    #   The locations of the numbers we want Halton numbers for.  First Halton number is 1, second is 2, etc.  This is the
    #   position in the Halton sequence we need.  After this, pos is a dim X n matrix where columns are the Halton numbers 
    #   we want. 
    pos <- t(sapply(start, FUN=function(x,k){ x:(x+k-1) }, k=n ))
    
    
    #   Find halton sequence numbers using the finite sum formula Blair cooked up
    if( (n == 1) & all(start==0) ){
      return( matrix(0,1,dim) )
    }
    
    n.sum.terms <- max(floor( log(start + n - 1)/ log(bases) ) + 1)  # number of digits in base b needed to represent maximum pos in each dim
    
    
    ans <- matrix( 0, nrow(pos), ncol(pos) )
    for( j in 0:n.sum.terms ){
      ans <- ans + ((pos %/% bases^j) %% bases) / bases^(j+1)
    }
    
    if( n == 1){
      return(matrix(ans,1))
    } else {
      return(t(ans))
    }
    
}

#integral function
f=function(x,y){
  tri_root=function(x){
    if (x>=0){
      return (x**(1/3))
    }
    else{
      return (-(-x)**(1/3))
    }
  }
  return (exp(-x*y)*(sin(6*pi*x)+tri_root(cos(2*pi*y))))
}

#Calculate for integral
N<-10000
b1<-2
b2<-3
halton<-halton_seq(N, dim=2)

value<-vector(l=N)
for (i in 1:N){
  value[i]<-f(halton[i,1],halton[i,2])
}
I<-mean(value)  
print(I)




