n<-100
p<-1600
R=0.5^abs(outer(1:p,1:p,"-"))
I=diag(rep(1,p))
T1=numeric(1000)
T2=numeric(1000)

theta1=rep(0,p) 
theta2=rep(0,p) #h0
#theta2=c(runif(p/2, min = 0, max = 1),rep(0,0.5*p))
#theta2=theta2/sqrt(crossprod(theta2))*sqrt(0.1*sqrt(sum(diag(R%*%R))+sum(diag(I%*%I))))

U=function(x) x/sqrt(crossprod(x))
for(k in 1:1000){
  mu=matrix(c(theta1,theta1),nrow=2,byrow=TRUE)
  Sigma=list(R,16*R)
  ww=c(0.9, 0.1)
  X1=rMVNmixture2(n, weight=ww, mean=mu, Sigma=Sigma)[[1]]
  
  mu=matrix(c(theta2,theta2),nrow=2,byrow=TRUE)
  Sigma=list(I,16*I)
  ww=c(0.9, 0.1)
  X2=rMVNmixture2(n, weight=ww, mean=mu, Sigma=Sigma)[[1]]
  #feng
  
  DD1=matrix(0,nrow=n,ncol=p)
  DD2=matrix(0,nrow=n,ncol=p)
  hat1=matrix(0,nrow=n,ncol=p)
  hat2=matrix(0,nrow=n,ncol=p)
  
  for(w in 1:n){
    X1rm=X1[-w,]
    hattheta1=apply(X1rm,2,mean)
    D1=apply(X1rm,2,var)
    for(i in 1:2){
      epsilon1=(diag(1/sqrt(D1)))%*%(t(X1rm)-hattheta1)
      deno1=sum(1/sqrt(apply(epsilon1,2,crossprod)))
      hattheta1=(hattheta1+(diag(sqrt(D1)))%*%apply(apply(epsilon1,2,U),1,sum)/deno1)[,1]
      D1=diag(p*diag(sqrt(D1))%*%diag(apply(apply(epsilon1,2,U),1,crossprod))%*%diag(sqrt(D1))/n)
    }
    DD1[w,]=1/sqrt(D1)
    hat1[w,]=hattheta1   
  }
  for(u in 1:n){
    X2rm=X2[-u,]
    hattheta2=apply(X2rm,2,mean)
    D2=apply(X2rm,2,var)
    for(j in 1:2){
      epsilon2=(diag(1/sqrt(D2)))%*%(t(X2rm)-hattheta2)
      deno2=sum(1/sqrt(apply(epsilon2,2,crossprod)))
      hattheta2=(hattheta2+diag(sqrt(D2))%*%apply(apply(epsilon2,2,U),1,sum)/deno2)[,1]
      D2=diag(p*diag(sqrt(D2))%*%diag(apply(apply(epsilon2,2,U),1,crossprod))%*%diag(sqrt(D2))/n)
    }
    DD2[u,]=1/sqrt(D2)
    hat2[u,]=hattheta2  
  }  
  #Rnn=0
  #for(w in 1:n){
  #  for(u in 1:n){
  #    Rnn=Rnn+U(as.vector((1/sqrt(DD1[[w]]))*t(t(X1[w,])-hat2[[u]])))%*%U(as.vector((1/sqrt(DD2[[u]]))*t(t(X2[u,])-hat1[[w]])))
  #  }
  #}
  Rnn=sum(apply((X1%x%rep(1,n)-rep(1,n)%x%hat2)*(DD1%x%rep(1,n)),1,U)*apply((rep(1,n)%x%X2-hat1%x%rep(1,n))*(rep(1,n)%x%DD2),1,U))
  
  #-----------
  hattheta1=apply(X1,2,mean)
  D1=apply(X1,2,var)
  for(i in 1:2){
    epsilon1=(diag(1/sqrt(D1)))%*%(t(X1)-hattheta1)
    deno1=sum(1/sqrt(apply(epsilon1,2,crossprod)))
    hattheta1=(hattheta1+(diag(sqrt(D1)))%*%apply(apply(epsilon1,2,U),1,sum)/deno1)[,1]
    D1=diag(p*diag(sqrt(D1))%*%diag(apply(apply(epsilon1,2,U),1,crossprod))%*%diag(sqrt(D1))/n)
  }
  
  hattheta2=apply(X2,2,mean)
  D2=apply(X2,2,var)
  for(j in 1:2){
    epsilon2=(diag(1/sqrt(D2)))%*%(t(X2)-hattheta2)
    deno2=sum(1/sqrt(apply(epsilon2,2,crossprod)))
    hattheta2=(hattheta2+diag(sqrt(D2))%*%apply(apply(epsilon2,2,U),1,sum)/deno2)[,1]
    D2=diag(p*diag(sqrt(D2))%*%diag(apply(apply(epsilon2,2,U),1,crossprod))%*%diag(sqrt(D2))/n)
  }
  
  hatc1=sum(1/sqrt(apply((diag(1/sqrt(D1)))%*%(t(X1)-hattheta1),2,crossprod)))/n
  hatc2=sum(1/sqrt(apply((diag(1/sqrt(D2)))%*%(t(X2)-hattheta2),2,crossprod)))/n
  
  UU1=apply((diag(1/sqrt(D1))%*%(t(X1)-hattheta2)),2,U)
  UU2=apply((diag(1/sqrt(D2))%*%(t(X2)-hattheta1)),2,U)
  
  tildeU1=apply((diag(1/sqrt(D1)))%*%(t(X1)-hattheta1),2,U)
  tildeU2=apply((diag(1/sqrt(D2)))%*%(t(X2)-hattheta2),2,U)
  
  A1=(t(tildeU1)%*%(diag(sqrt(D1)/sqrt(D2)))%*%tildeU1)^2
  A2=(t(tildeU2)%*%(diag(sqrt(D2)/sqrt(D1)))%*%tildeU2)^2  
  
  trA1=(sum(A1)-sum(diag(A1)))*p*p*(hatc2^2)/((hatc1^2)*n*(n-1))
  trA2=(sum(A2)-sum(diag(A2)))*p*p*(hatc1^2)/((hatc2^2)*n*(n-1))
  trA3=sum((t(tildeU1)%*%tildeU2)^2)*p*p/(n*n)
  
  bias1=hatc2/hatc1*sum(sqrt(D1/D2))/n/(p)
  bias2=hatc1/hatc2*sum(sqrt(D2/D1))/n/(p)
  sigman=2*trA1/(p*p*n*(n-1))+2*trA2/(p*p*n*(n-1))+4*trA3/(n*n*p*p)
  Rn=-sum(t(UU1)%*%(UU2))/(n*n)-bias1-bias2
  T1[k]=Rn/(sqrt(sigman))
  T2[k]=(-1)*Rnn/n/n/(sqrt(sigman))
  print(k)
  if(k%%100==0) {
    write(T1,"C:/Users/tp/Desktop/1.txt",ncol=1)
    write(T2,"C:/Users/tp/Desktop/2.txt",ncol=1)
  }
}

length(T1[T1>1.64])
length(T2[T2>1.64])
qqnorm(T1)
qqnorm(T2)


