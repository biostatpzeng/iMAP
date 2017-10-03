

library(MASS)
options(scipen=999)

ComProp=function(V1X,V2X,V3X,beta1,beta2,beta3,Sigma,pai,n)
{
invSigma=ginv(Sigma)
invV1=ginv(V1X)
invV2=ginv(V2X)
invV3=ginv(V3X)
Nn=diag(n)
ninvSigma=Nn%*%invSigma
M1=ninvSigma+invV1
M2=ninvSigma+invV2
M3=ninvSigma+invV3

c11=-0.5*log(det(M1     ))-0.5*log(det(V1X))
c21=-0.5*log(abs(M2[1,1]))-0.5*log(    V2X[1,1])
c31=-0.5*log(abs(M3[2,2]))-0.5*log(    V3X[2,2])

#c11=-0.5*log(det(M1))-0.5*log(det(V1X))
#c21=-0.5*log(det(M2))-0.5*log(det(V2X))
#c31=-0.5*log(det(M3))-0.5*log(det(V3X))

BB1=0.5*(beta1%*%M1)*beta1
BB2=0.5*(beta2%*%M2)*beta2
BB3=0.5*(beta3%*%M3)*beta3

BB1=BB1[,1]+BB1[,2]
BB2=BB2[,1]+BB2[,2]
BB3=BB3[,1]+BB3[,2]

c1=c11+log(pai[1])+BB1
c2=c21+log(pai[2])+BB2
c3=c31+log(pai[3])+BB3
c4=    log(pai[4])

propx=cbind(c1,c2,c3,c4)
Llike=sum(propx)
propx=exp(propx-apply(propx,1,max))
propx=propx/(propx[,1]+propx[,2]+propx[,3]+propx[,4])
object = list(propx=propx,Llike=Llike)
}



ComPropLD=function(V1X,V2X,V3X,beta1,beta2,beta3,Sigma,pai,n,LD)
{
invSigma=ginv(Sigma)
invV1=ginv(V1X)
invV2=ginv(V2X)
invV3=ginv(V3X)
Nn=diag(n)
ninvSigma=(Nn%*%invSigma)
m=dim(beta1)[1]
M1=array(matrix(NA,2,2),c(2,2,m))
M2=array(matrix(NA,2,2),c(2,2,m))
M3=array(matrix(NA,2,2),c(2,2,m))
c11=rep(NA,m)
c21=rep(NA,m)
c31=rep(NA,m)
BB1=matrix(NA,m,2)
BB2=matrix(NA,m,2)
BB3=matrix(NA,m,2)
for (j in 1:m) 
{
M1[,,j]=(1/LD[j])*ninvSigma+invV1
M2[,,j]=(1/LD[j])*ninvSigma+invV2
M3[,,j]=(1/LD[j])*ninvSigma+invV3
c11[j]=-0.5*log(det(M1[,,j]     ))-0.5*log(det(V1X))
c21[j]=-0.5*log(abs(M2[,,j][1,1]))-0.5*log(    V2X[1,1])
c31[j]=-0.5*log(abs(M3[,,j][2,2]))-0.5*log(    V3X[2,2])
BB1[j,]=0.5*(beta1[j,]%*%M1[,,j])*beta1[j,]
BB2[j,]=0.5*(beta2[j,]%*%M2[,,j])*beta2[j,]
BB3[j,]=0.5*(beta3[j,]%*%M3[,,j])*beta3[j,]
}

BB1=BB1[,1]+BB1[,2]
BB2=BB2[,1]+BB2[,2]
BB3=BB3[,1]+BB3[,2]

c1=c11+log(pai[1])+BB1
c2=c21+log(pai[2])+BB2
c3=c31+log(pai[3])+BB3
c4=    log(pai[4])

propx=cbind(c1,c2,c3,c4)
Llike=sum(propx)
propx=exp(propx-apply(propx,1,max))
propx=propx/(propx[,1]+propx[,2]+propx[,3]+propx[,4])
object = list(propx=propx,Llike=Llike)
}



ComProp_anno=function(V1X,V2X,V3X,beta1,beta2,beta3,Sigma,prop,n)
{
invSigma=ginv(Sigma)
invV1=ginv(V1X)
invV2=ginv(V2X)
invV3=ginv(V3X)
Nn=diag(n)
ninvSigma=Nn%*%invSigma
M1=ninvSigma+invV1
M2=ninvSigma+invV2
M3=ninvSigma+invV3

c11=-0.5*log(det(M1     ))-0.5*log(det(V1X))
c21=-0.5*log(abs(M2[1,1]))-0.5*log(    V2X[1,1])
c31=-0.5*log(abs(M3[2,2]))-0.5*log(    V3X[2,2])

#c11=-0.5*log(det(M1))-0.5*log(det(V1X))
#c21=-0.5*log(det(M2))-0.5*log(det(V2X))
#c31=-0.5*log(det(M3))-0.5*log(det(V3X))

BB1=0.5*(beta1%*%M1)*beta1
BB2=0.5*(beta2%*%M2)*beta2
BB3=0.5*(beta3%*%M3)*beta3

BB1=BB1[,1]+BB1[,2]
BB2=BB2[,1]+BB2[,2]
BB3=BB3[,1]+BB3[,2]

c1=c11+log(prop[,1])+BB1
c2=c21+log(prop[,2])+BB2
c3=c31+log(prop[,3])+BB3
c4=    log(prop[,4])

propx=cbind(c1,c2,c3,c4)
Llike=sum(propx)
propx=exp(propx-apply(propx,1,max))
propx=propx/(propx[,1]+propx[,2]+propx[,3]+propx[,4])
object = list(propx=propx,Llike=Llike)
}

ComProp_annoLD=function(V1X,V2X,V3X,beta1,beta2,beta3,Sigma,prop,n,LD)
{
invSigma=ginv(Sigma)
invV1=ginv(V1X)
invV2=ginv(V2X)
invV3=ginv(V3X)
Nn=diag(n)
ninvSigma=(Nn%*%invSigma)
m=dim(beta1)[1]
c11=rep(NA,m)
c21=rep(NA,m)
c31=rep(NA,m)
BB1=matrix(NA,m,2)
BB2=matrix(NA,m,2)
BB3=matrix(NA,m,2)
for (j in 1:m) 
{
M1=(1/LD[j])*ninvSigma+invV1
M2=(1/LD[j])*ninvSigma+invV2
M3=(1/LD[j])*ninvSigma+invV3
c11[j]=-0.5*log(det(M1     ))-0.5*log(det(V1X))
c21[j]=-0.5*log(abs(M2[1,1]))-0.5*log(    V2X[1,1])
c31[j]=-0.5*log(abs(M3[2,2]))-0.5*log(    V3X[2,2])
BB1[j,]=0.5*(beta1[j,]%*%M1)*beta1[j,]
BB2[j,]=0.5*(beta2[j,]%*%M2)*beta2[j,]
BB3[j,]=0.5*(beta3[j,]%*%M3)*beta3[j,]
}

BB1=BB1[,1]+BB1[,2]
BB2=BB2[,1]+BB2[,2]
BB3=BB3[,1]+BB3[,2]

c1=c11+log(prop[,1])+BB1
c2=c21+log(prop[,2])+BB2
c3=c31+log(prop[,3])+BB3
c4=    log(prop[,4])

propx=cbind(c1,c2,c3,c4)
Llike=sum(propx)
propx=exp(propx-apply(propx,1,max))
propx=propx/(propx[,1]+propx[,2]+propx[,3]+propx[,4])
object = list(propx=propx,Llike=Llike)
}


ComBeta=function(n,Sigma,V1X,V2X,V3X,Z)
{
invSigma=ginv(Sigma)
Nn=diag(n)
Z=as.matrix(Z)
ninvSigma=Nn%*%invSigma
yx=as.matrix(Z%*%sqrt(Nn))
ninvSigmaV1=ginv(ninvSigma+ginv(V1X))%*%invSigma
ninvSigmaV2=ginv(ninvSigma+ginv(V2X))%*%invSigma
ninvSigmaV3=ginv(ninvSigma+ginv(V3X))%*%invSigma
beta1  = yx%*%ninvSigmaV1
beta2  = yx%*%ninvSigmaV2
beta3  = yx%*%ninvSigmaV3
beta2[,2]  = 0
beta3[,1]  = 0
object = list(beta1=beta1,beta2=beta2,beta3=beta3)
}

ComBetaLD=function(n,Sigma,V1X,V2X,V3X,Z,LD)
{
invSigma=ginv(Sigma)
Nn=diag(n)
Z=as.matrix(Z)
ninvSigma=(Nn%*%invSigma)
m=length(LD)
yx=as.matrix((1/LD)*(Z%*%sqrt(Nn)))
beta1=matrix(NA,m,2)
beta2=matrix(NA,m,2)
beta3=matrix(NA,m,2)
for (j in 1:m)
{
ninvSigma_LD=ninvSigma/LD[j]
ninvSigmaV1=ginv(ninvSigma_LD+ginv(V1X))%*%invSigma
ninvSigmaV2=ginv(ninvSigma_LD+ginv(V2X))%*%invSigma
ninvSigmaV3=ginv(ninvSigma_LD+ginv(V3X))%*%invSigma
beta1[j,]= yx[j,]%*%ninvSigmaV1
beta2[j,]= yx[j,]%*%ninvSigmaV2
beta3[j,]= yx[j,]%*%ninvSigmaV3
}
beta2[,2]  = 0
beta3[,1]  = 0
object = list(beta1=beta1,beta2=beta2,beta3=beta3)
}


ComVk=function(beta1,beta2,beta3,propx,A0)
{
mk =apply(propx,2,sum)
beta1x=beta1*sqrt(propx[,1])
beta2x=beta2*sqrt(propx[,2])
beta3x=beta3*sqrt(propx[,3])
#beta1x=as.matrix(beta1x)
#beta2x=as.matrix(beta2x)
#beta3x=as.matrix(beta3x)
#m  =0.002*dim(beta1)[1]
m  =0.002*dim(beta1)[1]
V1=(t(beta1x)%*%beta1x+m*A0)/(mk[1]+m)
#m  =0.003*dim(beta1)[1]
m  =0.003*dim(beta1)[1]
V2=(t(beta2x)%*%beta2x+m*A0)/(mk[2]+m)
V3=(t(beta3x)%*%beta3x+m*A0)/(mk[3]+m)

V2[1,2]=V2[2,1]=V2[2,2]=0
V3[1,1]=V3[1,2]=V3[2,1]=0
object = list(V1=V1,V2=V2,V3=V3)
}



