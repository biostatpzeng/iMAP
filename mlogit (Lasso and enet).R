#*******************************************************
#*******************************************************
#Iteratively Re-Weighted Least Squares (IRWLS) for 
#binary logistic and multinomial logistic models 
#*******************************************************
#*******************************************************

options(scipen = 9)
#############################################
###############  multinomial logistic #######
#############################################

library(MASS)
library(glmnet)


#Lassomlogit=function(x0,propx,beta_ml,Var_ml)
Lassomlogit=function(x0,propx,beta_ml,Var_ml,alphax)
{
nL=50
m=dim(propx)[1]
coef0=as.vector(t(beta_ml))
covb=as.matrix(Var_ml)
sigma=ginv(covb)
spectrumd=eigen(sigma)
lam=spectrumd$values
V=spectrumd$vectors
L=length(lam)
qx=dim(x0)[2]
#penalty=rep(c(0,rep(1,qx)),3)*(1/abs(coef0))
penalty=rep(c(0,rep(1,qx)),3)
xstar=V%*%diag(sqrt(lam),nrow=L,ncol=L)%*%t(V)
ystar=xstar%*%as.vector(coef0)
options(warn=-1)

fit0=glmnet(xstar,ystar,alpha = alphax,family=c("gaussian"),
           standardize=F,intercept=F,penalty.factor=penalty)

lambdax=seq(from=range(fit0$lambda)[1]*0.2,to=range(fit0$lambda)[2]*0.8,length=nL)
fit0=glmnet(xstar,ystar,alpha = alphax,family=c("gaussian"),lambda=lambdax,
           standardize=F,intercept=F,penalty.factor=penalty)

BICx=matrix(NA,nL,3)
beta_path=matrix(NA,nL,L)

for (Lx in 1:nL)
{
beta_path[Lx,]=coef(fit0,s=lambdax[Lx])[-1,1]
dfx=sum(abs(beta_path[Lx,])>0)-3

#p1=exp(cbind(1,x0)%*%beta_path[Lx,(qx*0+1):(qx*1+1)])
#p2=exp(cbind(1,x0)%*%beta_path[Lx,(qx*1+2):(qx*2+2)])
#p3=exp(cbind(1,x0)%*%beta_path[Lx,(qx*2+3):(qx*3+3)])
#p1234=(p1+p2+p3+1)
#p1=p1/p1234
#p2=p2/p1234
#p3=p3/p1234
#p4= 1/p1234

BICx[Lx,1]=sum((ystar-xstar%*%beta_path[Lx,])^2)+log(length(ystar))*dfx
}
index=which.min(BICx[,1])
object = list(coef=matrix(beta_path[index,],3,qx+1,byrow=T),
	     index=index)
}

# in all the function, y can be 0-1 or a estimated values beteween 0-1;
# in mlogitfast y can be a vector (n X 1) or a matrix (n X k);
# mlogitfast==mlogitDA, but the mlogitfast is orders of magnitude 
# faster than mlogitDA


# in all the function, y can be 0-1 or a estimated values beteween 0-1;
# in mlogitfast y can be a vector (n X 1) or a matrix (n X k);
# mlogitfast==mlogitDA, but the mlogitfast is orders of magnitude 
# faster than mlogitDA


mlogitfast=function (x,y,intercept=TRUE,weights,reference,epsilon,maxit,...)
{

### control options
if (missing(reference))  reference=4
if (missing(maxit))      maxit    =200
if (missing(epsilon))    epsilon  =1e-6
if (missing(weights))    weights  =1

tx=1; delta=10e5

if (intercept==TRUE) x0 = cbind(1,x)
else                 x0 = x
p  = dim(x0)[2]
n  = dim(x0)[1]
#####################################################################
### process the data and sepcify the reference
### long y or wide y (here each row of y is 0 or 1)
qx=dim(as.matrix(y))[2] 
if (qx==1) {
k=length(sort(unique(y)))
##############  augment Y
catx=sort(unique(y))[-reference]
y0  =matrix(NA,n,k-1)
for (j in 1:(k-1))  { y0[,j]=(y==catx[j])*1 }
}
if (qx>1) {k=qx; y0=y[,-reference]}
############ process the data #################
#####################################################################


########  initial values
beta0=NULL
#system.time(
#for (j in 1:(k-1)) {
#beta0=c(beta0,glm(y0[,j]~x0-1)$coef)
#}
#)
#system.time(ginv(t(x0)%*%x0)%*%t(x0)%*%y0)
beta0=as.vector(ginv(t(x0)%*%x0)%*%t(x0)%*%y0)
betax=matrix(NA,maxit,p*(k-1))
betax[tx,]=beta0
########  initial values


########################################################################
##### logLike
logL=rep(0,maxit)
logL[tx]=-log(k)*n

logLike=function(x0,y0,betax)
{
p  = dim(x0)[2]
n  = dim(x0)[1]
k  = dim(y0)[2]+1
taux =matrix(NA,n,k-1)
T0=0
for (j in 1:(k-1))  {
taux[,j]=x0%*%betax[(p*(j-1)+1):(j*p)]
T0=T0+exp(taux[,j])
}

beta_x=0
y0_taux=taux*y0
for (j in 1:(k-1))  {
beta_x=beta_x+y0_taux[,j]
}
return (sum(beta_x-log(T0+1)))
}
##### logLike
########################################################################


while ((delta>epsilon) & (tx<maxit)) 
{
########################
taux=matrix(NA,n,k-1)
pai =matrix(NA,n,k-1)
for (j in 1:(k-1))  {
#taux[,j]=x0%*%betax[tx,(p*(j-1)+1):(j*p)]
taux[,j]=exp(x0%*%betax[tx,(p*(j-1)+1):(j*p)])
}

T0=0
#for (j in 1:(k-1))  {T0=T0+exp(taux[,j]) }
for (j in 1:(k-1))  {T0=T0+taux[,j] }
#T0=apply(taux,1,sum)
#for (j in 1:(k-1))  {
#pai[,j]=exp(taux[,j])/(T0+1)
#pai[,j]=taux[,j]/(T0+1)
#pai[,j]=ifelse(pai[,j]<1e-6,    1e-6,    pai[,j])
#pai[,j]=ifelse(pai[,j]>0.999999,0.999999,pai[,j])
#}

pai=taux/(T0+1)
pai=ifelse(pai<1e-6,    1e-6,    pai)
pai=ifelse(pai>0.999999,0.999999,pai)

##########   vector of U
#U=NULL
#for (j in 1:(k-1))  {U=c(U,t(x0)%*%(y0[,j]-pai[,j]))}
U=as.vector(t(x0)%*%(y0-pai))
U=as.vector(t(x0)%*%((y0-pai)*weights))# for weights
##########   matrix H 
#x_w=x0
H=matrix(0,p*(k-1),p*(k-1))
for (j in 1:(k-1)) {
	for (m in 1:j)  { 
		if (j==m) w=pai[,j]*(1-pai[,m])/2
		else      w=pai[,j]*(0-pai[,m])
		#for (i in 1:p) { x_w[,i]=x0[,i]*w} ### x_w=x*w
		w=weights*w # for weights
		H_jm=t(x0*w)%*%x0
		H[(p*(j-1)+1):(j*p),(p*(m-1)+1):(m*p)]=H_jm
	}
}

H=H+t(H)
b_tem=betax[tx,]+ginv(H)%*%U
delta1=max(abs(b_tem-betax[tx, ]))
tx=tx+1

#delta2=abs(logL[tx]-logL[tx-1])
#delta =min(delta2,delta1)
delta=delta1
betax[tx, ]=b_tem
#print(tx)
#print(delta)
}

betax=matrix(betax[tx, ],k-1,p,byrow=T)
if (intercept==TRUE) {colnames(betax)=(c("inter",colnames(x, do.NULL=F, prefix="x")))}
else colnames(betax)=colnames(x, do.NULL=F, prefix="x")

rownames(betax)=paste("cat_",seq(1:k)[-reference],sep="")
#,"_",rep(sort(unique(y))[-reference],each=p),sep="")

############## logLike
beta_x=0; y0_taux=taux*y0
for (j in 1:(k-1))  {beta_x=beta_x+y0_taux[,j]}
logL=sum(beta_x-log(T0+1))
############## logLike

#object = list(coef=betax,tx=tx,logL=logL[tx],Var=Var)
Var=ginv(H)
se=matrix(sqrt(diag(Var)),k-1,p,byrow=T)
zscore=betax/se
pvalue=(1-pnorm(abs(zscore),0,1))*2
colnames(se)=colnames(betax)
rownames(se)=rownames(betax)
colnames(zscore)=colnames(betax)
rownames(zscore)=rownames(betax)
colnames(pvalue)=colnames(betax)
rownames(pvalue)=rownames(betax)
colnames(Var)=rownames(Var)=rep(colnames(betax),k-1)
object = list(coef=betax,
	      se=se,zscore=zscore,
	      pvalue=pvalue,
	      Var=Var,
	      logL=logL,
	      tx=tx)
##end
}

