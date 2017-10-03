
source("multisummary.R")
source("MultiLogit_IRWLS.R")

iMAP=function(
	Zvalue, # z value; m by 2 matrix
	Var,    # Variance for z value; m by 2 matrix; it can be null
	anno,   # annotation for SNPs; m by q matrix; it can be null
	infor,  # basic information for SNPs;m by 4 matrix; 
	        # 1=rsid,2=chr,3=pos,4=maf (from reference) 
		# rsid must match with the order of the Zvalue matrix;
		# maf in infor is used to esmated the variance of the SNP effect sizes if no Var is given; importantly, either Var or infor should be given.
	#cutoff,
	n,       # sample sises for the two traits; 2 by 1 vector
	LSA=FALSE# if TRUE, perform Lasso for annotation selection;
	         # if FALSE, only estimate annotation coefficients.
	)  {
############ initial value for prop
m=dim(Zvalue)[1]
N=c(n[1],n[2])
z0    = 3
index1= (abs(Zvalue[,1])>=z0)*1
index2= (abs(Zvalue[,2])>=z0)*1
px1   = mean(index1)
px2   = mean(index2)
px12  = mean((index1==1) & (index2==1))
px1   = ifelse(px1==0,  1/m,     px1)
px2   = ifelse(px2==0,  1/m,     px2)
px12  = ifelse(px12==0, px1*px2, px12)
prop0 =c(px12,px1,px2,1-px12-px1-px2)
prop  =cbind(rep(prop0[1],m),rep(prop0[2],m),rep(prop0[3],m),rep(prop0[4],m))
if (missing(anno)) {prop=c(px12,px1,px2,1-px12-px1-px2)}
############ initial value for prop

############ initial value for beta
if (missing(Var)) {
se_est=1/(infor[,4]*(1-infor[,4]))
effect1=Zvalue[,1]*sqrt(se_est/n[1])
effect2=Zvalue[,2]*sqrt(se_est/n[2])
} 
else {
effect1=Zvalue[,1]*sqrt(Var[,1])
effect2=Zvalue[,2]*sqrt(Var[,2])
}
EFFECT=cbind(effect1,effect2)
############ initial value for beta

############ estimate sigma2
index5=which(apply(abs(Zvalue),1,max)<2)
V0=cor(Zvalue[index5,])
A0=V0/400
maxstep=100
delta = 100
V1=matrix(0,maxstep,4)
V2=matrix(0,maxstep,4)
V3=matrix(0,maxstep,4)
V00=V0
V1[1,] =as.vector(V00)
V2[1,1]=as.vector(V00)[1]
V3[1,4]=as.vector(V00)[4]
beta1=beta2=beta3=EFFECT
beta2[,2]=0
beta3[,1]=0
V1X=matrix(V1[1,],2,2)
V2X=matrix(V1[1,],2,2)
V3X=matrix(V1[1,],2,2)
Sigma=V0
beta1=as.matrix(beta1)
beta2=as.matrix(beta2)
beta3=as.matrix(beta3)
MS =500
mkx=matrix(0,MS+2,4)
pai=prop
S=1

########## start estimation
################################################
while ((delta>4) & (S<MS))
{

if  (missing(anno)) {
fit_prop=ComProp(V1X,V2X,V3X,beta1,beta2,beta3,Sigma,pai,N)
propx=fit_prop$propx
mk =apply(propx,2,sum)+prop0*m/100
pai=mk/sum(mk)
}

else {
fit_prop=ComProp_anno(V1X,V2X,V3X,beta1,beta2,beta3,Sigma,pai,N)
propx=fit_prop$propx
x0=as.matrix(anno)
fitmlogit=mlogitfast(x0,propx,intercept=T,reference=4)

if (LSA==TRUE) fitmlogit=Lassomlogit(x0,propx,fitmlogit$coef,fitmlogit$Var)

tau1=exp(cbind(1,x0)%*%(fitmlogit$coef)[1,])
tau2=exp(cbind(1,x0)%*%(fitmlogit$coef)[2,])
tau3=exp(cbind(1,x0)%*%(fitmlogit$coef)[3,])
taux=as.vector(1/(tau1+tau2+tau3+1))
pai=cbind(tau1,tau2,tau3,1)*taux
mk =apply(propx,2,sum)
}

fit_Vk=ComVk(beta1,beta2,beta3,propx,A0)
V1X=fit_Vk$V1
V2X=fit_Vk$V2
V3X=fit_Vk$V3

fit_beta=ComBeta(N,Sigma,V1=V1X,V2=V2X,V3=V3X,Zvalue)
beta1=fit_beta$beta1
beta2=fit_beta$beta2
beta3=fit_beta$beta3

mkx[S+1,]=mk
delta=max(abs(mkx[S+1,-4]-mkx[S,-4]))
print(round(mk),3)
print(delta)
S=S+1
}

pai=apply(propx,2,mean)
p11=pai[1]
p10=pai[1]+pai[2]
p01=pai[1]+pai[3]
#Ts=abs(sqrt(m)*(p11-p10*p01)/sqrt(p10*p01*(1-p10*p01)))
#pvalue=1-exp(-(Ts/sqrt(log(m)))^(-2))

#f1=sample(c(1,-1),m,replace=T)
#if (sum(Zvalue[,1]>=0)==m) Zvalue[,1]=Zvalue[,1]*f1
#if (sum(Zvalue[,2]>=0)==m) Zvalue[,2]=Zvalue[,2]*f1
#index=apply(abs(Zvalue),1,which.max)
#sig_fdr=locfdr(Zvalue[,index],plot=0)$fdr

#fdr1=locfdr((Zvalue[,1]),plot=0)$fdr
#fdr2=locfdr((Zvalue[,2]),plot=0)$fdr
#if (sum(Zvalue[,1]>=0)==m) fdr1=locfdr((Zvalue[,1])*f1,plot=0)$fdr
#if (sum(Zvalue[,2]>=0)==m) fdr2=locfdr((Zvalue[,2])*f1,plot=0)$fdr

fdr10=propx[,4]+propx[,3]
fdr01=propx[,4]+propx[,2]
fdr11=1-propx[,1]
fdr12=propx[,4]

##Lfdr=cbind(fdr10,fdr01,fdr11,fdr12,sig_fdr)
#Lfdr=cbind(fdr10,fdr01,fdr11,fdr12)
#if (missing(cutoff)) {
#cutoff0=c(0.3,0.2,0.1,0.05)
#cutoff=cutoff0
#fdr_est=matrix(NA,length(cutoff0),dim(Lfdr)[2])
#for (i in 1:length(cutoff0)) {
#fdr_est[i,]=apply(Lfdr<=cutoff0[i],2,sum)
#}
#}
#
#else {
#cutoff0=cutoff
#fdr_est=matrix(NA,length(cutoff0),dim(Lfdr)[2])
#for (i in 1:length(cutoff0)) {
#fdr_est[i,]=apply(Lfdr<=cutoff0[i],2,sum)
#}
#}

#fdr_est0=cbind(fdr_est[,1:3],fdr_est[,3]+fdr_est[,1],fdr_est[,3]+fdr_est[,2],fdr_est[,4:5])

colnames(Lfdr)    =c("fdr10","fdr01","fdr11","fdr12")
#colnames(fdr_est0)=c("fdr10","fdr01","fdr11","fdr10+fdr11","fdr01+fdr11")
#rownames(fdr_est0)=cutoff

if  (missing(anno)) {
object = list(
	mk  =rbind(apply(propx,2,sum),apply(propx,2,mean)),
	#fdr =fdr_est0,
	Lfdr=Lfdr,
	#Llike=fit_prop$Llike,
	#Tp  =c(Ts,pvalue),
	prop=cbind(infor,propx))}

else {
object = list(
	mk  =rbind(apply(propx,2,sum),apply(propx,2,mean)),
	#fdr =fdr_est0,
	Lfdr=Lfdr,
	#Llike=fit_prop$Llike,
	#Tp  =c(Ts,pvalue),
	prop=cbind(infor,propx),
	fitm=fitmlogit)}
}
