
## The function of iMAP
The input of iMAP includes z values (m by 2 matrix), annotation matrix (m by p matrix) and sample sizes for the two traits. An example data are ginven, including z values and annotation matrix. 

iMAP <- function(

	Zvalue, # z value; m by 2 matrix
	
	Var,    # Variance for z value; m by 2 matrix; it can be null
	
	anno,   # annotation for SNPs; m by q matrix; it can be null
	
	infor,  # basic information for SNPs;m by 4 matrix; 
	        # 1=rsid,2=chr,3=pos,4=maf (from reference) 
		# rsid must match with the order of the Zvalue matrix;
		# maf in infor is used to esmated the variance of the SNP effect sizes if no Var is given; importantly, either Var or infor should be given.
	n,       # sample sises for the two traits; 2 by 1 vector
	
	LSA=FALSE# if TRUE, perform Lasso for annotation selection;
	
	         # if FALSE, only estimate annotation coefficients.
	)  {

## No annotations are incorporated
n <-  c(10000,10000)
 
datx   <-  read.table("zvalue.txt", sep='\t', header=T)

Zvalue <-  datx[,c(4,6)]

Var    <-  datx[,c(5,7)]

inf1or  <-  datx[,1:3]

fit0   <-  iMAP(Zvalue,Var=Var,infor=infor,n=n,LSA=FALSE)

## Four annotations are incorporated

x1    <- read.table("annotation.txt", sep="\t", header=F)

annox <- apply(x1, 2, scale)

fit1  <- iMAP(Zvalue, anno=annox, Var=Var, infor=infor, n=n, LSA=FALSE)

> fit1$fitm$coef

 inter          V1          	 V2       	   V3     	     V4
 
-9.601150    -0.01445103	 -0.871568089  	0.59717230  	0.32469690
	
-3.483640     0.45298061 	 0.450329403 	 0.02177689  	0.04531874

-3.414173    0.02949652 	 0.003572322 	-0.46010296 	-0.48706708
 


## Annotation selection
fit2   <- iMAP(Zvalue, anno=annox, Var=Var, infor=infor, n=n, LSA=TRUE)
> fit2$fitm$coef

	 inter          V1     	 V2       	   V3       	   V4
 
 -11.285977	 0.0000000 	0.0000000  	0.0000000  	0.000000
 	
  -3.332366 	0.2883482 	0.2817387 	 0.0000000 	 0.000000
  
  -3.255403 	0.0000000	 0.0000000 	-0.2841287 	-0.321641
  


