
## The function of iMAP
The input of iMAP includes z values (m by 2 matrix), annotation matrix (m by p matrix) and sample sizes for the two traits. An example data are ginven, including z values and annotation matrix. 

iMAP=function(

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
