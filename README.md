 
iMAP: integrative MApping of Pleiotropic association
========================================================================================================
## Introduction
**iMAP** is a method which performs integrative mapping of pleiotropic association and functional annotations using penalized Gaussian mixture models. Specifically, iMAP directly models summary statistics from GWASs, uses a multivariate Gaussian distribution to account for phenotypic correlation between traits, simultaneously infers genome-wide SNP association pattern using mixture modeling, and has the potential to reveal causal relationship between traits. Importantly, iMAP can integrate a large number of binary and continuous SNP functional annotations to substantially improve association mapping power, and, with a sparsity-inducing penalty term, is capable of selecting informative annotations from a large, potentially noninformative set. To enable scalable inference of iMAP to association studies with hundreds of thousands of individuals and millions of SNPs, we further develop an efficient expectation maximization algorithm based on a recently proposed approximation algorithm for penalized regression.

**[iMAP](https://github.com/biostatpzeng/iMAP/blob/master/iMAP.R)** is implemented in R statistical environment.

## Cite
[Ping Zeng](https://github.com/biostatpzeng), Xingjie Hao and Xiang Zhou. Pleiotropic Mapping and Annotation Selection in Genome-wide Association Studies with Pe-nalized Gaussian Mixture Models. bioRxiv 2018. [Doi: 10.1101/256461](https://www.biorxiv.org/content/early/2018/01/31/256461).
## Contact
We are very grateful to any questions, comments, or bugs reports; and please contact [Ping Zeng](https://github.com/biostatpzeng) via zpstat@xzhmu.edu.cn or ~~pingzeng@umich.edu~~.

## Update
2018-02-01  iMAP was updated to include elastic enet for selecting correlated annotations.



