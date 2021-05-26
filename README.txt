README file for R package supporting the paper "Differential network analysis by simultaneously considering changes in gene interactions and gene expression".


Contents of this archive
------------------------
This archive contains 
(1) pkg: subdirectory that contains the R package.
(2) chNet-manual.pdf: reference manual.


The chNet package has the following R-package dependencies: MASS, Matrix, igraph, mvtnorm, glmnet. The dependents are automatically installed along with chNet. You can use the following commands to install chNet from GitHub. 


# Step 1. Install the devtools package. Invoke R and then type
install.packages("devtools") 

# Step 2. Load the devtools package.
library("devtools") 

# Step 3. Install the chNet package from GitHub.
install_github("Zhangxf-ccnu/chNet", subdir="pkg") 


Useage
Load the library chNet in R console, by running 
library(chNet)

Taking the TCGA breast cancer datasets as an example, run the following code:
data("TCGA.BRCA")
result = chNet(TCGA.BRCA$X,TCGA.BRCA$group, lambar = 2.85, parallel = FALSE, nCpus = 4)

Please do not hesitate to contact Dr. Xiao-Fei Zhang (zhangxf@mail.ccnu.edu.cn) or Miss Jia-Juan Tu (adaline_juan@mails.ccnu.edu.cn) to seek any clarifications regarding any  contents or operation of the archive.