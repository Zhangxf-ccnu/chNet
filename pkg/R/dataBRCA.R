
#' @title TCGA breast cancer data
#' @description The TCGA breast cancer gene expression dataset used in our case study. The data (level 3,
#' Agilent G450 microarray, version: May 6 2017) is diwnloaded from the TCGA database using the
#' TCGA2STAT R package. It includes gene expression measurements for 231 luminal A cancers
#' and 95 basal-like cancers. The data only includes expression measurement of genes that overlap
#' with the breast cancer pathway (hsa05224) collected from the Kyoto Encyclopedia of Genes and
#' Genomes database. It includes an expression matrix for which the rows represent the samples and
#' the columns represent the genes, and a vector denoting the group membership of each sample.
#' @usage data("TCGA.BRCA")
#' @source [1] The Cancer Genome Atlas Research Network (2012), Comprehensive molecular portraits of human breast tumors. Nature. 490 (7418), 61-70. (http://cancergenome.nih.gov/)
#' @references  Jia-Juan Tu, Le Ou-Yang, Hong Yan, Hong Qin and Xiao-Fei Zhang(2020), Differential network analysis
#'  by simultaneously considering  changes in  gene interactions and gene expression.
#' @author Jia-Juan Tu
#' @seealso { \code{\link{generate.data}}, \code{\link{TCGA.BRCA}}, \code{\link{GSE13159.AML}}}

#' @examples
#' data("TCGA.BRCA")

"TCGA.BRCA"



