#' diffanal2
#' @name diffanal2
#' @docType package
#' @import Biobase limma genefilter plyr
NULL

#' Include this data in my package. 
#' @title Lymphoma tissue samples from a study in 2010 
#' @docType data
#' @keywords datasets
#' @name lymphoma.2010
#' @usage data(lymphoma.2010)
#' @format An ExpressionSet with 116 samples and ~19000 genes from lymphoma patients
NULL

# Column names with aggregate meaning will have tokens joined by the JOIN.CHR
# e.g. cls1-cls2-p.value means the p.value comparing cls1 and cls2.
JOIN.CHR = ".."
SAMPLE.ID.COL <- '.rownames'
