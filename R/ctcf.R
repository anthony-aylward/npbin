#' @title Table of allele counts etc. from a ChIP-seq dataset.
#'
#' @format A data frame with 19782 rows and 12 variables:
#' \describe{
#'   \item{chr}{chromosome}
#'   \item{location}{genomic location based on hg19}
#'   \item{m}{total number of reads covering the SNP}
#'   \item{xm}{total number of reads at the SNP from the maternal allele}
#'   \item{winning.chip}{
#'     the allele with more ChIP-seq reads. "P" if xm < m/2 and "M" otherwise
#'   }
#'   \item{motif}{
#'     the ID and the transcription factor name of the motif in JASPAR database
#'     (Mathelier and others, 2013).
#'   }
#'   \item{pval.mat.atSNP}{
#'     the p-value of the motif on the maternal allele from R package
#'     \code{atSNP}
#'   }
#'   \item{winning.motif}{
#'     the allele with stronger motif, e.g. it is "M" if  pval.mat.atSNP <
#'     pval.pat.atSNP
#'   }
#'   \item{potential.TP}{
#'     whether it is a potential TP based on our criteria described in the
#'     Supplementary notes. The users can define it differently using different
#'     thresholds on the various p-values from atSNP.
#'   }
#'   \item{potential.FP}{
#'     whether it is a potential FP based on our criteria described in the
#'     Supplementary notes. The users can define it differently using different
#'     thresholds on the various p-values from atSNP.
#'   }
#' }
"ctcf"
