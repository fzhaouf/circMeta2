#' Example circRNA data from brain samples
#'
#' A GRanges object containing circRNA expression data from brain samples
#' including metadata for differential expression analysis.
#'
#' @format A GRanges object with circRNA coordinates and expression values
#' @source Processed from brain tissue samples
"BM10.circs"

#' Metadata for brain samples
#'
#' Sample metadata including condition, age, and sex information
#' for circRNA differential expression analysis.
#'
#' @format A data frame with sample information:
#' \describe{
#'   \item{sampleid}{Sample identifier}
#'   \item{condid}{Condition (1 = control, 2 = case)}
#'   \item{age}{Sample age}
#'   \item{sex}{Sample sex (M/F)}
#' }
#' @source Brain tissue sample collection
"metainfo"
