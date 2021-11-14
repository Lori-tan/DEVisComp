#' Scaled Gene Counts From Airway Experiment
#'
#' An RNAseq experiment characterize the human airway smooth muscle transcriptome at
#' baseline and treated with a glucocorticosteroid (i.e. Dexamethasone (Dex), 1μM
#' for 18h).
#'
#' @source University of Pennsylvania, United States.
#'
#' @format A data frame with columns:
#' \describe{
#'  \item{ensgene}{Gene codes}
#'  \item{SRR1039508}{Cell line N61311 first condition, the control.}
#'  \item{SRR1039509}{Cell line N61311 Second condition, treated with Dex.}
#'  \item{SRR1039512}{Cell line N052611 first condition, the control.}
#'  \item{SRR1039513}{Cell line N052611 Second condition, treated with Dex.}
#'  \item{SRR1039516}{Cell line N080611 first condition, the control.}
#'  \item{SRR1039517}{Cell line N080611 Second condition, treated with Dex.}
#'  \item{SRR1039520}{Cell line N061011 first condition, the control.}
#'  \item{SRR1039521}{Cell line N061011 Second condition, treated with Dex.}
#' }
#' @examples
#' \dontrun{
#'  airwayCounts
#' }
#' @references
#' Himes et al. (2014). RNA-Seq transcriptome profiling identifies CRISPLD2 as a
#' glucocorticoid responsive gene that modulates cytokine function in airway smooth
#' muscle cells. \emph{PloS one, 9}(6), e99625. \href{https://doi.org/10.1371/journal.pone.0099625}{Link}

"airwayCounts"

#' Metadata for Airway Experiment
#'
#' The metadata of an RNAseq experiment characterize the human airway smooth muscle
#' transcriptome at baseline and treated with a glucocorticosteroid (i.e.
#' Dexamethasone (Dex), 1μM for 18h).
#'
#' @source University of Pennsylvania, United States.
#'
#' @format A data frame with columns:
#' \describe{
#'  \item{id}{Sample ID.}
#'  \item{dex}{Whether or not the sample is treated with dex.}
#'  \item{celltype}{Cell line.}
#'  \item{geo_id}{GEO ID.}
#'  }
#' @examples
#' \dontrun{
#'  airwayMetadata
#' }
#' @references
#' Himes et al. (2014). RNA-Seq transcriptome profiling identifies CRISPLD2 as a
#' glucocorticoid responsive gene that modulates cytokine function in airway smooth
#' muscle cells. \emph{PloS one, 9}(6), e99625. \href{https://doi.org/10.1371/journal.pone.0099625}{Link}

"airwayMetadata"

# [END]
