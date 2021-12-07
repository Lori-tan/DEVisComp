#' Launch Shiny App for DEVisComp
#'
#' A function that launches the Shiny app for DEVisComp.
#' The shiny app permit to plot MA plots, and Volcano plots
#' comparing the DESeq2 and edgeR input.
#'
#' @return No return value but open up a Shiny page.
#'
#' @examples
#' \dontrun{
#'
#' DEVisComp::runDEVisComp()
#' }
#'
#' @references
#' Grolemund, G. (2015). Learn Shiny - Video Tutorials. \href{https://shiny.rstudio.com/tutorial/}{Link}
#'
#' @export
#' @importFrom shiny runApp

runDEVisComp <- function() {
  appDir <- system.file("shiny-scripts",
                        package = "DEVisComp")
  shiny::runApp(appDir, display.mode = "normal")
  return()
}
# [END]
