# support5860

#' @title Launch interactive shiny application
#' @description Launches interactive shiny session in the browser.
#' @param ... Arguments are passed to \code{\link{runApp}}.
#' @return This function does not return anything
#' @seealso \code{\link{runApp}}
#' @examples
#' \dontrun{
#' library(support5860)
#' runShiny()
#' }
#' @import shinythemes shiny shinyhelper data.table hdf5r reticulate ggplot2 
#' @importFrom DT DTOutput renderDT
#' @importFrom shinycssloaders withSpinner
#' @export
#' 
runShiny <- function(...) {
    appDir <- system.file("app", package="support5860")
    if (appDir == "") {
      stop("Could not find app directory. Try re-installing `support5860`.", call. = FALSE)
    }

  runApp(appDir,...)
}
