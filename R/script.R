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
#' @import shinythemes shiny shinyhelper
#' @importFrom DT DTOutput renderDT
#' @importFrom readxl read_xlsx
#' @export
#' 
runShiny <- function(appDir = NULL, ...) {
  if(is.null(appDir)) {
    appDir <- system.file("app-sc", package="support5860")
  } else {
    if (appDir == "") {
      stop("Could not find app directory. Try re-installing `support5860`.", call. = FALSE)
    }
  }

  runApp(appDir,...)
}
