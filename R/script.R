# support-5860

#' @title Launch interactive shiny application
#' @description Launches interactive shiny session in the browser.
#' @param display.mode Display mode. See \code{\link{runApp}}.
#' @param ... Further arguments are passed to \code{\link{runApp}}.
#' @return This function does not return anything
#' @seealso \code{\link{runApp}}
#' @examples
#' \dontrun{
#' library(support5860)
#' runShiny()
#' }
#' @import shinythemes shiny
#' @importFrom DT DTOutput renderDT
#' @importFrom readxl read_xlsx
#' @export
#' 
runShiny <- function(display.mode="normal",...) {
  appDir <- system.file("app", package="support5860")
  if (appDir == "") {
    stop("Could not find app directory. Try re-installing `support5860`.", call. = FALSE)
  }
  
  runApp(appDir,display.mode=display.mode,...)
}
