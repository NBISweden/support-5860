<div class="wrapper-logo"><img class="logo" src="assets/logo.svg"></div>

This repository contains project reports and R shiny web application.

***

## Old reports

- Integration of healthy BEC and LEC, Trajectory <span class="badge">13-Aug-2021</span>
- Integration of healthy and diseased samples <span class="badge">27-Sep-2021</span>
- Integration of healthy samples <span class="badge">27-Sep-2021</span>
- Integration of diseased samples <span class="badge">27-Sep-2021</span>

## Current reports

- [Integration combined](report-combined-integration-port.html) (healthy and diseased, healthy, diseased) <span class="badge">08-Nov-2021</span>
- [Comparing different integration methods](report-compare-integration-port.html) <span class="badge">09-Nov-2021</span>
- [Trajectory healthy-diseased LEC](report-healthy-diseased-lec-trajectory-port.html) <span class="badge">10-Jan-2021</span>
- [Trajectory healthy LEC](report-healthy-lec-trajectory-port.html) <span class="badge">12-Jan-2022</span>

***

## Shiny app

### Online

<img loading="lazy" src="https://www.scilifelab.se/wp-content/uploads/2021/03/scilifelab_logo_email.png" alt="scilifelab" height="40">

[`https://support5860.serve.scilifelab.se/`](https://support5860.serve.scilifelab.se/)

### Docker

<img loading="lazy" src="https://upload.wikimedia.org/wikipedia/commons/7/79/Docker_%28container_engine%29_logo.png" alt="docker" height="45">

The app can be run as a docker container if docker is installed on the system.

```
# run in a unix terminal
docker run --rm -p 8787:8787 royfrancis/support-5860:latest
```

### Local

<img loading="lazy" src="https://www.rstudio.com/assets/img/logo.svg" alt="rstudio" height="45">

The shiny app can be installed as an R package and run locally. Install R first.

```{r,eval=FALSE}
# run in R

# install dependencies (only installs packages that are not installed)
reqPkg = c("data.table", "DT", "ggdendro", "ggplot2", "ggplotify", "ggrepel", "glue",
	   "gridExtra", "hdf5r", "magrittr", "Matrix", "RColorBrewer", "readr",
	   "remotes", "reticulate", "R.utils", "Seurat", "shiny", "shinycssloaders",
	   "shinyhelper", "showtext", "shinythemes", "remotes")
newPkg = reqPkg[!(reqPkg %in% installed.packages()[,"Package"])]

if(length(newPkg)){
  install.packages(newPkg,repos="https://cloud.r-project.org/")
}

# install app package
remotes::install_github("NBISweden/support-5860")

# load package and run app
library(support5860)
runShiny()
```
