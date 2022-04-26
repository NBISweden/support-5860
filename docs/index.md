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

<img loading="lazy" src="https://www.scilifelab.se/wp-content/uploads/2021/03/scilifelab_logo_email.png" alt="scilifelab" width="194" height="43">

[SciLifeLab Serve](https://support5860.serve.scilifelab.se/)

### Docker

<img loading="lazy" src="https://marvel-b1-cdn.bc0a.com/f00000000152152/www.zend.com/sites/default/files/image/2019-09/logo-docker.jpg" alt="docker" height="50">

```
docker run --rm -p 8787:8787 royfrancis/support-5860:latest
```

### Local

<img loading="lazy" src="https://blog.desdelinux.net/wp-content/uploads/2019/02/rstudio-og.png.webp" alt="rstudio" height="50">

The shiny app can be installed as an R package and run locally.

```{r,eval=FALSE}
# install dependencies

reqPkg = c("data.table", "DT", "ggdendro", "ggplot2", "ggplotify", "ggrepel", "glue", "gridExtra", "hdf5r", "magrittr", "Matrix", "RColorBrewer", "readr", "remotes", "reticulate", "R.utils", "Seurat", "shiny", "shinycssloaders", "shinyhelper", "showtext", "shinythemes", "remotes")
newPkg = reqPkg[!(reqPkg %in% installed.packages()[,"Package"])]

if(length(newPkg)){
  install.packages(newPkg,repos="https://cloud.r-project.org/")
}

# install package
remotes::install_github("NBISweden/support-5860")

# load package
library(support5860)
runShiny()
```
