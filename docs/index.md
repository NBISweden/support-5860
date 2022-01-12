<div class="wrapper-logo"><img class="logo" src="assets/logo.svg"></div>

This repository contains project reports and R shiny web application.

***

## Old reports

**13-Aug-2021**

- Integration of healthy BEC and LEC, Trajectory

**27-Sep-2021**

- Integration of healthy and diseased samples
- Integration of healthy samples
- Integration of diseased samples

## Current reports

**08-Nov-2021**

- [Integration combined](report-combined-integration-port.html) (healthy and diseased, healthy, diseased)

**09-Nov-2021**

- [Comparing different integration methods](report-compare-integration-port.html)

**10-Jan-2021**

- [Trajectory healthy-diseased LEC](report-healthy-diseased-lec-trajectory-port.html)

***

## Shiny app

### Online

[NBIS RShiny Server](https://rshiny.nbis.se/shiny-server-apps/rshiny-support-5860/)

### Local

The shiny app can be installed as an R package and run locally.

```{r,eval=FALSE}
# install dependencies
install.packages(c("shiny","DT","readxl","shinythemes","remotes"),
repos="https://cloud.r-project.org/")

# install package
remotes::install_github("NBISweden/support5860")

# load package
library(support5860)
runShiny()
```
