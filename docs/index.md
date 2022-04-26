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

[SciLifeLab Serve](https://support5860.serve.scilifelab.se/)

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
