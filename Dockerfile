# support-5860

FROM rocker/shiny:4.0.5
LABEL authors="Roy Francis"
ARG REPO="NBISweden/support-5860"

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get clean && \
    apt-get install -y git libxml2-dev libudunits2-dev libhdf5-dev && \
    rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'reqPkg = c("data.table", "DT", "ggdendro", "ggplot2", "ggplotify", "ggrepel", "glue", "grid", "hdf5r", "magrittr", "Matrix", "patchwork", "RColorBrewer", "readr", "remotes", "reticulate", "R.utils", "Seurat", "shiny", "shinycssloaders", "shinyhelper", "showtext", "shinythemes"); newPkg = reqPkg[!(reqPkg %in% installed.packages()[,"Package"])]; if(length(newPkg)){install.packages(newPkg)}'

RUN cd /srv/shiny-server/ && \
    git clone https://github.com/${REPO}.git temp && \
    ln -s temp/inst/app app && \
    sudo chown -R shiny:shiny /srv/shiny-server/app

EXPOSE 3838

CMD ["R", "-e", "shiny::runApp('/srv/shiny-server/app/', host = '0.0.0.0', port = 8787)"]

# docker build -t royfrancis/support-5860:latest . 
# docker run --rm -p 8787:8787 royfrancis/support-5860:latest

