library(dplyr)
library(Seurat)

lec <- readRDS("~/proj-scendothelial/report/data/processed/healthy-diseased/lec/integrated/sct/16122021_human_LEC_SCT.rds")

lec@meta.data <- lec[[]][,c("nCount_RNA","nFeature_RNA","cell_type","study","tissue","condition","sample","percent_mito","percent_ribo","Phase","seurat_clusters")]
lec@assays$SCT <- NULL

library(ShinyCell)
scConf = createConfig(lec)
makeShinyApp(lec, scConf, gex.assay="RNA", gex.slot="data", gene.mapping = FALSE,
             shiny.title = "Support 5860", shiny.dir="inst/app1")
shiny::runApp("inst/app1/")


