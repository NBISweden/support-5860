library(shiny)
library(shinyhelper)
library(data.table)
library(Matrix)
library(DT)
library(magrittr)
library(ggplot2)
library(ggplotify)
library(ggrepel)
library(hdf5r)
library(ggdendro)
library(grid)
library(shinycssloaders)
library(patchwork)

# load font for plot
sysfonts::font_add_google(name = "Lato", family = "Lato")
showtext::showtext_auto()

if(!exists("lcconf")) lcconf = readRDS("lcconf.rds")
if(!exists("lcdef")) lcdef  = readRDS("lcdef.rds")
if(!exists("lcgene")) lcgene = readRDS("lcgene.rds")
if(!exists("lcmeta")) lcmeta = readRDS("lcmeta.rds")
if(!exists("lcmar")) lcmar = readRDS("lcmar.rds")
if(!exists("ldconf")) ldconf = readRDS("ldconf.rds")
if(!exists("lddef")) lddef  = readRDS("lddef.rds")
if(!exists("ldgene")) ldgene = readRDS("ldgene.rds")
if(!exists("ldmeta")) ldmeta = readRDS("ldmeta.rds")
if(!exists("ldmar")) ldmar = readRDS("ldmar.rds")
if(!exists("lhconf")) lhconf = readRDS("lhconf.rds")
if(!exists("lhdef")) lhdef  = readRDS("lhdef.rds")
if(!exists("lhgene")) lhgene = readRDS("lhgene.rds")
if(!exists("lhmeta")) lhmeta = readRDS("lhmeta.rds")
if(!exists("lhmar")) lhmar = readRDS("lhmar.rds")
if(!exists("vhconf")) vhconf = readRDS("vhconf.rds")
if(!exists("vhdef")) vhdef  = readRDS("vhdef.rds")
if(!exists("vhgene")) vhgene = readRDS("vhgene.rds")
if(!exists("vhmeta")) vhmeta = readRDS("vhmeta.rds")
if(!exists("vhmar")) vhmar = readRDS("vhmar.rds")

### Initialise variables and functions ----

# Colour palette
cList <- list(
  c(
    "grey85", "#FFF7EC", "#FEE8C8", "#FDD49E", "#FDBB84",
    "#FC8D59", "#EF6548", "#D7301F", "#B30000", "#7F0000"
  ),
  c(
    "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFBF",
    "#FEE090", "#FDAE61", "#F46D43", "#D73027"
  )[c(1, 1:9, 9)],
  c(
    "#FDE725", "#AADC32", "#5DC863", "#27AD81", "#21908C",
    "#2C728E", "#3B528B", "#472D7B", "#440154"
  )
)
names(cList) <- c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple")

# Panel sizes
pList <- c("400px", "600px", "800px")
names(pList) <- c("Small", "Medium", "Large")
pList2 <- c("500px", "700px", "900px")
names(pList2) <- c("Small", "Medium", "Large")
pList3 <- c("600px", "800px", "1000px")
names(pList3) <- c("Small", "Medium", "Large")
# baseplot font size
sList <- c(12, 14, 18, 22)
names(sList) <- c("Smaller", "Small", "Medium", "Large")
# ggrepel font size
lList <- c(5, 6, 7)
names(lList) <- c("Small", "Medium", "Large")

# Function to extract legend
g_legend <- function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

# progress indicator
show_progress <- function(...){
  return(shinycssloaders::withSpinner(..., type = 7, color = "#95a5a6"))
}

# Plot theme
# @description Custom ggplot theme
# @param base_size (Numeric) Base font size
# @param XYval (Logical) Show XY axes text?
# @param Xang (Numeric) X axis text angle
# @param XjusH (Numeric) X axis horizontal justification
# @param lpos (Character) Position of Legend
# @param font (Character) Google font
# @param col_text (Character) Text colour
# @param col_line (Character) Line colour
#
sctheme <- function(base_size = 24, XYval = TRUE, Xang = 0, XjusH = 0.5, lpos = "bottom", font = "Lato", col_text = "grey30", col_line = "grey60") {
  oupTheme <- theme(
    text = element_text(size = base_size, family = font),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", colour = "white"),
    axis.line = element_line(colour = col_line),
    axis.ticks = element_line(colour = col_line),
    axis.title = element_text(colour = col_text),
    axis.text = element_text(size = base_size, colour = col_text),
    axis.text.x = element_text(angle = Xang, hjust = XjusH),
    strip.background = element_rect(colour = "white"),
    legend.position = lpos,
    legend.key = element_rect(colour = NA, fill = NA)
  )
  if (!XYval) {
    oupTheme <- oupTheme +
      theme(
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank()
      )
  }
  
  return(oupTheme)
}

### Plotting functions ----

# @description DR scatterplot for gene expression
# @param dtab (Data.table) with columns X, Y, geneName and val
# @param bgCells (Logical) Background points
# @param inpdrX (Character) X axis variable for DR
# @param inpdrY (Character) Y axis variable for DR
# @param inpsiz (Numeric) Point size
# @param inpcol (Character) Custom colour label
# @param inpord (Character) Custom plotting order
# @param inpfsz (Character) Custom font size
# @param inppasp (Character) Custom aspect ratio
# @param inptxt (Logical) Show XY labels
#
scScatter <- function(dtab, bgCells = FALSE, inpdrX, inpdrY, inpsiz, inpcol, inpfsz, inpasp, inptxt){
  
  if(any(!c("X", "Y", "geneName", "val") %in% colnames(dtab))) "Input missing one or more columns: X, Y, geneName, val."
  
  ggOut <- ggplot(dtab, aes(X, Y, color = val))
  rat <- (max(dtab$X) - min(dtab$X)) / (max(dtab$Y) - min(dtab$Y))
  ltitle <- dtab$geneName[1]
  
  if (bgCells) {
    ggOut <- ggOut +
      geom_point(data = ggData2, color = "snow2", size = inpsiz, shape = 20)
  }
  
  ggOut <- ggOut +
    geom_point(size = inpsiz, shape = 20) + xlab(inpdrX) + ylab(inpdrY) +
    scale_color_gradientn(ltitle, colours = cList[[inpcol]]) +
    #guides(color = guide_colorbar(barwidth = 20)) + 
    sctheme(base_size = sList[inpfsz], XYval = inptxt, lpos = "right")
  
  if (inpasp == "Square") {
    ggOut <- ggOut + coord_fixed(ratio = rat)
  } else if (inpasp == "Fixed") {
    ggOut <- ggOut + coord_fixed()
  }
  
  return(ggOut)
}

# @description Plot gene expression on dimred for any number of input genes
# @param inpConf (data.frame) Configuration table
# @param inpMeta (data.frame) Metadata table
# @param inpdrX (Character) X axis variable for DR
# @param inpdrY (Character) Y axis variable for DR
# @param inp (Character) Gene name to use
# @param inpsub1 (Character) Name of metadata column for subsetting
# @param inpsub2 (Character/Vector) Levels under metadata column for subsetting
# @param inpH5 (Character) Path to expression h5
# @param inpGene (Numeric) Named gene expression vector
# @param inpsiz (Numeric) Point size
# @param inpcol (Character) Custom colour label
# @param inpord (Character) Custom plotting order
# @param inpfsz (Character) Custom font size
# @param inppasp (Character) Custom aspect ratio
# @param inptxt (Logical) Show XY labels
# @param inpncol (Integer) Number of rows of plots
# @details 
# Config table contains columns ID (Character, columns name in metadata), UI (Character, UI id), fID (Character, Levels for categorical data, | separated), fCL (Character, Colours for categorical data, | separated), fRow (Integer, number of rows), grp (Logical), dimred (Logical)
#
scFeature <- function(inpConf, inpMeta, inpdrX, inpdrY, inp, inpsub1, inpsub2, inpH5, inpGene, inpsiz, inpcol, inpord, inpfsz, inpasp, inptxt, inpncol = 0){
  
  if (is.null(inpsub1)) inpsub1 <- inpConf$UI[1]
  
  # Identify genes that are in our dataset
  geneList <- scGeneList(inp, inpGene)
  geneList <- geneList[present == TRUE]
  shiny::validate(need(nrow(geneList) <= 36, "More than 36 genes to plot! Please reduce the gene list!"))
  shiny::validate(need(nrow(geneList) > 0, "Please input at least 1 gene to plot!"))
  
  # Prepare ggData
  h5file <- H5File$new(inpH5, mode = "r")
  h5data <- h5file[["grp"]][["data"]]
  ggData <- data.table()
  for (iGene in geneList$gene) {
    tmp <- inpMeta[, c("sampleID",inpConf[UI == inpdrX]$ID, inpConf[UI == inpdrY]$ID, inpConf[UI == inpsub1]$ID),with = FALSE]
    colnames(tmp) <- c("sampleID", "X", "Y", "sub")
    tmp$geneName <- iGene
    tmp$val <- h5data$read(args = list(inpGene[iGene], quote(expr = )))
    ggData <- rbindlist(list(ggData, tmp))
  }
  h5file$close_all()
  if (length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)) ggData <- ggData[sub %in% inpsub2]
  
  bgCells <- FALSE
  if (length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)) {
    bgCells <- TRUE
    ggData2 <- ggData[!sub %in% inpsub2]
    ggData <- ggData[sub %in% inpsub2]
  }
  
  if (inpord == "Max") {
    ggData <- ggData[order(val)]
  } else if (inpord == "Min") {
    ggData <- ggData[order(-val)]
  } else if (inpord == "Random") {
    ggData <- ggData[sample(nrow(ggData))]
  }
  
  ggDataSplit <- split(ggData, by=c("geneName"), flatten=FALSE)
  plist <- vector("list", length = length(ggDataSplit))
  for(i in seq_along(ggDataSplit)){
    plist[[i]] <- scScatter(dtab = ggDataSplit[[i]], bgCells, inpdrX, inpdrY, inpsiz, inpcol, inpfsz, inpasp, inptxt)
  }

  if(inpncol < 1) inpncol <- floor(sqrt(length(plist)))
  ggOut <- wrap_plots(plist) + plot_layout(ncol = inpncol)
  return(ggOut)
}

# @description Plot cell information on dimred
# @param inpConf (data.frame) Configuration table
# @param inpMeta (data.frame) Metadata table
# @param inpdrX (Character) X axis variable for DR
# @param inpdrY (Character) Y axis variable for DR
# @param inp1 (Character) Gene name to use
# @param inpsub1 (Character) Name of metadata column for subsetting
# @param inpsub2 (Character/Vector) Levels under metadata column for subsetting
# @param inpsiz (Numeric) Point size
# @param inpcol (Character) Custom colour label
# @param inpord (Character) Custom plotting order
# @param inpfsz (Character) Custom font size
# @param inppasp (Character) Custom aspect ratio
# @param inptxt (Logical) Show XY labels
# @param inplab 
#
scDRcell <- function(inpConf, inpMeta, inpdrX, inpdrY, inp1, inpsub1, inpsub2, inpsiz, inpcol, inpord, inpfsz, inpasp, inptxt, inplab) {
  
  if (is.null(inpsub1)) { inpsub1 <- inpConf$UI[1] }

  # Prepare ggData
  ggData <- inpMeta[, c(
    inpConf[UI == inpdrX]$ID, inpConf[UI == inpdrY]$ID,
    inpConf[UI == inp1]$ID, inpConf[UI == inpsub1]$ID
  ),
  with = FALSE
  ]
  colnames(ggData) <- c("X", "Y", "val", "sub")
  rat <- (max(ggData$X) - min(ggData$X)) / (max(ggData$Y) - min(ggData$Y))
  bgCells <- FALSE
  if (length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)) {
    bgCells <- TRUE
    ggData2 <- ggData[!sub %in% inpsub2]
    ggData <- ggData[sub %in% inpsub2]
  }

  if (inpord == "Max") {
    ggData <- ggData[order(val)]
  } else if (inpord == "Min") {
    ggData <- ggData[order(-val)]
  } else if (inpord == "Random") {
    ggData <- ggData[sample(nrow(ggData))]
  }

  # Do factoring if required
  if (!is.na(inpConf[UI == inp1]$fCL)) {
    ggCol <- strsplit(inpConf[UI == inp1]$fCL, "\\|")[[1]]
    names(ggCol) <- levels(ggData$val)
    ggLvl <- levels(ggData$val)[levels(ggData$val) %in% unique(ggData$val)]
    ggData$val <- factor(ggData$val, levels = ggLvl)
    ggCol <- ggCol[ggLvl]
  }

  # Actual ggplot
  ggOut <- ggplot(ggData, aes(X, Y, color = val))
  if (bgCells) {
  ggOut <- ggOut +
    geom_point(data = ggData2, color = "snow2", size = inpsiz, shape = 20)
  }
  ggOut <- ggOut +
    geom_point(size = inpsiz, shape = 20) +
    xlab(inpdrX) +
    ylab(inpdrY) +
    sctheme(base_size = sList[inpfsz], XYval = inptxt)

  if (is.na(inpConf[UI == inp1]$fCL)) {
    ggOut <- ggOut +
    scale_color_gradientn("", colours = cList[[inpcol]]) +
    guides(color = guide_colorbar(barwidth = 20)) 
  } else {
    sListX <- min(nchar(paste0(levels(ggData$val), collapse = "")), 200)
    sListX <- 0.75 * (sList - (1.5 * floor(sListX / 50)))
    ggOut <- ggOut + scale_color_manual("", values = ggCol) +
      guides(color = guide_legend(
        override.aes = list(size = 5),
        nrow = inpConf[UI == inp1]$fRow
      )) +
      theme(legend.text = element_text(size = sListX[inpfsz]))

    if (inplab) { ggData3 <- ggData[, .(X = mean(X), Y = mean(Y)), by = "val"]
      lListX <- min(nchar(paste0(ggData3$val, collapse = "")), 200)
      lListX <- lList - (0.25 * floor(lListX / 50))
      ggOut <- ggOut +
        geom_text_repel(
          data = ggData3, aes(X, Y, label = val),
          color = "grey10", bg.color = "grey95", bg.r = 0.15,
          size = lListX[inpfsz], seed = 42
        ) }
    }

  if (inpasp == "Square") { ggOut <- ggOut + coord_fixed(ratio = rat) } else if (inpasp == "Fixed") { ggOut <- ggOut + coord_fixed() }

  if(is.numeric(inpMeta[[inp1]])) ggOut <- ggOut + guides(color = guide_colorbar(barwidth = 25))
  return(ggOut)
}

# @description Cell number, stats data
# @param inpConf (data.frame) Configuration table
# @param inpMeta (data.frame) Metadata table
# @param inp1 (Character) Name of metadata column for coloring points
# @param inp2 (Character) Gene name to use
# @param inpsub1 (Character) Name of metadata column for subsetting
# @param inpsub2 (Character/Vector) Levels under metadata column for subsetting
# @param inpH5 (Character) Name/path to h5 object
# @param inpGene (Character) Gene to use
# @param inpsplt (Character) Whether to split table by Quartile or Decile
#
scDRnum <- function(inpConf, inpMeta, inp1, inp2, inpsub1, inpsub2, inpH5, inpGene, inpsplt) {
  
  if (is.null(inpsub1)) { inpsub1 <- inpConf$UI[1] }

  # Prepare ggData
  ggData <- inpMeta[, c(inpConf[UI == inp1]$ID, inpConf[UI == inpsub1]$ID),
    with = FALSE
  ]

  colnames(ggData) <- c("group", "sub")
  h5file <- H5File$new(inpH5, mode = "r")
  h5data <- h5file[["grp"]][["data"]]
  ggData$val2 <- h5data$read(args = list(inpGene[inp2], quote(expr = T)))
  ggData[val2 < 0]$val2 <- 0
  h5file$close_all()
  if (length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)) { ggData <- ggData[sub %in% inpsub2] }

  # Split inp1 if necessary
  if (is.na(inpConf[UI == inp1]$fCL)) { if (inpsplt == "Quartile") { nBk <- 4 }
    if (inpsplt == "Decile") { nBk <- 10 }
    ggData$group <- cut(ggData$group, breaks = nBk) }

  # Actual data.table
  ggData$express <- FALSE
  ggData[val2 > 0]$express <- TRUE
  ggData1 <- ggData[express == TRUE, .(nExpress = .N), by = "group"]
  ggData <- ggData[, .(nCells = .N), by = "group"]
  ggData <- ggData1[ggData, on = "group"]
  ggData <- ggData[, c("group", "nCells", "nExpress"), with = FALSE]
  ggData[is.na(nExpress)]$nExpress <- 0
  ggData$pctExpress <- 100 * ggData$nExpress / ggData$nCells
  ggData <- ggData[order(group)]
  colnames(ggData)[3] <- paste0(colnames(ggData)[3], "_", inp2)

  return(ggData)
}

# @description Plot gene expression on dimred
# @param inpConf (data.frame) Configuration table
# @param inpMeta (data.frame) Metadata table
# @param inpdrX (Character) X axis variable for DR
# @param inpdrY (Character) Y axis variable for DR
# @param inp1 (Character) Gene name to use
# @param inpsub1 (Character) Name of metadata column for subsetting
# @param inpsub2 (Character/Vector) Levels under metadata column for subsetting
# @param inpH5 (Character) Path to gene expression h5 file (sc1gexpr.h5)
# @param inpGene (integer) Named integer vector of gene expression values (sc1gene.rds)
# @param inpsiz (Numeric) Point size
# @param inpcol (Character) Custom colour label
# @param inpord (Character) Custom plotting order
# @param inpfsz (Character) Custom font size
# @param inppasp (Character) Custom aspect ratio
# @param inptxt (Logical) Show XY labels
#
scDRgene <- function(inpConf, inpMeta, inpdrX, inpdrY, inp1, inpsub1, inpsub2, inpH5, inpGene, inpsiz, inpcol, inpord, inpfsz, inpasp, inptxt) {
  
  if (is.null(inpsub1)) { inpsub1 <- inpConf$UI[1] }

  # Prepare ggData
  ggData <- inpMeta[, c(
    inpConf[UI == inpdrX]$ID, inpConf[UI == inpdrY]$ID,
    inpConf[UI == inpsub1]$ID
  ),
  with = FALSE
  ]
  colnames(ggData) <- c("X", "Y", "sub")
  rat <- (max(ggData$X) - min(ggData$X)) / (max(ggData$Y) - min(ggData$Y))

  h5file <- H5File$new(inpH5, mode = "r")
  h5data <- h5file[["grp"]][["data"]]
  ggData$val <- h5data$read(args = list(inpGene[inp1], quote(expr = )))
  ggData[val < 0]$val <- 0
  h5file$close_all()
  bgCells <- FALSE

  if (length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)) { bgCells <- TRUE
    ggData2 <- ggData[!sub %in% inpsub2]
    ggData <- ggData[sub %in% inpsub2] }

  if (inpord == "Max") {
    ggData <- ggData[order(val)]
  } else if (inpord == "Min") { 
    ggData <- ggData[order(-val)] 
  } else if (inpord == "Random") { 
    ggData <- ggData[sample(nrow(ggData))] 
  }

  # Actual ggplot
  ggOut <- ggplot(ggData, aes(X, Y, color = val))
  if (bgCells) { 
    ggOut <- ggOut +
      geom_point(data = ggData2, color = "snow2", size = inpsiz, shape = 20)
  }

  ggOut <- ggOut +
    geom_point(size = inpsiz, shape = 20) + xlab(inpdrX) + ylab(inpdrY) +
    sctheme(base_size = sList[inpfsz], XYval = inptxt) +
    scale_color_gradientn(inp1, colours = cList[[inpcol]]) +
    guides(color = guide_colorbar(barwidth = 20))

  if (inpasp == "Square") { 
    ggOut <- ggOut + coord_fixed(ratio = rat) 
  } else if (inpasp == "Fixed") {
    ggOut <- ggOut + coord_fixed()
  }

  return(ggOut)
}

# Plot gene coexpression on dimred
bilinear <- function(x, y, xy, Q11, Q21, Q12, Q22) {
  oup <- (xy - x) * (xy - y) * Q11 + x * (xy - y) * Q21 + (xy - x) * y * Q12 + x * y * Q22
  oup <- oup / (xy * xy)
  return(oup)
}

# @description Gene Co-expression on dimred with legend
# @param inpConf (data.frame) Configuration table
# @param inpMeta (data.frame) Metadata table
# @param inpdrX (Character) X axis variable for DR
# @param inpdrY (Character) Y axis variable for DR
# @param inp1 (Character) Gene name to use
# @param inp2 (Character) Gene name to use
# @param inpsub1 (Character) Name of metadata column for subsetting
# @param inpsub2 (Character/Vector) Levels under metadata column for subsetting
# @param inpH5 (Character) Path to gene expression h5 file (sc1gexpr.h5)
# @param inpGene (integer) Named integer vector of gene expression values (sc1gene.rds)
# @param inpsiz (Numeric) Point size
# @param inpcol (Character) Custom colour label
# @param inpord (Character) Custom plotting order
# @param inpfsz (Character) Custom font size
# @param inppasp (Character) Custom aspect ratio
# @param inptxt (Logical) Show XY labels
#
scDRcoexFull <- function(inpConf, inpMeta, inpdrX, inpdrY, inp1, inp2, inpsub1,
                     inpsub2, inpH5, inpGene, inpsiz, inpcol, inpord, inpfsz,
                     inpasp, inptxt) {
  g1 <- scDRcoex(inpConf, inpMeta, inpdrX, inpdrY, inp1, inp2, inpsub1,
                       inpsub2, inpH5, inpGene, inpsiz, inpcol, inpord, inpfsz,
                       inpasp, inptxt)
  g2 <- scDRcoexLeg(inp1, inp2, inpcol, inpfsz)
  g <- wrap_plots(g1, g2) + plot_layout(ncol = 2, nrow = 1, widths = c(10, 2))
  return(g)
}

# @description Gene Co-expression on dimred
# @param inpConf (data.frame) Configuration table
# @param inpMeta (data.frame) Metadata table
# @param inpdrX (Character) X axis variable for DR
# @param inpdrY (Character) Y axis variable for DR
# @param inp1 (Character) Gene name to use
# @param inp2 (Character) Gene name to use
# @param inpsub1 (Character) Name of metadata column for subsetting
# @param inpsub2 (Character/Vector) Levels under metadata column for subsetting
# @param inpH5 (Character) Path to gene expression h5 file (sc1gexpr.h5)
# @param inpGene (integer) Named integer vector of gene expression values (sc1gene.rds)
# @param inpsiz (Numeric) Point size
# @param inpcol (Character) Custom colour label
# @param inpord (Character) Custom plotting order
# @param inpfsz (Character) Custom font size
# @param inppasp (Character) Custom aspect ratio
# @param inptxt (Logical) Show XY labels
#
scDRcoex <- function(inpConf, inpMeta, inpdrX, inpdrY, inp1, inp2, inpsub1,
                     inpsub2, inpH5, inpGene, inpsiz, inpcol, inpord, inpfsz,
                     inpasp, inptxt) {
  
  if (is.null(inpsub1)) { inpsub1 <- inpConf$UI[1] }

  # Prepare ggData
  ggData <- inpMeta[, c(
    inpConf[UI == inpdrX]$ID, inpConf[UI == inpdrY]$ID,
    inpConf[UI == inpsub1]$ID
  ),
  with = FALSE
  ]
  colnames(ggData) <- c("X", "Y", "sub")
  rat <- (max(ggData$X) - min(ggData$X)) / (max(ggData$Y) - min(ggData$Y))

  h5file <- H5File$new(inpH5, mode = "r")
  h5data <- h5file[["grp"]][["data"]]
  ggData$val1 <- h5data$read(args = list(inpGene[inp1], quote(expr = )))
  ggData[val1 < 0]$val1 <- 0
  ggData$val2 <- h5data$read(args = list(inpGene[inp2], quote(expr = )))
  ggData[val2 < 0]$val2 <- 0
  h5file$close_all()
  bgCells <- FALSE

  if (length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)) {
    bgCells <- TRUE
    ggData2 <- ggData[!sub %in% inpsub2]
    ggData <- ggData[sub %in% inpsub2]
  }

  # Generate coex color palette
  cInp <- strsplit(inpcol, "; ")[[1]]
  if (cInp[1] == "Red (Gene1)") {
    c10 <- c(255, 0, 0) 
  } else if (cInp[1] == "Orange (Gene1)") {
    c10 <- c(255, 140, 0)
  } else {
    c10 <- c(0, 255, 0)
  }

  if (cInp[2] == "Green (Gene2)") { c01 <- c(0, 255, 0) } else { c01 <- c(0, 0, 255) }

  c00 <- c(217, 217, 217)
  c11 <- c10 + c01
  nGrid <- 16
  nPad <- 2
  nTot <- nGrid + nPad * 2
  gg <- data.table(v1 = rep(0:nTot, nTot + 1), v2 = sort(rep(0:nTot, nTot + 1)))
  gg$vv1 <- gg$v1 - nPad
  gg[vv1 < 0]$vv1 <- 0
  gg[vv1 > nGrid]$vv1 <- nGrid
  gg$vv2 <- gg$v2 - nPad
  gg[vv2 < 0]$vv2 <- 0
  gg[vv2 > nGrid]$vv2 <- nGrid
  gg$cR <- bilinear(gg$vv1, gg$vv2, nGrid, c00[1], c10[1], c01[1], c11[1])
  gg$cG <- bilinear(gg$vv1, gg$vv2, nGrid, c00[2], c10[2], c01[2], c11[2])
  gg$cB <- bilinear(gg$vv1, gg$vv2, nGrid, c00[3], c10[3], c01[3], c11[3])
  gg$cMix <- rgb(gg$cR, gg$cG, gg$cB, maxColorValue = 255)
  gg <- gg[, c("v1", "v2", "cMix")]

  # Map colours
  ggData$v1 <- round(nTot * ggData$val1 / max(ggData$val1))
  ggData$v2 <- round(nTot * ggData$val2 / max(ggData$val2))
  ggData$v0 <- ggData$v1 + ggData$v2
  ggData <- gg[ggData, on = c("v1", "v2")]

  if (inpord == "Max") {
    ggData <- ggData[order(v0)]
  } else if (inpord == "Min") {
    ggData <- ggData[order(-v0)]
  } else if (inpord == "Random") {
    ggData <- ggData[sample(nrow(ggData))]
  }

  # Actual ggplot
  ggOut <- ggplot(ggData, aes(X, Y))
  if (bgCells) {
    ggOut <- ggOut +
    geom_point(data = ggData2, color = "snow2", size = inpsiz, shape = 20)
  }
  
  ggOut <- ggOut +
    geom_point(size = inpsiz, shape = 20, color = ggData$cMix) +
    xlab(inpdrX) + ylab(inpdrY) +
    sctheme(base_size = sList[inpfsz], XYval = inptxt) +
    scale_color_gradientn(inp1, colours = cList[[1]]) +
    guides(color = guide_colorbar(barwidth = 20))

  if (inpasp == "Square") {
    ggOut <- ggOut + coord_fixed(ratio = rat)
  } else if (inpasp == "Fixed") {
    ggOut <- ggOut + coord_fixed()
  }
  
  return(ggOut)
}
  
# Co-exp plot legend
# @param inp1
# @param inp2
# @param inpcol Colour
# @param inpfsz Font size
#
scDRcoexLeg <- function(inp1, inp2, inpcol, inpfsz) {
  
  # Generate coex color palette
  cInp <- strsplit(inpcol, "; ")[[1]]
  if (cInp[1] == "Red (Gene1)") {
    c10 <- c(255, 0, 0)
  } else if (cInp[1] == "Orange (Gene1)") { 
    c10 <- c(255, 140, 0) 
  } else { 
    c10 <- c(0, 255, 0) 
  }

  if (cInp[2] == "Green (Gene2)") { 
    c01 <- c(0, 255, 0) 
  } else { 
    c01 <- c(0, 0, 255) 
  }

  c00 <- c(217, 217, 217)
  c11 <- c10 + c01
  nGrid <- 16
  nPad <- 2
  nTot <- nGrid + nPad * 2
  gg <- data.table(v1 = rep(0:nTot, nTot + 1), v2 = sort(rep(0:nTot, nTot + 1)))
  gg$vv1 <- gg$v1 - nPad
  gg[vv1 < 0]$vv1 <- 0
  gg[vv1 > nGrid]$vv1 <- nGrid
  gg$vv2 <- gg$v2 - nPad
  gg[vv2 < 0]$vv2 <- 0
  gg[vv2 > nGrid]$vv2 <- nGrid
  gg$cR <- bilinear(gg$vv1, gg$vv2, nGrid, c00[1], c10[1], c01[1], c11[1])
  gg$cG <- bilinear(gg$vv1, gg$vv2, nGrid, c00[2], c10[2], c01[2], c11[2])
  gg$cB <- bilinear(gg$vv1, gg$vv2, nGrid, c00[3], c10[3], c01[3], c11[3])
  gg$cMix <- rgb(gg$cR, gg$cG, gg$cB, maxColorValue = 255)
  gg <- gg[, c("v1", "v2", "cMix")]

  # Actual ggplot
  ggOut <- ggplot(gg, aes(v1, v2)) +
    geom_tile(fill = gg$cMix) +
    xlab(inp1) +
    ylab(inp2) +
    coord_fixed(ratio = 1) +
    scale_x_continuous(breaks = c(0, nTot), label = c("low", "high")) +
    scale_y_continuous(breaks = c(0, nTot), label = c("low", "high")) +
    sctheme(base_size = sList[inpfsz], XYval = TRUE)
  return(ggOut)
}


scDRcoexNum <- function(inpConf, inpMeta, inp1, inp2, inpsub1, inpsub2, inpH5, inpGene) {
  
  if (is.null(inpsub1)) {
  inpsub1 <- inpConf$UI[1]
  }

  # Prepare ggData
  ggData <- inpMeta[, c(inpConf[UI == inpsub1]$ID), with = FALSE]
  colnames(ggData) <- c("sub")
  h5file <- H5File$new(inpH5, mode = "r")
  h5data <- h5file[["grp"]][["data"]]
  ggData$val1 <- h5data$read(args = list(inpGene[inp1], quote(expr = )))
  ggData[val1 < 0]$val1 <- 0
  ggData$val2 <- h5data$read(args = list(inpGene[inp2], quote(expr = )))
  ggData[val2 < 0]$val2 <- 0
  h5file$close_all()
  if (length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)) { ggData <- ggData[sub %in% inpsub2] }

  # Actual data.table
  ggData$express <- "none"
  ggData[val1 > 0]$express <- inp1
  ggData[val2 > 0]$express <- inp2
  ggData[val1 > 0 & val2 > 0]$express <- "both"
  ggData$express <- factor(ggData$express, levels = unique(c("both", inp1, inp2, "none")))
  ggData <- ggData[, .(nCells = .N), by = "express"]
  ggData$percent <- 100 * ggData$nCells / sum(ggData$nCells)
  ggData <- ggData[order(express)]
  colnames(ggData)[1] <- "expression > 0"
  
  return(ggData)
}

# Plot violin / boxplot / lineplot
# @description Violin plot gene expression
# @param inpConf (data.frame) Configuration table
# @param inpMeta (data.frame) Metadata table
# @param inp1 (Character) X axis cell info
# @param inp2 (Character) Y axis cell info / gene
# @param inpsub1 (Character) Name of metadata column for subsetting
# @param inpsub2 (Character/Vector) Levels under metadata column for subsetting
# @param inpH5 (Character) Path to gene expression h5 file (sc1gexpr.h5)
# @param inpGene (integer) Named integer vector of gene expression values (sc1gene.rds)
# @param inptyp (Character) Plot type. "violin", "boxplot" or "lineplot".
# @param inppts (Logical) Should points be displayed?
# @param inpsiz (Numeric) Point size
# @param inpfsz (Character) Custom font size
# @param inpbarsz (Numeric) Bar size for lineplot
#
scVioBox <- function(inpConf, inpMeta, inp1, inp2, inpsub1, inpsub2, inpH5, inpGene, inptyp, inppts, inpsiz, inpfsz, inpbarsz) {
  
  if (is.null(inpsub1)) { inpsub1 <- inpConf$UI[1] }
  if (is.null(inpbarsz)) inpbarsz <- 0.3

  # Prepare ggData
  ggData <- inpMeta[, c(inpConf[UI == inp1]$ID, inpConf[UI == inpsub1]$ID),
    with = FALSE
  ]
  colnames(ggData) <- c("X", "sub")

  # Load in either cell meta or gene expr
  if (inp2 %in% inpConf$UI) {
    ggData$val <- inpMeta[[inpConf[UI == inp2]$ID]]
  } else {
    h5file <- H5File$new(inpH5, mode = "r")
    h5data <- h5file[["grp"]][["data"]]
    ggData$val <- h5data$read(args = list(inpGene[inp2], quote(expr = )))
    ggData[val < 0]$val <- 0
    set.seed(42)
    tmpNoise <- rnorm(length(ggData$val)) * diff(range(ggData$val)) / 1000
    ggData$val <- ggData$val + tmpNoise
    h5file$close_all()
  }

  if (length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)) {
    ggData <- ggData[sub %in% inpsub2]
  }

  # Do factoring
  ggCol <- strsplit(inpConf[UI == inp1]$fCL, "\\|")[[1]]
  names(ggCol) <- levels(ggData$X)
  ggLvl <- levels(ggData$X)[levels(ggData$X) %in% unique(ggData$X)]
  ggData$X <- factor(ggData$X, levels = ggLvl)
  ggCol <- ggCol[ggLvl]

  # Actual ggplot
  ggOut <- ggplot(ggData, aes(X, val, fill = X))
  
  if (inptyp == "violin") {
    ggOut <- ggOut + geom_violin(scale = "width")
  } else if (inptyp == "boxplot")  {
    ggOut <- ggOut + geom_boxplot()
  } else if (inptyp == "lineplot") {
    ggOut <- ggOut + geom_col(aes(col = X), position = position_dodge2(preserve = "single"), size = inpbarsz)
  }
  
  if (inppts) ggOut <- ggOut + geom_jitter(size = inpsiz, shape = 20, alpha = 0.4)

  ggOut <- ggOut + xlab(inp1) + ylab(inp2) +
    sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +
    scale_colour_manual("", values = ggCol) +
    scale_fill_manual("", values = ggCol) +
    theme(legend.position = "none")
    
  return(ggOut)
  }

# Plot proportion plot
# @description Proportion barplot
# @param inpConf (data.frame) Configuration table
# @param inpMeta (data.frame) Metadata table
# @param inp1 (Character) Gene name to use
# @param inp2 (Character) Gene name to use
# @param inpsub1 (Character) Name of metadata column for subsetting
# @param inpsub2 (Character/Vector) Levels under metadata column for subsetting
# @param inptyp (Character) Plot type. "violin" else boxplot.
# @param inpflp (Logical) Flip coordinates?
# @param inpfsz (Character) Custom font size
#
scProp <- function(inpConf, inpMeta, inp1, inp2, inpsub1, inpsub2, inptyp, inpflp, inpfsz) {
  
  if (is.null(inpsub1)) { inpsub1 <- inpConf$UI[1] }

  # Prepare ggData
  ggData <- inpMeta[, c(
    inpConf[UI == inp1]$ID, inpConf[UI == inp2]$ID,
    inpConf[UI == inpsub1]$ID
  ),
  with = FALSE
  ]
  colnames(ggData) <- c("X", "grp", "sub")

  if (length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)) { ggData <- ggData[sub %in% inpsub2] }
  ggData <- ggData[, .(nCells = .N), by = c("X", "grp")]
  ggData <- ggData[, { tot <- sum(nCells)
    .SD[, .(
      pctCells = 100 * sum(nCells) / tot,
      nCells = nCells
    ), by = "grp"] }, by = "X"]

  # Do factoring
  ggCol <- strsplit(inpConf[UI == inp2]$fCL, "\\|")[[1]]
  names(ggCol) <- levels(ggData$grp)
  ggLvl <- levels(ggData$grp)[levels(ggData$grp) %in% unique(ggData$grp)]
  ggData$grp <- factor(ggData$grp, levels = ggLvl)
  ggCol <- ggCol[ggLvl]

  # Actual ggplot
  if (inptyp == "Proportion") {
    ggOut <- ggplot(ggData, aes(X, pctCells, fill = grp)) +
      geom_col() +
      ylab("Cell Proportion (%)") 
  } else { 
    ggOut <- ggplot(ggData, aes(X, nCells, fill = grp)) +
      geom_col() +
      ylab("Number of Cells")
  }

  if (inpflp) ggOut <- ggOut + coord_flip()

  ggOut <- ggOut + xlab(inp1) +
    sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +
    scale_fill_manual("", values = ggCol) +
    theme(legend.position = "right")

  return(ggOut)
}

# Get gene list
scGeneList <- function(inp, inpGene) {
  geneList <- data.table(
  gene = inp,
  present = TRUE
  )
  geneList[!gene %in% names(inpGene)]$present <- FALSE
  return(geneList)
}

# Plot gene expression bubbleplot / heatmap
# @description dotplot / heatmap gene expression
# @param inpConf (data.frame) Configuration table
# @param inpMeta (data.frame) Metadata table
# @param inp (Character) Gene names
# @param inpGrp
# @param inpPlt
# @param inpsub1 (Character) Name of metadata column for subsetting
# @param inpsub2 (Character/Vector) Levels under metadata column for subsetting
# @param inpH5 (Character) Path to gene expression h5 file (sc1gexpr.h5)
# @param inpGene (integer) Named integer vector of gene expression values (sc1gene.rds)
# @param inpScl
# @param inpRow
# @param inpCol
# @param inpcols
# @param inpfsz (Character) Custom font size
# @param col_line (Character) Line colour
#
scBubbHeat <- function(inpConf, inpMeta, inp, inpGrp, inpPlt, inpsub1, inpsub2, inpH5, inpGene, inpScl, inpRow, inpCol, inpcols, inpfsz, col_line = "grey60") {
  
  if (is.null(inpsub1)) { inpsub1 <- inpConf$UI[1] }

  # Identify genes that are in our dataset
  geneList <- scGeneList(inp, inpGene)
  geneList <- geneList[present == TRUE]
  shiny::validate(need(nrow(geneList) <= 50, "More than 50 genes to plot! Please reduce the gene list!"))
  shiny::validate(need(nrow(geneList) > 1, "Please input at least 2 genes to plot!"))

  # Prepare ggData
  h5file <- H5File$new(inpH5, mode = "r")
  h5data <- h5file[["grp"]][["data"]]
  ggData <- data.table()
  for (iGene in geneList$gene) {
    tmp <- inpMeta[, c("sampleID", inpConf[UI == inpsub1]$ID), with = FALSE]
    colnames(tmp) <- c("sampleID", "sub")
    tmp$grpBy <- inpMeta[[inpConf[UI == inpGrp]$ID]]
    tmp$geneName <- iGene
    tmp$val <- h5data$read(args = list(inpGene[iGene], quote(expr = )))
    ggData <- rbindlist(list(ggData, tmp)) 
  }
  h5file$close_all()
  
  if (length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)) {
    ggData <- ggData[sub %in% inpsub2] 
  }
  shiny::validate(need(uniqueN(ggData$grpBy) > 1, "Only 1 group present, unable to plot!"))

  # Aggregate
  ggData$val <- expm1(ggData$val)
  ggData <- ggData[, .(val = mean(val), prop = sum(val > 0) / length(sampleID)),
    by = c("geneName", "grpBy")
  ]
  ggData$val <- log1p(ggData$val)

  # Scale if required
  colRange <- range(ggData$val)
  if (inpScl) {
    ggData[, val := scale(val), keyby = "geneName"]
    colRange <- c(-max(abs(range(ggData$val))), max(abs(range(ggData$val)))) 
  }

  # hclust row/col if necessary
  ggMat <- dcast.data.table(ggData, geneName ~ grpBy, value.var = "val")
  tmp <- ggMat$geneName
  ggMat <- as.matrix(ggMat[, -1])
  rownames(ggMat) <- tmp
  if (inpRow) {
    hcRow <- dendro_data(as.dendrogram(hclust(dist(ggMat))))
    ggRow <- ggplot() +
      coord_flip() +
      geom_segment(data = hcRow$segments, aes(x = x, y = y, xend = xend, yend = yend), col = col_line) +
      scale_y_continuous(
        breaks = rep(0, uniqueN(ggData$grpBy)),
        labels = unique(ggData$grpBy), expand = c(0, 0)
      ) +
      scale_x_continuous(
        breaks = seq_along(hcRow$labels$label),
        labels = hcRow$labels$label, expand = c(0, 0.5)
      ) +
      sctheme(base_size = sList[inpfsz]) +
      theme(
        axis.title = element_blank(), axis.line = element_blank(),
        axis.ticks = element_blank(), axis.text.y = element_blank(),
        axis.text.x = element_text(color = "white", angle = 45, hjust = 1)
      )
    ggData$geneName <- factor(ggData$geneName, levels = hcRow$labels$label) 
  } else { 
    ggData$geneName <- factor(ggData$geneName, levels = rev(geneList$gene)) 
  }

  if (inpCol) {
    hcCol <- dendro_data(as.dendrogram(hclust(dist(t(ggMat)))))
    ggCol <- ggplot() +
      geom_segment(data = hcCol$segments, aes(x = x, y = y, xend = xend, yend = yend), col = col_line) +
      scale_x_continuous(
        breaks = seq_along(hcCol$labels$label),
        labels = hcCol$labels$label, expand = c(0.05, 0)
      ) +
      scale_y_continuous(
        breaks = rep(0, uniqueN(ggData$geneName)),
        labels = unique(ggData$geneName), expand = c(0, 0)
      ) +
      sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +
      theme(
        axis.title = element_blank(), axis.line = element_blank(),
        axis.ticks = element_blank(), axis.text.x = element_blank(),
        axis.text.y = element_text(color = "white")
      )
    ggData$grpBy <- factor(ggData$grpBy, levels = hcCol$labels$label)
    }

  # Actual plot according to plottype
  if (inpPlt == "Bubbleplot") {
    # Bubbleplot
    ggOut <- ggplot(ggData, aes(grpBy, geneName, color = val, size = prop)) +
      geom_point() +
      sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +
      scale_x_discrete(expand = c(0.05, 0)) +
      scale_y_discrete(expand = c(0, 0.5)) +
      scale_size_continuous("proportion",
        range = c(0, 8),
        limits = c(0, 1), breaks = c(0.00, 0.25, 0.50, 0.75, 1.00)
      ) +
      scale_color_gradientn("expression", limits = colRange, colours = cList[[inpcols]]) +
      guides(color = guide_colorbar(barwidth = 20)) +
      theme(axis.title = element_blank(), legend.box = "vertical") 
  } else {
    # Heatmap
    ggOut <- ggplot(ggData, aes(grpBy, geneName, fill = val)) +
      geom_tile() +
      sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +
      scale_x_discrete(expand = c(0.05, 0)) +
      scale_y_discrete(expand = c(0, 0.5)) +
      scale_fill_gradientn("expression", limits = colRange, colours = cList[[inpcols]]) +
      guides(fill = guide_colorbar(barwidth = 20)) +
      theme(axis.title = element_blank()) 
  }

  # Final tidy
  ggLeg <- g_legend(ggOut)
  ggOut <- ggOut + theme(legend.position = "none")
  
  if (inpRow & inpCol) {
    ggOut <- wrap_plots(ggCol, plot_spacer(), ggOut, ggRow, ggplotify::as.ggplot(ggLeg))+
              plot_layout(ncol = 2, widths = c(7,1), heights = c(1,7,2))
  } else if (inpRow) {
    ggOut <- wrap_plots(ggOut, ggRow, ggplotify::as.ggplot(ggLeg))+
              plot_layout(ncol = 2, widths = c(7,1), heights = c(7,2))
  } else if (inpCol) {
    ggOut <- wrap_plots(ggCol, ggOut, ggplotify::as.ggplot(ggLeg))+
               plot_layout(ncol = 1, heights = c(1,7,2))
  } else {
    ggOut <- wrap_plots(ggOut, ggplotify::as.ggplot(ggLeg))+
      plot_layout(ncol = 1, heights = c(7,2))
  }
  
  return(ggOut)
}


### Server code ----
shinyServer(function(input, output, session) {

### For all tags and Server-side selectize
observe_helpers()
optCrt="{ option_create: function(data,escape) {return('<div class=\"create\"><strong>' + '</strong></div>');} }"
updateSelectizeInput(session, "lc_civge_inp2", choices = names(lcgene), server = TRUE,
                     selected = lcdef$gene1, options = list(
                       maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
updateSelectizeInput(session, "lc_gevge_inp1", choices = names(lcgene), server = TRUE,
                     selected = lcdef$gene1, options = list(
                       maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
updateSelectizeInput(session, "lc_gevge_inp2", choices = names(lcgene), server = TRUE,
                     selected = lcdef$gene2, options = list(
                       maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
updateSelectizeInput(session, "lc_gec_inp1", choices = names(lcgene), server = TRUE,
                     selected = lcdef$gene1, options = list(
                       maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
updateSelectizeInput(session, "lc_gec_inp2", choices = names(lcgene), server = TRUE,
                     selected = lcdef$gene2, options = list(
                       maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
updateSelectizeInput(session, "lc_gem_inp", choices = names(lcgene), server = TRUE,
                     selected = lcdef$genes[1:9], options = list(
                     create = TRUE, persist = TRUE, render = I(optCrt)))
updateSelectizeInput(session, "lc_hea_inp", choices = names(lcgene), server = TRUE,
                     selected = lcdef$genes[1:12], options = list(
                     create = TRUE, persist = TRUE, render = I(optCrt)))
updateSelectizeInput(session, "lc_vio_inp2", server = TRUE,
                     choices = c(lcconf[is.na(fID)]$UI,names(lcgene)),
                     selected = lcconf[is.na(fID)]$UI[1], options = list(
                       maxOptions = length(lcconf[is.na(fID)]$UI) + 3,
                       create = TRUE, persist = TRUE, render = I(optCrt)))  
### Tab civge cell info vs gene exp ----

  output$lc_civge_sub1.ui <- renderUI({
    sub = strsplit(lcconf[UI == input$lc_civge_sub1]$fID, "\\|")[[1]]
    checkboxGroupInput("lc_civge_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
  })
  observeEvent(input$lc_civge_sub1non, {
    sub = strsplit(lcconf[UI == input$lc_civge_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "lc_civge_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
  })
  observeEvent(input$lc_civge_sub1all, {
    sub = strsplit(lcconf[UI == input$lc_civge_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "lc_civge_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
  })

output$lc_civge_oup1 <- renderPlot({
  req(input$lc_civge_inp1)
  scDRcell(lcconf, lcmeta, input$lc_civge_drX, input$lc_civge_drY, input$lc_civge_inp1, input$lc_civge_sub1, input$lc_civge_sub2, input$lc_civge_siz, input$lc_civge_col1, input$lc_civge_ord1, input$lc_civge_fsz, input$lc_civge_asp, input$lc_civge_txt, input$lc_civge_lab1)
})

output$lc_civge_oup1.ui <- renderUI({
  show_progress(imageOutput("lc_civge_oup1", height = pList[input$lc_civge_psz]))
})

output$lc_civge_oup1.png <- downloadHandler(
 filename = function() { tolower(paste0("lc", "_", input$lc_civge_drX, "_", input$lc_civge_drY, "_", input$lc_civge_inp1, ".png")) },
 content = function(file) {
   ggsave(
   file, device = "png", dpi = input$lc_civge_oup1.res, bg = "white",
   plot = scDRcell(lcconf, lcmeta, input$lc_civge_drX, input$lc_civge_drY, input$lc_civge_inp1,   input$lc_civge_sub1, input$lc_civge_sub2, input$lc_civge_siz, input$lc_civge_col1, input$lc_civge_ord1,  input$lc_civge_fsz, input$lc_civge_asp, input$lc_civge_txt, input$lc_civge_lab1)
   )
})

output$lc_civge_oup1.pdf <- downloadHandler(
 filename = function() { tolower(paste0("lc", "_", input$lc_civge_drX,"_", input$lc_civge_drY,"_", input$lc_civge_inp1,".pdf")) },
 content = function(file) {
   ggsave(
   file, device = "pdf", useDingbats = FALSE, bg = "white", onefile = TRUE,
   plot = scDRcell(lcconf, lcmeta, input$lc_civge_drX, input$lc_civge_drY, input$lc_civge_inp1,   input$lc_civge_sub1, input$lc_civge_sub2, input$lc_civge_siz, input$lc_civge_col1, input$lc_civge_ord1,  input$lc_civge_fsz, input$lc_civge_asp, input$lc_civge_txt, input$lc_civge_lab1)
   )
})

output$lc_civge_oup1.svg <- downloadHandler(
 filename = function() { tolower(paste0("lc", "_", input$lc_civge_drX,"_", input$lc_civge_drY,"_", input$lc_civge_inp1,".svg")) },
 content = function(file) {
   ggsave(
   file, device = "svg", bg = "white",
   plot = scDRcell(lcconf, lcmeta, input$lc_civge_drX, input$lc_civge_drY, input$lc_civge_inp1,   input$lc_civge_sub1, input$lc_civge_sub2, input$lc_civge_siz, input$lc_civge_col1, input$lc_civge_ord1,  input$lc_civge_fsz, input$lc_civge_asp, input$lc_civge_txt, input$lc_civge_lab1)
   )
})

output$lc_civge_.dt <- renderDataTable({
 req(input$lc_civge_inp2)
 ggData = scDRnum(lcconf, lcmeta, input$lc_civge_inp1, input$lc_civge_inp2, input$lc_civge_sub1, input$lc_civge_sub2, "lcgexpr.h5", lcgene, input$lc_civge_splt)
 datatable(ggData, rownames = FALSE, extensions = "Buttons", options = list(pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel"))) %>%
   formatRound(columns = c("pctExpress"), digits = 2)
})

output$lc_civge_oup2 <- renderPlot({
 req(input$lc_civge_inp2)
 scDRgene(lcconf, lcmeta, input$lc_civge_drX, input$lc_civge_drY, input$lc_civge_inp2, input$lc_civge_sub1, input$lc_civge_sub2, "lcgexpr.h5", lcgene, input$lc_civge_siz, input$lc_civge_col2, input$lc_civge_ord2, input$lc_civge_fsz, input$lc_civge_asp, input$lc_civge_txt)
})

output$lc_civge_oup2.ui <- renderUI({
 show_progress(imageOutput("lc_civge_oup2", height = pList[input$lc_civge_psz]))
})

output$lc_civge_oup2.png <- downloadHandler(
 filename = function() { tolower(paste0("lc", "_", input$lc_civge_drX,"_",input$lc_civge_drY,"_", input$lc_civge_inp2,".png")) },
 content = function(file) {
   ggsave(
   file, device = "png", dpi = input$lc_civge_oup2.res, bg = "white",
   plot = scDRgene(lcconf, lcmeta, input$lc_civge_drX, input$lc_civge_drY, input$lc_civge_inp2, input$lc_civge_sub1, input$lc_civge_sub2, "lcgexpr.h5", lcgene, input$lc_civge_siz, input$lc_civge_col2, input$lc_civge_ord2, input$lc_civge_fsz, input$lc_civge_asp, input$lc_civge_txt)
   )
})

output$lc_civge_oup2.pdf <- downloadHandler(
 filename = function() { tolower(paste0("lc", "_", input$lc_civge_drX,"_",input$lc_civge_drY,"_", input$lc_civge_inp2,".pdf")) },
 content = function(file) {
   ggsave(
   file, device = "pdf", useDingbats = FALSE, bg = "white", onefile = TRUE,
   plot = scDRgene(lcconf, lcmeta, input$lc_civge_drX, input$lc_civge_drY, input$lc_civge_inp2,  input$lc_civge_sub1, input$lc_civge_sub2, "lcgexpr.h5", lcgene, input$lc_civge_siz, input$lc_civge_col2, input$lc_civge_ord2, input$lc_civge_fsz, input$lc_civge_asp, input$lc_civge_txt)
   )
}) 

output$lc_civge_oup2.svg <- downloadHandler(
 filename = function() { tolower(paste0("lc", "_", input$lc_civge_drX,"_",input$lc_civge_drY,"_", input$lc_civge_inp2,".svg")) },
 content = function(file) {
   ggsave(
   file, device = "svg", bg = "white",
   plot = scDRgene(lcconf, lcmeta, input$lc_civge_drX, input$lc_civge_drY, input$lc_civge_inp2,  input$lc_civge_sub1, input$lc_civge_sub2, "lcgexpr.h5", lcgene, input$lc_civge_siz, input$lc_civge_col2, input$lc_civge_ord2, input$lc_civge_fsz, input$lc_civge_asp, input$lc_civge_txt)
   )
}) # End of tab civge



### Tab civci cell info vs cell info ----
  
  output$lc_civci_sub1.ui <- renderUI({
    sub = strsplit(lcconf[UI == input$lc_civci_sub1]$fID, "\\|")[[1]]
    checkboxGroupInput("lc_civci_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
  })
  observeEvent(input$lc_civci_sub1non, {
    sub = strsplit(lcconf[UI == input$lc_civci_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "lc_civci_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
  })
  observeEvent(input$lc_civci_sub1all, {
    sub = strsplit(lcconf[UI == input$lc_civci_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "lc_civci_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
  })

output$lc_civci_oup1 <- renderPlot({
  req(input$lc_civci_inp1)
  scDRcell(lcconf, lcmeta, input$lc_civci_drX, input$lc_civci_drY, input$lc_civci_inp1, input$lc_civci_sub1, input$lc_civci_sub2, input$lc_civci_siz, input$lc_civci_col1, input$lc_civci_ord1, input$lc_civci_fsz, input$lc_civci_asp, input$lc_civci_txt, input$lc_civci_lab1)
})

output$lc_civci_oup1.ui <- renderUI({
  show_progress(imageOutput("lc_civci_oup1", height = pList[input$lc_civci_psz]))
})

output$lc_civci_oup1.png <- downloadHandler(
  filename = function() { tolower(paste0("lc", "_", input$lc_civci_drX, "_", input$lc_civci_drY, "_", input$lc_civci_inp1, ".png")) },
  content = function(file) {
    ggsave(
    file, device = "png", dpi = input$lc_civci_oup1.res, bg = "white",
    plot = scDRcell(lcconf, lcmeta, input$lc_civci_drX, input$lc_civci_drY, input$lc_civci_inp1, input$lc_civci_sub1, input$lc_civci_sub2, input$lc_civci_siz, input$lc_civci_col1, input$lc_civci_ord1, input$lc_civci_fsz, input$lc_civci_asp, input$lc_civci_txt, input$lc_civci_lab1)
    )
})

output$lc_civci_oup1.pdf <- downloadHandler(
  filename = function() { tolower(paste0("lc", "_", input$lc_civci_drX, "_", input$lc_civci_drY, "_", input$lc_civci_inp1, ".pdf")) },
  content = function(file) { ggsave(
    file, device = "pdf", useDingbats = FALSE, bg = "white", onefile = TRUE,
    plot = scDRcell(lcconf, lcmeta, input$lc_civci_drX, input$lc_civci_drY, input$lc_civci_inp1, input$lc_civci_sub1, input$lc_civci_sub2, input$lc_civci_siz, input$lc_civci_col1, input$lc_civci_ord1, input$lc_civci_fsz, input$lc_civci_asp, input$lc_civci_txt, input$lc_civci_lab1) )
})

output$lc_civci_oup1.svg <- downloadHandler(
  filename = function() { tolower(paste0("lc", "_", input$lc_civci_drX, "_", input$lc_civci_drY, "_", input$lc_civci_inp1, ".svg")) },
  content = function(file) { ggsave(
    file, device = "svg", bg = "white",
    plot = scDRcell(lcconf, lcmeta, input$lc_civci_drX, input$lc_civci_drY, input$lc_civci_inp1, input$lc_civci_sub1, input$lc_civci_sub2, input$lc_civci_siz, input$lc_civci_col1, input$lc_civci_ord1, input$lc_civci_fsz, input$lc_civci_asp, input$lc_civci_txt, input$lc_civci_lab1) )
})

output$lc_civci_oup2 <- renderPlot({
  req(input$lc_civci_inp2)
  scDRcell(lcconf, lcmeta, input$lc_civci_drX, input$lc_civci_drY, input$lc_civci_inp2, input$lc_civci_sub1, input$lc_civci_sub2, input$lc_civci_siz, input$lc_civci_col2, input$lc_civci_ord2, input$lc_civci_fsz, input$lc_civci_asp, input$lc_civci_txt, input$lc_civci_lab2)
})

output$lc_civci_oup2.ui <- renderUI({
  show_progress(imageOutput("lc_civci_oup2", height = pList[input$lc_civci_psz]))
})

output$lc_civci_oup2.png <- downloadHandler(
  filename = function() { tolower(paste0("lc", "_", input$lc_civci_drX,"_",input$lc_civci_drY,"_", input$lc_civci_inp2,".png")) },
  content = function(file) {
    ggsave(
    file, device = "png", bg = "white", dpi = input$lc_civci_oup2.res,
    plot = scDRcell(lcconf, lcmeta, input$lc_civci_drX, input$lc_civci_drY, input$lc_civci_inp2, input$lc_civci_sub1, input$lc_civci_sub2, input$lc_civci_siz, input$lc_civci_col2, input$lc_civci_ord2, input$lc_civci_fsz, input$lc_civci_asp, input$lc_civci_txt, input$lc_civci_lab2)
    )
}) 

output$lc_civci_oup2.pdf <- downloadHandler(
  filename = function() { tolower(paste0("lc", "_", input$lc_civci_drX,"_",input$lc_civci_drY,"_", input$lc_civci_inp2,".pdf")) },
  content = function(file) {
    ggsave(
    file, device = "pdf", useDingbats = FALSE, bg = "white", onefile = TRUE,
    plot = scDRcell(lcconf, lcmeta, input$lc_civci_drX, input$lc_civci_drY, input$lc_civci_inp2, input$lc_civci_sub1, input$lc_civci_sub2, input$lc_civci_siz, input$lc_civci_col2, input$lc_civci_ord2, input$lc_civci_fsz, input$lc_civci_asp, input$lc_civci_txt, input$lc_civci_lab2) 
    )
})

output$lc_civci_oup2.svg <- downloadHandler(
  filename = function() { tolower(paste0("lc", "_", input$lc_civci_drX,"_",input$lc_civci_drY,"_", input$lc_civci_inp2,".svg")) },
  content = function(file) {
    ggsave(
    file, device = "svg", bg = "white",
    plot = scDRcell(lcconf, lcmeta, input$lc_civci_drX, input$lc_civci_drY, input$lc_civci_inp2, input$lc_civci_sub1, input$lc_civci_sub2, input$lc_civci_siz, input$lc_civci_col2, input$lc_civci_ord2, input$lc_civci_fsz, input$lc_civci_asp, input$lc_civci_txt, input$lc_civci_lab2) 
    )
}) # End of tab civci



### Tab gevge gene exp vs gene exp ----

  output$lc_gevge_sub1.ui <- renderUI({
    sub = strsplit(lcconf[UI == input$lc_gevge_sub1]$fID, "\\|")[[1]]
    checkboxGroupInput("lc_gevge_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
  })
  observeEvent(input$lc_gevge_sub1non, {
    sub = strsplit(lcconf[UI == input$lc_gevge_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "lc_gevge_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
  })
  observeEvent(input$lc_gevge_sub1all, {
    sub = strsplit(lcconf[UI == input$lc_gevge_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "lc_gevge_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
  })

output$lc_gevge_oup1 <- renderPlot({
  req(input$lc_gevge_inp1)
  scDRgene(lcconf, lcmeta, input$lc_gevge_drX, input$lc_gevge_drY, input$lc_gevge_inp1, input$lc_gevge_sub1, input$lc_gevge_sub2, "lcgexpr.h5", lcgene, input$lc_gevge_siz, input$lc_gevge_col1, input$lc_gevge_ord1, input$lc_gevge_fsz, input$lc_gevge_asp, input$lc_gevge_txt)
})

output$lc_gevge_oup1.ui <- renderUI({
  show_progress(imageOutput("lc_gevge_oup1", height = pList[input$lc_gevge_psz]))
})

output$lc_gevge_oup1.png <- downloadHandler(
  filename = function() { tolower(paste0("lc", "_", input$lc_gevge_drX, "_", input$lc_gevge_drY, "_", input$lc_gevge_inp1,".png")) },
  content = function(file) {
    ggsave(
    file, device = "png", bg = "white", dpi = input$lc_gevge_oup1.res,
    plot = scDRgene(lcconf, lcmeta, input$lc_gevge_drX, input$lc_gevge_drY, input$lc_gevge_inp1, 
                    input$lc_gevge_sub1, input$lc_gevge_sub2,
                    "lcgexpr.h5", lcgene,
                    input$lc_gevge_siz, input$lc_gevge_col1, input$lc_gevge_ord1,
                    input$lc_gevge_fsz, input$lc_gevge_asp, input$lc_gevge_txt)
    )
})

output$lc_gevge_oup1.pdf <- downloadHandler(
  filename = function() { tolower(paste0("lc", "_", input$lc_gevge_drX, "_", input$lc_gevge_drY, "_", input$lc_gevge_inp1, ".pdf")) },
  content = function(file) {
    ggsave(
    file, device = "pdf", useDingbats = FALSE, bg = "white", onefile = TRUE,
    plot = scDRgene(lcconf, lcmeta, input$lc_gevge_drX, input$lc_gevge_drY, input$lc_gevge_inp1, input$lc_gevge_sub1, input$lc_gevge_sub2, "lcgexpr.h5", lcgene, input$lc_gevge_siz, input$lc_gevge_col1, input$lc_gevge_ord1, input$lc_gevge_fsz, input$lc_gevge_asp, input$lc_gevge_txt)
    )
})

output$lc_gevge_oup1.svg <- downloadHandler(
  filename = function() { tolower(paste0("lc", "_", input$lc_gevge_drX, "_", input$lc_gevge_drY, "_", input$lc_gevge_inp1, ".svg")) },
  content = function(file) {
    ggsave(
    file, device = "svg", bg = "white",
    plot = scDRgene(lcconf, lcmeta, input$lc_gevge_drX, input$lc_gevge_drY, input$lc_gevge_inp1, input$lc_gevge_sub1, input$lc_gevge_sub2, "lcgexpr.h5", lcgene, input$lc_gevge_siz, input$lc_gevge_col1, input$lc_gevge_ord1, input$lc_gevge_fsz, input$lc_gevge_asp, input$lc_gevge_txt)
    )
})

output$lc_gevge_oup2 <- renderPlot({
  req(input$lc_gevge_inp2)
  scDRgene(lcconf, lcmeta, input$lc_gevge_drX, input$lc_gevge_drY, input$lc_gevge_inp2, input$lc_gevge_sub1, input$lc_gevge_sub2, "lcgexpr.h5", lcgene, input$lc_gevge_siz, input$lc_gevge_col2, input$lc_gevge_ord2, input$lc_gevge_fsz, input$lc_gevge_asp, input$lc_gevge_txt)
})

output$lc_gevge_oup2.ui <- renderUI({
  show_progress(imageOutput("lc_gevge_oup2", height = pList[input$lc_gevge_psz]))
})

output$lc_gevge_oup2.png <- downloadHandler(
  filename = function() { tolower(paste0("lc", "_", input$lc_gevge_drX, "_", input$lc_gevge_drY, "_", input$lc_gevge_inp2,".png")) },
  content = function(file) {
    ggsave(
    file, device = "png", bg = "white", dpi = input$lc_gevge_oup2.res,
    plot = scDRgene(lcconf, lcmeta, input$lc_gevge_drX, input$lc_gevge_drY, input$lc_gevge_inp2, input$lc_gevge_sub1, input$lc_gevge_sub2, "lcgexpr.h5", lcgene, input$lc_gevge_siz, input$lc_gevge_col2, input$lc_gevge_ord2, input$lc_gevge_fsz, input$lc_gevge_asp, input$lc_gevge_txt)
    )
}) 

output$lc_gevge_oup2.pdf <- downloadHandler(
  filename = function() { tolower(paste0("lc", "_", input$lc_gevge_drX, "_", input$lc_gevge_drY, "_", input$lc_gevge_inp2,".pdf")) },
  content = function(file) {
    ggsave(
    file, device = "pdf", useDingbats = FALSE, bg = "white", onefile = TRUE,
    plot = scDRgene(lcconf, lcmeta, input$lc_gevge_drX, input$lc_gevge_drY, input$lc_gevge_inp2, 
                    input$lc_gevge_sub1, input$lc_gevge_sub2,
                    "lcgexpr.h5", lcgene,
                    input$lc_gevge_siz, input$lc_gevge_col2, input$lc_gevge_ord2,
                    input$lc_gevge_fsz, input$lc_gevge_asp, input$lc_gevge_txt)
    )
})

output$lc_gevge_oup2.svg <- downloadHandler(
  filename = function() { tolower(paste0("lc", "_", input$lc_gevge_drX, "_", input$lc_gevge_drY, "_", input$lc_gevge_inp2,".svg")) },
  content = function(file) {
    ggsave(
    file, device = "svg", bg = "white",
    plot = scDRgene(lcconf, lcmeta, input$lc_gevge_drX, input$lc_gevge_drY, input$lc_gevge_inp2, 
                    input$lc_gevge_sub1, input$lc_gevge_sub2,
                    "lcgexpr.h5", lcgene,
                    input$lc_gevge_siz, input$lc_gevge_col2, input$lc_gevge_ord2,
                    input$lc_gevge_fsz, input$lc_gevge_asp, input$lc_gevge_txt)
    )
}) # End of tab gevge




### Tab gem gene expression multi ----

output$lc_gem_sub1.ui <- renderUI({
  sub = strsplit(lcconf[UI == input$lc_gem_sub1]$fID, "\\|")[[1]]
  checkboxGroupInput("lc_gem_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
})
observeEvent(input$lc_gem_sub1non, {
  sub = strsplit(lcconf[UI == input$lc_gem_sub1]$fID, "\\|")[[1]]
  updateCheckboxGroupInput(session, inputId = "lc_gem_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
})
observeEvent(input$lc_gem_sub1all, {
  sub = strsplit(lcconf[UI == input$lc_gem_sub1]$fID, "\\|")[[1]]
  updateCheckboxGroupInput(session, inputId = "lc_gem_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
})

output$lc_gem_oup1 <- renderPlot({
  req(input$lc_gem_inp)
  
  scFeature(lcconf, lcmeta, input$lc_gem_drX, input$lc_gem_drY, input$lc_gem_inp, input$lc_gem_sub1, input$lc_gem_sub2, "lcgexpr.h5", lcgene, input$lc_gem_siz, input$lc_gem_col, input$lc_gem_ord, input$lc_gem_fsz, input$lc_gem_asp, input$lc_gem_txt, input$lc_gem_ncol)
})

output$lc_gem_oup1.ui <- renderUI({
  show_progress(imageOutput("lc_gem_oup1", height = pList[input$lc_gem_psz]))
})

output$lc_gem_oup1.png <- downloadHandler(
  filename = function() { tolower(paste0("lc", "_", input$lc_gem_drX, "_", input$lc_gem_drY, "_expression.png")) },
  content = function(file) {
    ggsave(
      file, device = "png", height = input$lc_gem_oup1.height, width = input$lc_gem_oup1.width, dpi = input$lc_gem_oup1.res, units = "cm", bg = "white",
      plot = scFeature(lcconf, lcmeta, input$lc_gem_drX, input$lc_gem_drY, input$lc_gem_inp, input$lc_gem_sub1, input$lc_gem_sub2, "lcgexpr.h5", lcgene, input$lc_gem_siz, input$lc_gem_col, input$lc_gem_ord, input$lc_gem_fsz, input$lc_gem_asp, input$lc_gem_txt, input$lc_gem_ncol)
    )
}) 

output$lc_gem_oup1.pdf <- downloadHandler(
  filename = function() { tolower(paste0("lc", "_", input$lc_gem_drX, "_", input$lc_gem_drY, "_expression.pdf")) },
  content = function(file) {
  pdf(file, useDingbats = FALSE, height = input$lc_gem_oup1.height/2.54, width = input$lc_gem_oup1.width/2.54, bg = "white", onefile = TRUE)
  showtext::showtext_begin()
  print(scFeature(lcconf, lcmeta, input$lc_gem_drX, input$lc_gem_drY, input$lc_gem_inp, input$lc_gem_sub1, input$lc_gem_sub2, "lcgexpr.h5", lcgene, input$lc_gem_siz, input$lc_gem_col, input$lc_gem_ord, input$lc_gem_fsz, input$lc_gem_asp, input$lc_gem_txt, input$lc_gem_ncol))
  showtext::showtext_end()
  dev.off()
})

output$lc_gem_oup1.svg <- downloadHandler(
  filename = function() { tolower(paste0("lc", "_", input$lc_gem_drX, "_", input$lc_gem_drY, "_expression.svg")) },
  content = function(file) { ggsave(
    file, device = "svg", height = input$lc_gem_oup1.height, width = input$lc_gem_oup1.width, units = "cm", bg = "white",
    plot = scFeature(lcconf, lcmeta, input$lc_gem_drX, input$lc_gem_drY, input$lc_gem_inp, input$lc_gem_sub1, input$lc_gem_sub2, "lcgexpr.h5", lcgene, input$lc_gem_siz, input$lc_gem_col, input$lc_gem_ord, input$lc_gem_fsz, input$lc_gem_asp, input$lc_gem_txt, input$lc_gem_ncol))
}) # End of tab gem



### Tab gec gene co-expression ----

  output$lc_gec_sub1.ui <- renderUI({
    sub = strsplit(lcconf[UI == input$lc_gec_sub1]$fID, "\\|")[[1]]
    checkboxGroupInput("lc_gec_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
  })
  observeEvent(input$lc_gec_sub1non, {
    sub = strsplit(lcconf[UI == input$lc_gec_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "lc_gec_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
  })
  observeEvent(input$lc_gec_sub1all, {
    sub = strsplit(lcconf[UI == input$lc_gec_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "lc_gec_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
  })

output$lc_gec_oup1 <- renderPlot({
  scDRcoexFull(lcconf, lcmeta, input$lc_gec_drX, input$lc_gec_drY, input$lc_gec_inp1, input$lc_gec_inp2, input$lc_gec_sub1, input$lc_gec_sub2, "lcgexpr.h5", lcgene, input$lc_gec_siz, input$lc_gec_col1, input$lc_gec_ord1, input$lc_gec_fsz, input$lc_gec_asp, input$lc_gec_txt)
})

output$lc_gec_oup1.ui <- renderUI({
  show_progress(imageOutput("lc_gec_oup1", height = pList2[input$lc_gec_psz]))
})

output$lc_gec_oup1.png <- downloadHandler(
  filename = function() { tolower(paste0("lc", "_", input$lc_gec_drX, "_", input$lc_gec_drY, "_", input$lc_gec_inp1, "_", input$lc_gec_inp2, ".png")) },
  content = function(file) { ggsave(
    file, device = "png", bg = "white", dpi = input$lc_gec_oup1.res,
    plot = scDRcoexFull(lcconf, lcmeta, input$lc_gec_drX, input$lc_gec_drY, input$lc_gec_inp1, input$lc_gec_inp2, input$lc_gec_sub1, input$lc_gec_sub2, "lcgexpr.h5", lcgene, input$lc_gec_siz, input$lc_gec_col1, input$lc_gec_ord1, input$lc_gec_fsz, input$lc_gec_asp, input$lc_gec_txt) )
})

output$lc_gec_oup1.pdf <- downloadHandler(
  filename = function() { tolower(paste0("lc", "_", input$lc_gec_drX, "_", input$lc_gec_drY, "_", input$lc_gec_inp1, "_", input$lc_gec_inp2, ".pdf")) },
  content = function(file) {
    ggsave(
    file, device = "pdf", useDingbats = FALSE, bg = "white", onefile = TRUE,
    plot = scDRcoexFull(lcconf, lcmeta, input$lc_gec_drX, input$lc_gec_drY, input$lc_gec_inp1, input$lc_gec_inp2, input$lc_gec_sub1, input$lc_gec_sub2, "lcgexpr.h5", lcgene, input$lc_gec_siz, input$lc_gec_col1, input$lc_gec_ord1, input$lc_gec_fsz, input$lc_gec_asp, input$lc_gec_txt)
    )
})

output$lc_gec_oup1.svg <- downloadHandler(
  filename = function() { tolower(paste0("lc", "_", input$lc_gec_drX, "_", input$lc_gec_drY, "_", input$lc_gec_inp1, "_", input$lc_gec_inp2, ".svg")) },
  content = function(file) {
    ggsave(
    file, device = "svg", bg = "white",
    plot = scDRcoexFull(lcconf, lcmeta, input$lc_gec_drX, input$lc_gec_drY, input$lc_gec_inp1, input$lc_gec_inp2, input$lc_gec_sub1, input$lc_gec_sub2, "lcgexpr.h5", lcgene, input$lc_gec_siz, input$lc_gec_col1, input$lc_gec_ord1, input$lc_gec_fsz, input$lc_gec_asp, input$lc_gec_txt)
    )
})

output$lc_gec_.dt <- renderDataTable({
  ggData = scDRcoexNum(lcconf, lcmeta, input$lc_gec_inp1, input$lc_gec_inp2, input$lc_gec_sub1, input$lc_gec_sub2, "lcgexpr.h5", lcgene)
  datatable(ggData, rownames = FALSE, extensions = "Buttons", options = list(pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel"))) %>%
            formatRound(columns = c("percent"), digits = 2)
}) # End of tab gec



### Tab vio violinplot / boxplot / lineplot ----

  output$lc_vio_sub1.ui <- renderUI({
    sub = strsplit(lcconf[UI == input$lc_vio_sub1]$fID, "\\|")[[1]]
    checkboxGroupInput("lc_vio_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
  })
  observeEvent(input$lc_vio_sub1non, {
    sub = strsplit(lcconf[UI == input$lc_vio_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "lc_vio_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
  })
  observeEvent(input$lc_vio_sub1all, {
    sub = strsplit(lcconf[UI == input$lc_vio_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "lc_vio_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
  })

output$lc_vio_oup <- renderPlot({
  gh5 <- ifelse(input$lc_vio_datatype == "normalised","lcgexpr.h5","lcgexpr2.h5")
  scVioBox(lcconf, lcmeta, input$lc_vio_inp1, input$lc_vio_inp2, input$lc_vio_sub1, input$lc_vio_sub2, gh5, lcgene, input$lc_vio_typ, input$lc_vio_pts, input$lc_vio_siz, input$lc_vio_fsz, input$lc_vio_barsz)
})

output$lc_vio_oup.ui <- renderUI({
  show_progress(imageOutput("lc_vio_oup", height = pList2[input$lc_vio_psz]))
})

output$lc_vio_oup.png <- downloadHandler(
  filename = function() { tolower(paste0("lc", "_", input$lc_vio_typ, "_", input$lc_vio_datatype, "_", input$lc_vio_inp1, "_", input$lc_vio_inp2,".png")) },
  content = function(file) {
    gh5 <- ifelse(input$lc_vio_datatype == "normalised","lcgexpr.h5","lcgexpr2.h5")
    ggsave(
    file, device = "png", bg = "white", dpi = input$lc_vio_oup.res,
    plot = scVioBox(lcconf, lcmeta, input$lc_vio_inp1, input$lc_vio_inp2, input$lc_vio_sub1, input$lc_vio_sub2, gh5, lcgene, input$lc_vio_typ, input$lc_vio_pts, input$lc_vio_siz, input$lc_vio_fsz, input$lc_vio_barsz)
    )
}) 

output$lc_vio_oup.pdf <- downloadHandler(
  filename = function() { tolower(paste0("lc", "_", input$lc_vio_typ, "_", input$lc_vio_datatype, "_", input$lc_vio_inp1, "_", input$lc_vio_inp2, ".pdf")) },
  content = function(file) {
    gh5 <- ifelse(input$lc_vio_datatype == "normalised","lcgexpr.h5","lcgexpr2.h5")
    ggsave(
    file, device = "pdf", useDingbats = FALSE, bg = "white", onefile = TRUE,
    plot = scVioBox(lcconf, lcmeta, input$lc_vio_inp1, input$lc_vio_inp2, input$lc_vio_sub1, input$lc_vio_sub2, gh5, lcgene, input$lc_vio_typ, input$lc_vio_pts, input$lc_vio_siz, input$lc_vio_fsz, input$lc_vio_barsz)
    )
})

output$lc_vio_oup.svg <- downloadHandler(
  filename = function() { tolower(paste0("lc", "_", input$lc_vio_typ, "_", input$lc_vio_datatype, "_", input$lc_vio_inp1, "_", input$lc_vio_inp2, ".svg")) },
  content = function(file) {
    gh5 <- ifelse(input$lc_vio_datatype == "normalised","lcgexpr.h5","lcgexpr2.h5")
    ggsave(
    file, device = "svg", bg = "white",
    plot = scVioBox(lcconf, lcmeta, input$lc_vio_inp1, input$lc_vio_inp2, input$lc_vio_sub1, input$lc_vio_sub2, gh5, lcgene, input$lc_vio_typ, input$lc_vio_pts, input$lc_vio_siz, input$lc_vio_fsz, input$lc_vio_barsz)
    )
}) # End of tab vio




### Tab pro proportion plot ----

  output$lc_pro_sub1.ui <- renderUI({
    sub = strsplit(lcconf[UI == input$lc_pro_sub1]$fID, "\\|")[[1]]
    checkboxGroupInput("lc_pro_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
  })
  observeEvent(input$lc_pro_sub1non, {
    sub = strsplit(lcconf[UI == input$lc_pro_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "lc_pro_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
  })
  observeEvent(input$lc_pro_sub1all, {
    sub = strsplit(lcconf[UI == input$lc_pro_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "lc_pro_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
  })

output$lc_pro_oup <- renderPlot({
  scProp(lcconf, lcmeta, input$lc_pro_inp1, input$lc_pro_inp2, input$lc_pro_sub1, input$lc_pro_sub2, input$lc_pro_typ, input$lc_pro_flp, input$lc_pro_fsz)
})

output$lc_pro_oup.ui <- renderUI({
  show_progress(imageOutput("lc_pro_oup", height = pList2[input$lc_pro_psz]))
})

output$lc_pro_oup.png <- downloadHandler(
  filename = function() { tolower(paste0("lc", "_", input$lc_pro_typ, "_", input$lc_pro_inp1, "_", input$lc_pro_inp2, ".png")) },
  content = function(file) {
    ggsave(
    file, device = "png", bg = "white", dpi = input$lc_pro_oup.res,
    plot = scProp(lcconf, lcmeta, input$lc_pro_inp1, input$lc_pro_inp2, input$lc_pro_sub1, input$lc_pro_sub2, input$lc_pro_typ, input$lc_pro_flp, input$lc_pro_fsz)
    )
  }) 
  
output$lc_pro_oup.pdf <- downloadHandler(
  filename = function() { tolower(paste0("lc", "_", input$lc_pro_typ, "_", input$lc_pro_inp1, "_", input$lc_pro_inp2, ".pdf")) },
  content = function(file) {
    ggsave(
    file, device = "pdf", useDingbats = FALSE, bg = "white", onefile = TRUE,
    plot = scProp(lcconf, lcmeta, input$lc_pro_inp1, input$lc_pro_inp2, input$lc_pro_sub1, input$lc_pro_sub2, input$lc_pro_typ, input$lc_pro_flp, input$lc_pro_fsz)
    )
  })
  
output$lc_pro_oup.svg <- downloadHandler(
  filename = function() { tolower(paste0("lc", "_", input$lc_pro_typ, "_", input$lc_pro_inp1, "_", input$lc_pro_inp2, ".svg")) },
  content = function(file) {
    ggsave(
    file, device = "svg", bg = "white",
    plot = scProp(lcconf, lcmeta, input$lc_pro_inp1, input$lc_pro_inp2, input$lc_pro_sub1, input$lc_pro_sub2, input$lc_pro_typ, input$lc_pro_flp, input$lc_pro_fsz)
    )
  }) # End of tab pro



### Tab hea heatmap / dotplot ----

  output$lc_hea_sub1.ui <- renderUI({
    sub = strsplit(lcconf[UI == input$lc_hea_sub1]$fID, "\\|")[[1]]
    checkboxGroupInput("lc_hea_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
  })
  observeEvent(input$lc_hea_sub1non, {
    sub = strsplit(lcconf[UI == input$lc_hea_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "lc_hea_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
  })
  observeEvent(input$lc_hea_sub1all, {
    sub = strsplit(lcconf[UI == input$lc_hea_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "lc_hea_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
  })

output$lc_hea_oupTxt <- renderUI({
  geneList = scGeneList(input$lc_hea_inp, lcgene)
  if(nrow(geneList) > 50){
    HTML("More than 50 input genes! Please reduce the gene list!")
  }
})

output$lc_hea_oup <- renderPlot({
  scBubbHeat(lcconf, lcmeta, input$lc_hea_inp, input$lc_hea_grp, input$lc_hea_plt, input$lc_hea_sub1, input$lc_hea_sub2, "lcgexpr.h5", lcgene, input$lc_hea_scl, input$lc_hea_row, input$lc_hea_col, input$lc_hea_cols, input$lc_hea_fsz)
})

output$lc_hea_oup.ui <- renderUI({
  show_progress(imageOutput("lc_hea_oup", height = pList3[input$lc_hea_psz]))
})

output$lc_hea_oup.png <- downloadHandler(
  filename = function() { tolower(paste0("lc", "_", input$lc_hea_plt,"_",input$lc_hea_grp,".png")) },
  content = function(file) {
    ggsave(
    file, device = "png", bg = "white", height = input$lc_hea_oup.height, width = input$lc_hea_oup.width, units = "cm", dpi = input$lc_hea_oup.res, plot = scBubbHeat(lcconf, lcmeta, input$lc_hea_inp, input$lc_hea_grp, input$lc_hea_plt, input$lc_hea_sub1, input$lc_hea_sub2, "lcgexpr.h5", lcgene, input$lc_hea_scl, input$lc_hea_row, input$lc_hea_col, input$lc_hea_cols, input$lc_hea_fsz)
    )
})

output$lc_hea_oup.pdf <- downloadHandler(
  filename = function() { tolower(paste0("lc", "_", input$lc_hea_plt,"_",input$lc_hea_grp,".pdf")) },
  content = function(file) {
    pdf(file, useDingbats = FALSE, bg = "white", height = input$lc_hea_oup.height/2.54, width = input$lc_hea_oup.width/2.54, onefile = TRUE)
    showtext::showtext_begin()
    print(scBubbHeat(lcconf, lcmeta, input$lc_hea_inp, input$lc_hea_grp, input$lc_hea_plt, input$lc_hea_sub1, input$lc_hea_sub2, "lcgexpr.h5", lcgene, input$lc_hea_scl, input$lc_hea_row, input$lc_hea_col, input$lc_hea_cols, input$lc_hea_fsz))
    showtext::showtext_end()
    dev.off()
})

output$lc_hea_oup.svg <- downloadHandler(
  filename = function() { tolower(paste0("lc", "_", input$lc_hea_plt,"_",input$lc_hea_grp,".svg")) },
  content = function(file) {
    ggsave(
    file, device = "svg", bg = "white", height = input$lc_hea_oup.height, width = input$lc_hea_oup.width, units = "cm", 
    plot = scBubbHeat(lcconf, lcmeta, input$lc_hea_inp, input$lc_hea_grp, input$lc_hea_plt, input$lc_hea_sub1, input$lc_hea_sub2, "lcgexpr.h5", lcgene, input$lc_hea_scl, input$lc_hea_row, input$lc_hea_col, input$lc_hea_cols, input$lc_hea_fsz)
    )
}) # End of tab hea      
       


### Tab markers ----

output$lc_mar_table <- renderDataTable({
  req(input$lc_mar_cls)
  datatable(lcmar[[input$lc_mar_cls]], rownames = FALSE, extensions = "Buttons", options = list(dom = "lftiprB", buttons = c("copy", "csv", "excel")))
}) # End of tab mar
optCrt="{ option_create: function(data,escape) {return('<div class=\"create\"><strong>' + '</strong></div>');} }"
updateSelectizeInput(session, "ld_civge_inp2", choices = names(ldgene), server = TRUE,
                     selected = lddef$gene1, options = list(
                       maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
updateSelectizeInput(session, "ld_gevge_inp1", choices = names(ldgene), server = TRUE,
                     selected = lddef$gene1, options = list(
                       maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
updateSelectizeInput(session, "ld_gevge_inp2", choices = names(ldgene), server = TRUE,
                     selected = lddef$gene2, options = list(
                       maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
updateSelectizeInput(session, "ld_gec_inp1", choices = names(ldgene), server = TRUE,
                     selected = lddef$gene1, options = list(
                       maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
updateSelectizeInput(session, "ld_gec_inp2", choices = names(ldgene), server = TRUE,
                     selected = lddef$gene2, options = list(
                       maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
updateSelectizeInput(session, "ld_gem_inp", choices = names(ldgene), server = TRUE,
                     selected = lddef$genes[1:9], options = list(
                     create = TRUE, persist = TRUE, render = I(optCrt)))
updateSelectizeInput(session, "ld_hea_inp", choices = names(ldgene), server = TRUE,
                     selected = lddef$genes[1:12], options = list(
                     create = TRUE, persist = TRUE, render = I(optCrt)))
updateSelectizeInput(session, "ld_vio_inp2", server = TRUE,
                     choices = c(ldconf[is.na(fID)]$UI,names(ldgene)),
                     selected = ldconf[is.na(fID)]$UI[1], options = list(
                       maxOptions = length(ldconf[is.na(fID)]$UI) + 3,
                       create = TRUE, persist = TRUE, render = I(optCrt)))  
### Tab civge cell info vs gene exp ----

  output$ld_civge_sub1.ui <- renderUI({
    sub = strsplit(ldconf[UI == input$ld_civge_sub1]$fID, "\\|")[[1]]
    checkboxGroupInput("ld_civge_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
  })
  observeEvent(input$ld_civge_sub1non, {
    sub = strsplit(ldconf[UI == input$ld_civge_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "ld_civge_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
  })
  observeEvent(input$ld_civge_sub1all, {
    sub = strsplit(ldconf[UI == input$ld_civge_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "ld_civge_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
  })

output$ld_civge_oup1 <- renderPlot({
  req(input$ld_civge_inp1)
  scDRcell(ldconf, ldmeta, input$ld_civge_drX, input$ld_civge_drY, input$ld_civge_inp1, input$ld_civge_sub1, input$ld_civge_sub2, input$ld_civge_siz, input$ld_civge_col1, input$ld_civge_ord1, input$ld_civge_fsz, input$ld_civge_asp, input$ld_civge_txt, input$ld_civge_lab1)
})

output$ld_civge_oup1.ui <- renderUI({
  show_progress(imageOutput("ld_civge_oup1", height = pList[input$ld_civge_psz]))
})

output$ld_civge_oup1.png <- downloadHandler(
 filename = function() { tolower(paste0("ld", "_", input$ld_civge_drX, "_", input$ld_civge_drY, "_", input$ld_civge_inp1, ".png")) },
 content = function(file) {
   ggsave(
   file, device = "png", dpi = input$ld_civge_oup1.res, bg = "white",
   plot = scDRcell(ldconf, ldmeta, input$ld_civge_drX, input$ld_civge_drY, input$ld_civge_inp1,   input$ld_civge_sub1, input$ld_civge_sub2, input$ld_civge_siz, input$ld_civge_col1, input$ld_civge_ord1,  input$ld_civge_fsz, input$ld_civge_asp, input$ld_civge_txt, input$ld_civge_lab1)
   )
})

output$ld_civge_oup1.pdf <- downloadHandler(
 filename = function() { tolower(paste0("ld", "_", input$ld_civge_drX,"_", input$ld_civge_drY,"_", input$ld_civge_inp1,".pdf")) },
 content = function(file) {
   ggsave(
   file, device = "pdf", useDingbats = FALSE, bg = "white", onefile = TRUE,
   plot = scDRcell(ldconf, ldmeta, input$ld_civge_drX, input$ld_civge_drY, input$ld_civge_inp1,   input$ld_civge_sub1, input$ld_civge_sub2, input$ld_civge_siz, input$ld_civge_col1, input$ld_civge_ord1,  input$ld_civge_fsz, input$ld_civge_asp, input$ld_civge_txt, input$ld_civge_lab1)
   )
})

output$ld_civge_oup1.svg <- downloadHandler(
 filename = function() { tolower(paste0("ld", "_", input$ld_civge_drX,"_", input$ld_civge_drY,"_", input$ld_civge_inp1,".svg")) },
 content = function(file) {
   ggsave(
   file, device = "svg", bg = "white",
   plot = scDRcell(ldconf, ldmeta, input$ld_civge_drX, input$ld_civge_drY, input$ld_civge_inp1,   input$ld_civge_sub1, input$ld_civge_sub2, input$ld_civge_siz, input$ld_civge_col1, input$ld_civge_ord1,  input$ld_civge_fsz, input$ld_civge_asp, input$ld_civge_txt, input$ld_civge_lab1)
   )
})

output$ld_civge_.dt <- renderDataTable({
 req(input$ld_civge_inp2)
 ggData = scDRnum(ldconf, ldmeta, input$ld_civge_inp1, input$ld_civge_inp2, input$ld_civge_sub1, input$ld_civge_sub2, "ldgexpr.h5", ldgene, input$ld_civge_splt)
 datatable(ggData, rownames = FALSE, extensions = "Buttons", options = list(pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel"))) %>%
   formatRound(columns = c("pctExpress"), digits = 2)
})

output$ld_civge_oup2 <- renderPlot({
 req(input$ld_civge_inp2)
 scDRgene(ldconf, ldmeta, input$ld_civge_drX, input$ld_civge_drY, input$ld_civge_inp2, input$ld_civge_sub1, input$ld_civge_sub2, "ldgexpr.h5", ldgene, input$ld_civge_siz, input$ld_civge_col2, input$ld_civge_ord2, input$ld_civge_fsz, input$ld_civge_asp, input$ld_civge_txt)
})

output$ld_civge_oup2.ui <- renderUI({
 show_progress(imageOutput("ld_civge_oup2", height = pList[input$ld_civge_psz]))
})

output$ld_civge_oup2.png <- downloadHandler(
 filename = function() { tolower(paste0("ld", "_", input$ld_civge_drX,"_",input$ld_civge_drY,"_", input$ld_civge_inp2,".png")) },
 content = function(file) {
   ggsave(
   file, device = "png", dpi = input$ld_civge_oup2.res, bg = "white",
   plot = scDRgene(ldconf, ldmeta, input$ld_civge_drX, input$ld_civge_drY, input$ld_civge_inp2, input$ld_civge_sub1, input$ld_civge_sub2, "ldgexpr.h5", ldgene, input$ld_civge_siz, input$ld_civge_col2, input$ld_civge_ord2, input$ld_civge_fsz, input$ld_civge_asp, input$ld_civge_txt)
   )
})

output$ld_civge_oup2.pdf <- downloadHandler(
 filename = function() { tolower(paste0("ld", "_", input$ld_civge_drX,"_",input$ld_civge_drY,"_", input$ld_civge_inp2,".pdf")) },
 content = function(file) {
   ggsave(
   file, device = "pdf", useDingbats = FALSE, bg = "white", onefile = TRUE,
   plot = scDRgene(ldconf, ldmeta, input$ld_civge_drX, input$ld_civge_drY, input$ld_civge_inp2,  input$ld_civge_sub1, input$ld_civge_sub2, "ldgexpr.h5", ldgene, input$ld_civge_siz, input$ld_civge_col2, input$ld_civge_ord2, input$ld_civge_fsz, input$ld_civge_asp, input$ld_civge_txt)
   )
}) 

output$ld_civge_oup2.svg <- downloadHandler(
 filename = function() { tolower(paste0("ld", "_", input$ld_civge_drX,"_",input$ld_civge_drY,"_", input$ld_civge_inp2,".svg")) },
 content = function(file) {
   ggsave(
   file, device = "svg", bg = "white",
   plot = scDRgene(ldconf, ldmeta, input$ld_civge_drX, input$ld_civge_drY, input$ld_civge_inp2,  input$ld_civge_sub1, input$ld_civge_sub2, "ldgexpr.h5", ldgene, input$ld_civge_siz, input$ld_civge_col2, input$ld_civge_ord2, input$ld_civge_fsz, input$ld_civge_asp, input$ld_civge_txt)
   )
}) # End of tab civge



### Tab civci cell info vs cell info ----
  
  output$ld_civci_sub1.ui <- renderUI({
    sub = strsplit(ldconf[UI == input$ld_civci_sub1]$fID, "\\|")[[1]]
    checkboxGroupInput("ld_civci_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
  })
  observeEvent(input$ld_civci_sub1non, {
    sub = strsplit(ldconf[UI == input$ld_civci_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "ld_civci_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
  })
  observeEvent(input$ld_civci_sub1all, {
    sub = strsplit(ldconf[UI == input$ld_civci_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "ld_civci_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
  })

output$ld_civci_oup1 <- renderPlot({
  req(input$ld_civci_inp1)
  scDRcell(ldconf, ldmeta, input$ld_civci_drX, input$ld_civci_drY, input$ld_civci_inp1, input$ld_civci_sub1, input$ld_civci_sub2, input$ld_civci_siz, input$ld_civci_col1, input$ld_civci_ord1, input$ld_civci_fsz, input$ld_civci_asp, input$ld_civci_txt, input$ld_civci_lab1)
})

output$ld_civci_oup1.ui <- renderUI({
  show_progress(imageOutput("ld_civci_oup1", height = pList[input$ld_civci_psz]))
})

output$ld_civci_oup1.png <- downloadHandler(
  filename = function() { tolower(paste0("ld", "_", input$ld_civci_drX, "_", input$ld_civci_drY, "_", input$ld_civci_inp1, ".png")) },
  content = function(file) {
    ggsave(
    file, device = "png", dpi = input$ld_civci_oup1.res, bg = "white",
    plot = scDRcell(ldconf, ldmeta, input$ld_civci_drX, input$ld_civci_drY, input$ld_civci_inp1, input$ld_civci_sub1, input$ld_civci_sub2, input$ld_civci_siz, input$ld_civci_col1, input$ld_civci_ord1, input$ld_civci_fsz, input$ld_civci_asp, input$ld_civci_txt, input$ld_civci_lab1)
    )
})

output$ld_civci_oup1.pdf <- downloadHandler(
  filename = function() { tolower(paste0("ld", "_", input$ld_civci_drX, "_", input$ld_civci_drY, "_", input$ld_civci_inp1, ".pdf")) },
  content = function(file) { ggsave(
    file, device = "pdf", useDingbats = FALSE, bg = "white", onefile = TRUE,
    plot = scDRcell(ldconf, ldmeta, input$ld_civci_drX, input$ld_civci_drY, input$ld_civci_inp1, input$ld_civci_sub1, input$ld_civci_sub2, input$ld_civci_siz, input$ld_civci_col1, input$ld_civci_ord1, input$ld_civci_fsz, input$ld_civci_asp, input$ld_civci_txt, input$ld_civci_lab1) )
})

output$ld_civci_oup1.svg <- downloadHandler(
  filename = function() { tolower(paste0("ld", "_", input$ld_civci_drX, "_", input$ld_civci_drY, "_", input$ld_civci_inp1, ".svg")) },
  content = function(file) { ggsave(
    file, device = "svg", bg = "white",
    plot = scDRcell(ldconf, ldmeta, input$ld_civci_drX, input$ld_civci_drY, input$ld_civci_inp1, input$ld_civci_sub1, input$ld_civci_sub2, input$ld_civci_siz, input$ld_civci_col1, input$ld_civci_ord1, input$ld_civci_fsz, input$ld_civci_asp, input$ld_civci_txt, input$ld_civci_lab1) )
})

output$ld_civci_oup2 <- renderPlot({
  req(input$ld_civci_inp2)
  scDRcell(ldconf, ldmeta, input$ld_civci_drX, input$ld_civci_drY, input$ld_civci_inp2, input$ld_civci_sub1, input$ld_civci_sub2, input$ld_civci_siz, input$ld_civci_col2, input$ld_civci_ord2, input$ld_civci_fsz, input$ld_civci_asp, input$ld_civci_txt, input$ld_civci_lab2)
})

output$ld_civci_oup2.ui <- renderUI({
  show_progress(imageOutput("ld_civci_oup2", height = pList[input$ld_civci_psz]))
})

output$ld_civci_oup2.png <- downloadHandler(
  filename = function() { tolower(paste0("ld", "_", input$ld_civci_drX,"_",input$ld_civci_drY,"_", input$ld_civci_inp2,".png")) },
  content = function(file) {
    ggsave(
    file, device = "png", bg = "white", dpi = input$ld_civci_oup2.res,
    plot = scDRcell(ldconf, ldmeta, input$ld_civci_drX, input$ld_civci_drY, input$ld_civci_inp2, input$ld_civci_sub1, input$ld_civci_sub2, input$ld_civci_siz, input$ld_civci_col2, input$ld_civci_ord2, input$ld_civci_fsz, input$ld_civci_asp, input$ld_civci_txt, input$ld_civci_lab2)
    )
}) 

output$ld_civci_oup2.pdf <- downloadHandler(
  filename = function() { tolower(paste0("ld", "_", input$ld_civci_drX,"_",input$ld_civci_drY,"_", input$ld_civci_inp2,".pdf")) },
  content = function(file) {
    ggsave(
    file, device = "pdf", useDingbats = FALSE, bg = "white", onefile = TRUE,
    plot = scDRcell(ldconf, ldmeta, input$ld_civci_drX, input$ld_civci_drY, input$ld_civci_inp2, input$ld_civci_sub1, input$ld_civci_sub2, input$ld_civci_siz, input$ld_civci_col2, input$ld_civci_ord2, input$ld_civci_fsz, input$ld_civci_asp, input$ld_civci_txt, input$ld_civci_lab2) 
    )
})

output$ld_civci_oup2.svg <- downloadHandler(
  filename = function() { tolower(paste0("ld", "_", input$ld_civci_drX,"_",input$ld_civci_drY,"_", input$ld_civci_inp2,".svg")) },
  content = function(file) {
    ggsave(
    file, device = "svg", bg = "white",
    plot = scDRcell(ldconf, ldmeta, input$ld_civci_drX, input$ld_civci_drY, input$ld_civci_inp2, input$ld_civci_sub1, input$ld_civci_sub2, input$ld_civci_siz, input$ld_civci_col2, input$ld_civci_ord2, input$ld_civci_fsz, input$ld_civci_asp, input$ld_civci_txt, input$ld_civci_lab2) 
    )
}) # End of tab civci



### Tab gevge gene exp vs gene exp ----

  output$ld_gevge_sub1.ui <- renderUI({
    sub = strsplit(ldconf[UI == input$ld_gevge_sub1]$fID, "\\|")[[1]]
    checkboxGroupInput("ld_gevge_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
  })
  observeEvent(input$ld_gevge_sub1non, {
    sub = strsplit(ldconf[UI == input$ld_gevge_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "ld_gevge_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
  })
  observeEvent(input$ld_gevge_sub1all, {
    sub = strsplit(ldconf[UI == input$ld_gevge_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "ld_gevge_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
  })

output$ld_gevge_oup1 <- renderPlot({
  req(input$ld_gevge_inp1)
  scDRgene(ldconf, ldmeta, input$ld_gevge_drX, input$ld_gevge_drY, input$ld_gevge_inp1, input$ld_gevge_sub1, input$ld_gevge_sub2, "ldgexpr.h5", ldgene, input$ld_gevge_siz, input$ld_gevge_col1, input$ld_gevge_ord1, input$ld_gevge_fsz, input$ld_gevge_asp, input$ld_gevge_txt)
})

output$ld_gevge_oup1.ui <- renderUI({
  show_progress(imageOutput("ld_gevge_oup1", height = pList[input$ld_gevge_psz]))
})

output$ld_gevge_oup1.png <- downloadHandler(
  filename = function() { tolower(paste0("ld", "_", input$ld_gevge_drX, "_", input$ld_gevge_drY, "_", input$ld_gevge_inp1,".png")) },
  content = function(file) {
    ggsave(
    file, device = "png", bg = "white", dpi = input$ld_gevge_oup1.res,
    plot = scDRgene(ldconf, ldmeta, input$ld_gevge_drX, input$ld_gevge_drY, input$ld_gevge_inp1, 
                    input$ld_gevge_sub1, input$ld_gevge_sub2,
                    "ldgexpr.h5", ldgene,
                    input$ld_gevge_siz, input$ld_gevge_col1, input$ld_gevge_ord1,
                    input$ld_gevge_fsz, input$ld_gevge_asp, input$ld_gevge_txt)
    )
})

output$ld_gevge_oup1.pdf <- downloadHandler(
  filename = function() { tolower(paste0("ld", "_", input$ld_gevge_drX, "_", input$ld_gevge_drY, "_", input$ld_gevge_inp1, ".pdf")) },
  content = function(file) {
    ggsave(
    file, device = "pdf", useDingbats = FALSE, bg = "white", onefile = TRUE,
    plot = scDRgene(ldconf, ldmeta, input$ld_gevge_drX, input$ld_gevge_drY, input$ld_gevge_inp1, input$ld_gevge_sub1, input$ld_gevge_sub2, "ldgexpr.h5", ldgene, input$ld_gevge_siz, input$ld_gevge_col1, input$ld_gevge_ord1, input$ld_gevge_fsz, input$ld_gevge_asp, input$ld_gevge_txt)
    )
})

output$ld_gevge_oup1.svg <- downloadHandler(
  filename = function() { tolower(paste0("ld", "_", input$ld_gevge_drX, "_", input$ld_gevge_drY, "_", input$ld_gevge_inp1, ".svg")) },
  content = function(file) {
    ggsave(
    file, device = "svg", bg = "white",
    plot = scDRgene(ldconf, ldmeta, input$ld_gevge_drX, input$ld_gevge_drY, input$ld_gevge_inp1, input$ld_gevge_sub1, input$ld_gevge_sub2, "ldgexpr.h5", ldgene, input$ld_gevge_siz, input$ld_gevge_col1, input$ld_gevge_ord1, input$ld_gevge_fsz, input$ld_gevge_asp, input$ld_gevge_txt)
    )
})

output$ld_gevge_oup2 <- renderPlot({
  req(input$ld_gevge_inp2)
  scDRgene(ldconf, ldmeta, input$ld_gevge_drX, input$ld_gevge_drY, input$ld_gevge_inp2, input$ld_gevge_sub1, input$ld_gevge_sub2, "ldgexpr.h5", ldgene, input$ld_gevge_siz, input$ld_gevge_col2, input$ld_gevge_ord2, input$ld_gevge_fsz, input$ld_gevge_asp, input$ld_gevge_txt)
})

output$ld_gevge_oup2.ui <- renderUI({
  show_progress(imageOutput("ld_gevge_oup2", height = pList[input$ld_gevge_psz]))
})

output$ld_gevge_oup2.png <- downloadHandler(
  filename = function() { tolower(paste0("ld", "_", input$ld_gevge_drX, "_", input$ld_gevge_drY, "_", input$ld_gevge_inp2,".png")) },
  content = function(file) {
    ggsave(
    file, device = "png", bg = "white", dpi = input$ld_gevge_oup2.res,
    plot = scDRgene(ldconf, ldmeta, input$ld_gevge_drX, input$ld_gevge_drY, input$ld_gevge_inp2, input$ld_gevge_sub1, input$ld_gevge_sub2, "ldgexpr.h5", ldgene, input$ld_gevge_siz, input$ld_gevge_col2, input$ld_gevge_ord2, input$ld_gevge_fsz, input$ld_gevge_asp, input$ld_gevge_txt)
    )
}) 

output$ld_gevge_oup2.pdf <- downloadHandler(
  filename = function() { tolower(paste0("ld", "_", input$ld_gevge_drX, "_", input$ld_gevge_drY, "_", input$ld_gevge_inp2,".pdf")) },
  content = function(file) {
    ggsave(
    file, device = "pdf", useDingbats = FALSE, bg = "white", onefile = TRUE,
    plot = scDRgene(ldconf, ldmeta, input$ld_gevge_drX, input$ld_gevge_drY, input$ld_gevge_inp2, 
                    input$ld_gevge_sub1, input$ld_gevge_sub2,
                    "ldgexpr.h5", ldgene,
                    input$ld_gevge_siz, input$ld_gevge_col2, input$ld_gevge_ord2,
                    input$ld_gevge_fsz, input$ld_gevge_asp, input$ld_gevge_txt)
    )
})

output$ld_gevge_oup2.svg <- downloadHandler(
  filename = function() { tolower(paste0("ld", "_", input$ld_gevge_drX, "_", input$ld_gevge_drY, "_", input$ld_gevge_inp2,".svg")) },
  content = function(file) {
    ggsave(
    file, device = "svg", bg = "white",
    plot = scDRgene(ldconf, ldmeta, input$ld_gevge_drX, input$ld_gevge_drY, input$ld_gevge_inp2, 
                    input$ld_gevge_sub1, input$ld_gevge_sub2,
                    "ldgexpr.h5", ldgene,
                    input$ld_gevge_siz, input$ld_gevge_col2, input$ld_gevge_ord2,
                    input$ld_gevge_fsz, input$ld_gevge_asp, input$ld_gevge_txt)
    )
}) # End of tab gevge




### Tab gem gene expression multi ----

output$ld_gem_sub1.ui <- renderUI({
  sub = strsplit(ldconf[UI == input$ld_gem_sub1]$fID, "\\|")[[1]]
  checkboxGroupInput("ld_gem_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
})
observeEvent(input$ld_gem_sub1non, {
  sub = strsplit(ldconf[UI == input$ld_gem_sub1]$fID, "\\|")[[1]]
  updateCheckboxGroupInput(session, inputId = "ld_gem_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
})
observeEvent(input$ld_gem_sub1all, {
  sub = strsplit(ldconf[UI == input$ld_gem_sub1]$fID, "\\|")[[1]]
  updateCheckboxGroupInput(session, inputId = "ld_gem_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
})

output$ld_gem_oup1 <- renderPlot({
  req(input$ld_gem_inp)
  
  scFeature(ldconf, ldmeta, input$ld_gem_drX, input$ld_gem_drY, input$ld_gem_inp, input$ld_gem_sub1, input$ld_gem_sub2, "ldgexpr.h5", ldgene, input$ld_gem_siz, input$ld_gem_col, input$ld_gem_ord, input$ld_gem_fsz, input$ld_gem_asp, input$ld_gem_txt, input$ld_gem_ncol)
})

output$ld_gem_oup1.ui <- renderUI({
  show_progress(imageOutput("ld_gem_oup1", height = pList[input$ld_gem_psz]))
})

output$ld_gem_oup1.png <- downloadHandler(
  filename = function() { tolower(paste0("ld", "_", input$ld_gem_drX, "_", input$ld_gem_drY, "_expression.png")) },
  content = function(file) {
    ggsave(
      file, device = "png", height = input$ld_gem_oup1.height, width = input$ld_gem_oup1.width, dpi = input$ld_gem_oup1.res, units = "cm", bg = "white",
      plot = scFeature(ldconf, ldmeta, input$ld_gem_drX, input$ld_gem_drY, input$ld_gem_inp, input$ld_gem_sub1, input$ld_gem_sub2, "ldgexpr.h5", ldgene, input$ld_gem_siz, input$ld_gem_col, input$ld_gem_ord, input$ld_gem_fsz, input$ld_gem_asp, input$ld_gem_txt, input$ld_gem_ncol)
    )
}) 

output$ld_gem_oup1.pdf <- downloadHandler(
  filename = function() { tolower(paste0("ld", "_", input$ld_gem_drX, "_", input$ld_gem_drY, "_expression.pdf")) },
  content = function(file) {
  pdf(file, useDingbats = FALSE, height = input$ld_gem_oup1.height/2.54, width = input$ld_gem_oup1.width/2.54, bg = "white", onefile = TRUE)
  showtext::showtext_begin()
  print(scFeature(ldconf, ldmeta, input$ld_gem_drX, input$ld_gem_drY, input$ld_gem_inp, input$ld_gem_sub1, input$ld_gem_sub2, "ldgexpr.h5", ldgene, input$ld_gem_siz, input$ld_gem_col, input$ld_gem_ord, input$ld_gem_fsz, input$ld_gem_asp, input$ld_gem_txt, input$ld_gem_ncol))
  showtext::showtext_end()
  dev.off()
})

output$ld_gem_oup1.svg <- downloadHandler(
  filename = function() { tolower(paste0("ld", "_", input$ld_gem_drX, "_", input$ld_gem_drY, "_expression.svg")) },
  content = function(file) { ggsave(
    file, device = "svg", height = input$ld_gem_oup1.height, width = input$ld_gem_oup1.width, units = "cm", bg = "white",
    plot = scFeature(ldconf, ldmeta, input$ld_gem_drX, input$ld_gem_drY, input$ld_gem_inp, input$ld_gem_sub1, input$ld_gem_sub2, "ldgexpr.h5", ldgene, input$ld_gem_siz, input$ld_gem_col, input$ld_gem_ord, input$ld_gem_fsz, input$ld_gem_asp, input$ld_gem_txt, input$ld_gem_ncol))
}) # End of tab gem



### Tab gec gene co-expression ----

  output$ld_gec_sub1.ui <- renderUI({
    sub = strsplit(ldconf[UI == input$ld_gec_sub1]$fID, "\\|")[[1]]
    checkboxGroupInput("ld_gec_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
  })
  observeEvent(input$ld_gec_sub1non, {
    sub = strsplit(ldconf[UI == input$ld_gec_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "ld_gec_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
  })
  observeEvent(input$ld_gec_sub1all, {
    sub = strsplit(ldconf[UI == input$ld_gec_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "ld_gec_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
  })

output$ld_gec_oup1 <- renderPlot({
  scDRcoexFull(ldconf, ldmeta, input$ld_gec_drX, input$ld_gec_drY, input$ld_gec_inp1, input$ld_gec_inp2, input$ld_gec_sub1, input$ld_gec_sub2, "ldgexpr.h5", ldgene, input$ld_gec_siz, input$ld_gec_col1, input$ld_gec_ord1, input$ld_gec_fsz, input$ld_gec_asp, input$ld_gec_txt)
})

output$ld_gec_oup1.ui <- renderUI({
  show_progress(imageOutput("ld_gec_oup1", height = pList2[input$ld_gec_psz]))
})

output$ld_gec_oup1.png <- downloadHandler(
  filename = function() { tolower(paste0("ld", "_", input$ld_gec_drX, "_", input$ld_gec_drY, "_", input$ld_gec_inp1, "_", input$ld_gec_inp2, ".png")) },
  content = function(file) { ggsave(
    file, device = "png", bg = "white", dpi = input$ld_gec_oup1.res,
    plot = scDRcoexFull(ldconf, ldmeta, input$ld_gec_drX, input$ld_gec_drY, input$ld_gec_inp1, input$ld_gec_inp2, input$ld_gec_sub1, input$ld_gec_sub2, "ldgexpr.h5", ldgene, input$ld_gec_siz, input$ld_gec_col1, input$ld_gec_ord1, input$ld_gec_fsz, input$ld_gec_asp, input$ld_gec_txt) )
})

output$ld_gec_oup1.pdf <- downloadHandler(
  filename = function() { tolower(paste0("ld", "_", input$ld_gec_drX, "_", input$ld_gec_drY, "_", input$ld_gec_inp1, "_", input$ld_gec_inp2, ".pdf")) },
  content = function(file) {
    ggsave(
    file, device = "pdf", useDingbats = FALSE, bg = "white", onefile = TRUE,
    plot = scDRcoexFull(ldconf, ldmeta, input$ld_gec_drX, input$ld_gec_drY, input$ld_gec_inp1, input$ld_gec_inp2, input$ld_gec_sub1, input$ld_gec_sub2, "ldgexpr.h5", ldgene, input$ld_gec_siz, input$ld_gec_col1, input$ld_gec_ord1, input$ld_gec_fsz, input$ld_gec_asp, input$ld_gec_txt)
    )
})

output$ld_gec_oup1.svg <- downloadHandler(
  filename = function() { tolower(paste0("ld", "_", input$ld_gec_drX, "_", input$ld_gec_drY, "_", input$ld_gec_inp1, "_", input$ld_gec_inp2, ".svg")) },
  content = function(file) {
    ggsave(
    file, device = "svg", bg = "white",
    plot = scDRcoexFull(ldconf, ldmeta, input$ld_gec_drX, input$ld_gec_drY, input$ld_gec_inp1, input$ld_gec_inp2, input$ld_gec_sub1, input$ld_gec_sub2, "ldgexpr.h5", ldgene, input$ld_gec_siz, input$ld_gec_col1, input$ld_gec_ord1, input$ld_gec_fsz, input$ld_gec_asp, input$ld_gec_txt)
    )
})

output$ld_gec_.dt <- renderDataTable({
  ggData = scDRcoexNum(ldconf, ldmeta, input$ld_gec_inp1, input$ld_gec_inp2, input$ld_gec_sub1, input$ld_gec_sub2, "ldgexpr.h5", ldgene)
  datatable(ggData, rownames = FALSE, extensions = "Buttons", options = list(pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel"))) %>%
            formatRound(columns = c("percent"), digits = 2)
}) # End of tab gec



### Tab vio violinplot / boxplot / lineplot ----

  output$ld_vio_sub1.ui <- renderUI({
    sub = strsplit(ldconf[UI == input$ld_vio_sub1]$fID, "\\|")[[1]]
    checkboxGroupInput("ld_vio_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
  })
  observeEvent(input$ld_vio_sub1non, {
    sub = strsplit(ldconf[UI == input$ld_vio_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "ld_vio_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
  })
  observeEvent(input$ld_vio_sub1all, {
    sub = strsplit(ldconf[UI == input$ld_vio_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "ld_vio_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
  })

output$ld_vio_oup <- renderPlot({
  gh5 <- ifelse(input$ld_vio_datatype == "normalised","ldgexpr.h5","ldgexpr2.h5")
  scVioBox(ldconf, ldmeta, input$ld_vio_inp1, input$ld_vio_inp2, input$ld_vio_sub1, input$ld_vio_sub2, gh5, ldgene, input$ld_vio_typ, input$ld_vio_pts, input$ld_vio_siz, input$ld_vio_fsz, input$ld_vio_barsz)
})

output$ld_vio_oup.ui <- renderUI({
  show_progress(imageOutput("ld_vio_oup", height = pList2[input$ld_vio_psz]))
})

output$ld_vio_oup.png <- downloadHandler(
  filename = function() { tolower(paste0("ld", "_", input$ld_vio_typ, "_", input$ld_vio_datatype, "_", input$ld_vio_inp1, "_", input$ld_vio_inp2,".png")) },
  content = function(file) {
    gh5 <- ifelse(input$ld_vio_datatype == "normalised","ldgexpr.h5","ldgexpr2.h5")
    ggsave(
    file, device = "png", bg = "white", dpi = input$ld_vio_oup.res,
    plot = scVioBox(ldconf, ldmeta, input$ld_vio_inp1, input$ld_vio_inp2, input$ld_vio_sub1, input$ld_vio_sub2, gh5, ldgene, input$ld_vio_typ, input$ld_vio_pts, input$ld_vio_siz, input$ld_vio_fsz, input$ld_vio_barsz)
    )
}) 

output$ld_vio_oup.pdf <- downloadHandler(
  filename = function() { tolower(paste0("ld", "_", input$ld_vio_typ, "_", input$ld_vio_datatype, "_", input$ld_vio_inp1, "_", input$ld_vio_inp2, ".pdf")) },
  content = function(file) {
    gh5 <- ifelse(input$ld_vio_datatype == "normalised","ldgexpr.h5","ldgexpr2.h5")
    ggsave(
    file, device = "pdf", useDingbats = FALSE, bg = "white", onefile = TRUE,
    plot = scVioBox(ldconf, ldmeta, input$ld_vio_inp1, input$ld_vio_inp2, input$ld_vio_sub1, input$ld_vio_sub2, gh5, ldgene, input$ld_vio_typ, input$ld_vio_pts, input$ld_vio_siz, input$ld_vio_fsz, input$ld_vio_barsz)
    )
})

output$ld_vio_oup.svg <- downloadHandler(
  filename = function() { tolower(paste0("ld", "_", input$ld_vio_typ, "_", input$ld_vio_datatype, "_", input$ld_vio_inp1, "_", input$ld_vio_inp2, ".svg")) },
  content = function(file) {
    gh5 <- ifelse(input$ld_vio_datatype == "normalised","ldgexpr.h5","ldgexpr2.h5")
    ggsave(
    file, device = "svg", bg = "white",
    plot = scVioBox(ldconf, ldmeta, input$ld_vio_inp1, input$ld_vio_inp2, input$ld_vio_sub1, input$ld_vio_sub2, gh5, ldgene, input$ld_vio_typ, input$ld_vio_pts, input$ld_vio_siz, input$ld_vio_fsz, input$ld_vio_barsz)
    )
}) # End of tab vio




### Tab pro proportion plot ----

  output$ld_pro_sub1.ui <- renderUI({
    sub = strsplit(ldconf[UI == input$ld_pro_sub1]$fID, "\\|")[[1]]
    checkboxGroupInput("ld_pro_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
  })
  observeEvent(input$ld_pro_sub1non, {
    sub = strsplit(ldconf[UI == input$ld_pro_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "ld_pro_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
  })
  observeEvent(input$ld_pro_sub1all, {
    sub = strsplit(ldconf[UI == input$ld_pro_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "ld_pro_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
  })

output$ld_pro_oup <- renderPlot({
  scProp(ldconf, ldmeta, input$ld_pro_inp1, input$ld_pro_inp2, input$ld_pro_sub1, input$ld_pro_sub2, input$ld_pro_typ, input$ld_pro_flp, input$ld_pro_fsz)
})

output$ld_pro_oup.ui <- renderUI({
  show_progress(imageOutput("ld_pro_oup", height = pList2[input$ld_pro_psz]))
})

output$ld_pro_oup.png <- downloadHandler(
  filename = function() { tolower(paste0("ld", "_", input$ld_pro_typ, "_", input$ld_pro_inp1, "_", input$ld_pro_inp2, ".png")) },
  content = function(file) {
    ggsave(
    file, device = "png", bg = "white", dpi = input$ld_pro_oup.res,
    plot = scProp(ldconf, ldmeta, input$ld_pro_inp1, input$ld_pro_inp2, input$ld_pro_sub1, input$ld_pro_sub2, input$ld_pro_typ, input$ld_pro_flp, input$ld_pro_fsz)
    )
  }) 
  
output$ld_pro_oup.pdf <- downloadHandler(
  filename = function() { tolower(paste0("ld", "_", input$ld_pro_typ, "_", input$ld_pro_inp1, "_", input$ld_pro_inp2, ".pdf")) },
  content = function(file) {
    ggsave(
    file, device = "pdf", useDingbats = FALSE, bg = "white", onefile = TRUE,
    plot = scProp(ldconf, ldmeta, input$ld_pro_inp1, input$ld_pro_inp2, input$ld_pro_sub1, input$ld_pro_sub2, input$ld_pro_typ, input$ld_pro_flp, input$ld_pro_fsz)
    )
  })
  
output$ld_pro_oup.svg <- downloadHandler(
  filename = function() { tolower(paste0("ld", "_", input$ld_pro_typ, "_", input$ld_pro_inp1, "_", input$ld_pro_inp2, ".svg")) },
  content = function(file) {
    ggsave(
    file, device = "svg", bg = "white",
    plot = scProp(ldconf, ldmeta, input$ld_pro_inp1, input$ld_pro_inp2, input$ld_pro_sub1, input$ld_pro_sub2, input$ld_pro_typ, input$ld_pro_flp, input$ld_pro_fsz)
    )
  }) # End of tab pro



### Tab hea heatmap / dotplot ----

  output$ld_hea_sub1.ui <- renderUI({
    sub = strsplit(ldconf[UI == input$ld_hea_sub1]$fID, "\\|")[[1]]
    checkboxGroupInput("ld_hea_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
  })
  observeEvent(input$ld_hea_sub1non, {
    sub = strsplit(ldconf[UI == input$ld_hea_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "ld_hea_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
  })
  observeEvent(input$ld_hea_sub1all, {
    sub = strsplit(ldconf[UI == input$ld_hea_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "ld_hea_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
  })

output$ld_hea_oupTxt <- renderUI({
  geneList = scGeneList(input$ld_hea_inp, ldgene)
  if(nrow(geneList) > 50){
    HTML("More than 50 input genes! Please reduce the gene list!")
  }
})

output$ld_hea_oup <- renderPlot({
  scBubbHeat(ldconf, ldmeta, input$ld_hea_inp, input$ld_hea_grp, input$ld_hea_plt, input$ld_hea_sub1, input$ld_hea_sub2, "ldgexpr.h5", ldgene, input$ld_hea_scl, input$ld_hea_row, input$ld_hea_col, input$ld_hea_cols, input$ld_hea_fsz)
})

output$ld_hea_oup.ui <- renderUI({
  show_progress(imageOutput("ld_hea_oup", height = pList3[input$ld_hea_psz]))
})

output$ld_hea_oup.png <- downloadHandler(
  filename = function() { tolower(paste0("ld", "_", input$ld_hea_plt,"_",input$ld_hea_grp,".png")) },
  content = function(file) {
    ggsave(
    file, device = "png", bg = "white", height = input$ld_hea_oup.height, width = input$ld_hea_oup.width, units = "cm", dpi = input$ld_hea_oup.res, plot = scBubbHeat(ldconf, ldmeta, input$ld_hea_inp, input$ld_hea_grp, input$ld_hea_plt, input$ld_hea_sub1, input$ld_hea_sub2, "ldgexpr.h5", ldgene, input$ld_hea_scl, input$ld_hea_row, input$ld_hea_col, input$ld_hea_cols, input$ld_hea_fsz)
    )
})

output$ld_hea_oup.pdf <- downloadHandler(
  filename = function() { tolower(paste0("ld", "_", input$ld_hea_plt,"_",input$ld_hea_grp,".pdf")) },
  content = function(file) {
    pdf(file, useDingbats = FALSE, bg = "white", height = input$ld_hea_oup.height/2.54, width = input$ld_hea_oup.width/2.54, onefile = TRUE)
    showtext::showtext_begin()
    print(scBubbHeat(ldconf, ldmeta, input$ld_hea_inp, input$ld_hea_grp, input$ld_hea_plt, input$ld_hea_sub1, input$ld_hea_sub2, "ldgexpr.h5", ldgene, input$ld_hea_scl, input$ld_hea_row, input$ld_hea_col, input$ld_hea_cols, input$ld_hea_fsz))
    showtext::showtext_end()
    dev.off()
})

output$ld_hea_oup.svg <- downloadHandler(
  filename = function() { tolower(paste0("ld", "_", input$ld_hea_plt,"_",input$ld_hea_grp,".svg")) },
  content = function(file) {
    ggsave(
    file, device = "svg", bg = "white", height = input$ld_hea_oup.height, width = input$ld_hea_oup.width, units = "cm", 
    plot = scBubbHeat(ldconf, ldmeta, input$ld_hea_inp, input$ld_hea_grp, input$ld_hea_plt, input$ld_hea_sub1, input$ld_hea_sub2, "ldgexpr.h5", ldgene, input$ld_hea_scl, input$ld_hea_row, input$ld_hea_col, input$ld_hea_cols, input$ld_hea_fsz)
    )
}) # End of tab hea      
       


### Tab markers ----

output$ld_mar_table <- renderDataTable({
  req(input$ld_mar_cls)
  datatable(ldmar[[input$ld_mar_cls]], rownames = FALSE, extensions = "Buttons", options = list(dom = "lftiprB", buttons = c("copy", "csv", "excel")))
}) # End of tab mar
optCrt="{ option_create: function(data,escape) {return('<div class=\"create\"><strong>' + '</strong></div>');} }"
updateSelectizeInput(session, "lh_civge_inp2", choices = names(lhgene), server = TRUE,
                     selected = lhdef$gene1, options = list(
                       maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
updateSelectizeInput(session, "lh_gevge_inp1", choices = names(lhgene), server = TRUE,
                     selected = lhdef$gene1, options = list(
                       maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
updateSelectizeInput(session, "lh_gevge_inp2", choices = names(lhgene), server = TRUE,
                     selected = lhdef$gene2, options = list(
                       maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
updateSelectizeInput(session, "lh_gec_inp1", choices = names(lhgene), server = TRUE,
                     selected = lhdef$gene1, options = list(
                       maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
updateSelectizeInput(session, "lh_gec_inp2", choices = names(lhgene), server = TRUE,
                     selected = lhdef$gene2, options = list(
                       maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
updateSelectizeInput(session, "lh_gem_inp", choices = names(lhgene), server = TRUE,
                     selected = lhdef$genes[1:9], options = list(
                     create = TRUE, persist = TRUE, render = I(optCrt)))
updateSelectizeInput(session, "lh_hea_inp", choices = names(lhgene), server = TRUE,
                     selected = lhdef$genes[1:12], options = list(
                     create = TRUE, persist = TRUE, render = I(optCrt)))
updateSelectizeInput(session, "lh_vio_inp2", server = TRUE,
                     choices = c(lhconf[is.na(fID)]$UI,names(lhgene)),
                     selected = lhconf[is.na(fID)]$UI[1], options = list(
                       maxOptions = length(lhconf[is.na(fID)]$UI) + 3,
                       create = TRUE, persist = TRUE, render = I(optCrt)))  
### Tab civge cell info vs gene exp ----

  output$lh_civge_sub1.ui <- renderUI({
    sub = strsplit(lhconf[UI == input$lh_civge_sub1]$fID, "\\|")[[1]]
    checkboxGroupInput("lh_civge_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
  })
  observeEvent(input$lh_civge_sub1non, {
    sub = strsplit(lhconf[UI == input$lh_civge_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "lh_civge_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
  })
  observeEvent(input$lh_civge_sub1all, {
    sub = strsplit(lhconf[UI == input$lh_civge_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "lh_civge_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
  })

output$lh_civge_oup1 <- renderPlot({
  req(input$lh_civge_inp1)
  scDRcell(lhconf, lhmeta, input$lh_civge_drX, input$lh_civge_drY, input$lh_civge_inp1, input$lh_civge_sub1, input$lh_civge_sub2, input$lh_civge_siz, input$lh_civge_col1, input$lh_civge_ord1, input$lh_civge_fsz, input$lh_civge_asp, input$lh_civge_txt, input$lh_civge_lab1)
})

output$lh_civge_oup1.ui <- renderUI({
  show_progress(imageOutput("lh_civge_oup1", height = pList[input$lh_civge_psz]))
})

output$lh_civge_oup1.png <- downloadHandler(
 filename = function() { tolower(paste0("lh", "_", input$lh_civge_drX, "_", input$lh_civge_drY, "_", input$lh_civge_inp1, ".png")) },
 content = function(file) {
   ggsave(
   file, device = "png", dpi = input$lh_civge_oup1.res, bg = "white",
   plot = scDRcell(lhconf, lhmeta, input$lh_civge_drX, input$lh_civge_drY, input$lh_civge_inp1,   input$lh_civge_sub1, input$lh_civge_sub2, input$lh_civge_siz, input$lh_civge_col1, input$lh_civge_ord1,  input$lh_civge_fsz, input$lh_civge_asp, input$lh_civge_txt, input$lh_civge_lab1)
   )
})

output$lh_civge_oup1.pdf <- downloadHandler(
 filename = function() { tolower(paste0("lh", "_", input$lh_civge_drX,"_", input$lh_civge_drY,"_", input$lh_civge_inp1,".pdf")) },
 content = function(file) {
   ggsave(
   file, device = "pdf", useDingbats = FALSE, bg = "white", onefile = TRUE,
   plot = scDRcell(lhconf, lhmeta, input$lh_civge_drX, input$lh_civge_drY, input$lh_civge_inp1,   input$lh_civge_sub1, input$lh_civge_sub2, input$lh_civge_siz, input$lh_civge_col1, input$lh_civge_ord1,  input$lh_civge_fsz, input$lh_civge_asp, input$lh_civge_txt, input$lh_civge_lab1)
   )
})

output$lh_civge_oup1.svg <- downloadHandler(
 filename = function() { tolower(paste0("lh", "_", input$lh_civge_drX,"_", input$lh_civge_drY,"_", input$lh_civge_inp1,".svg")) },
 content = function(file) {
   ggsave(
   file, device = "svg", bg = "white",
   plot = scDRcell(lhconf, lhmeta, input$lh_civge_drX, input$lh_civge_drY, input$lh_civge_inp1,   input$lh_civge_sub1, input$lh_civge_sub2, input$lh_civge_siz, input$lh_civge_col1, input$lh_civge_ord1,  input$lh_civge_fsz, input$lh_civge_asp, input$lh_civge_txt, input$lh_civge_lab1)
   )
})

output$lh_civge_.dt <- renderDataTable({
 req(input$lh_civge_inp2)
 ggData = scDRnum(lhconf, lhmeta, input$lh_civge_inp1, input$lh_civge_inp2, input$lh_civge_sub1, input$lh_civge_sub2, "lhgexpr.h5", lhgene, input$lh_civge_splt)
 datatable(ggData, rownames = FALSE, extensions = "Buttons", options = list(pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel"))) %>%
   formatRound(columns = c("pctExpress"), digits = 2)
})

output$lh_civge_oup2 <- renderPlot({
 req(input$lh_civge_inp2)
 scDRgene(lhconf, lhmeta, input$lh_civge_drX, input$lh_civge_drY, input$lh_civge_inp2, input$lh_civge_sub1, input$lh_civge_sub2, "lhgexpr.h5", lhgene, input$lh_civge_siz, input$lh_civge_col2, input$lh_civge_ord2, input$lh_civge_fsz, input$lh_civge_asp, input$lh_civge_txt)
})

output$lh_civge_oup2.ui <- renderUI({
 show_progress(imageOutput("lh_civge_oup2", height = pList[input$lh_civge_psz]))
})

output$lh_civge_oup2.png <- downloadHandler(
 filename = function() { tolower(paste0("lh", "_", input$lh_civge_drX,"_",input$lh_civge_drY,"_", input$lh_civge_inp2,".png")) },
 content = function(file) {
   ggsave(
   file, device = "png", dpi = input$lh_civge_oup2.res, bg = "white",
   plot = scDRgene(lhconf, lhmeta, input$lh_civge_drX, input$lh_civge_drY, input$lh_civge_inp2, input$lh_civge_sub1, input$lh_civge_sub2, "lhgexpr.h5", lhgene, input$lh_civge_siz, input$lh_civge_col2, input$lh_civge_ord2, input$lh_civge_fsz, input$lh_civge_asp, input$lh_civge_txt)
   )
})

output$lh_civge_oup2.pdf <- downloadHandler(
 filename = function() { tolower(paste0("lh", "_", input$lh_civge_drX,"_",input$lh_civge_drY,"_", input$lh_civge_inp2,".pdf")) },
 content = function(file) {
   ggsave(
   file, device = "pdf", useDingbats = FALSE, bg = "white", onefile = TRUE,
   plot = scDRgene(lhconf, lhmeta, input$lh_civge_drX, input$lh_civge_drY, input$lh_civge_inp2,  input$lh_civge_sub1, input$lh_civge_sub2, "lhgexpr.h5", lhgene, input$lh_civge_siz, input$lh_civge_col2, input$lh_civge_ord2, input$lh_civge_fsz, input$lh_civge_asp, input$lh_civge_txt)
   )
}) 

output$lh_civge_oup2.svg <- downloadHandler(
 filename = function() { tolower(paste0("lh", "_", input$lh_civge_drX,"_",input$lh_civge_drY,"_", input$lh_civge_inp2,".svg")) },
 content = function(file) {
   ggsave(
   file, device = "svg", bg = "white",
   plot = scDRgene(lhconf, lhmeta, input$lh_civge_drX, input$lh_civge_drY, input$lh_civge_inp2,  input$lh_civge_sub1, input$lh_civge_sub2, "lhgexpr.h5", lhgene, input$lh_civge_siz, input$lh_civge_col2, input$lh_civge_ord2, input$lh_civge_fsz, input$lh_civge_asp, input$lh_civge_txt)
   )
}) # End of tab civge



### Tab civci cell info vs cell info ----
  
  output$lh_civci_sub1.ui <- renderUI({
    sub = strsplit(lhconf[UI == input$lh_civci_sub1]$fID, "\\|")[[1]]
    checkboxGroupInput("lh_civci_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
  })
  observeEvent(input$lh_civci_sub1non, {
    sub = strsplit(lhconf[UI == input$lh_civci_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "lh_civci_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
  })
  observeEvent(input$lh_civci_sub1all, {
    sub = strsplit(lhconf[UI == input$lh_civci_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "lh_civci_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
  })

output$lh_civci_oup1 <- renderPlot({
  req(input$lh_civci_inp1)
  scDRcell(lhconf, lhmeta, input$lh_civci_drX, input$lh_civci_drY, input$lh_civci_inp1, input$lh_civci_sub1, input$lh_civci_sub2, input$lh_civci_siz, input$lh_civci_col1, input$lh_civci_ord1, input$lh_civci_fsz, input$lh_civci_asp, input$lh_civci_txt, input$lh_civci_lab1)
})

output$lh_civci_oup1.ui <- renderUI({
  show_progress(imageOutput("lh_civci_oup1", height = pList[input$lh_civci_psz]))
})

output$lh_civci_oup1.png <- downloadHandler(
  filename = function() { tolower(paste0("lh", "_", input$lh_civci_drX, "_", input$lh_civci_drY, "_", input$lh_civci_inp1, ".png")) },
  content = function(file) {
    ggsave(
    file, device = "png", dpi = input$lh_civci_oup1.res, bg = "white",
    plot = scDRcell(lhconf, lhmeta, input$lh_civci_drX, input$lh_civci_drY, input$lh_civci_inp1, input$lh_civci_sub1, input$lh_civci_sub2, input$lh_civci_siz, input$lh_civci_col1, input$lh_civci_ord1, input$lh_civci_fsz, input$lh_civci_asp, input$lh_civci_txt, input$lh_civci_lab1)
    )
})

output$lh_civci_oup1.pdf <- downloadHandler(
  filename = function() { tolower(paste0("lh", "_", input$lh_civci_drX, "_", input$lh_civci_drY, "_", input$lh_civci_inp1, ".pdf")) },
  content = function(file) { ggsave(
    file, device = "pdf", useDingbats = FALSE, bg = "white", onefile = TRUE,
    plot = scDRcell(lhconf, lhmeta, input$lh_civci_drX, input$lh_civci_drY, input$lh_civci_inp1, input$lh_civci_sub1, input$lh_civci_sub2, input$lh_civci_siz, input$lh_civci_col1, input$lh_civci_ord1, input$lh_civci_fsz, input$lh_civci_asp, input$lh_civci_txt, input$lh_civci_lab1) )
})

output$lh_civci_oup1.svg <- downloadHandler(
  filename = function() { tolower(paste0("lh", "_", input$lh_civci_drX, "_", input$lh_civci_drY, "_", input$lh_civci_inp1, ".svg")) },
  content = function(file) { ggsave(
    file, device = "svg", bg = "white",
    plot = scDRcell(lhconf, lhmeta, input$lh_civci_drX, input$lh_civci_drY, input$lh_civci_inp1, input$lh_civci_sub1, input$lh_civci_sub2, input$lh_civci_siz, input$lh_civci_col1, input$lh_civci_ord1, input$lh_civci_fsz, input$lh_civci_asp, input$lh_civci_txt, input$lh_civci_lab1) )
})

output$lh_civci_oup2 <- renderPlot({
  req(input$lh_civci_inp2)
  scDRcell(lhconf, lhmeta, input$lh_civci_drX, input$lh_civci_drY, input$lh_civci_inp2, input$lh_civci_sub1, input$lh_civci_sub2, input$lh_civci_siz, input$lh_civci_col2, input$lh_civci_ord2, input$lh_civci_fsz, input$lh_civci_asp, input$lh_civci_txt, input$lh_civci_lab2)
})

output$lh_civci_oup2.ui <- renderUI({
  show_progress(imageOutput("lh_civci_oup2", height = pList[input$lh_civci_psz]))
})

output$lh_civci_oup2.png <- downloadHandler(
  filename = function() { tolower(paste0("lh", "_", input$lh_civci_drX,"_",input$lh_civci_drY,"_", input$lh_civci_inp2,".png")) },
  content = function(file) {
    ggsave(
    file, device = "png", bg = "white", dpi = input$lh_civci_oup2.res,
    plot = scDRcell(lhconf, lhmeta, input$lh_civci_drX, input$lh_civci_drY, input$lh_civci_inp2, input$lh_civci_sub1, input$lh_civci_sub2, input$lh_civci_siz, input$lh_civci_col2, input$lh_civci_ord2, input$lh_civci_fsz, input$lh_civci_asp, input$lh_civci_txt, input$lh_civci_lab2)
    )
}) 

output$lh_civci_oup2.pdf <- downloadHandler(
  filename = function() { tolower(paste0("lh", "_", input$lh_civci_drX,"_",input$lh_civci_drY,"_", input$lh_civci_inp2,".pdf")) },
  content = function(file) {
    ggsave(
    file, device = "pdf", useDingbats = FALSE, bg = "white", onefile = TRUE,
    plot = scDRcell(lhconf, lhmeta, input$lh_civci_drX, input$lh_civci_drY, input$lh_civci_inp2, input$lh_civci_sub1, input$lh_civci_sub2, input$lh_civci_siz, input$lh_civci_col2, input$lh_civci_ord2, input$lh_civci_fsz, input$lh_civci_asp, input$lh_civci_txt, input$lh_civci_lab2) 
    )
})

output$lh_civci_oup2.svg <- downloadHandler(
  filename = function() { tolower(paste0("lh", "_", input$lh_civci_drX,"_",input$lh_civci_drY,"_", input$lh_civci_inp2,".svg")) },
  content = function(file) {
    ggsave(
    file, device = "svg", bg = "white",
    plot = scDRcell(lhconf, lhmeta, input$lh_civci_drX, input$lh_civci_drY, input$lh_civci_inp2, input$lh_civci_sub1, input$lh_civci_sub2, input$lh_civci_siz, input$lh_civci_col2, input$lh_civci_ord2, input$lh_civci_fsz, input$lh_civci_asp, input$lh_civci_txt, input$lh_civci_lab2) 
    )
}) # End of tab civci



### Tab gevge gene exp vs gene exp ----

  output$lh_gevge_sub1.ui <- renderUI({
    sub = strsplit(lhconf[UI == input$lh_gevge_sub1]$fID, "\\|")[[1]]
    checkboxGroupInput("lh_gevge_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
  })
  observeEvent(input$lh_gevge_sub1non, {
    sub = strsplit(lhconf[UI == input$lh_gevge_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "lh_gevge_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
  })
  observeEvent(input$lh_gevge_sub1all, {
    sub = strsplit(lhconf[UI == input$lh_gevge_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "lh_gevge_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
  })

output$lh_gevge_oup1 <- renderPlot({
  req(input$lh_gevge_inp1)
  scDRgene(lhconf, lhmeta, input$lh_gevge_drX, input$lh_gevge_drY, input$lh_gevge_inp1, input$lh_gevge_sub1, input$lh_gevge_sub2, "lhgexpr.h5", lhgene, input$lh_gevge_siz, input$lh_gevge_col1, input$lh_gevge_ord1, input$lh_gevge_fsz, input$lh_gevge_asp, input$lh_gevge_txt)
})

output$lh_gevge_oup1.ui <- renderUI({
  show_progress(imageOutput("lh_gevge_oup1", height = pList[input$lh_gevge_psz]))
})

output$lh_gevge_oup1.png <- downloadHandler(
  filename = function() { tolower(paste0("lh", "_", input$lh_gevge_drX, "_", input$lh_gevge_drY, "_", input$lh_gevge_inp1,".png")) },
  content = function(file) {
    ggsave(
    file, device = "png", bg = "white", dpi = input$lh_gevge_oup1.res,
    plot = scDRgene(lhconf, lhmeta, input$lh_gevge_drX, input$lh_gevge_drY, input$lh_gevge_inp1, 
                    input$lh_gevge_sub1, input$lh_gevge_sub2,
                    "lhgexpr.h5", lhgene,
                    input$lh_gevge_siz, input$lh_gevge_col1, input$lh_gevge_ord1,
                    input$lh_gevge_fsz, input$lh_gevge_asp, input$lh_gevge_txt)
    )
})

output$lh_gevge_oup1.pdf <- downloadHandler(
  filename = function() { tolower(paste0("lh", "_", input$lh_gevge_drX, "_", input$lh_gevge_drY, "_", input$lh_gevge_inp1, ".pdf")) },
  content = function(file) {
    ggsave(
    file, device = "pdf", useDingbats = FALSE, bg = "white", onefile = TRUE,
    plot = scDRgene(lhconf, lhmeta, input$lh_gevge_drX, input$lh_gevge_drY, input$lh_gevge_inp1, input$lh_gevge_sub1, input$lh_gevge_sub2, "lhgexpr.h5", lhgene, input$lh_gevge_siz, input$lh_gevge_col1, input$lh_gevge_ord1, input$lh_gevge_fsz, input$lh_gevge_asp, input$lh_gevge_txt)
    )
})

output$lh_gevge_oup1.svg <- downloadHandler(
  filename = function() { tolower(paste0("lh", "_", input$lh_gevge_drX, "_", input$lh_gevge_drY, "_", input$lh_gevge_inp1, ".svg")) },
  content = function(file) {
    ggsave(
    file, device = "svg", bg = "white",
    plot = scDRgene(lhconf, lhmeta, input$lh_gevge_drX, input$lh_gevge_drY, input$lh_gevge_inp1, input$lh_gevge_sub1, input$lh_gevge_sub2, "lhgexpr.h5", lhgene, input$lh_gevge_siz, input$lh_gevge_col1, input$lh_gevge_ord1, input$lh_gevge_fsz, input$lh_gevge_asp, input$lh_gevge_txt)
    )
})

output$lh_gevge_oup2 <- renderPlot({
  req(input$lh_gevge_inp2)
  scDRgene(lhconf, lhmeta, input$lh_gevge_drX, input$lh_gevge_drY, input$lh_gevge_inp2, input$lh_gevge_sub1, input$lh_gevge_sub2, "lhgexpr.h5", lhgene, input$lh_gevge_siz, input$lh_gevge_col2, input$lh_gevge_ord2, input$lh_gevge_fsz, input$lh_gevge_asp, input$lh_gevge_txt)
})

output$lh_gevge_oup2.ui <- renderUI({
  show_progress(imageOutput("lh_gevge_oup2", height = pList[input$lh_gevge_psz]))
})

output$lh_gevge_oup2.png <- downloadHandler(
  filename = function() { tolower(paste0("lh", "_", input$lh_gevge_drX, "_", input$lh_gevge_drY, "_", input$lh_gevge_inp2,".png")) },
  content = function(file) {
    ggsave(
    file, device = "png", bg = "white", dpi = input$lh_gevge_oup2.res,
    plot = scDRgene(lhconf, lhmeta, input$lh_gevge_drX, input$lh_gevge_drY, input$lh_gevge_inp2, input$lh_gevge_sub1, input$lh_gevge_sub2, "lhgexpr.h5", lhgene, input$lh_gevge_siz, input$lh_gevge_col2, input$lh_gevge_ord2, input$lh_gevge_fsz, input$lh_gevge_asp, input$lh_gevge_txt)
    )
}) 

output$lh_gevge_oup2.pdf <- downloadHandler(
  filename = function() { tolower(paste0("lh", "_", input$lh_gevge_drX, "_", input$lh_gevge_drY, "_", input$lh_gevge_inp2,".pdf")) },
  content = function(file) {
    ggsave(
    file, device = "pdf", useDingbats = FALSE, bg = "white", onefile = TRUE,
    plot = scDRgene(lhconf, lhmeta, input$lh_gevge_drX, input$lh_gevge_drY, input$lh_gevge_inp2, 
                    input$lh_gevge_sub1, input$lh_gevge_sub2,
                    "lhgexpr.h5", lhgene,
                    input$lh_gevge_siz, input$lh_gevge_col2, input$lh_gevge_ord2,
                    input$lh_gevge_fsz, input$lh_gevge_asp, input$lh_gevge_txt)
    )
})

output$lh_gevge_oup2.svg <- downloadHandler(
  filename = function() { tolower(paste0("lh", "_", input$lh_gevge_drX, "_", input$lh_gevge_drY, "_", input$lh_gevge_inp2,".svg")) },
  content = function(file) {
    ggsave(
    file, device = "svg", bg = "white",
    plot = scDRgene(lhconf, lhmeta, input$lh_gevge_drX, input$lh_gevge_drY, input$lh_gevge_inp2, 
                    input$lh_gevge_sub1, input$lh_gevge_sub2,
                    "lhgexpr.h5", lhgene,
                    input$lh_gevge_siz, input$lh_gevge_col2, input$lh_gevge_ord2,
                    input$lh_gevge_fsz, input$lh_gevge_asp, input$lh_gevge_txt)
    )
}) # End of tab gevge




### Tab gem gene expression multi ----

output$lh_gem_sub1.ui <- renderUI({
  sub = strsplit(lhconf[UI == input$lh_gem_sub1]$fID, "\\|")[[1]]
  checkboxGroupInput("lh_gem_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
})
observeEvent(input$lh_gem_sub1non, {
  sub = strsplit(lhconf[UI == input$lh_gem_sub1]$fID, "\\|")[[1]]
  updateCheckboxGroupInput(session, inputId = "lh_gem_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
})
observeEvent(input$lh_gem_sub1all, {
  sub = strsplit(lhconf[UI == input$lh_gem_sub1]$fID, "\\|")[[1]]
  updateCheckboxGroupInput(session, inputId = "lh_gem_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
})

output$lh_gem_oup1 <- renderPlot({
  req(input$lh_gem_inp)
  
  scFeature(lhconf, lhmeta, input$lh_gem_drX, input$lh_gem_drY, input$lh_gem_inp, input$lh_gem_sub1, input$lh_gem_sub2, "lhgexpr.h5", lhgene, input$lh_gem_siz, input$lh_gem_col, input$lh_gem_ord, input$lh_gem_fsz, input$lh_gem_asp, input$lh_gem_txt, input$lh_gem_ncol)
})

output$lh_gem_oup1.ui <- renderUI({
  show_progress(imageOutput("lh_gem_oup1", height = pList[input$lh_gem_psz]))
})

output$lh_gem_oup1.png <- downloadHandler(
  filename = function() { tolower(paste0("lh", "_", input$lh_gem_drX, "_", input$lh_gem_drY, "_expression.png")) },
  content = function(file) {
    ggsave(
      file, device = "png", height = input$lh_gem_oup1.height, width = input$lh_gem_oup1.width, dpi = input$lh_gem_oup1.res, units = "cm", bg = "white",
      plot = scFeature(lhconf, lhmeta, input$lh_gem_drX, input$lh_gem_drY, input$lh_gem_inp, input$lh_gem_sub1, input$lh_gem_sub2, "lhgexpr.h5", lhgene, input$lh_gem_siz, input$lh_gem_col, input$lh_gem_ord, input$lh_gem_fsz, input$lh_gem_asp, input$lh_gem_txt, input$lh_gem_ncol)
    )
}) 

output$lh_gem_oup1.pdf <- downloadHandler(
  filename = function() { tolower(paste0("lh", "_", input$lh_gem_drX, "_", input$lh_gem_drY, "_expression.pdf")) },
  content = function(file) {
  pdf(file, useDingbats = FALSE, height = input$lh_gem_oup1.height/2.54, width = input$lh_gem_oup1.width/2.54, bg = "white", onefile = TRUE)
  showtext::showtext_begin()
  print(scFeature(lhconf, lhmeta, input$lh_gem_drX, input$lh_gem_drY, input$lh_gem_inp, input$lh_gem_sub1, input$lh_gem_sub2, "lhgexpr.h5", lhgene, input$lh_gem_siz, input$lh_gem_col, input$lh_gem_ord, input$lh_gem_fsz, input$lh_gem_asp, input$lh_gem_txt, input$lh_gem_ncol))
  showtext::showtext_end()
  dev.off()
})

output$lh_gem_oup1.svg <- downloadHandler(
  filename = function() { tolower(paste0("lh", "_", input$lh_gem_drX, "_", input$lh_gem_drY, "_expression.svg")) },
  content = function(file) { ggsave(
    file, device = "svg", height = input$lh_gem_oup1.height, width = input$lh_gem_oup1.width, units = "cm", bg = "white",
    plot = scFeature(lhconf, lhmeta, input$lh_gem_drX, input$lh_gem_drY, input$lh_gem_inp, input$lh_gem_sub1, input$lh_gem_sub2, "lhgexpr.h5", lhgene, input$lh_gem_siz, input$lh_gem_col, input$lh_gem_ord, input$lh_gem_fsz, input$lh_gem_asp, input$lh_gem_txt, input$lh_gem_ncol))
}) # End of tab gem



### Tab gec gene co-expression ----

  output$lh_gec_sub1.ui <- renderUI({
    sub = strsplit(lhconf[UI == input$lh_gec_sub1]$fID, "\\|")[[1]]
    checkboxGroupInput("lh_gec_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
  })
  observeEvent(input$lh_gec_sub1non, {
    sub = strsplit(lhconf[UI == input$lh_gec_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "lh_gec_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
  })
  observeEvent(input$lh_gec_sub1all, {
    sub = strsplit(lhconf[UI == input$lh_gec_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "lh_gec_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
  })

output$lh_gec_oup1 <- renderPlot({
  scDRcoexFull(lhconf, lhmeta, input$lh_gec_drX, input$lh_gec_drY, input$lh_gec_inp1, input$lh_gec_inp2, input$lh_gec_sub1, input$lh_gec_sub2, "lhgexpr.h5", lhgene, input$lh_gec_siz, input$lh_gec_col1, input$lh_gec_ord1, input$lh_gec_fsz, input$lh_gec_asp, input$lh_gec_txt)
})

output$lh_gec_oup1.ui <- renderUI({
  show_progress(imageOutput("lh_gec_oup1", height = pList2[input$lh_gec_psz]))
})

output$lh_gec_oup1.png <- downloadHandler(
  filename = function() { tolower(paste0("lh", "_", input$lh_gec_drX, "_", input$lh_gec_drY, "_", input$lh_gec_inp1, "_", input$lh_gec_inp2, ".png")) },
  content = function(file) { ggsave(
    file, device = "png", bg = "white", dpi = input$lh_gec_oup1.res,
    plot = scDRcoexFull(lhconf, lhmeta, input$lh_gec_drX, input$lh_gec_drY, input$lh_gec_inp1, input$lh_gec_inp2, input$lh_gec_sub1, input$lh_gec_sub2, "lhgexpr.h5", lhgene, input$lh_gec_siz, input$lh_gec_col1, input$lh_gec_ord1, input$lh_gec_fsz, input$lh_gec_asp, input$lh_gec_txt) )
})

output$lh_gec_oup1.pdf <- downloadHandler(
  filename = function() { tolower(paste0("lh", "_", input$lh_gec_drX, "_", input$lh_gec_drY, "_", input$lh_gec_inp1, "_", input$lh_gec_inp2, ".pdf")) },
  content = function(file) {
    ggsave(
    file, device = "pdf", useDingbats = FALSE, bg = "white", onefile = TRUE,
    plot = scDRcoexFull(lhconf, lhmeta, input$lh_gec_drX, input$lh_gec_drY, input$lh_gec_inp1, input$lh_gec_inp2, input$lh_gec_sub1, input$lh_gec_sub2, "lhgexpr.h5", lhgene, input$lh_gec_siz, input$lh_gec_col1, input$lh_gec_ord1, input$lh_gec_fsz, input$lh_gec_asp, input$lh_gec_txt)
    )
})

output$lh_gec_oup1.svg <- downloadHandler(
  filename = function() { tolower(paste0("lh", "_", input$lh_gec_drX, "_", input$lh_gec_drY, "_", input$lh_gec_inp1, "_", input$lh_gec_inp2, ".svg")) },
  content = function(file) {
    ggsave(
    file, device = "svg", bg = "white",
    plot = scDRcoexFull(lhconf, lhmeta, input$lh_gec_drX, input$lh_gec_drY, input$lh_gec_inp1, input$lh_gec_inp2, input$lh_gec_sub1, input$lh_gec_sub2, "lhgexpr.h5", lhgene, input$lh_gec_siz, input$lh_gec_col1, input$lh_gec_ord1, input$lh_gec_fsz, input$lh_gec_asp, input$lh_gec_txt)
    )
})

output$lh_gec_.dt <- renderDataTable({
  ggData = scDRcoexNum(lhconf, lhmeta, input$lh_gec_inp1, input$lh_gec_inp2, input$lh_gec_sub1, input$lh_gec_sub2, "lhgexpr.h5", lhgene)
  datatable(ggData, rownames = FALSE, extensions = "Buttons", options = list(pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel"))) %>%
            formatRound(columns = c("percent"), digits = 2)
}) # End of tab gec



### Tab vio violinplot / boxplot / lineplot ----

  output$lh_vio_sub1.ui <- renderUI({
    sub = strsplit(lhconf[UI == input$lh_vio_sub1]$fID, "\\|")[[1]]
    checkboxGroupInput("lh_vio_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
  })
  observeEvent(input$lh_vio_sub1non, {
    sub = strsplit(lhconf[UI == input$lh_vio_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "lh_vio_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
  })
  observeEvent(input$lh_vio_sub1all, {
    sub = strsplit(lhconf[UI == input$lh_vio_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "lh_vio_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
  })

output$lh_vio_oup <- renderPlot({
  gh5 <- ifelse(input$lh_vio_datatype == "normalised","lhgexpr.h5","lhgexpr2.h5")
  scVioBox(lhconf, lhmeta, input$lh_vio_inp1, input$lh_vio_inp2, input$lh_vio_sub1, input$lh_vio_sub2, gh5, lhgene, input$lh_vio_typ, input$lh_vio_pts, input$lh_vio_siz, input$lh_vio_fsz, input$lh_vio_barsz)
})

output$lh_vio_oup.ui <- renderUI({
  show_progress(imageOutput("lh_vio_oup", height = pList2[input$lh_vio_psz]))
})

output$lh_vio_oup.png <- downloadHandler(
  filename = function() { tolower(paste0("lh", "_", input$lh_vio_typ, "_", input$lh_vio_datatype, "_", input$lh_vio_inp1, "_", input$lh_vio_inp2,".png")) },
  content = function(file) {
    gh5 <- ifelse(input$lh_vio_datatype == "normalised","lhgexpr.h5","lhgexpr2.h5")
    ggsave(
    file, device = "png", bg = "white", dpi = input$lh_vio_oup.res,
    plot = scVioBox(lhconf, lhmeta, input$lh_vio_inp1, input$lh_vio_inp2, input$lh_vio_sub1, input$lh_vio_sub2, gh5, lhgene, input$lh_vio_typ, input$lh_vio_pts, input$lh_vio_siz, input$lh_vio_fsz, input$lh_vio_barsz)
    )
}) 

output$lh_vio_oup.pdf <- downloadHandler(
  filename = function() { tolower(paste0("lh", "_", input$lh_vio_typ, "_", input$lh_vio_datatype, "_", input$lh_vio_inp1, "_", input$lh_vio_inp2, ".pdf")) },
  content = function(file) {
    gh5 <- ifelse(input$lh_vio_datatype == "normalised","lhgexpr.h5","lhgexpr2.h5")
    ggsave(
    file, device = "pdf", useDingbats = FALSE, bg = "white", onefile = TRUE,
    plot = scVioBox(lhconf, lhmeta, input$lh_vio_inp1, input$lh_vio_inp2, input$lh_vio_sub1, input$lh_vio_sub2, gh5, lhgene, input$lh_vio_typ, input$lh_vio_pts, input$lh_vio_siz, input$lh_vio_fsz, input$lh_vio_barsz)
    )
})

output$lh_vio_oup.svg <- downloadHandler(
  filename = function() { tolower(paste0("lh", "_", input$lh_vio_typ, "_", input$lh_vio_datatype, "_", input$lh_vio_inp1, "_", input$lh_vio_inp2, ".svg")) },
  content = function(file) {
    gh5 <- ifelse(input$lh_vio_datatype == "normalised","lhgexpr.h5","lhgexpr2.h5")
    ggsave(
    file, device = "svg", bg = "white",
    plot = scVioBox(lhconf, lhmeta, input$lh_vio_inp1, input$lh_vio_inp2, input$lh_vio_sub1, input$lh_vio_sub2, gh5, lhgene, input$lh_vio_typ, input$lh_vio_pts, input$lh_vio_siz, input$lh_vio_fsz, input$lh_vio_barsz)
    )
}) # End of tab vio




### Tab pro proportion plot ----

  output$lh_pro_sub1.ui <- renderUI({
    sub = strsplit(lhconf[UI == input$lh_pro_sub1]$fID, "\\|")[[1]]
    checkboxGroupInput("lh_pro_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
  })
  observeEvent(input$lh_pro_sub1non, {
    sub = strsplit(lhconf[UI == input$lh_pro_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "lh_pro_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
  })
  observeEvent(input$lh_pro_sub1all, {
    sub = strsplit(lhconf[UI == input$lh_pro_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "lh_pro_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
  })

output$lh_pro_oup <- renderPlot({
  scProp(lhconf, lhmeta, input$lh_pro_inp1, input$lh_pro_inp2, input$lh_pro_sub1, input$lh_pro_sub2, input$lh_pro_typ, input$lh_pro_flp, input$lh_pro_fsz)
})

output$lh_pro_oup.ui <- renderUI({
  show_progress(imageOutput("lh_pro_oup", height = pList2[input$lh_pro_psz]))
})

output$lh_pro_oup.png <- downloadHandler(
  filename = function() { tolower(paste0("lh", "_", input$lh_pro_typ, "_", input$lh_pro_inp1, "_", input$lh_pro_inp2, ".png")) },
  content = function(file) {
    ggsave(
    file, device = "png", bg = "white", dpi = input$lh_pro_oup.res,
    plot = scProp(lhconf, lhmeta, input$lh_pro_inp1, input$lh_pro_inp2, input$lh_pro_sub1, input$lh_pro_sub2, input$lh_pro_typ, input$lh_pro_flp, input$lh_pro_fsz)
    )
  }) 
  
output$lh_pro_oup.pdf <- downloadHandler(
  filename = function() { tolower(paste0("lh", "_", input$lh_pro_typ, "_", input$lh_pro_inp1, "_", input$lh_pro_inp2, ".pdf")) },
  content = function(file) {
    ggsave(
    file, device = "pdf", useDingbats = FALSE, bg = "white", onefile = TRUE,
    plot = scProp(lhconf, lhmeta, input$lh_pro_inp1, input$lh_pro_inp2, input$lh_pro_sub1, input$lh_pro_sub2, input$lh_pro_typ, input$lh_pro_flp, input$lh_pro_fsz)
    )
  })
  
output$lh_pro_oup.svg <- downloadHandler(
  filename = function() { tolower(paste0("lh", "_", input$lh_pro_typ, "_", input$lh_pro_inp1, "_", input$lh_pro_inp2, ".svg")) },
  content = function(file) {
    ggsave(
    file, device = "svg", bg = "white",
    plot = scProp(lhconf, lhmeta, input$lh_pro_inp1, input$lh_pro_inp2, input$lh_pro_sub1, input$lh_pro_sub2, input$lh_pro_typ, input$lh_pro_flp, input$lh_pro_fsz)
    )
  }) # End of tab pro



### Tab hea heatmap / dotplot ----

  output$lh_hea_sub1.ui <- renderUI({
    sub = strsplit(lhconf[UI == input$lh_hea_sub1]$fID, "\\|")[[1]]
    checkboxGroupInput("lh_hea_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
  })
  observeEvent(input$lh_hea_sub1non, {
    sub = strsplit(lhconf[UI == input$lh_hea_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "lh_hea_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
  })
  observeEvent(input$lh_hea_sub1all, {
    sub = strsplit(lhconf[UI == input$lh_hea_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "lh_hea_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
  })

output$lh_hea_oupTxt <- renderUI({
  geneList = scGeneList(input$lh_hea_inp, lhgene)
  if(nrow(geneList) > 50){
    HTML("More than 50 input genes! Please reduce the gene list!")
  }
})

output$lh_hea_oup <- renderPlot({
  scBubbHeat(lhconf, lhmeta, input$lh_hea_inp, input$lh_hea_grp, input$lh_hea_plt, input$lh_hea_sub1, input$lh_hea_sub2, "lhgexpr.h5", lhgene, input$lh_hea_scl, input$lh_hea_row, input$lh_hea_col, input$lh_hea_cols, input$lh_hea_fsz)
})

output$lh_hea_oup.ui <- renderUI({
  show_progress(imageOutput("lh_hea_oup", height = pList3[input$lh_hea_psz]))
})

output$lh_hea_oup.png <- downloadHandler(
  filename = function() { tolower(paste0("lh", "_", input$lh_hea_plt,"_",input$lh_hea_grp,".png")) },
  content = function(file) {
    ggsave(
    file, device = "png", bg = "white", height = input$lh_hea_oup.height, width = input$lh_hea_oup.width, units = "cm", dpi = input$lh_hea_oup.res, plot = scBubbHeat(lhconf, lhmeta, input$lh_hea_inp, input$lh_hea_grp, input$lh_hea_plt, input$lh_hea_sub1, input$lh_hea_sub2, "lhgexpr.h5", lhgene, input$lh_hea_scl, input$lh_hea_row, input$lh_hea_col, input$lh_hea_cols, input$lh_hea_fsz)
    )
})

output$lh_hea_oup.pdf <- downloadHandler(
  filename = function() { tolower(paste0("lh", "_", input$lh_hea_plt,"_",input$lh_hea_grp,".pdf")) },
  content = function(file) {
    pdf(file, useDingbats = FALSE, bg = "white", height = input$lh_hea_oup.height/2.54, width = input$lh_hea_oup.width/2.54, onefile = TRUE)
    showtext::showtext_begin()
    print(scBubbHeat(lhconf, lhmeta, input$lh_hea_inp, input$lh_hea_grp, input$lh_hea_plt, input$lh_hea_sub1, input$lh_hea_sub2, "lhgexpr.h5", lhgene, input$lh_hea_scl, input$lh_hea_row, input$lh_hea_col, input$lh_hea_cols, input$lh_hea_fsz))
    showtext::showtext_end()
    dev.off()
})

output$lh_hea_oup.svg <- downloadHandler(
  filename = function() { tolower(paste0("lh", "_", input$lh_hea_plt,"_",input$lh_hea_grp,".svg")) },
  content = function(file) {
    ggsave(
    file, device = "svg", bg = "white", height = input$lh_hea_oup.height, width = input$lh_hea_oup.width, units = "cm", 
    plot = scBubbHeat(lhconf, lhmeta, input$lh_hea_inp, input$lh_hea_grp, input$lh_hea_plt, input$lh_hea_sub1, input$lh_hea_sub2, "lhgexpr.h5", lhgene, input$lh_hea_scl, input$lh_hea_row, input$lh_hea_col, input$lh_hea_cols, input$lh_hea_fsz)
    )
}) # End of tab hea      
       


### Tab markers ----

output$lh_mar_table <- renderDataTable({
  req(input$lh_mar_cls)
  datatable(lhmar[[input$lh_mar_cls]], rownames = FALSE, extensions = "Buttons", options = list(dom = "lftiprB", buttons = c("copy", "csv", "excel")))
}) # End of tab mar
optCrt="{ option_create: function(data,escape) {return('<div class=\"create\"><strong>' + '</strong></div>');} }"
updateSelectizeInput(session, "vh_civge_inp2", choices = names(vhgene), server = TRUE,
                     selected = vhdef$gene1, options = list(
                       maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
updateSelectizeInput(session, "vh_gevge_inp1", choices = names(vhgene), server = TRUE,
                     selected = vhdef$gene1, options = list(
                       maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
updateSelectizeInput(session, "vh_gevge_inp2", choices = names(vhgene), server = TRUE,
                     selected = vhdef$gene2, options = list(
                       maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
updateSelectizeInput(session, "vh_gec_inp1", choices = names(vhgene), server = TRUE,
                     selected = vhdef$gene1, options = list(
                       maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
updateSelectizeInput(session, "vh_gec_inp2", choices = names(vhgene), server = TRUE,
                     selected = vhdef$gene2, options = list(
                       maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
updateSelectizeInput(session, "vh_gem_inp", choices = names(vhgene), server = TRUE,
                     selected = vhdef$genes[1:9], options = list(
                     create = TRUE, persist = TRUE, render = I(optCrt)))
updateSelectizeInput(session, "vh_hea_inp", choices = names(vhgene), server = TRUE,
                     selected = vhdef$genes[1:12], options = list(
                     create = TRUE, persist = TRUE, render = I(optCrt)))
updateSelectizeInput(session, "vh_vio_inp2", server = TRUE,
                     choices = c(vhconf[is.na(fID)]$UI,names(vhgene)),
                     selected = vhconf[is.na(fID)]$UI[1], options = list(
                       maxOptions = length(vhconf[is.na(fID)]$UI) + 3,
                       create = TRUE, persist = TRUE, render = I(optCrt)))  
### Tab civge cell info vs gene exp ----

  output$vh_civge_sub1.ui <- renderUI({
    sub = strsplit(vhconf[UI == input$vh_civge_sub1]$fID, "\\|")[[1]]
    checkboxGroupInput("vh_civge_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
  })
  observeEvent(input$vh_civge_sub1non, {
    sub = strsplit(vhconf[UI == input$vh_civge_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "vh_civge_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
  })
  observeEvent(input$vh_civge_sub1all, {
    sub = strsplit(vhconf[UI == input$vh_civge_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "vh_civge_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
  })

output$vh_civge_oup1 <- renderPlot({
  req(input$vh_civge_inp1)
  scDRcell(vhconf, vhmeta, input$vh_civge_drX, input$vh_civge_drY, input$vh_civge_inp1, input$vh_civge_sub1, input$vh_civge_sub2, input$vh_civge_siz, input$vh_civge_col1, input$vh_civge_ord1, input$vh_civge_fsz, input$vh_civge_asp, input$vh_civge_txt, input$vh_civge_lab1)
})

output$vh_civge_oup1.ui <- renderUI({
  show_progress(imageOutput("vh_civge_oup1", height = pList[input$vh_civge_psz]))
})

output$vh_civge_oup1.png <- downloadHandler(
 filename = function() { tolower(paste0("vh", "_", input$vh_civge_drX, "_", input$vh_civge_drY, "_", input$vh_civge_inp1, ".png")) },
 content = function(file) {
   ggsave(
   file, device = "png", dpi = input$vh_civge_oup1.res, bg = "white",
   plot = scDRcell(vhconf, vhmeta, input$vh_civge_drX, input$vh_civge_drY, input$vh_civge_inp1,   input$vh_civge_sub1, input$vh_civge_sub2, input$vh_civge_siz, input$vh_civge_col1, input$vh_civge_ord1,  input$vh_civge_fsz, input$vh_civge_asp, input$vh_civge_txt, input$vh_civge_lab1)
   )
})

output$vh_civge_oup1.pdf <- downloadHandler(
 filename = function() { tolower(paste0("vh", "_", input$vh_civge_drX,"_", input$vh_civge_drY,"_", input$vh_civge_inp1,".pdf")) },
 content = function(file) {
   ggsave(
   file, device = "pdf", useDingbats = FALSE, bg = "white", onefile = TRUE,
   plot = scDRcell(vhconf, vhmeta, input$vh_civge_drX, input$vh_civge_drY, input$vh_civge_inp1,   input$vh_civge_sub1, input$vh_civge_sub2, input$vh_civge_siz, input$vh_civge_col1, input$vh_civge_ord1,  input$vh_civge_fsz, input$vh_civge_asp, input$vh_civge_txt, input$vh_civge_lab1)
   )
})

output$vh_civge_oup1.svg <- downloadHandler(
 filename = function() { tolower(paste0("vh", "_", input$vh_civge_drX,"_", input$vh_civge_drY,"_", input$vh_civge_inp1,".svg")) },
 content = function(file) {
   ggsave(
   file, device = "svg", bg = "white",
   plot = scDRcell(vhconf, vhmeta, input$vh_civge_drX, input$vh_civge_drY, input$vh_civge_inp1,   input$vh_civge_sub1, input$vh_civge_sub2, input$vh_civge_siz, input$vh_civge_col1, input$vh_civge_ord1,  input$vh_civge_fsz, input$vh_civge_asp, input$vh_civge_txt, input$vh_civge_lab1)
   )
})

output$vh_civge_.dt <- renderDataTable({
 req(input$vh_civge_inp2)
 ggData = scDRnum(vhconf, vhmeta, input$vh_civge_inp1, input$vh_civge_inp2, input$vh_civge_sub1, input$vh_civge_sub2, "vhgexpr.h5", vhgene, input$vh_civge_splt)
 datatable(ggData, rownames = FALSE, extensions = "Buttons", options = list(pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel"))) %>%
   formatRound(columns = c("pctExpress"), digits = 2)
})

output$vh_civge_oup2 <- renderPlot({
 req(input$vh_civge_inp2)
 scDRgene(vhconf, vhmeta, input$vh_civge_drX, input$vh_civge_drY, input$vh_civge_inp2, input$vh_civge_sub1, input$vh_civge_sub2, "vhgexpr.h5", vhgene, input$vh_civge_siz, input$vh_civge_col2, input$vh_civge_ord2, input$vh_civge_fsz, input$vh_civge_asp, input$vh_civge_txt)
})

output$vh_civge_oup2.ui <- renderUI({
 show_progress(imageOutput("vh_civge_oup2", height = pList[input$vh_civge_psz]))
})

output$vh_civge_oup2.png <- downloadHandler(
 filename = function() { tolower(paste0("vh", "_", input$vh_civge_drX,"_",input$vh_civge_drY,"_", input$vh_civge_inp2,".png")) },
 content = function(file) {
   ggsave(
   file, device = "png", dpi = input$vh_civge_oup2.res, bg = "white",
   plot = scDRgene(vhconf, vhmeta, input$vh_civge_drX, input$vh_civge_drY, input$vh_civge_inp2, input$vh_civge_sub1, input$vh_civge_sub2, "vhgexpr.h5", vhgene, input$vh_civge_siz, input$vh_civge_col2, input$vh_civge_ord2, input$vh_civge_fsz, input$vh_civge_asp, input$vh_civge_txt)
   )
})

output$vh_civge_oup2.pdf <- downloadHandler(
 filename = function() { tolower(paste0("vh", "_", input$vh_civge_drX,"_",input$vh_civge_drY,"_", input$vh_civge_inp2,".pdf")) },
 content = function(file) {
   ggsave(
   file, device = "pdf", useDingbats = FALSE, bg = "white", onefile = TRUE,
   plot = scDRgene(vhconf, vhmeta, input$vh_civge_drX, input$vh_civge_drY, input$vh_civge_inp2,  input$vh_civge_sub1, input$vh_civge_sub2, "vhgexpr.h5", vhgene, input$vh_civge_siz, input$vh_civge_col2, input$vh_civge_ord2, input$vh_civge_fsz, input$vh_civge_asp, input$vh_civge_txt)
   )
}) 

output$vh_civge_oup2.svg <- downloadHandler(
 filename = function() { tolower(paste0("vh", "_", input$vh_civge_drX,"_",input$vh_civge_drY,"_", input$vh_civge_inp2,".svg")) },
 content = function(file) {
   ggsave(
   file, device = "svg", bg = "white",
   plot = scDRgene(vhconf, vhmeta, input$vh_civge_drX, input$vh_civge_drY, input$vh_civge_inp2,  input$vh_civge_sub1, input$vh_civge_sub2, "vhgexpr.h5", vhgene, input$vh_civge_siz, input$vh_civge_col2, input$vh_civge_ord2, input$vh_civge_fsz, input$vh_civge_asp, input$vh_civge_txt)
   )
}) # End of tab civge



### Tab civci cell info vs cell info ----
  
  output$vh_civci_sub1.ui <- renderUI({
    sub = strsplit(vhconf[UI == input$vh_civci_sub1]$fID, "\\|")[[1]]
    checkboxGroupInput("vh_civci_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
  })
  observeEvent(input$vh_civci_sub1non, {
    sub = strsplit(vhconf[UI == input$vh_civci_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "vh_civci_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
  })
  observeEvent(input$vh_civci_sub1all, {
    sub = strsplit(vhconf[UI == input$vh_civci_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "vh_civci_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
  })

output$vh_civci_oup1 <- renderPlot({
  req(input$vh_civci_inp1)
  scDRcell(vhconf, vhmeta, input$vh_civci_drX, input$vh_civci_drY, input$vh_civci_inp1, input$vh_civci_sub1, input$vh_civci_sub2, input$vh_civci_siz, input$vh_civci_col1, input$vh_civci_ord1, input$vh_civci_fsz, input$vh_civci_asp, input$vh_civci_txt, input$vh_civci_lab1)
})

output$vh_civci_oup1.ui <- renderUI({
  show_progress(imageOutput("vh_civci_oup1", height = pList[input$vh_civci_psz]))
})

output$vh_civci_oup1.png <- downloadHandler(
  filename = function() { tolower(paste0("vh", "_", input$vh_civci_drX, "_", input$vh_civci_drY, "_", input$vh_civci_inp1, ".png")) },
  content = function(file) {
    ggsave(
    file, device = "png", dpi = input$vh_civci_oup1.res, bg = "white",
    plot = scDRcell(vhconf, vhmeta, input$vh_civci_drX, input$vh_civci_drY, input$vh_civci_inp1, input$vh_civci_sub1, input$vh_civci_sub2, input$vh_civci_siz, input$vh_civci_col1, input$vh_civci_ord1, input$vh_civci_fsz, input$vh_civci_asp, input$vh_civci_txt, input$vh_civci_lab1)
    )
})

output$vh_civci_oup1.pdf <- downloadHandler(
  filename = function() { tolower(paste0("vh", "_", input$vh_civci_drX, "_", input$vh_civci_drY, "_", input$vh_civci_inp1, ".pdf")) },
  content = function(file) { ggsave(
    file, device = "pdf", useDingbats = FALSE, bg = "white", onefile = TRUE,
    plot = scDRcell(vhconf, vhmeta, input$vh_civci_drX, input$vh_civci_drY, input$vh_civci_inp1, input$vh_civci_sub1, input$vh_civci_sub2, input$vh_civci_siz, input$vh_civci_col1, input$vh_civci_ord1, input$vh_civci_fsz, input$vh_civci_asp, input$vh_civci_txt, input$vh_civci_lab1) )
})

output$vh_civci_oup1.svg <- downloadHandler(
  filename = function() { tolower(paste0("vh", "_", input$vh_civci_drX, "_", input$vh_civci_drY, "_", input$vh_civci_inp1, ".svg")) },
  content = function(file) { ggsave(
    file, device = "svg", bg = "white",
    plot = scDRcell(vhconf, vhmeta, input$vh_civci_drX, input$vh_civci_drY, input$vh_civci_inp1, input$vh_civci_sub1, input$vh_civci_sub2, input$vh_civci_siz, input$vh_civci_col1, input$vh_civci_ord1, input$vh_civci_fsz, input$vh_civci_asp, input$vh_civci_txt, input$vh_civci_lab1) )
})

output$vh_civci_oup2 <- renderPlot({
  req(input$vh_civci_inp2)
  scDRcell(vhconf, vhmeta, input$vh_civci_drX, input$vh_civci_drY, input$vh_civci_inp2, input$vh_civci_sub1, input$vh_civci_sub2, input$vh_civci_siz, input$vh_civci_col2, input$vh_civci_ord2, input$vh_civci_fsz, input$vh_civci_asp, input$vh_civci_txt, input$vh_civci_lab2)
})

output$vh_civci_oup2.ui <- renderUI({
  show_progress(imageOutput("vh_civci_oup2", height = pList[input$vh_civci_psz]))
})

output$vh_civci_oup2.png <- downloadHandler(
  filename = function() { tolower(paste0("vh", "_", input$vh_civci_drX,"_",input$vh_civci_drY,"_", input$vh_civci_inp2,".png")) },
  content = function(file) {
    ggsave(
    file, device = "png", bg = "white", dpi = input$vh_civci_oup2.res,
    plot = scDRcell(vhconf, vhmeta, input$vh_civci_drX, input$vh_civci_drY, input$vh_civci_inp2, input$vh_civci_sub1, input$vh_civci_sub2, input$vh_civci_siz, input$vh_civci_col2, input$vh_civci_ord2, input$vh_civci_fsz, input$vh_civci_asp, input$vh_civci_txt, input$vh_civci_lab2)
    )
}) 

output$vh_civci_oup2.pdf <- downloadHandler(
  filename = function() { tolower(paste0("vh", "_", input$vh_civci_drX,"_",input$vh_civci_drY,"_", input$vh_civci_inp2,".pdf")) },
  content = function(file) {
    ggsave(
    file, device = "pdf", useDingbats = FALSE, bg = "white", onefile = TRUE,
    plot = scDRcell(vhconf, vhmeta, input$vh_civci_drX, input$vh_civci_drY, input$vh_civci_inp2, input$vh_civci_sub1, input$vh_civci_sub2, input$vh_civci_siz, input$vh_civci_col2, input$vh_civci_ord2, input$vh_civci_fsz, input$vh_civci_asp, input$vh_civci_txt, input$vh_civci_lab2) 
    )
})

output$vh_civci_oup2.svg <- downloadHandler(
  filename = function() { tolower(paste0("vh", "_", input$vh_civci_drX,"_",input$vh_civci_drY,"_", input$vh_civci_inp2,".svg")) },
  content = function(file) {
    ggsave(
    file, device = "svg", bg = "white",
    plot = scDRcell(vhconf, vhmeta, input$vh_civci_drX, input$vh_civci_drY, input$vh_civci_inp2, input$vh_civci_sub1, input$vh_civci_sub2, input$vh_civci_siz, input$vh_civci_col2, input$vh_civci_ord2, input$vh_civci_fsz, input$vh_civci_asp, input$vh_civci_txt, input$vh_civci_lab2) 
    )
}) # End of tab civci



### Tab gevge gene exp vs gene exp ----

  output$vh_gevge_sub1.ui <- renderUI({
    sub = strsplit(vhconf[UI == input$vh_gevge_sub1]$fID, "\\|")[[1]]
    checkboxGroupInput("vh_gevge_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
  })
  observeEvent(input$vh_gevge_sub1non, {
    sub = strsplit(vhconf[UI == input$vh_gevge_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "vh_gevge_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
  })
  observeEvent(input$vh_gevge_sub1all, {
    sub = strsplit(vhconf[UI == input$vh_gevge_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "vh_gevge_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
  })

output$vh_gevge_oup1 <- renderPlot({
  req(input$vh_gevge_inp1)
  scDRgene(vhconf, vhmeta, input$vh_gevge_drX, input$vh_gevge_drY, input$vh_gevge_inp1, input$vh_gevge_sub1, input$vh_gevge_sub2, "vhgexpr.h5", vhgene, input$vh_gevge_siz, input$vh_gevge_col1, input$vh_gevge_ord1, input$vh_gevge_fsz, input$vh_gevge_asp, input$vh_gevge_txt)
})

output$vh_gevge_oup1.ui <- renderUI({
  show_progress(imageOutput("vh_gevge_oup1", height = pList[input$vh_gevge_psz]))
})

output$vh_gevge_oup1.png <- downloadHandler(
  filename = function() { tolower(paste0("vh", "_", input$vh_gevge_drX, "_", input$vh_gevge_drY, "_", input$vh_gevge_inp1,".png")) },
  content = function(file) {
    ggsave(
    file, device = "png", bg = "white", dpi = input$vh_gevge_oup1.res,
    plot = scDRgene(vhconf, vhmeta, input$vh_gevge_drX, input$vh_gevge_drY, input$vh_gevge_inp1, 
                    input$vh_gevge_sub1, input$vh_gevge_sub2,
                    "vhgexpr.h5", vhgene,
                    input$vh_gevge_siz, input$vh_gevge_col1, input$vh_gevge_ord1,
                    input$vh_gevge_fsz, input$vh_gevge_asp, input$vh_gevge_txt)
    )
})

output$vh_gevge_oup1.pdf <- downloadHandler(
  filename = function() { tolower(paste0("vh", "_", input$vh_gevge_drX, "_", input$vh_gevge_drY, "_", input$vh_gevge_inp1, ".pdf")) },
  content = function(file) {
    ggsave(
    file, device = "pdf", useDingbats = FALSE, bg = "white", onefile = TRUE,
    plot = scDRgene(vhconf, vhmeta, input$vh_gevge_drX, input$vh_gevge_drY, input$vh_gevge_inp1, input$vh_gevge_sub1, input$vh_gevge_sub2, "vhgexpr.h5", vhgene, input$vh_gevge_siz, input$vh_gevge_col1, input$vh_gevge_ord1, input$vh_gevge_fsz, input$vh_gevge_asp, input$vh_gevge_txt)
    )
})

output$vh_gevge_oup1.svg <- downloadHandler(
  filename = function() { tolower(paste0("vh", "_", input$vh_gevge_drX, "_", input$vh_gevge_drY, "_", input$vh_gevge_inp1, ".svg")) },
  content = function(file) {
    ggsave(
    file, device = "svg", bg = "white",
    plot = scDRgene(vhconf, vhmeta, input$vh_gevge_drX, input$vh_gevge_drY, input$vh_gevge_inp1, input$vh_gevge_sub1, input$vh_gevge_sub2, "vhgexpr.h5", vhgene, input$vh_gevge_siz, input$vh_gevge_col1, input$vh_gevge_ord1, input$vh_gevge_fsz, input$vh_gevge_asp, input$vh_gevge_txt)
    )
})

output$vh_gevge_oup2 <- renderPlot({
  req(input$vh_gevge_inp2)
  scDRgene(vhconf, vhmeta, input$vh_gevge_drX, input$vh_gevge_drY, input$vh_gevge_inp2, input$vh_gevge_sub1, input$vh_gevge_sub2, "vhgexpr.h5", vhgene, input$vh_gevge_siz, input$vh_gevge_col2, input$vh_gevge_ord2, input$vh_gevge_fsz, input$vh_gevge_asp, input$vh_gevge_txt)
})

output$vh_gevge_oup2.ui <- renderUI({
  show_progress(imageOutput("vh_gevge_oup2", height = pList[input$vh_gevge_psz]))
})

output$vh_gevge_oup2.png <- downloadHandler(
  filename = function() { tolower(paste0("vh", "_", input$vh_gevge_drX, "_", input$vh_gevge_drY, "_", input$vh_gevge_inp2,".png")) },
  content = function(file) {
    ggsave(
    file, device = "png", bg = "white", dpi = input$vh_gevge_oup2.res,
    plot = scDRgene(vhconf, vhmeta, input$vh_gevge_drX, input$vh_gevge_drY, input$vh_gevge_inp2, input$vh_gevge_sub1, input$vh_gevge_sub2, "vhgexpr.h5", vhgene, input$vh_gevge_siz, input$vh_gevge_col2, input$vh_gevge_ord2, input$vh_gevge_fsz, input$vh_gevge_asp, input$vh_gevge_txt)
    )
}) 

output$vh_gevge_oup2.pdf <- downloadHandler(
  filename = function() { tolower(paste0("vh", "_", input$vh_gevge_drX, "_", input$vh_gevge_drY, "_", input$vh_gevge_inp2,".pdf")) },
  content = function(file) {
    ggsave(
    file, device = "pdf", useDingbats = FALSE, bg = "white", onefile = TRUE,
    plot = scDRgene(vhconf, vhmeta, input$vh_gevge_drX, input$vh_gevge_drY, input$vh_gevge_inp2, 
                    input$vh_gevge_sub1, input$vh_gevge_sub2,
                    "vhgexpr.h5", vhgene,
                    input$vh_gevge_siz, input$vh_gevge_col2, input$vh_gevge_ord2,
                    input$vh_gevge_fsz, input$vh_gevge_asp, input$vh_gevge_txt)
    )
})

output$vh_gevge_oup2.svg <- downloadHandler(
  filename = function() { tolower(paste0("vh", "_", input$vh_gevge_drX, "_", input$vh_gevge_drY, "_", input$vh_gevge_inp2,".svg")) },
  content = function(file) {
    ggsave(
    file, device = "svg", bg = "white",
    plot = scDRgene(vhconf, vhmeta, input$vh_gevge_drX, input$vh_gevge_drY, input$vh_gevge_inp2, 
                    input$vh_gevge_sub1, input$vh_gevge_sub2,
                    "vhgexpr.h5", vhgene,
                    input$vh_gevge_siz, input$vh_gevge_col2, input$vh_gevge_ord2,
                    input$vh_gevge_fsz, input$vh_gevge_asp, input$vh_gevge_txt)
    )
}) # End of tab gevge




### Tab gem gene expression multi ----

output$vh_gem_sub1.ui <- renderUI({
  sub = strsplit(vhconf[UI == input$vh_gem_sub1]$fID, "\\|")[[1]]
  checkboxGroupInput("vh_gem_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
})
observeEvent(input$vh_gem_sub1non, {
  sub = strsplit(vhconf[UI == input$vh_gem_sub1]$fID, "\\|")[[1]]
  updateCheckboxGroupInput(session, inputId = "vh_gem_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
})
observeEvent(input$vh_gem_sub1all, {
  sub = strsplit(vhconf[UI == input$vh_gem_sub1]$fID, "\\|")[[1]]
  updateCheckboxGroupInput(session, inputId = "vh_gem_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
})

output$vh_gem_oup1 <- renderPlot({
  req(input$vh_gem_inp)
  
  scFeature(vhconf, vhmeta, input$vh_gem_drX, input$vh_gem_drY, input$vh_gem_inp, input$vh_gem_sub1, input$vh_gem_sub2, "vhgexpr.h5", vhgene, input$vh_gem_siz, input$vh_gem_col, input$vh_gem_ord, input$vh_gem_fsz, input$vh_gem_asp, input$vh_gem_txt, input$vh_gem_ncol)
})

output$vh_gem_oup1.ui <- renderUI({
  show_progress(imageOutput("vh_gem_oup1", height = pList[input$vh_gem_psz]))
})

output$vh_gem_oup1.png <- downloadHandler(
  filename = function() { tolower(paste0("vh", "_", input$vh_gem_drX, "_", input$vh_gem_drY, "_expression.png")) },
  content = function(file) {
    ggsave(
      file, device = "png", height = input$vh_gem_oup1.height, width = input$vh_gem_oup1.width, dpi = input$vh_gem_oup1.res, units = "cm", bg = "white",
      plot = scFeature(vhconf, vhmeta, input$vh_gem_drX, input$vh_gem_drY, input$vh_gem_inp, input$vh_gem_sub1, input$vh_gem_sub2, "vhgexpr.h5", vhgene, input$vh_gem_siz, input$vh_gem_col, input$vh_gem_ord, input$vh_gem_fsz, input$vh_gem_asp, input$vh_gem_txt, input$vh_gem_ncol)
    )
}) 

output$vh_gem_oup1.pdf <- downloadHandler(
  filename = function() { tolower(paste0("vh", "_", input$vh_gem_drX, "_", input$vh_gem_drY, "_expression.pdf")) },
  content = function(file) {
  pdf(file, useDingbats = FALSE, height = input$vh_gem_oup1.height/2.54, width = input$vh_gem_oup1.width/2.54, bg = "white", onefile = TRUE)
  showtext::showtext_begin()
  print(scFeature(vhconf, vhmeta, input$vh_gem_drX, input$vh_gem_drY, input$vh_gem_inp, input$vh_gem_sub1, input$vh_gem_sub2, "vhgexpr.h5", vhgene, input$vh_gem_siz, input$vh_gem_col, input$vh_gem_ord, input$vh_gem_fsz, input$vh_gem_asp, input$vh_gem_txt, input$vh_gem_ncol))
  showtext::showtext_end()
  dev.off()
})

output$vh_gem_oup1.svg <- downloadHandler(
  filename = function() { tolower(paste0("vh", "_", input$vh_gem_drX, "_", input$vh_gem_drY, "_expression.svg")) },
  content = function(file) { ggsave(
    file, device = "svg", height = input$vh_gem_oup1.height, width = input$vh_gem_oup1.width, units = "cm", bg = "white",
    plot = scFeature(vhconf, vhmeta, input$vh_gem_drX, input$vh_gem_drY, input$vh_gem_inp, input$vh_gem_sub1, input$vh_gem_sub2, "vhgexpr.h5", vhgene, input$vh_gem_siz, input$vh_gem_col, input$vh_gem_ord, input$vh_gem_fsz, input$vh_gem_asp, input$vh_gem_txt, input$vh_gem_ncol))
}) # End of tab gem



### Tab gec gene co-expression ----

  output$vh_gec_sub1.ui <- renderUI({
    sub = strsplit(vhconf[UI == input$vh_gec_sub1]$fID, "\\|")[[1]]
    checkboxGroupInput("vh_gec_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
  })
  observeEvent(input$vh_gec_sub1non, {
    sub = strsplit(vhconf[UI == input$vh_gec_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "vh_gec_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
  })
  observeEvent(input$vh_gec_sub1all, {
    sub = strsplit(vhconf[UI == input$vh_gec_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "vh_gec_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
  })

output$vh_gec_oup1 <- renderPlot({
  scDRcoexFull(vhconf, vhmeta, input$vh_gec_drX, input$vh_gec_drY, input$vh_gec_inp1, input$vh_gec_inp2, input$vh_gec_sub1, input$vh_gec_sub2, "vhgexpr.h5", vhgene, input$vh_gec_siz, input$vh_gec_col1, input$vh_gec_ord1, input$vh_gec_fsz, input$vh_gec_asp, input$vh_gec_txt)
})

output$vh_gec_oup1.ui <- renderUI({
  show_progress(imageOutput("vh_gec_oup1", height = pList2[input$vh_gec_psz]))
})

output$vh_gec_oup1.png <- downloadHandler(
  filename = function() { tolower(paste0("vh", "_", input$vh_gec_drX, "_", input$vh_gec_drY, "_", input$vh_gec_inp1, "_", input$vh_gec_inp2, ".png")) },
  content = function(file) { ggsave(
    file, device = "png", bg = "white", dpi = input$vh_gec_oup1.res,
    plot = scDRcoexFull(vhconf, vhmeta, input$vh_gec_drX, input$vh_gec_drY, input$vh_gec_inp1, input$vh_gec_inp2, input$vh_gec_sub1, input$vh_gec_sub2, "vhgexpr.h5", vhgene, input$vh_gec_siz, input$vh_gec_col1, input$vh_gec_ord1, input$vh_gec_fsz, input$vh_gec_asp, input$vh_gec_txt) )
})

output$vh_gec_oup1.pdf <- downloadHandler(
  filename = function() { tolower(paste0("vh", "_", input$vh_gec_drX, "_", input$vh_gec_drY, "_", input$vh_gec_inp1, "_", input$vh_gec_inp2, ".pdf")) },
  content = function(file) {
    ggsave(
    file, device = "pdf", useDingbats = FALSE, bg = "white", onefile = TRUE,
    plot = scDRcoexFull(vhconf, vhmeta, input$vh_gec_drX, input$vh_gec_drY, input$vh_gec_inp1, input$vh_gec_inp2, input$vh_gec_sub1, input$vh_gec_sub2, "vhgexpr.h5", vhgene, input$vh_gec_siz, input$vh_gec_col1, input$vh_gec_ord1, input$vh_gec_fsz, input$vh_gec_asp, input$vh_gec_txt)
    )
})

output$vh_gec_oup1.svg <- downloadHandler(
  filename = function() { tolower(paste0("vh", "_", input$vh_gec_drX, "_", input$vh_gec_drY, "_", input$vh_gec_inp1, "_", input$vh_gec_inp2, ".svg")) },
  content = function(file) {
    ggsave(
    file, device = "svg", bg = "white",
    plot = scDRcoexFull(vhconf, vhmeta, input$vh_gec_drX, input$vh_gec_drY, input$vh_gec_inp1, input$vh_gec_inp2, input$vh_gec_sub1, input$vh_gec_sub2, "vhgexpr.h5", vhgene, input$vh_gec_siz, input$vh_gec_col1, input$vh_gec_ord1, input$vh_gec_fsz, input$vh_gec_asp, input$vh_gec_txt)
    )
})

output$vh_gec_.dt <- renderDataTable({
  ggData = scDRcoexNum(vhconf, vhmeta, input$vh_gec_inp1, input$vh_gec_inp2, input$vh_gec_sub1, input$vh_gec_sub2, "vhgexpr.h5", vhgene)
  datatable(ggData, rownames = FALSE, extensions = "Buttons", options = list(pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel"))) %>%
            formatRound(columns = c("percent"), digits = 2)
}) # End of tab gec



### Tab vio violinplot / boxplot / lineplot ----

  output$vh_vio_sub1.ui <- renderUI({
    sub = strsplit(vhconf[UI == input$vh_vio_sub1]$fID, "\\|")[[1]]
    checkboxGroupInput("vh_vio_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
  })
  observeEvent(input$vh_vio_sub1non, {
    sub = strsplit(vhconf[UI == input$vh_vio_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "vh_vio_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
  })
  observeEvent(input$vh_vio_sub1all, {
    sub = strsplit(vhconf[UI == input$vh_vio_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "vh_vio_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
  })

output$vh_vio_oup <- renderPlot({
  gh5 <- ifelse(input$vh_vio_datatype == "normalised","vhgexpr.h5","vhgexpr2.h5")
  scVioBox(vhconf, vhmeta, input$vh_vio_inp1, input$vh_vio_inp2, input$vh_vio_sub1, input$vh_vio_sub2, gh5, vhgene, input$vh_vio_typ, input$vh_vio_pts, input$vh_vio_siz, input$vh_vio_fsz, input$vh_vio_barsz)
})

output$vh_vio_oup.ui <- renderUI({
  show_progress(imageOutput("vh_vio_oup", height = pList2[input$vh_vio_psz]))
})

output$vh_vio_oup.png <- downloadHandler(
  filename = function() { tolower(paste0("vh", "_", input$vh_vio_typ, "_", input$vh_vio_datatype, "_", input$vh_vio_inp1, "_", input$vh_vio_inp2,".png")) },
  content = function(file) {
    gh5 <- ifelse(input$vh_vio_datatype == "normalised","vhgexpr.h5","vhgexpr2.h5")
    ggsave(
    file, device = "png", bg = "white", dpi = input$vh_vio_oup.res,
    plot = scVioBox(vhconf, vhmeta, input$vh_vio_inp1, input$vh_vio_inp2, input$vh_vio_sub1, input$vh_vio_sub2, gh5, vhgene, input$vh_vio_typ, input$vh_vio_pts, input$vh_vio_siz, input$vh_vio_fsz, input$vh_vio_barsz)
    )
}) 

output$vh_vio_oup.pdf <- downloadHandler(
  filename = function() { tolower(paste0("vh", "_", input$vh_vio_typ, "_", input$vh_vio_datatype, "_", input$vh_vio_inp1, "_", input$vh_vio_inp2, ".pdf")) },
  content = function(file) {
    gh5 <- ifelse(input$vh_vio_datatype == "normalised","vhgexpr.h5","vhgexpr2.h5")
    ggsave(
    file, device = "pdf", useDingbats = FALSE, bg = "white", onefile = TRUE,
    plot = scVioBox(vhconf, vhmeta, input$vh_vio_inp1, input$vh_vio_inp2, input$vh_vio_sub1, input$vh_vio_sub2, gh5, vhgene, input$vh_vio_typ, input$vh_vio_pts, input$vh_vio_siz, input$vh_vio_fsz, input$vh_vio_barsz)
    )
})

output$vh_vio_oup.svg <- downloadHandler(
  filename = function() { tolower(paste0("vh", "_", input$vh_vio_typ, "_", input$vh_vio_datatype, "_", input$vh_vio_inp1, "_", input$vh_vio_inp2, ".svg")) },
  content = function(file) {
    gh5 <- ifelse(input$vh_vio_datatype == "normalised","vhgexpr.h5","vhgexpr2.h5")
    ggsave(
    file, device = "svg", bg = "white",
    plot = scVioBox(vhconf, vhmeta, input$vh_vio_inp1, input$vh_vio_inp2, input$vh_vio_sub1, input$vh_vio_sub2, gh5, vhgene, input$vh_vio_typ, input$vh_vio_pts, input$vh_vio_siz, input$vh_vio_fsz, input$vh_vio_barsz)
    )
}) # End of tab vio




### Tab pro proportion plot ----

  output$vh_pro_sub1.ui <- renderUI({
    sub = strsplit(vhconf[UI == input$vh_pro_sub1]$fID, "\\|")[[1]]
    checkboxGroupInput("vh_pro_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
  })
  observeEvent(input$vh_pro_sub1non, {
    sub = strsplit(vhconf[UI == input$vh_pro_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "vh_pro_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
  })
  observeEvent(input$vh_pro_sub1all, {
    sub = strsplit(vhconf[UI == input$vh_pro_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "vh_pro_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
  })

output$vh_pro_oup <- renderPlot({
  scProp(vhconf, vhmeta, input$vh_pro_inp1, input$vh_pro_inp2, input$vh_pro_sub1, input$vh_pro_sub2, input$vh_pro_typ, input$vh_pro_flp, input$vh_pro_fsz)
})

output$vh_pro_oup.ui <- renderUI({
  show_progress(imageOutput("vh_pro_oup", height = pList2[input$vh_pro_psz]))
})

output$vh_pro_oup.png <- downloadHandler(
  filename = function() { tolower(paste0("vh", "_", input$vh_pro_typ, "_", input$vh_pro_inp1, "_", input$vh_pro_inp2, ".png")) },
  content = function(file) {
    ggsave(
    file, device = "png", bg = "white", dpi = input$vh_pro_oup.res,
    plot = scProp(vhconf, vhmeta, input$vh_pro_inp1, input$vh_pro_inp2, input$vh_pro_sub1, input$vh_pro_sub2, input$vh_pro_typ, input$vh_pro_flp, input$vh_pro_fsz)
    )
  }) 
  
output$vh_pro_oup.pdf <- downloadHandler(
  filename = function() { tolower(paste0("vh", "_", input$vh_pro_typ, "_", input$vh_pro_inp1, "_", input$vh_pro_inp2, ".pdf")) },
  content = function(file) {
    ggsave(
    file, device = "pdf", useDingbats = FALSE, bg = "white", onefile = TRUE,
    plot = scProp(vhconf, vhmeta, input$vh_pro_inp1, input$vh_pro_inp2, input$vh_pro_sub1, input$vh_pro_sub2, input$vh_pro_typ, input$vh_pro_flp, input$vh_pro_fsz)
    )
  })
  
output$vh_pro_oup.svg <- downloadHandler(
  filename = function() { tolower(paste0("vh", "_", input$vh_pro_typ, "_", input$vh_pro_inp1, "_", input$vh_pro_inp2, ".svg")) },
  content = function(file) {
    ggsave(
    file, device = "svg", bg = "white",
    plot = scProp(vhconf, vhmeta, input$vh_pro_inp1, input$vh_pro_inp2, input$vh_pro_sub1, input$vh_pro_sub2, input$vh_pro_typ, input$vh_pro_flp, input$vh_pro_fsz)
    )
  }) # End of tab pro



### Tab hea heatmap / dotplot ----

  output$vh_hea_sub1.ui <- renderUI({
    sub = strsplit(vhconf[UI == input$vh_hea_sub1]$fID, "\\|")[[1]]
    checkboxGroupInput("vh_hea_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
  })
  observeEvent(input$vh_hea_sub1non, {
    sub = strsplit(vhconf[UI == input$vh_hea_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "vh_hea_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
  })
  observeEvent(input$vh_hea_sub1all, {
    sub = strsplit(vhconf[UI == input$vh_hea_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "vh_hea_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
  })

output$vh_hea_oupTxt <- renderUI({
  geneList = scGeneList(input$vh_hea_inp, vhgene)
  if(nrow(geneList) > 50){
    HTML("More than 50 input genes! Please reduce the gene list!")
  }
})

output$vh_hea_oup <- renderPlot({
  scBubbHeat(vhconf, vhmeta, input$vh_hea_inp, input$vh_hea_grp, input$vh_hea_plt, input$vh_hea_sub1, input$vh_hea_sub2, "vhgexpr.h5", vhgene, input$vh_hea_scl, input$vh_hea_row, input$vh_hea_col, input$vh_hea_cols, input$vh_hea_fsz)
})

output$vh_hea_oup.ui <- renderUI({
  show_progress(imageOutput("vh_hea_oup", height = pList3[input$vh_hea_psz]))
})

output$vh_hea_oup.png <- downloadHandler(
  filename = function() { tolower(paste0("vh", "_", input$vh_hea_plt,"_",input$vh_hea_grp,".png")) },
  content = function(file) {
    ggsave(
    file, device = "png", bg = "white", height = input$vh_hea_oup.height, width = input$vh_hea_oup.width, units = "cm", dpi = input$vh_hea_oup.res, plot = scBubbHeat(vhconf, vhmeta, input$vh_hea_inp, input$vh_hea_grp, input$vh_hea_plt, input$vh_hea_sub1, input$vh_hea_sub2, "vhgexpr.h5", vhgene, input$vh_hea_scl, input$vh_hea_row, input$vh_hea_col, input$vh_hea_cols, input$vh_hea_fsz)
    )
})

output$vh_hea_oup.pdf <- downloadHandler(
  filename = function() { tolower(paste0("vh", "_", input$vh_hea_plt,"_",input$vh_hea_grp,".pdf")) },
  content = function(file) {
    pdf(file, useDingbats = FALSE, bg = "white", height = input$vh_hea_oup.height/2.54, width = input$vh_hea_oup.width/2.54, onefile = TRUE)
    showtext::showtext_begin()
    print(scBubbHeat(vhconf, vhmeta, input$vh_hea_inp, input$vh_hea_grp, input$vh_hea_plt, input$vh_hea_sub1, input$vh_hea_sub2, "vhgexpr.h5", vhgene, input$vh_hea_scl, input$vh_hea_row, input$vh_hea_col, input$vh_hea_cols, input$vh_hea_fsz))
    showtext::showtext_end()
    dev.off()
})

output$vh_hea_oup.svg <- downloadHandler(
  filename = function() { tolower(paste0("vh", "_", input$vh_hea_plt,"_",input$vh_hea_grp,".svg")) },
  content = function(file) {
    ggsave(
    file, device = "svg", bg = "white", height = input$vh_hea_oup.height, width = input$vh_hea_oup.width, units = "cm", 
    plot = scBubbHeat(vhconf, vhmeta, input$vh_hea_inp, input$vh_hea_grp, input$vh_hea_plt, input$vh_hea_sub1, input$vh_hea_sub2, "vhgexpr.h5", vhgene, input$vh_hea_scl, input$vh_hea_row, input$vh_hea_col, input$vh_hea_cols, input$vh_hea_fsz)
    )
}) # End of tab hea      
       


### Tab markers ----

output$vh_mar_table <- renderDataTable({
  req(input$vh_mar_cls)
  datatable(vhmar[[input$vh_mar_cls]], rownames = FALSE, extensions = "Buttons", options = list(dom = "lftiprB", buttons = c("copy", "csv", "excel")))
}) # End of tab mar

})


