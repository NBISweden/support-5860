library(shiny)
library(shinyhelper)
library(shinythemes)
library(data.table)
library(Matrix)
library(DT)
library(magrittr)

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

### UI code
shinyUI(
fluidPage(style="margin:0;padding:0;",
tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")),

theme = shinythemes::shinytheme("flatly"),
navbarPage(
"Support 8560"
,navbarMenu("LEC Combined",
# tab civge ----
tabPanel(
  "CellInfo vs GeneExpr",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Cell information vs Gene expression"),
          p("Cell information and gene expression side-by-side on low-dimensional represention.")
        )
      ),
      # row 2 ----
      fluidRow(
        class = "tab-section",
        # row 2 col 1
        column(
          4,
          div(
            class = "input-panel",
            h4("Dimension Reduction"),
            selectInput("lc_civge_drX", "X-axis:",
              choices = lcconf[dimred == TRUE]$UI,
              selected = lcdef$dimred[1]
            ),
            selectInput("lc_civge_drY", "Y-axis:",
              choices = lcconf[dimred == TRUE]$UI,
              selected = lcdef$dimred[2]
            )
          )
        ), # row 1 col 1
        # row 2 col 2
        column(
          4,
          div(
            class = "input-panel",
            checkboxInput("lc_civge_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.lc_civge_togL == true",
              selectInput("lc_civge_sub1", "Cell info to subset:",
                choices = lcconf[grp == TRUE]$UI,
                selected = lcdef$grp1
              ),
              uiOutput("lc_civge_sub1.ui"),
              actionButton("lc_civge_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("lc_civge_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            )
          )
        ), # End of column
        # row 2 col 3
        column(
          4,
          div(
            class = "input-panel",
            checkboxInput("lc_civge_tog0", "Adjust graphics"),
            conditionalPanel(
              condition = "input.lc_civge_tog0 == true",
              sliderInput("lc_civge_siz", "Point size:",
                min = 0, max = 4, value = 1.25, step = 0.25
              ),
              radioButtons("lc_civge_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("lc_civge_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              radioButtons("lc_civge_asp", "Aspect ratio:",
                choices = c("Square", "Fixed", "Free"),
                selected = "Square", inline = TRUE
              ),
              checkboxInput("lc_civge_txt", "Show axis text", value = FALSE)
            )
          )
        ) # row 2 col 3
      ), # row 2
      # row 3 ----
      fluidRow(
        class = "tab-section",
        # row 3 col 1
        column(6,
          style = "border-right: 2px solid #f3f6f4",
          h4("Cell information"),
          # row 3 col 1 row 1
          fluidRow(
            class = "tab-section",
            column(
              6,
              div(
                class = "input-panel",
                selectInput("lc_civge_inp1", "Cell info:",
                  choices = lcconf$UI,
                  selected = lcdef$meta1
                ) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Cell info to colour cells",
                    content = c(
                      "- Select cell info to colour cells",
                      "- Categorical covariates have a fixed colour palette",
                      "- Continuous covariates are coloured in a Blue-Yellow-Red colour scheme, which can be changed in the plot controls"
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("lc_civge_tog1", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.lc_civge_tog1 == true",
                  radioButtons("lc_civge_col1", "Colour (Continuous data):",
                    choices = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple"),
                    selected = "Blue-Yellow-Red"
                  ),
                  radioButtons("lc_civge_ord1", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Original", inline = TRUE
                  ),
                  checkboxInput("lc_civge_lab1", "Show cell info labels", value = TRUE)
                )
              )
            )
          ),
          # row 3 col 1 row 2
          fluidRow(
            class = "tab-section",
            column(
              12,
              uiOutput("lc_civge_oup1.ui")
            )
          ),
          # row 3 col 1 row 3
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("lc_civge_oup1.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("lc_civge_oup1.png", "Download PNG", class = "btn-sm"),
                downloadButton("lc_civge_oup1.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("lc_civge_oup1.svg", "Download SVG", class = "btn-sm"),
                checkboxInput("lc_civge_tog9", "Show cell numbers / statistics")
              )
            )
          ),
          conditionalPanel(
            condition = "input.lc_civge_tog9 == true",
            h4("Cell numbers / statistics"),
            radioButtons("lc_civge_splt", "Split continuous cell info into:",
              choices = c("Quartile", "Decile"),
              selected = "Decile", inline = TRUE
            ),
            dataTableOutput("lc_civge_.dt")
          )
        ), # row 3 col 1
        # row 3 col 2
        column(
          6,
          h4("Gene expression"),
          fluidRow(
            class = "tab-section",
            column(
              6,
              div(
                class = "input-panel",
                selectInput("lc_civge_inp2", "Gene name:", choices = NULL) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Gene expression to colour cells",
                    content = c(
                      "- Select gene to colour cells by gene expression",
                      "- Type in gene names for unlisted genes",
                      "- Gene expression are coloured in a White-Red colour scheme which can be changed in the plot controls"
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("lc_civge_tog2", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.lc_civge_tog2 == true",
                  radioButtons("lc_civge_col2", "Colour:",
                    choices = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple"),
                    selected = "White-Red"
                  ),
                  radioButtons("lc_civge_ord2", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Max", inline = TRUE
                  )
                )
              )
            )
          ),
          fluidRow(
            class = "tab-section",
            column(
              12,
              uiOutput("lc_civge_oup2.ui")
            )
          ),
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("lc_civge_oup2.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("lc_civge_oup2.png", "Download PNG", class = "btn-sm"),
                downloadButton("lc_civge_oup2.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("lc_civge_oup2.svg", "Download SVG", class = "btn-sm")
              )
            )
          )
        ) # row 3 col 2
      ), # row 3
      hr()
    )
  )
) # End of tab civge

,
# tab civci ----
tabPanel(
  "CellInfo vs CellInfo",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Cell info vs cell info"),
          p("Two cell infos side-by-side on low-dimensional represention.")
        ) # row 1 col 1
      ), # row 1
      # row 2 ----
      fluidRow(
        class = "tab-section",
        column(
          4,
          div(
            class = "input-panel",
            h4("Dimension Reduction"),
            selectInput("lc_civci_drX", "X-axis:",
              choices = lcconf[dimred == TRUE]$UI,
              selected = lcdef$dimred[1]
            ),
            selectInput("lc_civci_drY", "Y-axis:",
              choices = lcconf[dimred == TRUE]$UI,
              selected = lcdef$dimred[2]
            )
          )
        ), # row 2 col 2
        # row 2 col 2
        column(4,
          div(
            class = "input-panel",
            checkboxInput("lc_civci_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.lc_civci_togL == true",
              selectInput("lc_civci_sub1", "Cell info to subset:",
                choices = lcconf[grp == TRUE]$UI,
                selected = lcdef$grp1
              ),
              uiOutput("lc_civci_sub1.ui"),
              actionButton("lc_civci_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("lc_civci_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            )
          )
        ), # row 2 col 2
        # row 2 col 3
        column(
          4,
          div(
            class = "input-panel",
            checkboxInput("lc_civci_tog0", "Adjust graphics"),
            conditionalPanel(
              condition = "input.lc_civci_tog0 == true",
              sliderInput("lc_civci_siz", "Point size:",
                min = 0, max = 4, value = 1.25, step = 0.25
              ),
              radioButtons("lc_civci_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("lc_civci_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              radioButtons("lc_civci_asp", "Aspect ratio:",
                choices = c("Square", "Fixed", "Free"),
                selected = "Square", inline = TRUE
              ),
              checkboxInput("lc_civci_txt", "Show axis text", value = FALSE)
            )
          )
        ) # row 2 col 3
      ), # row 2
      # row 3 ----
      fluidRow(
        class = "tab-section",
        # row 3 col 1
        column(
          6,
          style = "border-right: 2px solid #f3f6f4",
          h4("Cell info 1"),
          # row 3 col 1 row 1
          fluidRow(
            class = "tab-section",
            column(
              6,
              div(
                class = "input-panel",
                selectInput("lc_civci_inp1", "Cell info:",
                  choices = lcconf$UI,
                  selected = lcdef$meta1
                ) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Cell info to colour cells",
                    content = c(
                      "- Select cell info to colour cells",
                      "- Categorical covariates have a fixed colour palette",
                      "- Continuous covariates are coloured in a Blue-Yellow-Red colour scheme, which can be changed in the plot controls"
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("lc_civci_tog1", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.lc_civci_tog1 == true",
                  radioButtons("lc_civci_col1", "Colour (Continuous data):",
                    choices = c(
                      "White-Red", "Blue-Yellow-Red",
                      "Yellow-Green-Purple"
                    ),
                    selected = "Blue-Yellow-Red"
                  ),
                  radioButtons("lc_civci_ord1", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Original", inline = TRUE
                  ),
                  checkboxInput("lc_civci_lab1", "Show cell info labels", value = TRUE)
                )
              )
            )
          ),
          # row 3 col 1 row 2
          fluidRow(column(12, uiOutput("lc_civci_oup1.ui"))),
          # row 3 col 1 row 3
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("lc_civci_oup1.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("lc_civci_oup1.png", "Download PNG", class = "btn-sm"),
                downloadButton("lc_civci_oup1.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("lc_civci_oup1.svg", "Download svg", class = "btn-sm")
              )
            )
          )
        ), # row 3 col 1
        # row 3 col 2
        column(
          6, h4("Cell info 2"),
          fluidRow(
            class = "tab-section",
            column(
              6,
              div(
                class = "input-panel",
                selectInput("lc_civci_inp2", "Cell info:",
                  choices = lcconf$UI,
                  selected = lcdef$meta2
                ) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Cell info to colour cells",
                    content = c(
                      "- Select cell info to colour cells",
                      "- Categorical covariates have a fixed colour palette",
                      "- Continuous covariates are coloured in a Blue-Yellow-Red colour scheme, which can be changed in the plot controls"
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("lc_civci_tog2", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.lc_civci_tog2 == true",
                  radioButtons("lc_civci_col2", "Colour (Continuous data):",
                    choices = c(
                      "White-Red", "Blue-Yellow-Red",
                      "Yellow-Green-Purple"
                    ),
                    selected = "Blue-Yellow-Red"
                  ),
                  radioButtons("lc_civci_ord2", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Original", inline = TRUE
                  ),
                  checkboxInput("lc_civci_lab2", "Show cell info labels", value = TRUE)
                )
              )
            )
          ),
          fluidRow(column(12, uiOutput("lc_civci_oup2.ui"))),
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("lc_civci_oup2.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("lc_civci_oup2.png", "Download PNG", class = "btn-sm"),
                downloadButton("lc_civci_oup2.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("lc_civci_oup2.svg", "Download SVG", class = "btn-sm")
              )
            )
          )
        ) # row 3 col 2
      ), # row 3
      hr()
    )
  )
) # End of tab civci

,
# tab gevge ----
tabPanel(
  "GeneExpr vs GeneExpr",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Gene expression vs Gene expression"),
          p("Visualise two gene expressions side-by-side on low-dimensional representions.")
        )
      ),
      # row 2 ----
      fluidRow(
        class = "tab-section",
        # row 2 col 1
        column(
          4,
          div(
            class = "input-panel",
            h4("Dimension Reduction"),
            selectInput("lc_gevge_drX", "X-axis:",
              choices = lcconf[dimred == TRUE]$UI,
              selected = lcdef$dimred[1]
            ),
            selectInput("lc_gevge_drY", "Y-axis:",
              choices = lcconf[dimred == TRUE]$UI,
              selected = lcdef$dimred[2]
            )
          )
        ), # row 1 col 1
        # row 2 col 2
        column(4,
          div(
            class = "input-panel",
            checkboxInput("lc_gevge_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.lc_gevge_togL == true",
              selectInput("lc_gevge_sub1", "Cell info to subset:",
                choices = lcconf[grp == TRUE]$UI,
                selected = lcdef$grp1
              ),
              uiOutput("lc_gevge_sub1.ui"),
              actionButton("lc_gevge_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("lc_gevge_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            )
          )
        ), # End of column
        # row 2 col 3
        column(
          4,
          div(
            class = "input-panel",
            checkboxInput("lc_gevge_tog0", "Adjust graphics"),
            conditionalPanel(
              condition = "input.lc_gevge_tog0 == true",
              sliderInput("lc_gevge_siz", "Point size:",
                min = 0, max = 4, value = 1.25, step = 0.25
              ),
              radioButtons("lc_gevge_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("lc_gevge_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              radioButtons("lc_gevge_asp", "Aspect ratio:",
                choices = c("Square", "Fixed", "Free"),
                selected = "Square", inline = TRUE
              ),
              checkboxInput("lc_gevge_txt", "Show axis text", value = FALSE)
            )
          )
        ) # row 2 col 3
      ), # row 2
      # row 3 ----
      fluidRow(
        class = "tab-section",
        column(
          6,
          style = "border-right: 2px solid #f3f6f4",
          h4("Gene expression 1"),
          fluidRow(
            class = "tab-section",
            column(
              6,
              div(
                class = "input-panel",
                selectInput("lc_gevge_inp1", "Gene name:", choices = NULL) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Gene expression to colour cells",
                    content = c(
                      "- Select gene to colour cells by gene expression",
                      "- Type in gene names for unlisted genes",
                      "- Gene expression are coloured in a White-Red colour scheme which can be changed in the plot controls"
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("lc_gevge_tog1", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.lc_gevge_tog1 == true",
                  radioButtons("lc_gevge_col1", "Colour:",
                    choices = c(
                      "White-Red", "Blue-Yellow-Red",
                      "Yellow-Green-Purple"
                    ),
                    selected = "White-Red"
                  ),
                  radioButtons("lc_gevge_ord1", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Max", inline = TRUE
                  )
                )
              )
            )
          ),
          # row 3 col 1 row 2
          fluidRow(
            class = "tab-section",
            column(12, uiOutput("lc_gevge_oup1.ui"))
          ),
          # row 3 col 1 row 3
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("lc_gevge_oup1.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("lc_gevge_oup1.png", "Download PNG", class = "btn-sm"),
                downloadButton("lc_gevge_oup1.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("lc_gevge_oup1.svg", "Download SVG", class = "btn-sm")
              )
            )
          )
        ), # row 3 col 1
        # row 3 col 2
        column(
          6, h4("Gene expression 2"),
          fluidRow(
            class = "tab-section",
            column(
              6,
              div(
                class = "input-panel",
                selectInput("lc_gevge_inp2", "Gene name:", choices = NULL) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Gene expression to colour cells",
                    content = c(
                      "- Select gene to colour cells by gene expression",
                      "- Type in gene names for unlisted genes",
                      "- Gene expression are coloured in a White-Red colour scheme which can be changed in the plot controls"
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("lc_gevge_tog2", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.lc_gevge_tog2 == true",
                  radioButtons("lc_gevge_col2", "Colour:",
                    choices = c(
                      "White-Red", "Blue-Yellow-Red",
                      "Yellow-Green-Purple"
                    ),
                    selected = "White-Red"
                  ),
                  radioButtons("lc_gevge_ord2", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Max", inline = TRUE
                  )
                )
              )
            )
          ),
          fluidRow(
            class = "tab-section",
            column(12, uiOutput("lc_gevge_oup2.ui"))
          ),
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                  numericInput("lc_gevge_oup2.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                  downloadButton("lc_gevge_oup2.png", "Download PNG", class = "btn-sm"),
                  downloadButton("lc_gevge_oup2.pdf", "Download PDF", class = "btn-sm"),
                  downloadButton("lc_gevge_oup2.svg", "Download SVG", class = "btn-sm")
              )
            )
          )
        ) # row 3 col 2
      ), # row 3
      hr()
    )
  )
) # End of tab gevge

,
# tab gem ----
tabPanel(
  "Expression",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Gene expression"),
          p("Explore gene expression on low-dimensional represention.")
        ) # row 1 col 1
      ), # row 1
      # row 2 ----
      fluidRow(
        class = "tab-section",
        column(4,
               fluidRow(
        column(
          12,
          div(
            class = "input-panel input-panel-section",
            selectInput("lc_gem_inp", "Genes:", choices = NULL, multiple = TRUE) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "List of genes to plot on bubbleplot / heatmap",
                content = c(
                  "- Input genes to plot",
                  "- Type in gene names for unlisted genes"
                )
              ),
            selectInput("lc_gem_drX", "X-axis:",
                        choices = lcconf[dimred == TRUE]$UI,
                        selected = lcdef$dimred[1]
            ),
            selectInput("lc_gem_drY", "Y-axis:",
                        choices = lcconf[dimred == TRUE]$UI,
                        selected = lcdef$dimred[2]
            ),
            checkboxInput("lc_gem_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.lc_gem_togL == true",
              selectInput("lc_gem_sub1", "Cell info to subset:",
                          choices = lcconf[grp == TRUE]$UI,
                          selected = lcdef$grp1
              ),
              uiOutput("lc_gem_sub1.ui"),
              actionButton("lc_gem_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("lc_gem_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            ),
            checkboxInput("lc_gem_tog0", "Adjust graphics"),
            conditionalPanel(
              condition = "input.lc_gem_tog0 == true",
              sliderInput("lc_gem_siz", "Point size:",
                          min = 0, max = 3, value = 0.5, step = 0.1
              ),
              radioButtons("lc_gem_psz", "Plot size:",
                           choices = c("Small", "Medium", "Large"),
                           selected = "Medium", inline = TRUE
              ),
              radioButtons("lc_gem_fsz", "Font size:",
                           choices = c("Smaller", "Small", "Medium", "Large"),
                           selected = "Small", inline = TRUE
              ),
              radioButtons("lc_gem_asp", "Aspect ratio:",
                           choices = c("Square", "Fixed", "Free"),
                           selected = "Square", inline = TRUE
              ),
              checkboxInput("lc_gem_txt", "Show axis text", value = FALSE),
              radioButtons("lc_gem_col", "Colour (Continuous data):",
                           choices = c(
                             "White-Red", "Blue-Yellow-Red",
                             "Yellow-Green-Purple"
                           ),
                           selected = "Blue-Yellow-Red"
              ),
              radioButtons("lc_gem_ord", "Plot order:",
                           choices = c("Max", "Min", "Original", "Random"),
                           selected = "Max", inline = TRUE
              ),
              numericInput("lc_gem_ncol", "Number of columns", value = 0, min = 0, step = 1)
            )
          ),
          div(
            class = "input-panel",
            fluidRow(
            column(4,
              numericInput("lc_gem_oup1.height", "Height:", min = 1, max = 50, value = 25, step = 2)
            ),
            column(4,
              numericInput("lc_gem_oup1.width", "Width:", min = 1, max = 50, value = 25, step = 2)
            ),
            column(4,
              numericInput("lc_gem_oup1.res", "Res:", min = 72, max = 600, value = 150, step = 5)
            )
            ),
            downloadButton("lc_gem_oup1.png", "Download PNG", class = "btn-sm"),
            downloadButton("lc_gem_oup1.pdf", "Download PDF", class = "btn-sm"),
            downloadButton("lc_gem_oup1.svg", "Download SVG", class = "btn-sm")
          )
        )
               )
      ),
      column(8,
             uiOutput("lc_gem_oup1.ui")
      )
      ),
      hr()
    )
  )
) # End of tab gem

,
# tab gec ----
tabPanel(
  "Gene coexpression",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Coexpression of two genes on reduced dimensions"),
          p("Visualise the coexpression of two genes on low-dimensional representions.")
        )
      ),
      # row 2 ----
      fluidRow(
        class = "tab-section",
        # row 2 col 1
        column(4,
               column(
                 12,
                 div(
                   class = "input-panel input-panel-section",
                   h4("Dimension Reduction"),
                   selectInput("lc_gec_drX", "X-axis:",
                               choices = lcconf[dimred == TRUE]$UI,
                               selected = lcdef$dimred[1]
                   ),
                   selectInput("lc_gec_drY", "Y-axis:",
                               choices = lcconf[dimred == TRUE]$UI,
                               selected = lcdef$dimred[2]
                   ),
                   selectInput("lc_gec_inp1", "Gene 1:", choices = NULL) %>%
                     helper(
                       type = "inline", size = "m", fade = TRUE,
                       title = "Gene expression to colour cells",
                       content = c(
                         "- Select gene to colour cells by gene expression",
                         "- Type in gene names for unlisted genes",
                         "- Gene expression are coloured in a White-Red colour scheme which can be changed in the plot controls"
                       )
                     ),
                   selectInput("lc_gec_inp2", "Gene 2:", choices = NULL) %>%
                     helper(
                       type = "inline", size = "m", fade = TRUE,
                       title = "Gene expression to colour cells",
                       content = c(
                         "- Select gene to colour cells by gene expression",
                         "- Type in gene names for unlisted genes",
                         "- Gene expression are coloured in a White-Blue colour scheme which can be changed in the plot controls"
                       )
                     ),
                   checkboxInput("lc_gec_togL", "Subset cells"),
                   conditionalPanel(
                    condition = "input.lc_gec_togL == true",
                    selectInput("lc_gec_sub1", "Cell info to subset:", choices = lcconf[grp == TRUE]$UI, selected = lcdef$grp1),
                     uiOutput("lc_gec_sub1.ui"),
                     actionButton("lc_gec_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
                     actionButton("lc_gec_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
                   ),
                   checkboxInput("lc_gec_tog0", "Adjust graphics"),
                   conditionalPanel(
                     condition = "input.lc_gec_tog0 == true",
                     radioButtons("lc_gec_col1", "Colour:",
                                  choices = c(
                                    "Red (Gene1); Blue (Gene2)",
                                    "Orange (Gene1); Blue (Gene2)",
                                    "Red (Gene1); Green (Gene2)",
                                    "Green (Gene1); Blue (Gene2)"
                                  ),
                                  selected = "Red (Gene1); Blue (Gene2)"
                     ),
                     radioButtons("lc_gec_ord1", "Plot order:",
                                  choices = c("Max", "Min", "Original", "Random"),
                                  selected = "Max", inline = TRUE
                     ),
                     sliderInput("lc_gec_siz", "Point size:",
                                 min = 0, max = 4, value = 1.25, step = 0.25
                     ),
                     radioButtons("lc_gec_psz", "Plot size:",
                                  choices = c("Small", "Medium", "Large"),
                                  selected = "Medium", inline = TRUE
                     ),
                     radioButtons("lc_gec_fsz", "Font size:",
                                  choices = c("Small", "Medium", "Large"),
                                  selected = "Small", inline = TRUE
                     ),
                     radioButtons("lc_gec_asp", "Aspect ratio:",
                                  choices = c("Square", "Fixed", "Free"),
                                  selected = "Square", inline = TRUE
                     ),
                     checkboxInput("lc_gec_txt", "Show axis text", value = FALSE)
                   )
                 ),
                 div(class="input-panel input-panel-section",
                     numericInput("lc_gec_oup1.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                     downloadButton("lc_gec_oup1.png", "Download PNG", class = "btn-sm"),
                     downloadButton("lc_gec_oup1.pdf", "Download PDF", class = "btn-sm"),
                     downloadButton("lc_gec_oup1.svg", "Download SVG", class = "btn-sm"),
                 ),
                 div(class="input-panel-section",
                     h4("Cell numbers"),
                     dataTableOutput("lc_gec_.dt")
                 )
               )
        ), # row 2 col 1
        # row 2 col 2
        column(
          8,
          uiOutput("lc_gec_oup1.ui"),
        )
      ), # end of row 2
      hr()
    ) # col
  ) # row
) # End of tab gec

,
# tab vio ----
tabPanel(
  "Violinplot / Boxplot",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Cell information / gene expression violin plot / box plot"),
          p("Visualise the gene expression or continuous cell information (e.g. Number of UMIs / module score) across groups of cells (e.g. libary / clusters).")
        )
      ),
      # row 2 ----
      fluidRow(
        class = "tab-section",
        # row 2 col 1
        column(
          3,
          div(
            class = "input-panel",
            style = "border-right: 2px solid #f3f6f4",
            selectInput("lc_vio_inp1", "Cell info (X-axis):",
              choices = lcconf[grp == TRUE]$UI,
              selected = lcdef$grp1
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell info to group cells",
                content = c(
                  "- Select categorical cell info to group cells by",
                  "- Single cells are grouped by this categorical covariate",
                  "- Plotted as the X-axis of the violin plot / box plot"
                )
              ),
            selectInput("lc_vio_inp2", "Cell info / Gene (Y-axis):", choices = NULL) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell Info / Gene to plot",
                content = c(
                  "- Select cell info / gene to plot on Y-axis",
                  "- Can be continuous cell info (e.g. nUMIs / scores)",
                  "- Can also be gene expression",
                  "- Type in gene names for unlisted genes"
                )
              ),
            radioButtons("lc_vio_typ", "Plot type:",
              choices = c("violin", "boxplot", "lineplot"),
              selected = "violin", inline = TRUE
            ),
            checkboxInput("lc_vio_pts", "Show data points", value = FALSE),
            checkboxInput("lc_vio_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.lc_vio_togL == true",
              selectInput("lc_vio_sub1", "Cell info to subset:",
                choices = lcconf[grp == TRUE]$UI,
                selected = lcdef$grp1
              ),
              uiOutput("lc_vio_sub1.ui"),
              actionButton("lc_vio_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("lc_vio_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            ),
            checkboxInput("lc_vio_tog", "Adjust graphics"),
            conditionalPanel(
              condition = "input.lc_vio_tog == true",
              sliderInput("lc_vio_siz", "Data point size:",
                min = 0, max = 4, value = 1.25, step = 0.25
              ),
              radioButtons("lc_vio_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              radioButtons("lc_vio_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              conditionalPanel(
              condition = "input.lc_vio_typ == 'lineplot'",
              sliderInput("lc_vio_barsz", "Line size", min = 0.05, max = 0.5, step = 0.01, value = 0.3)
              )
            ),
            selectInput("lc_vio_datatype", "Data type", choices = c("normalised", "raw"), selected = "normalised"),
          ),
          div(
            class = "input-panel",
            numericInput("lc_vio_oup.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
            downloadButton("lc_vio_oup.png", "Download PNG", class = "btn-sm"),
            downloadButton("lc_vio_oup.pdf", "Download PDF", class = "btn-sm"),
            downloadButton("lc_vio_oup.svg", "Download SVG", class = "btn-sm")
          )
        ), # row 2 col 1
        # row 2 col 2
        column(
          9, uiOutput("lc_vio_oup.ui")
        ) # row 2 col 2
      ), # row 2
      hr()
    )
  )
) # End of tab vio

,
# tab pro ----
tabPanel(
  "Proportion plot",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Proportion / cell numbers across different cell information"),
          p("Visualise the composition of single cells based on one discrete cell information across another discrete cell information.")
        )
      ),
      # row 2 ----
      fluidRow(
        class = "tab-section",
        # row 2 col 1
        column(
          3,
          div(
            class = "input-panel",
            selectInput("lc_pro_inp1", "Cell info to plot (X-axis):",
              choices = lcconf[grp == TRUE]$UI,
              selected = lcdef$grp2
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell info to group cells",
                content = c(
                  "- Select categorical cell info to group cells",
                  "- Plotted as the X-axis of the proportion plot"
                )
              ),
            selectInput("lc_pro_inp2", "Cell info to group / colour by:",
              choices = lcconf[grp == TRUE]$UI,
              selected = lcdef$grp1
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell info to group / colour cells",
                content = c(
                  "- Select categorical cell info to group / colour cells",
                  "- Proportion / cell numbers are shown in different colours"
                )
              ),
            radioButtons("lc_pro_typ", "Plot value:",
              choices = c("Proportion", "CellNumbers"),
              selected = "Proportion", inline = TRUE
            ),
            checkboxInput("lc_pro_flp", "Flip X/Y", value = FALSE),
            checkboxInput("lc_pro_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.lc_pro_togL == true",
              selectInput("lc_pro_sub1", "Cell info to subset:",
                choices = lcconf[grp == TRUE]$UI,
                selected = lcdef$grp1
              ),
              uiOutput("lc_pro_sub1.ui"),
              actionButton("lc_pro_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("lc_pro_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            ),
            checkboxInput("lc_pro_tog", "Adjust graphics"),
            conditionalPanel(
              condition = "input.lc_pro_tog == true",
              radioButtons("lc_pro_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("lc_pro_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              )
            )
          ),
          div(
            class = "input-panel",
            numericInput("lc_pro_oup.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
            downloadButton("lc_pro_oup.png", "Download PNG", class = "btn-sm"),
            downloadButton("lc_pro_oup.pdf", "Download PDF", class = "btn-sm"),
            downloadButton("lc_pro_oup.svg", "Download SVG", class = "btn-sm")
          )
        ), # row 2 col 1
        # row 2 col 2
        column(
          9, uiOutput("lc_pro_oup.ui")
        ) # row 2 col 2
      ), # row 2
      hr()
    )
  )
) # End of tab pro

,
# tab hea ----
tabPanel(
  "Bubbleplot / Heatmap",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Gene expression bubbleplot / heatmap"),
          p("Visualise the gene expression patterns of multiple genes grouped by categorical cell information (e.g. library / cluster). The normalised expression are averaged, log-transformed and then plotted.")
        )
      ),
      # row 2 ----
      fluidRow(
        class = "tab-section",
        column(
          4,
          div(
            class = "input-panel",
            selectInput("lc_hea_inp", "Genes:", choices = NULL, multiple = TRUE) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "List of genes to plot on bubbleplot / heatmap",
                content = c(
                  "- Input genes to plot",
                  "- Type in gene names for unlisted genes"
                )
              ),
            selectInput("lc_hea_grp", "Group by:",
              choices = lcconf[grp == TRUE]$UI,
              selected = lcconf[grp == TRUE]$UI[1]
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell info to group cells",
                content = c(
                  "- Select categorical cell info to group cells",
                  "- Single cells are grouped by this categorical covariate",
                  "- Plotted as the X-axis of the bubbleplot / heatmap"
                )
              ),
            radioButtons("lc_hea_plt", "Plot type:",
              choices = c("Bubbleplot", "Heatmap"),
              selected = "Bubbleplot", inline = TRUE
            ),
            checkboxInput("lc_hea_scl", "Scale gene expression", value = TRUE),
            checkboxInput("lc_hea_row", "Cluster rows (genes)", value = TRUE),
            checkboxInput("lc_hea_col", "Cluster columns (samples)", value = FALSE),
            checkboxInput("lc_hea_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.lc_hea_togL == true",
              selectInput("lc_hea_sub1", "Cell info to subset:",
                choices = lcconf[grp == TRUE]$UI,
                selected = lcdef$grp1
              ),
              uiOutput("lc_hea_sub1.ui"),
              actionButton("lc_hea_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("lc_hea_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            ),
            checkboxInput("lc_hea_tog", "Adjust graphics"),
            conditionalPanel(
              condition = "input.lc_hea_tog == true",
              radioButtons("lc_hea_cols", "Colour scheme:",
                choices = c(
                  "White-Red", "Blue-Yellow-Red",
                  "Yellow-Green-Purple"
                ),
                selected = "Blue-Yellow-Red"
              ),
              radioButtons("lc_hea_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("lc_hea_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              )
            )
          ),
          div(
            class = "input-panel",
            fluidRow(
              column(4,
                     numericInput("lc_hea_oup.height", "Height:", min = 5, max = 100, value = 18, step = 2)
              ),
              column(4,
                     numericInput("lc_hea_oup.width", "Width:", min = 5, max = 100, value = 18, step = 2)
              ),
              column(4,
                     numericInput("lc_hea_oup.res", "Res:", min = 72, max = 600, value = 150, step = 5)
              )
            ),
            downloadButton("lc_hea_oup.png", "Download PNG", class = "btn-sm"),
            downloadButton("lc_hea_oup.pdf", "Download PDF", class = "btn-sm"),
            downloadButton("lc_hea_oup.svg", "Download SVG", class = "btn-sm")
          )
        ), # row 2 col 1
        # row 2 col 2
        column(
          8, h4(htmlOutput("lc_hea_oupTxt")),
          uiOutput("lc_hea_oup.ui")
        ) # row 2 col 2
      ), # row 2
      hr()
    )
  )
) # End of tab hea

,
# tab mar ----
tabPanel(
  "Markers",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Markers"),
          p("Explore markers for different clustering.")
        ) # row 1 col 1
      ), # row 1
      # row 2 ----
      fluidRow(
        class = "tab-section",
        column(3,
               fluidRow(
                 column(
                   12,
                   div(
                     class = "input-panel input-panel-section",
                     selectInput("lc_mar_cls","Select clustering:", choices = names(lcmar),selected = 1)
                   )
                 )
               )
        )
      ), # end of row 2
      # row 3 ----
      fluidRow(
      column(12,
        DTOutput("lc_mar_table")
      )
      ), # end of row 3
      hr()
    )
  )
) # End of tab mar

)
,navbarMenu("LEC Diseased",
# tab civge ----
tabPanel(
  "CellInfo vs GeneExpr",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Cell information vs Gene expression"),
          p("Cell information and gene expression side-by-side on low-dimensional represention.")
        )
      ),
      # row 2 ----
      fluidRow(
        class = "tab-section",
        # row 2 col 1
        column(
          4,
          div(
            class = "input-panel",
            h4("Dimension Reduction"),
            selectInput("ld_civge_drX", "X-axis:",
              choices = ldconf[dimred == TRUE]$UI,
              selected = lddef$dimred[1]
            ),
            selectInput("ld_civge_drY", "Y-axis:",
              choices = ldconf[dimred == TRUE]$UI,
              selected = lddef$dimred[2]
            )
          )
        ), # row 1 col 1
        # row 2 col 2
        column(
          4,
          div(
            class = "input-panel",
            checkboxInput("ld_civge_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.ld_civge_togL == true",
              selectInput("ld_civge_sub1", "Cell info to subset:",
                choices = ldconf[grp == TRUE]$UI,
                selected = lddef$grp1
              ),
              uiOutput("ld_civge_sub1.ui"),
              actionButton("ld_civge_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("ld_civge_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            )
          )
        ), # End of column
        # row 2 col 3
        column(
          4,
          div(
            class = "input-panel",
            checkboxInput("ld_civge_tog0", "Adjust graphics"),
            conditionalPanel(
              condition = "input.ld_civge_tog0 == true",
              sliderInput("ld_civge_siz", "Point size:",
                min = 0, max = 4, value = 1.25, step = 0.25
              ),
              radioButtons("ld_civge_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("ld_civge_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              radioButtons("ld_civge_asp", "Aspect ratio:",
                choices = c("Square", "Fixed", "Free"),
                selected = "Square", inline = TRUE
              ),
              checkboxInput("ld_civge_txt", "Show axis text", value = FALSE)
            )
          )
        ) # row 2 col 3
      ), # row 2
      # row 3 ----
      fluidRow(
        class = "tab-section",
        # row 3 col 1
        column(6,
          style = "border-right: 2px solid #f3f6f4",
          h4("Cell information"),
          # row 3 col 1 row 1
          fluidRow(
            class = "tab-section",
            column(
              6,
              div(
                class = "input-panel",
                selectInput("ld_civge_inp1", "Cell info:",
                  choices = ldconf$UI,
                  selected = lddef$meta1
                ) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Cell info to colour cells",
                    content = c(
                      "- Select cell info to colour cells",
                      "- Categorical covariates have a fixed colour palette",
                      "- Continuous covariates are coloured in a Blue-Yellow-Red colour scheme, which can be changed in the plot controls"
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("ld_civge_tog1", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.ld_civge_tog1 == true",
                  radioButtons("ld_civge_col1", "Colour (Continuous data):",
                    choices = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple"),
                    selected = "Blue-Yellow-Red"
                  ),
                  radioButtons("ld_civge_ord1", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Original", inline = TRUE
                  ),
                  checkboxInput("ld_civge_lab1", "Show cell info labels", value = TRUE)
                )
              )
            )
          ),
          # row 3 col 1 row 2
          fluidRow(
            class = "tab-section",
            column(
              12,
              uiOutput("ld_civge_oup1.ui")
            )
          ),
          # row 3 col 1 row 3
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("ld_civge_oup1.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("ld_civge_oup1.png", "Download PNG", class = "btn-sm"),
                downloadButton("ld_civge_oup1.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("ld_civge_oup1.svg", "Download SVG", class = "btn-sm"),
                checkboxInput("ld_civge_tog9", "Show cell numbers / statistics")
              )
            )
          ),
          conditionalPanel(
            condition = "input.ld_civge_tog9 == true",
            h4("Cell numbers / statistics"),
            radioButtons("ld_civge_splt", "Split continuous cell info into:",
              choices = c("Quartile", "Decile"),
              selected = "Decile", inline = TRUE
            ),
            dataTableOutput("ld_civge_.dt")
          )
        ), # row 3 col 1
        # row 3 col 2
        column(
          6,
          h4("Gene expression"),
          fluidRow(
            class = "tab-section",
            column(
              6,
              div(
                class = "input-panel",
                selectInput("ld_civge_inp2", "Gene name:", choices = NULL) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Gene expression to colour cells",
                    content = c(
                      "- Select gene to colour cells by gene expression",
                      "- Type in gene names for unlisted genes",
                      "- Gene expression are coloured in a White-Red colour scheme which can be changed in the plot controls"
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("ld_civge_tog2", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.ld_civge_tog2 == true",
                  radioButtons("ld_civge_col2", "Colour:",
                    choices = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple"),
                    selected = "White-Red"
                  ),
                  radioButtons("ld_civge_ord2", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Max", inline = TRUE
                  )
                )
              )
            )
          ),
          fluidRow(
            class = "tab-section",
            column(
              12,
              uiOutput("ld_civge_oup2.ui")
            )
          ),
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("ld_civge_oup2.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("ld_civge_oup2.png", "Download PNG", class = "btn-sm"),
                downloadButton("ld_civge_oup2.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("ld_civge_oup2.svg", "Download SVG", class = "btn-sm")
              )
            )
          )
        ) # row 3 col 2
      ), # row 3
      hr()
    )
  )
) # End of tab civge

,
# tab civci ----
tabPanel(
  "CellInfo vs CellInfo",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Cell info vs cell info"),
          p("Two cell infos side-by-side on low-dimensional represention.")
        ) # row 1 col 1
      ), # row 1
      # row 2 ----
      fluidRow(
        class = "tab-section",
        column(
          4,
          div(
            class = "input-panel",
            h4("Dimension Reduction"),
            selectInput("ld_civci_drX", "X-axis:",
              choices = ldconf[dimred == TRUE]$UI,
              selected = lddef$dimred[1]
            ),
            selectInput("ld_civci_drY", "Y-axis:",
              choices = ldconf[dimred == TRUE]$UI,
              selected = lddef$dimred[2]
            )
          )
        ), # row 2 col 2
        # row 2 col 2
        column(4,
          div(
            class = "input-panel",
            checkboxInput("ld_civci_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.ld_civci_togL == true",
              selectInput("ld_civci_sub1", "Cell info to subset:",
                choices = ldconf[grp == TRUE]$UI,
                selected = lddef$grp1
              ),
              uiOutput("ld_civci_sub1.ui"),
              actionButton("ld_civci_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("ld_civci_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            )
          )
        ), # row 2 col 2
        # row 2 col 3
        column(
          4,
          div(
            class = "input-panel",
            checkboxInput("ld_civci_tog0", "Adjust graphics"),
            conditionalPanel(
              condition = "input.ld_civci_tog0 == true",
              sliderInput("ld_civci_siz", "Point size:",
                min = 0, max = 4, value = 1.25, step = 0.25
              ),
              radioButtons("ld_civci_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("ld_civci_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              radioButtons("ld_civci_asp", "Aspect ratio:",
                choices = c("Square", "Fixed", "Free"),
                selected = "Square", inline = TRUE
              ),
              checkboxInput("ld_civci_txt", "Show axis text", value = FALSE)
            )
          )
        ) # row 2 col 3
      ), # row 2
      # row 3 ----
      fluidRow(
        class = "tab-section",
        # row 3 col 1
        column(
          6,
          style = "border-right: 2px solid #f3f6f4",
          h4("Cell info 1"),
          # row 3 col 1 row 1
          fluidRow(
            class = "tab-section",
            column(
              6,
              div(
                class = "input-panel",
                selectInput("ld_civci_inp1", "Cell info:",
                  choices = ldconf$UI,
                  selected = lddef$meta1
                ) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Cell info to colour cells",
                    content = c(
                      "- Select cell info to colour cells",
                      "- Categorical covariates have a fixed colour palette",
                      "- Continuous covariates are coloured in a Blue-Yellow-Red colour scheme, which can be changed in the plot controls"
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("ld_civci_tog1", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.ld_civci_tog1 == true",
                  radioButtons("ld_civci_col1", "Colour (Continuous data):",
                    choices = c(
                      "White-Red", "Blue-Yellow-Red",
                      "Yellow-Green-Purple"
                    ),
                    selected = "Blue-Yellow-Red"
                  ),
                  radioButtons("ld_civci_ord1", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Original", inline = TRUE
                  ),
                  checkboxInput("ld_civci_lab1", "Show cell info labels", value = TRUE)
                )
              )
            )
          ),
          # row 3 col 1 row 2
          fluidRow(column(12, uiOutput("ld_civci_oup1.ui"))),
          # row 3 col 1 row 3
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("ld_civci_oup1.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("ld_civci_oup1.png", "Download PNG", class = "btn-sm"),
                downloadButton("ld_civci_oup1.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("ld_civci_oup1.svg", "Download svg", class = "btn-sm")
              )
            )
          )
        ), # row 3 col 1
        # row 3 col 2
        column(
          6, h4("Cell info 2"),
          fluidRow(
            class = "tab-section",
            column(
              6,
              div(
                class = "input-panel",
                selectInput("ld_civci_inp2", "Cell info:",
                  choices = ldconf$UI,
                  selected = lddef$meta2
                ) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Cell info to colour cells",
                    content = c(
                      "- Select cell info to colour cells",
                      "- Categorical covariates have a fixed colour palette",
                      "- Continuous covariates are coloured in a Blue-Yellow-Red colour scheme, which can be changed in the plot controls"
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("ld_civci_tog2", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.ld_civci_tog2 == true",
                  radioButtons("ld_civci_col2", "Colour (Continuous data):",
                    choices = c(
                      "White-Red", "Blue-Yellow-Red",
                      "Yellow-Green-Purple"
                    ),
                    selected = "Blue-Yellow-Red"
                  ),
                  radioButtons("ld_civci_ord2", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Original", inline = TRUE
                  ),
                  checkboxInput("ld_civci_lab2", "Show cell info labels", value = TRUE)
                )
              )
            )
          ),
          fluidRow(column(12, uiOutput("ld_civci_oup2.ui"))),
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("ld_civci_oup2.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("ld_civci_oup2.png", "Download PNG", class = "btn-sm"),
                downloadButton("ld_civci_oup2.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("ld_civci_oup2.svg", "Download SVG", class = "btn-sm")
              )
            )
          )
        ) # row 3 col 2
      ), # row 3
      hr()
    )
  )
) # End of tab civci

,
# tab gevge ----
tabPanel(
  "GeneExpr vs GeneExpr",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Gene expression vs Gene expression"),
          p("Visualise two gene expressions side-by-side on low-dimensional representions.")
        )
      ),
      # row 2 ----
      fluidRow(
        class = "tab-section",
        # row 2 col 1
        column(
          4,
          div(
            class = "input-panel",
            h4("Dimension Reduction"),
            selectInput("ld_gevge_drX", "X-axis:",
              choices = ldconf[dimred == TRUE]$UI,
              selected = lddef$dimred[1]
            ),
            selectInput("ld_gevge_drY", "Y-axis:",
              choices = ldconf[dimred == TRUE]$UI,
              selected = lddef$dimred[2]
            )
          )
        ), # row 1 col 1
        # row 2 col 2
        column(4,
          div(
            class = "input-panel",
            checkboxInput("ld_gevge_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.ld_gevge_togL == true",
              selectInput("ld_gevge_sub1", "Cell info to subset:",
                choices = ldconf[grp == TRUE]$UI,
                selected = lddef$grp1
              ),
              uiOutput("ld_gevge_sub1.ui"),
              actionButton("ld_gevge_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("ld_gevge_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            )
          )
        ), # End of column
        # row 2 col 3
        column(
          4,
          div(
            class = "input-panel",
            checkboxInput("ld_gevge_tog0", "Adjust graphics"),
            conditionalPanel(
              condition = "input.ld_gevge_tog0 == true",
              sliderInput("ld_gevge_siz", "Point size:",
                min = 0, max = 4, value = 1.25, step = 0.25
              ),
              radioButtons("ld_gevge_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("ld_gevge_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              radioButtons("ld_gevge_asp", "Aspect ratio:",
                choices = c("Square", "Fixed", "Free"),
                selected = "Square", inline = TRUE
              ),
              checkboxInput("ld_gevge_txt", "Show axis text", value = FALSE)
            )
          )
        ) # row 2 col 3
      ), # row 2
      # row 3 ----
      fluidRow(
        class = "tab-section",
        column(
          6,
          style = "border-right: 2px solid #f3f6f4",
          h4("Gene expression 1"),
          fluidRow(
            class = "tab-section",
            column(
              6,
              div(
                class = "input-panel",
                selectInput("ld_gevge_inp1", "Gene name:", choices = NULL) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Gene expression to colour cells",
                    content = c(
                      "- Select gene to colour cells by gene expression",
                      "- Type in gene names for unlisted genes",
                      "- Gene expression are coloured in a White-Red colour scheme which can be changed in the plot controls"
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("ld_gevge_tog1", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.ld_gevge_tog1 == true",
                  radioButtons("ld_gevge_col1", "Colour:",
                    choices = c(
                      "White-Red", "Blue-Yellow-Red",
                      "Yellow-Green-Purple"
                    ),
                    selected = "White-Red"
                  ),
                  radioButtons("ld_gevge_ord1", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Max", inline = TRUE
                  )
                )
              )
            )
          ),
          # row 3 col 1 row 2
          fluidRow(
            class = "tab-section",
            column(12, uiOutput("ld_gevge_oup1.ui"))
          ),
          # row 3 col 1 row 3
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("ld_gevge_oup1.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("ld_gevge_oup1.png", "Download PNG", class = "btn-sm"),
                downloadButton("ld_gevge_oup1.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("ld_gevge_oup1.svg", "Download SVG", class = "btn-sm")
              )
            )
          )
        ), # row 3 col 1
        # row 3 col 2
        column(
          6, h4("Gene expression 2"),
          fluidRow(
            class = "tab-section",
            column(
              6,
              div(
                class = "input-panel",
                selectInput("ld_gevge_inp2", "Gene name:", choices = NULL) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Gene expression to colour cells",
                    content = c(
                      "- Select gene to colour cells by gene expression",
                      "- Type in gene names for unlisted genes",
                      "- Gene expression are coloured in a White-Red colour scheme which can be changed in the plot controls"
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("ld_gevge_tog2", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.ld_gevge_tog2 == true",
                  radioButtons("ld_gevge_col2", "Colour:",
                    choices = c(
                      "White-Red", "Blue-Yellow-Red",
                      "Yellow-Green-Purple"
                    ),
                    selected = "White-Red"
                  ),
                  radioButtons("ld_gevge_ord2", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Max", inline = TRUE
                  )
                )
              )
            )
          ),
          fluidRow(
            class = "tab-section",
            column(12, uiOutput("ld_gevge_oup2.ui"))
          ),
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                  numericInput("ld_gevge_oup2.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                  downloadButton("ld_gevge_oup2.png", "Download PNG", class = "btn-sm"),
                  downloadButton("ld_gevge_oup2.pdf", "Download PDF", class = "btn-sm"),
                  downloadButton("ld_gevge_oup2.svg", "Download SVG", class = "btn-sm")
              )
            )
          )
        ) # row 3 col 2
      ), # row 3
      hr()
    )
  )
) # End of tab gevge

,
# tab gem ----
tabPanel(
  "Expression",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Gene expression"),
          p("Explore gene expression on low-dimensional represention.")
        ) # row 1 col 1
      ), # row 1
      # row 2 ----
      fluidRow(
        class = "tab-section",
        column(4,
               fluidRow(
        column(
          12,
          div(
            class = "input-panel input-panel-section",
            selectInput("ld_gem_inp", "Genes:", choices = NULL, multiple = TRUE) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "List of genes to plot on bubbleplot / heatmap",
                content = c(
                  "- Input genes to plot",
                  "- Type in gene names for unlisted genes"
                )
              ),
            selectInput("ld_gem_drX", "X-axis:",
                        choices = ldconf[dimred == TRUE]$UI,
                        selected = lddef$dimred[1]
            ),
            selectInput("ld_gem_drY", "Y-axis:",
                        choices = ldconf[dimred == TRUE]$UI,
                        selected = lddef$dimred[2]
            ),
            checkboxInput("ld_gem_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.ld_gem_togL == true",
              selectInput("ld_gem_sub1", "Cell info to subset:",
                          choices = ldconf[grp == TRUE]$UI,
                          selected = lddef$grp1
              ),
              uiOutput("ld_gem_sub1.ui"),
              actionButton("ld_gem_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("ld_gem_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            ),
            checkboxInput("ld_gem_tog0", "Adjust graphics"),
            conditionalPanel(
              condition = "input.ld_gem_tog0 == true",
              sliderInput("ld_gem_siz", "Point size:",
                          min = 0, max = 3, value = 0.5, step = 0.1
              ),
              radioButtons("ld_gem_psz", "Plot size:",
                           choices = c("Small", "Medium", "Large"),
                           selected = "Medium", inline = TRUE
              ),
              radioButtons("ld_gem_fsz", "Font size:",
                           choices = c("Smaller", "Small", "Medium", "Large"),
                           selected = "Small", inline = TRUE
              ),
              radioButtons("ld_gem_asp", "Aspect ratio:",
                           choices = c("Square", "Fixed", "Free"),
                           selected = "Square", inline = TRUE
              ),
              checkboxInput("ld_gem_txt", "Show axis text", value = FALSE),
              radioButtons("ld_gem_col", "Colour (Continuous data):",
                           choices = c(
                             "White-Red", "Blue-Yellow-Red",
                             "Yellow-Green-Purple"
                           ),
                           selected = "Blue-Yellow-Red"
              ),
              radioButtons("ld_gem_ord", "Plot order:",
                           choices = c("Max", "Min", "Original", "Random"),
                           selected = "Max", inline = TRUE
              ),
              numericInput("ld_gem_ncol", "Number of columns", value = 0, min = 0, step = 1)
            )
          ),
          div(
            class = "input-panel",
            fluidRow(
            column(4,
              numericInput("ld_gem_oup1.height", "Height:", min = 1, max = 50, value = 25, step = 2)
            ),
            column(4,
              numericInput("ld_gem_oup1.width", "Width:", min = 1, max = 50, value = 25, step = 2)
            ),
            column(4,
              numericInput("ld_gem_oup1.res", "Res:", min = 72, max = 600, value = 150, step = 5)
            )
            ),
            downloadButton("ld_gem_oup1.png", "Download PNG", class = "btn-sm"),
            downloadButton("ld_gem_oup1.pdf", "Download PDF", class = "btn-sm"),
            downloadButton("ld_gem_oup1.svg", "Download SVG", class = "btn-sm")
          )
        )
               )
      ),
      column(8,
             uiOutput("ld_gem_oup1.ui")
      )
      ),
      hr()
    )
  )
) # End of tab gem

,
# tab gec ----
tabPanel(
  "Gene coexpression",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Coexpression of two genes on reduced dimensions"),
          p("Visualise the coexpression of two genes on low-dimensional representions.")
        )
      ),
      # row 2 ----
      fluidRow(
        class = "tab-section",
        # row 2 col 1
        column(4,
               column(
                 12,
                 div(
                   class = "input-panel input-panel-section",
                   h4("Dimension Reduction"),
                   selectInput("ld_gec_drX", "X-axis:",
                               choices = ldconf[dimred == TRUE]$UI,
                               selected = lddef$dimred[1]
                   ),
                   selectInput("ld_gec_drY", "Y-axis:",
                               choices = ldconf[dimred == TRUE]$UI,
                               selected = lddef$dimred[2]
                   ),
                   selectInput("ld_gec_inp1", "Gene 1:", choices = NULL) %>%
                     helper(
                       type = "inline", size = "m", fade = TRUE,
                       title = "Gene expression to colour cells",
                       content = c(
                         "- Select gene to colour cells by gene expression",
                         "- Type in gene names for unlisted genes",
                         "- Gene expression are coloured in a White-Red colour scheme which can be changed in the plot controls"
                       )
                     ),
                   selectInput("ld_gec_inp2", "Gene 2:", choices = NULL) %>%
                     helper(
                       type = "inline", size = "m", fade = TRUE,
                       title = "Gene expression to colour cells",
                       content = c(
                         "- Select gene to colour cells by gene expression",
                         "- Type in gene names for unlisted genes",
                         "- Gene expression are coloured in a White-Blue colour scheme which can be changed in the plot controls"
                       )
                     ),
                   checkboxInput("ld_gec_togL", "Subset cells"),
                   conditionalPanel(
                    condition = "input.ld_gec_togL == true",
                    selectInput("ld_gec_sub1", "Cell info to subset:", choices = ldconf[grp == TRUE]$UI, selected = lddef$grp1),
                     uiOutput("ld_gec_sub1.ui"),
                     actionButton("ld_gec_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
                     actionButton("ld_gec_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
                   ),
                   checkboxInput("ld_gec_tog0", "Adjust graphics"),
                   conditionalPanel(
                     condition = "input.ld_gec_tog0 == true",
                     radioButtons("ld_gec_col1", "Colour:",
                                  choices = c(
                                    "Red (Gene1); Blue (Gene2)",
                                    "Orange (Gene1); Blue (Gene2)",
                                    "Red (Gene1); Green (Gene2)",
                                    "Green (Gene1); Blue (Gene2)"
                                  ),
                                  selected = "Red (Gene1); Blue (Gene2)"
                     ),
                     radioButtons("ld_gec_ord1", "Plot order:",
                                  choices = c("Max", "Min", "Original", "Random"),
                                  selected = "Max", inline = TRUE
                     ),
                     sliderInput("ld_gec_siz", "Point size:",
                                 min = 0, max = 4, value = 1.25, step = 0.25
                     ),
                     radioButtons("ld_gec_psz", "Plot size:",
                                  choices = c("Small", "Medium", "Large"),
                                  selected = "Medium", inline = TRUE
                     ),
                     radioButtons("ld_gec_fsz", "Font size:",
                                  choices = c("Small", "Medium", "Large"),
                                  selected = "Small", inline = TRUE
                     ),
                     radioButtons("ld_gec_asp", "Aspect ratio:",
                                  choices = c("Square", "Fixed", "Free"),
                                  selected = "Square", inline = TRUE
                     ),
                     checkboxInput("ld_gec_txt", "Show axis text", value = FALSE)
                   )
                 ),
                 div(class="input-panel input-panel-section",
                     numericInput("ld_gec_oup1.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                     downloadButton("ld_gec_oup1.png", "Download PNG", class = "btn-sm"),
                     downloadButton("ld_gec_oup1.pdf", "Download PDF", class = "btn-sm"),
                     downloadButton("ld_gec_oup1.svg", "Download SVG", class = "btn-sm"),
                 ),
                 div(class="input-panel-section",
                     h4("Cell numbers"),
                     dataTableOutput("ld_gec_.dt")
                 )
               )
        ), # row 2 col 1
        # row 2 col 2
        column(
          8,
          uiOutput("ld_gec_oup1.ui"),
        )
      ), # end of row 2
      hr()
    ) # col
  ) # row
) # End of tab gec

,
# tab vio ----
tabPanel(
  "Violinplot / Boxplot",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Cell information / gene expression violin plot / box plot"),
          p("Visualise the gene expression or continuous cell information (e.g. Number of UMIs / module score) across groups of cells (e.g. libary / clusters).")
        )
      ),
      # row 2 ----
      fluidRow(
        class = "tab-section",
        # row 2 col 1
        column(
          3,
          div(
            class = "input-panel",
            style = "border-right: 2px solid #f3f6f4",
            selectInput("ld_vio_inp1", "Cell info (X-axis):",
              choices = ldconf[grp == TRUE]$UI,
              selected = lddef$grp1
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell info to group cells",
                content = c(
                  "- Select categorical cell info to group cells by",
                  "- Single cells are grouped by this categorical covariate",
                  "- Plotted as the X-axis of the violin plot / box plot"
                )
              ),
            selectInput("ld_vio_inp2", "Cell info / Gene (Y-axis):", choices = NULL) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell Info / Gene to plot",
                content = c(
                  "- Select cell info / gene to plot on Y-axis",
                  "- Can be continuous cell info (e.g. nUMIs / scores)",
                  "- Can also be gene expression",
                  "- Type in gene names for unlisted genes"
                )
              ),
            radioButtons("ld_vio_typ", "Plot type:",
              choices = c("violin", "boxplot", "lineplot"),
              selected = "violin", inline = TRUE
            ),
            checkboxInput("ld_vio_pts", "Show data points", value = FALSE),
            checkboxInput("ld_vio_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.ld_vio_togL == true",
              selectInput("ld_vio_sub1", "Cell info to subset:",
                choices = ldconf[grp == TRUE]$UI,
                selected = lddef$grp1
              ),
              uiOutput("ld_vio_sub1.ui"),
              actionButton("ld_vio_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("ld_vio_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            ),
            checkboxInput("ld_vio_tog", "Adjust graphics"),
            conditionalPanel(
              condition = "input.ld_vio_tog == true",
              sliderInput("ld_vio_siz", "Data point size:",
                min = 0, max = 4, value = 1.25, step = 0.25
              ),
              radioButtons("ld_vio_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              radioButtons("ld_vio_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              conditionalPanel(
              condition = "input.ld_vio_typ == 'lineplot'",
              sliderInput("ld_vio_barsz", "Line size", min = 0.05, max = 0.5, step = 0.01, value = 0.3)
              )
            ),
            selectInput("ld_vio_datatype", "Data type", choices = c("normalised", "raw"), selected = "normalised"),
          ),
          div(
            class = "input-panel",
            numericInput("ld_vio_oup.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
            downloadButton("ld_vio_oup.png", "Download PNG", class = "btn-sm"),
            downloadButton("ld_vio_oup.pdf", "Download PDF", class = "btn-sm"),
            downloadButton("ld_vio_oup.svg", "Download SVG", class = "btn-sm")
          )
        ), # row 2 col 1
        # row 2 col 2
        column(
          9, uiOutput("ld_vio_oup.ui")
        ) # row 2 col 2
      ), # row 2
      hr()
    )
  )
) # End of tab vio

,
# tab pro ----
tabPanel(
  "Proportion plot",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Proportion / cell numbers across different cell information"),
          p("Visualise the composition of single cells based on one discrete cell information across another discrete cell information.")
        )
      ),
      # row 2 ----
      fluidRow(
        class = "tab-section",
        # row 2 col 1
        column(
          3,
          div(
            class = "input-panel",
            selectInput("ld_pro_inp1", "Cell info to plot (X-axis):",
              choices = ldconf[grp == TRUE]$UI,
              selected = lddef$grp2
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell info to group cells",
                content = c(
                  "- Select categorical cell info to group cells",
                  "- Plotted as the X-axis of the proportion plot"
                )
              ),
            selectInput("ld_pro_inp2", "Cell info to group / colour by:",
              choices = ldconf[grp == TRUE]$UI,
              selected = lddef$grp1
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell info to group / colour cells",
                content = c(
                  "- Select categorical cell info to group / colour cells",
                  "- Proportion / cell numbers are shown in different colours"
                )
              ),
            radioButtons("ld_pro_typ", "Plot value:",
              choices = c("Proportion", "CellNumbers"),
              selected = "Proportion", inline = TRUE
            ),
            checkboxInput("ld_pro_flp", "Flip X/Y", value = FALSE),
            checkboxInput("ld_pro_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.ld_pro_togL == true",
              selectInput("ld_pro_sub1", "Cell info to subset:",
                choices = ldconf[grp == TRUE]$UI,
                selected = lddef$grp1
              ),
              uiOutput("ld_pro_sub1.ui"),
              actionButton("ld_pro_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("ld_pro_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            ),
            checkboxInput("ld_pro_tog", "Adjust graphics"),
            conditionalPanel(
              condition = "input.ld_pro_tog == true",
              radioButtons("ld_pro_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("ld_pro_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              )
            )
          ),
          div(
            class = "input-panel",
            numericInput("ld_pro_oup.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
            downloadButton("ld_pro_oup.png", "Download PNG", class = "btn-sm"),
            downloadButton("ld_pro_oup.pdf", "Download PDF", class = "btn-sm"),
            downloadButton("ld_pro_oup.svg", "Download SVG", class = "btn-sm")
          )
        ), # row 2 col 1
        # row 2 col 2
        column(
          9, uiOutput("ld_pro_oup.ui")
        ) # row 2 col 2
      ), # row 2
      hr()
    )
  )
) # End of tab pro

,
# tab hea ----
tabPanel(
  "Bubbleplot / Heatmap",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Gene expression bubbleplot / heatmap"),
          p("Visualise the gene expression patterns of multiple genes grouped by categorical cell information (e.g. library / cluster). The normalised expression are averaged, log-transformed and then plotted.")
        )
      ),
      # row 2 ----
      fluidRow(
        class = "tab-section",
        column(
          4,
          div(
            class = "input-panel",
            selectInput("ld_hea_inp", "Genes:", choices = NULL, multiple = TRUE) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "List of genes to plot on bubbleplot / heatmap",
                content = c(
                  "- Input genes to plot",
                  "- Type in gene names for unlisted genes"
                )
              ),
            selectInput("ld_hea_grp", "Group by:",
              choices = ldconf[grp == TRUE]$UI,
              selected = ldconf[grp == TRUE]$UI[1]
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell info to group cells",
                content = c(
                  "- Select categorical cell info to group cells",
                  "- Single cells are grouped by this categorical covariate",
                  "- Plotted as the X-axis of the bubbleplot / heatmap"
                )
              ),
            radioButtons("ld_hea_plt", "Plot type:",
              choices = c("Bubbleplot", "Heatmap"),
              selected = "Bubbleplot", inline = TRUE
            ),
            checkboxInput("ld_hea_scl", "Scale gene expression", value = TRUE),
            checkboxInput("ld_hea_row", "Cluster rows (genes)", value = TRUE),
            checkboxInput("ld_hea_col", "Cluster columns (samples)", value = FALSE),
            checkboxInput("ld_hea_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.ld_hea_togL == true",
              selectInput("ld_hea_sub1", "Cell info to subset:",
                choices = ldconf[grp == TRUE]$UI,
                selected = lddef$grp1
              ),
              uiOutput("ld_hea_sub1.ui"),
              actionButton("ld_hea_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("ld_hea_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            ),
            checkboxInput("ld_hea_tog", "Adjust graphics"),
            conditionalPanel(
              condition = "input.ld_hea_tog == true",
              radioButtons("ld_hea_cols", "Colour scheme:",
                choices = c(
                  "White-Red", "Blue-Yellow-Red",
                  "Yellow-Green-Purple"
                ),
                selected = "Blue-Yellow-Red"
              ),
              radioButtons("ld_hea_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("ld_hea_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              )
            )
          ),
          div(
            class = "input-panel",
            fluidRow(
              column(4,
                     numericInput("ld_hea_oup.height", "Height:", min = 5, max = 100, value = 18, step = 2)
              ),
              column(4,
                     numericInput("ld_hea_oup.width", "Width:", min = 5, max = 100, value = 18, step = 2)
              ),
              column(4,
                     numericInput("ld_hea_oup.res", "Res:", min = 72, max = 600, value = 150, step = 5)
              )
            ),
            downloadButton("ld_hea_oup.png", "Download PNG", class = "btn-sm"),
            downloadButton("ld_hea_oup.pdf", "Download PDF", class = "btn-sm"),
            downloadButton("ld_hea_oup.svg", "Download SVG", class = "btn-sm")
          )
        ), # row 2 col 1
        # row 2 col 2
        column(
          8, h4(htmlOutput("ld_hea_oupTxt")),
          uiOutput("ld_hea_oup.ui")
        ) # row 2 col 2
      ), # row 2
      hr()
    )
  )
) # End of tab hea

,
# tab mar ----
tabPanel(
  "Markers",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Markers"),
          p("Explore markers for different clustering.")
        ) # row 1 col 1
      ), # row 1
      # row 2 ----
      fluidRow(
        class = "tab-section",
        column(3,
               fluidRow(
                 column(
                   12,
                   div(
                     class = "input-panel input-panel-section",
                     selectInput("ld_mar_cls","Select clustering:", choices = names(ldmar),selected = 1)
                   )
                 )
               )
        )
      ), # end of row 2
      # row 3 ----
      fluidRow(
      column(12,
        DTOutput("ld_mar_table")
      )
      ), # end of row 3
      hr()
    )
  )
) # End of tab mar

)
,navbarMenu("LEC Healthy",
# tab civge ----
tabPanel(
  "CellInfo vs GeneExpr",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Cell information vs Gene expression"),
          p("Cell information and gene expression side-by-side on low-dimensional represention.")
        )
      ),
      # row 2 ----
      fluidRow(
        class = "tab-section",
        # row 2 col 1
        column(
          4,
          div(
            class = "input-panel",
            h4("Dimension Reduction"),
            selectInput("lh_civge_drX", "X-axis:",
              choices = lhconf[dimred == TRUE]$UI,
              selected = lhdef$dimred[1]
            ),
            selectInput("lh_civge_drY", "Y-axis:",
              choices = lhconf[dimred == TRUE]$UI,
              selected = lhdef$dimred[2]
            )
          )
        ), # row 1 col 1
        # row 2 col 2
        column(
          4,
          div(
            class = "input-panel",
            checkboxInput("lh_civge_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.lh_civge_togL == true",
              selectInput("lh_civge_sub1", "Cell info to subset:",
                choices = lhconf[grp == TRUE]$UI,
                selected = lhdef$grp1
              ),
              uiOutput("lh_civge_sub1.ui"),
              actionButton("lh_civge_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("lh_civge_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            )
          )
        ), # End of column
        # row 2 col 3
        column(
          4,
          div(
            class = "input-panel",
            checkboxInput("lh_civge_tog0", "Adjust graphics"),
            conditionalPanel(
              condition = "input.lh_civge_tog0 == true",
              sliderInput("lh_civge_siz", "Point size:",
                min = 0, max = 4, value = 1.25, step = 0.25
              ),
              radioButtons("lh_civge_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("lh_civge_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              radioButtons("lh_civge_asp", "Aspect ratio:",
                choices = c("Square", "Fixed", "Free"),
                selected = "Square", inline = TRUE
              ),
              checkboxInput("lh_civge_txt", "Show axis text", value = FALSE)
            )
          )
        ) # row 2 col 3
      ), # row 2
      # row 3 ----
      fluidRow(
        class = "tab-section",
        # row 3 col 1
        column(6,
          style = "border-right: 2px solid #f3f6f4",
          h4("Cell information"),
          # row 3 col 1 row 1
          fluidRow(
            class = "tab-section",
            column(
              6,
              div(
                class = "input-panel",
                selectInput("lh_civge_inp1", "Cell info:",
                  choices = lhconf$UI,
                  selected = lhdef$meta1
                ) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Cell info to colour cells",
                    content = c(
                      "- Select cell info to colour cells",
                      "- Categorical covariates have a fixed colour palette",
                      "- Continuous covariates are coloured in a Blue-Yellow-Red colour scheme, which can be changed in the plot controls"
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("lh_civge_tog1", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.lh_civge_tog1 == true",
                  radioButtons("lh_civge_col1", "Colour (Continuous data):",
                    choices = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple"),
                    selected = "Blue-Yellow-Red"
                  ),
                  radioButtons("lh_civge_ord1", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Original", inline = TRUE
                  ),
                  checkboxInput("lh_civge_lab1", "Show cell info labels", value = TRUE)
                )
              )
            )
          ),
          # row 3 col 1 row 2
          fluidRow(
            class = "tab-section",
            column(
              12,
              uiOutput("lh_civge_oup1.ui")
            )
          ),
          # row 3 col 1 row 3
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("lh_civge_oup1.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("lh_civge_oup1.png", "Download PNG", class = "btn-sm"),
                downloadButton("lh_civge_oup1.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("lh_civge_oup1.svg", "Download SVG", class = "btn-sm"),
                checkboxInput("lh_civge_tog9", "Show cell numbers / statistics")
              )
            )
          ),
          conditionalPanel(
            condition = "input.lh_civge_tog9 == true",
            h4("Cell numbers / statistics"),
            radioButtons("lh_civge_splt", "Split continuous cell info into:",
              choices = c("Quartile", "Decile"),
              selected = "Decile", inline = TRUE
            ),
            dataTableOutput("lh_civge_.dt")
          )
        ), # row 3 col 1
        # row 3 col 2
        column(
          6,
          h4("Gene expression"),
          fluidRow(
            class = "tab-section",
            column(
              6,
              div(
                class = "input-panel",
                selectInput("lh_civge_inp2", "Gene name:", choices = NULL) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Gene expression to colour cells",
                    content = c(
                      "- Select gene to colour cells by gene expression",
                      "- Type in gene names for unlisted genes",
                      "- Gene expression are coloured in a White-Red colour scheme which can be changed in the plot controls"
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("lh_civge_tog2", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.lh_civge_tog2 == true",
                  radioButtons("lh_civge_col2", "Colour:",
                    choices = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple"),
                    selected = "White-Red"
                  ),
                  radioButtons("lh_civge_ord2", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Max", inline = TRUE
                  )
                )
              )
            )
          ),
          fluidRow(
            class = "tab-section",
            column(
              12,
              uiOutput("lh_civge_oup2.ui")
            )
          ),
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("lh_civge_oup2.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("lh_civge_oup2.png", "Download PNG", class = "btn-sm"),
                downloadButton("lh_civge_oup2.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("lh_civge_oup2.svg", "Download SVG", class = "btn-sm")
              )
            )
          )
        ) # row 3 col 2
      ), # row 3
      hr()
    )
  )
) # End of tab civge

,
# tab civci ----
tabPanel(
  "CellInfo vs CellInfo",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Cell info vs cell info"),
          p("Two cell infos side-by-side on low-dimensional represention.")
        ) # row 1 col 1
      ), # row 1
      # row 2 ----
      fluidRow(
        class = "tab-section",
        column(
          4,
          div(
            class = "input-panel",
            h4("Dimension Reduction"),
            selectInput("lh_civci_drX", "X-axis:",
              choices = lhconf[dimred == TRUE]$UI,
              selected = lhdef$dimred[1]
            ),
            selectInput("lh_civci_drY", "Y-axis:",
              choices = lhconf[dimred == TRUE]$UI,
              selected = lhdef$dimred[2]
            )
          )
        ), # row 2 col 2
        # row 2 col 2
        column(4,
          div(
            class = "input-panel",
            checkboxInput("lh_civci_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.lh_civci_togL == true",
              selectInput("lh_civci_sub1", "Cell info to subset:",
                choices = lhconf[grp == TRUE]$UI,
                selected = lhdef$grp1
              ),
              uiOutput("lh_civci_sub1.ui"),
              actionButton("lh_civci_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("lh_civci_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            )
          )
        ), # row 2 col 2
        # row 2 col 3
        column(
          4,
          div(
            class = "input-panel",
            checkboxInput("lh_civci_tog0", "Adjust graphics"),
            conditionalPanel(
              condition = "input.lh_civci_tog0 == true",
              sliderInput("lh_civci_siz", "Point size:",
                min = 0, max = 4, value = 1.25, step = 0.25
              ),
              radioButtons("lh_civci_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("lh_civci_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              radioButtons("lh_civci_asp", "Aspect ratio:",
                choices = c("Square", "Fixed", "Free"),
                selected = "Square", inline = TRUE
              ),
              checkboxInput("lh_civci_txt", "Show axis text", value = FALSE)
            )
          )
        ) # row 2 col 3
      ), # row 2
      # row 3 ----
      fluidRow(
        class = "tab-section",
        # row 3 col 1
        column(
          6,
          style = "border-right: 2px solid #f3f6f4",
          h4("Cell info 1"),
          # row 3 col 1 row 1
          fluidRow(
            class = "tab-section",
            column(
              6,
              div(
                class = "input-panel",
                selectInput("lh_civci_inp1", "Cell info:",
                  choices = lhconf$UI,
                  selected = lhdef$meta1
                ) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Cell info to colour cells",
                    content = c(
                      "- Select cell info to colour cells",
                      "- Categorical covariates have a fixed colour palette",
                      "- Continuous covariates are coloured in a Blue-Yellow-Red colour scheme, which can be changed in the plot controls"
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("lh_civci_tog1", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.lh_civci_tog1 == true",
                  radioButtons("lh_civci_col1", "Colour (Continuous data):",
                    choices = c(
                      "White-Red", "Blue-Yellow-Red",
                      "Yellow-Green-Purple"
                    ),
                    selected = "Blue-Yellow-Red"
                  ),
                  radioButtons("lh_civci_ord1", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Original", inline = TRUE
                  ),
                  checkboxInput("lh_civci_lab1", "Show cell info labels", value = TRUE)
                )
              )
            )
          ),
          # row 3 col 1 row 2
          fluidRow(column(12, uiOutput("lh_civci_oup1.ui"))),
          # row 3 col 1 row 3
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("lh_civci_oup1.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("lh_civci_oup1.png", "Download PNG", class = "btn-sm"),
                downloadButton("lh_civci_oup1.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("lh_civci_oup1.svg", "Download svg", class = "btn-sm")
              )
            )
          )
        ), # row 3 col 1
        # row 3 col 2
        column(
          6, h4("Cell info 2"),
          fluidRow(
            class = "tab-section",
            column(
              6,
              div(
                class = "input-panel",
                selectInput("lh_civci_inp2", "Cell info:",
                  choices = lhconf$UI,
                  selected = lhdef$meta2
                ) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Cell info to colour cells",
                    content = c(
                      "- Select cell info to colour cells",
                      "- Categorical covariates have a fixed colour palette",
                      "- Continuous covariates are coloured in a Blue-Yellow-Red colour scheme, which can be changed in the plot controls"
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("lh_civci_tog2", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.lh_civci_tog2 == true",
                  radioButtons("lh_civci_col2", "Colour (Continuous data):",
                    choices = c(
                      "White-Red", "Blue-Yellow-Red",
                      "Yellow-Green-Purple"
                    ),
                    selected = "Blue-Yellow-Red"
                  ),
                  radioButtons("lh_civci_ord2", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Original", inline = TRUE
                  ),
                  checkboxInput("lh_civci_lab2", "Show cell info labels", value = TRUE)
                )
              )
            )
          ),
          fluidRow(column(12, uiOutput("lh_civci_oup2.ui"))),
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("lh_civci_oup2.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("lh_civci_oup2.png", "Download PNG", class = "btn-sm"),
                downloadButton("lh_civci_oup2.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("lh_civci_oup2.svg", "Download SVG", class = "btn-sm")
              )
            )
          )
        ) # row 3 col 2
      ), # row 3
      hr()
    )
  )
) # End of tab civci

,
# tab gevge ----
tabPanel(
  "GeneExpr vs GeneExpr",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Gene expression vs Gene expression"),
          p("Visualise two gene expressions side-by-side on low-dimensional representions.")
        )
      ),
      # row 2 ----
      fluidRow(
        class = "tab-section",
        # row 2 col 1
        column(
          4,
          div(
            class = "input-panel",
            h4("Dimension Reduction"),
            selectInput("lh_gevge_drX", "X-axis:",
              choices = lhconf[dimred == TRUE]$UI,
              selected = lhdef$dimred[1]
            ),
            selectInput("lh_gevge_drY", "Y-axis:",
              choices = lhconf[dimred == TRUE]$UI,
              selected = lhdef$dimred[2]
            )
          )
        ), # row 1 col 1
        # row 2 col 2
        column(4,
          div(
            class = "input-panel",
            checkboxInput("lh_gevge_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.lh_gevge_togL == true",
              selectInput("lh_gevge_sub1", "Cell info to subset:",
                choices = lhconf[grp == TRUE]$UI,
                selected = lhdef$grp1
              ),
              uiOutput("lh_gevge_sub1.ui"),
              actionButton("lh_gevge_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("lh_gevge_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            )
          )
        ), # End of column
        # row 2 col 3
        column(
          4,
          div(
            class = "input-panel",
            checkboxInput("lh_gevge_tog0", "Adjust graphics"),
            conditionalPanel(
              condition = "input.lh_gevge_tog0 == true",
              sliderInput("lh_gevge_siz", "Point size:",
                min = 0, max = 4, value = 1.25, step = 0.25
              ),
              radioButtons("lh_gevge_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("lh_gevge_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              radioButtons("lh_gevge_asp", "Aspect ratio:",
                choices = c("Square", "Fixed", "Free"),
                selected = "Square", inline = TRUE
              ),
              checkboxInput("lh_gevge_txt", "Show axis text", value = FALSE)
            )
          )
        ) # row 2 col 3
      ), # row 2
      # row 3 ----
      fluidRow(
        class = "tab-section",
        column(
          6,
          style = "border-right: 2px solid #f3f6f4",
          h4("Gene expression 1"),
          fluidRow(
            class = "tab-section",
            column(
              6,
              div(
                class = "input-panel",
                selectInput("lh_gevge_inp1", "Gene name:", choices = NULL) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Gene expression to colour cells",
                    content = c(
                      "- Select gene to colour cells by gene expression",
                      "- Type in gene names for unlisted genes",
                      "- Gene expression are coloured in a White-Red colour scheme which can be changed in the plot controls"
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("lh_gevge_tog1", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.lh_gevge_tog1 == true",
                  radioButtons("lh_gevge_col1", "Colour:",
                    choices = c(
                      "White-Red", "Blue-Yellow-Red",
                      "Yellow-Green-Purple"
                    ),
                    selected = "White-Red"
                  ),
                  radioButtons("lh_gevge_ord1", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Max", inline = TRUE
                  )
                )
              )
            )
          ),
          # row 3 col 1 row 2
          fluidRow(
            class = "tab-section",
            column(12, uiOutput("lh_gevge_oup1.ui"))
          ),
          # row 3 col 1 row 3
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("lh_gevge_oup1.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("lh_gevge_oup1.png", "Download PNG", class = "btn-sm"),
                downloadButton("lh_gevge_oup1.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("lh_gevge_oup1.svg", "Download SVG", class = "btn-sm")
              )
            )
          )
        ), # row 3 col 1
        # row 3 col 2
        column(
          6, h4("Gene expression 2"),
          fluidRow(
            class = "tab-section",
            column(
              6,
              div(
                class = "input-panel",
                selectInput("lh_gevge_inp2", "Gene name:", choices = NULL) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Gene expression to colour cells",
                    content = c(
                      "- Select gene to colour cells by gene expression",
                      "- Type in gene names for unlisted genes",
                      "- Gene expression are coloured in a White-Red colour scheme which can be changed in the plot controls"
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("lh_gevge_tog2", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.lh_gevge_tog2 == true",
                  radioButtons("lh_gevge_col2", "Colour:",
                    choices = c(
                      "White-Red", "Blue-Yellow-Red",
                      "Yellow-Green-Purple"
                    ),
                    selected = "White-Red"
                  ),
                  radioButtons("lh_gevge_ord2", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Max", inline = TRUE
                  )
                )
              )
            )
          ),
          fluidRow(
            class = "tab-section",
            column(12, uiOutput("lh_gevge_oup2.ui"))
          ),
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                  numericInput("lh_gevge_oup2.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                  downloadButton("lh_gevge_oup2.png", "Download PNG", class = "btn-sm"),
                  downloadButton("lh_gevge_oup2.pdf", "Download PDF", class = "btn-sm"),
                  downloadButton("lh_gevge_oup2.svg", "Download SVG", class = "btn-sm")
              )
            )
          )
        ) # row 3 col 2
      ), # row 3
      hr()
    )
  )
) # End of tab gevge

,
# tab gem ----
tabPanel(
  "Expression",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Gene expression"),
          p("Explore gene expression on low-dimensional represention.")
        ) # row 1 col 1
      ), # row 1
      # row 2 ----
      fluidRow(
        class = "tab-section",
        column(4,
               fluidRow(
        column(
          12,
          div(
            class = "input-panel input-panel-section",
            selectInput("lh_gem_inp", "Genes:", choices = NULL, multiple = TRUE) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "List of genes to plot on bubbleplot / heatmap",
                content = c(
                  "- Input genes to plot",
                  "- Type in gene names for unlisted genes"
                )
              ),
            selectInput("lh_gem_drX", "X-axis:",
                        choices = lhconf[dimred == TRUE]$UI,
                        selected = lhdef$dimred[1]
            ),
            selectInput("lh_gem_drY", "Y-axis:",
                        choices = lhconf[dimred == TRUE]$UI,
                        selected = lhdef$dimred[2]
            ),
            checkboxInput("lh_gem_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.lh_gem_togL == true",
              selectInput("lh_gem_sub1", "Cell info to subset:",
                          choices = lhconf[grp == TRUE]$UI,
                          selected = lhdef$grp1
              ),
              uiOutput("lh_gem_sub1.ui"),
              actionButton("lh_gem_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("lh_gem_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            ),
            checkboxInput("lh_gem_tog0", "Adjust graphics"),
            conditionalPanel(
              condition = "input.lh_gem_tog0 == true",
              sliderInput("lh_gem_siz", "Point size:",
                          min = 0, max = 3, value = 0.5, step = 0.1
              ),
              radioButtons("lh_gem_psz", "Plot size:",
                           choices = c("Small", "Medium", "Large"),
                           selected = "Medium", inline = TRUE
              ),
              radioButtons("lh_gem_fsz", "Font size:",
                           choices = c("Smaller", "Small", "Medium", "Large"),
                           selected = "Small", inline = TRUE
              ),
              radioButtons("lh_gem_asp", "Aspect ratio:",
                           choices = c("Square", "Fixed", "Free"),
                           selected = "Square", inline = TRUE
              ),
              checkboxInput("lh_gem_txt", "Show axis text", value = FALSE),
              radioButtons("lh_gem_col", "Colour (Continuous data):",
                           choices = c(
                             "White-Red", "Blue-Yellow-Red",
                             "Yellow-Green-Purple"
                           ),
                           selected = "Blue-Yellow-Red"
              ),
              radioButtons("lh_gem_ord", "Plot order:",
                           choices = c("Max", "Min", "Original", "Random"),
                           selected = "Max", inline = TRUE
              ),
              numericInput("lh_gem_ncol", "Number of columns", value = 0, min = 0, step = 1)
            )
          ),
          div(
            class = "input-panel",
            fluidRow(
            column(4,
              numericInput("lh_gem_oup1.height", "Height:", min = 1, max = 50, value = 25, step = 2)
            ),
            column(4,
              numericInput("lh_gem_oup1.width", "Width:", min = 1, max = 50, value = 25, step = 2)
            ),
            column(4,
              numericInput("lh_gem_oup1.res", "Res:", min = 72, max = 600, value = 150, step = 5)
            )
            ),
            downloadButton("lh_gem_oup1.png", "Download PNG", class = "btn-sm"),
            downloadButton("lh_gem_oup1.pdf", "Download PDF", class = "btn-sm"),
            downloadButton("lh_gem_oup1.svg", "Download SVG", class = "btn-sm")
          )
        )
               )
      ),
      column(8,
             uiOutput("lh_gem_oup1.ui")
      )
      ),
      hr()
    )
  )
) # End of tab gem

,
# tab gec ----
tabPanel(
  "Gene coexpression",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Coexpression of two genes on reduced dimensions"),
          p("Visualise the coexpression of two genes on low-dimensional representions.")
        )
      ),
      # row 2 ----
      fluidRow(
        class = "tab-section",
        # row 2 col 1
        column(4,
               column(
                 12,
                 div(
                   class = "input-panel input-panel-section",
                   h4("Dimension Reduction"),
                   selectInput("lh_gec_drX", "X-axis:",
                               choices = lhconf[dimred == TRUE]$UI,
                               selected = lhdef$dimred[1]
                   ),
                   selectInput("lh_gec_drY", "Y-axis:",
                               choices = lhconf[dimred == TRUE]$UI,
                               selected = lhdef$dimred[2]
                   ),
                   selectInput("lh_gec_inp1", "Gene 1:", choices = NULL) %>%
                     helper(
                       type = "inline", size = "m", fade = TRUE,
                       title = "Gene expression to colour cells",
                       content = c(
                         "- Select gene to colour cells by gene expression",
                         "- Type in gene names for unlisted genes",
                         "- Gene expression are coloured in a White-Red colour scheme which can be changed in the plot controls"
                       )
                     ),
                   selectInput("lh_gec_inp2", "Gene 2:", choices = NULL) %>%
                     helper(
                       type = "inline", size = "m", fade = TRUE,
                       title = "Gene expression to colour cells",
                       content = c(
                         "- Select gene to colour cells by gene expression",
                         "- Type in gene names for unlisted genes",
                         "- Gene expression are coloured in a White-Blue colour scheme which can be changed in the plot controls"
                       )
                     ),
                   checkboxInput("lh_gec_togL", "Subset cells"),
                   conditionalPanel(
                    condition = "input.lh_gec_togL == true",
                    selectInput("lh_gec_sub1", "Cell info to subset:", choices = lhconf[grp == TRUE]$UI, selected = lhdef$grp1),
                     uiOutput("lh_gec_sub1.ui"),
                     actionButton("lh_gec_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
                     actionButton("lh_gec_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
                   ),
                   checkboxInput("lh_gec_tog0", "Adjust graphics"),
                   conditionalPanel(
                     condition = "input.lh_gec_tog0 == true",
                     radioButtons("lh_gec_col1", "Colour:",
                                  choices = c(
                                    "Red (Gene1); Blue (Gene2)",
                                    "Orange (Gene1); Blue (Gene2)",
                                    "Red (Gene1); Green (Gene2)",
                                    "Green (Gene1); Blue (Gene2)"
                                  ),
                                  selected = "Red (Gene1); Blue (Gene2)"
                     ),
                     radioButtons("lh_gec_ord1", "Plot order:",
                                  choices = c("Max", "Min", "Original", "Random"),
                                  selected = "Max", inline = TRUE
                     ),
                     sliderInput("lh_gec_siz", "Point size:",
                                 min = 0, max = 4, value = 1.25, step = 0.25
                     ),
                     radioButtons("lh_gec_psz", "Plot size:",
                                  choices = c("Small", "Medium", "Large"),
                                  selected = "Medium", inline = TRUE
                     ),
                     radioButtons("lh_gec_fsz", "Font size:",
                                  choices = c("Small", "Medium", "Large"),
                                  selected = "Small", inline = TRUE
                     ),
                     radioButtons("lh_gec_asp", "Aspect ratio:",
                                  choices = c("Square", "Fixed", "Free"),
                                  selected = "Square", inline = TRUE
                     ),
                     checkboxInput("lh_gec_txt", "Show axis text", value = FALSE)
                   )
                 ),
                 div(class="input-panel input-panel-section",
                     numericInput("lh_gec_oup1.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                     downloadButton("lh_gec_oup1.png", "Download PNG", class = "btn-sm"),
                     downloadButton("lh_gec_oup1.pdf", "Download PDF", class = "btn-sm"),
                     downloadButton("lh_gec_oup1.svg", "Download SVG", class = "btn-sm"),
                 ),
                 div(class="input-panel-section",
                     h4("Cell numbers"),
                     dataTableOutput("lh_gec_.dt")
                 )
               )
        ), # row 2 col 1
        # row 2 col 2
        column(
          8,
          uiOutput("lh_gec_oup1.ui"),
        )
      ), # end of row 2
      hr()
    ) # col
  ) # row
) # End of tab gec

,
# tab vio ----
tabPanel(
  "Violinplot / Boxplot",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Cell information / gene expression violin plot / box plot"),
          p("Visualise the gene expression or continuous cell information (e.g. Number of UMIs / module score) across groups of cells (e.g. libary / clusters).")
        )
      ),
      # row 2 ----
      fluidRow(
        class = "tab-section",
        # row 2 col 1
        column(
          3,
          div(
            class = "input-panel",
            style = "border-right: 2px solid #f3f6f4",
            selectInput("lh_vio_inp1", "Cell info (X-axis):",
              choices = lhconf[grp == TRUE]$UI,
              selected = lhdef$grp1
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell info to group cells",
                content = c(
                  "- Select categorical cell info to group cells by",
                  "- Single cells are grouped by this categorical covariate",
                  "- Plotted as the X-axis of the violin plot / box plot"
                )
              ),
            selectInput("lh_vio_inp2", "Cell info / Gene (Y-axis):", choices = NULL) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell Info / Gene to plot",
                content = c(
                  "- Select cell info / gene to plot on Y-axis",
                  "- Can be continuous cell info (e.g. nUMIs / scores)",
                  "- Can also be gene expression",
                  "- Type in gene names for unlisted genes"
                )
              ),
            radioButtons("lh_vio_typ", "Plot type:",
              choices = c("violin", "boxplot", "lineplot"),
              selected = "violin", inline = TRUE
            ),
            checkboxInput("lh_vio_pts", "Show data points", value = FALSE),
            checkboxInput("lh_vio_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.lh_vio_togL == true",
              selectInput("lh_vio_sub1", "Cell info to subset:",
                choices = lhconf[grp == TRUE]$UI,
                selected = lhdef$grp1
              ),
              uiOutput("lh_vio_sub1.ui"),
              actionButton("lh_vio_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("lh_vio_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            ),
            checkboxInput("lh_vio_tog", "Adjust graphics"),
            conditionalPanel(
              condition = "input.lh_vio_tog == true",
              sliderInput("lh_vio_siz", "Data point size:",
                min = 0, max = 4, value = 1.25, step = 0.25
              ),
              radioButtons("lh_vio_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              radioButtons("lh_vio_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              conditionalPanel(
              condition = "input.lh_vio_typ == 'lineplot'",
              sliderInput("lh_vio_barsz", "Line size", min = 0.05, max = 0.5, step = 0.01, value = 0.3)
              )
            ),
            selectInput("lh_vio_datatype", "Data type", choices = c("normalised", "raw"), selected = "normalised"),
          ),
          div(
            class = "input-panel",
            numericInput("lh_vio_oup.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
            downloadButton("lh_vio_oup.png", "Download PNG", class = "btn-sm"),
            downloadButton("lh_vio_oup.pdf", "Download PDF", class = "btn-sm"),
            downloadButton("lh_vio_oup.svg", "Download SVG", class = "btn-sm")
          )
        ), # row 2 col 1
        # row 2 col 2
        column(
          9, uiOutput("lh_vio_oup.ui")
        ) # row 2 col 2
      ), # row 2
      hr()
    )
  )
) # End of tab vio

,
# tab pro ----
tabPanel(
  "Proportion plot",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Proportion / cell numbers across different cell information"),
          p("Visualise the composition of single cells based on one discrete cell information across another discrete cell information.")
        )
      ),
      # row 2 ----
      fluidRow(
        class = "tab-section",
        # row 2 col 1
        column(
          3,
          div(
            class = "input-panel",
            selectInput("lh_pro_inp1", "Cell info to plot (X-axis):",
              choices = lhconf[grp == TRUE]$UI,
              selected = lhdef$grp2
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell info to group cells",
                content = c(
                  "- Select categorical cell info to group cells",
                  "- Plotted as the X-axis of the proportion plot"
                )
              ),
            selectInput("lh_pro_inp2", "Cell info to group / colour by:",
              choices = lhconf[grp == TRUE]$UI,
              selected = lhdef$grp1
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell info to group / colour cells",
                content = c(
                  "- Select categorical cell info to group / colour cells",
                  "- Proportion / cell numbers are shown in different colours"
                )
              ),
            radioButtons("lh_pro_typ", "Plot value:",
              choices = c("Proportion", "CellNumbers"),
              selected = "Proportion", inline = TRUE
            ),
            checkboxInput("lh_pro_flp", "Flip X/Y", value = FALSE),
            checkboxInput("lh_pro_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.lh_pro_togL == true",
              selectInput("lh_pro_sub1", "Cell info to subset:",
                choices = lhconf[grp == TRUE]$UI,
                selected = lhdef$grp1
              ),
              uiOutput("lh_pro_sub1.ui"),
              actionButton("lh_pro_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("lh_pro_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            ),
            checkboxInput("lh_pro_tog", "Adjust graphics"),
            conditionalPanel(
              condition = "input.lh_pro_tog == true",
              radioButtons("lh_pro_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("lh_pro_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              )
            )
          ),
          div(
            class = "input-panel",
            numericInput("lh_pro_oup.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
            downloadButton("lh_pro_oup.png", "Download PNG", class = "btn-sm"),
            downloadButton("lh_pro_oup.pdf", "Download PDF", class = "btn-sm"),
            downloadButton("lh_pro_oup.svg", "Download SVG", class = "btn-sm")
          )
        ), # row 2 col 1
        # row 2 col 2
        column(
          9, uiOutput("lh_pro_oup.ui")
        ) # row 2 col 2
      ), # row 2
      hr()
    )
  )
) # End of tab pro

,
# tab hea ----
tabPanel(
  "Bubbleplot / Heatmap",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Gene expression bubbleplot / heatmap"),
          p("Visualise the gene expression patterns of multiple genes grouped by categorical cell information (e.g. library / cluster). The normalised expression are averaged, log-transformed and then plotted.")
        )
      ),
      # row 2 ----
      fluidRow(
        class = "tab-section",
        column(
          4,
          div(
            class = "input-panel",
            selectInput("lh_hea_inp", "Genes:", choices = NULL, multiple = TRUE) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "List of genes to plot on bubbleplot / heatmap",
                content = c(
                  "- Input genes to plot",
                  "- Type in gene names for unlisted genes"
                )
              ),
            selectInput("lh_hea_grp", "Group by:",
              choices = lhconf[grp == TRUE]$UI,
              selected = lhconf[grp == TRUE]$UI[1]
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell info to group cells",
                content = c(
                  "- Select categorical cell info to group cells",
                  "- Single cells are grouped by this categorical covariate",
                  "- Plotted as the X-axis of the bubbleplot / heatmap"
                )
              ),
            radioButtons("lh_hea_plt", "Plot type:",
              choices = c("Bubbleplot", "Heatmap"),
              selected = "Bubbleplot", inline = TRUE
            ),
            checkboxInput("lh_hea_scl", "Scale gene expression", value = TRUE),
            checkboxInput("lh_hea_row", "Cluster rows (genes)", value = TRUE),
            checkboxInput("lh_hea_col", "Cluster columns (samples)", value = FALSE),
            checkboxInput("lh_hea_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.lh_hea_togL == true",
              selectInput("lh_hea_sub1", "Cell info to subset:",
                choices = lhconf[grp == TRUE]$UI,
                selected = lhdef$grp1
              ),
              uiOutput("lh_hea_sub1.ui"),
              actionButton("lh_hea_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("lh_hea_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            ),
            checkboxInput("lh_hea_tog", "Adjust graphics"),
            conditionalPanel(
              condition = "input.lh_hea_tog == true",
              radioButtons("lh_hea_cols", "Colour scheme:",
                choices = c(
                  "White-Red", "Blue-Yellow-Red",
                  "Yellow-Green-Purple"
                ),
                selected = "Blue-Yellow-Red"
              ),
              radioButtons("lh_hea_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("lh_hea_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              )
            )
          ),
          div(
            class = "input-panel",
            fluidRow(
              column(4,
                     numericInput("lh_hea_oup.height", "Height:", min = 5, max = 100, value = 18, step = 2)
              ),
              column(4,
                     numericInput("lh_hea_oup.width", "Width:", min = 5, max = 100, value = 18, step = 2)
              ),
              column(4,
                     numericInput("lh_hea_oup.res", "Res:", min = 72, max = 600, value = 150, step = 5)
              )
            ),
            downloadButton("lh_hea_oup.png", "Download PNG", class = "btn-sm"),
            downloadButton("lh_hea_oup.pdf", "Download PDF", class = "btn-sm"),
            downloadButton("lh_hea_oup.svg", "Download SVG", class = "btn-sm")
          )
        ), # row 2 col 1
        # row 2 col 2
        column(
          8, h4(htmlOutput("lh_hea_oupTxt")),
          uiOutput("lh_hea_oup.ui")
        ) # row 2 col 2
      ), # row 2
      hr()
    )
  )
) # End of tab hea

,
# tab mar ----
tabPanel(
  "Markers",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Markers"),
          p("Explore markers for different clustering.")
        ) # row 1 col 1
      ), # row 1
      # row 2 ----
      fluidRow(
        class = "tab-section",
        column(3,
               fluidRow(
                 column(
                   12,
                   div(
                     class = "input-panel input-panel-section",
                     selectInput("lh_mar_cls","Select clustering:", choices = names(lhmar),selected = 1)
                   )
                 )
               )
        )
      ), # end of row 2
      # row 3 ----
      fluidRow(
      column(12,
        DTOutput("lh_mar_table")
      )
      ), # end of row 3
      hr()
    )
  )
) # End of tab mar

)
,navbarMenu("Valve Healthy",
# tab civge ----
tabPanel(
  "CellInfo vs GeneExpr",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Cell information vs Gene expression"),
          p("Cell information and gene expression side-by-side on low-dimensional represention.")
        )
      ),
      # row 2 ----
      fluidRow(
        class = "tab-section",
        # row 2 col 1
        column(
          4,
          div(
            class = "input-panel",
            h4("Dimension Reduction"),
            selectInput("vh_civge_drX", "X-axis:",
              choices = vhconf[dimred == TRUE]$UI,
              selected = vhdef$dimred[1]
            ),
            selectInput("vh_civge_drY", "Y-axis:",
              choices = vhconf[dimred == TRUE]$UI,
              selected = vhdef$dimred[2]
            )
          )
        ), # row 1 col 1
        # row 2 col 2
        column(
          4,
          div(
            class = "input-panel",
            checkboxInput("vh_civge_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.vh_civge_togL == true",
              selectInput("vh_civge_sub1", "Cell info to subset:",
                choices = vhconf[grp == TRUE]$UI,
                selected = vhdef$grp1
              ),
              uiOutput("vh_civge_sub1.ui"),
              actionButton("vh_civge_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("vh_civge_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            )
          )
        ), # End of column
        # row 2 col 3
        column(
          4,
          div(
            class = "input-panel",
            checkboxInput("vh_civge_tog0", "Adjust graphics"),
            conditionalPanel(
              condition = "input.vh_civge_tog0 == true",
              sliderInput("vh_civge_siz", "Point size:",
                min = 0, max = 4, value = 1.25, step = 0.25
              ),
              radioButtons("vh_civge_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("vh_civge_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              radioButtons("vh_civge_asp", "Aspect ratio:",
                choices = c("Square", "Fixed", "Free"),
                selected = "Square", inline = TRUE
              ),
              checkboxInput("vh_civge_txt", "Show axis text", value = FALSE)
            )
          )
        ) # row 2 col 3
      ), # row 2
      # row 3 ----
      fluidRow(
        class = "tab-section",
        # row 3 col 1
        column(6,
          style = "border-right: 2px solid #f3f6f4",
          h4("Cell information"),
          # row 3 col 1 row 1
          fluidRow(
            class = "tab-section",
            column(
              6,
              div(
                class = "input-panel",
                selectInput("vh_civge_inp1", "Cell info:",
                  choices = vhconf$UI,
                  selected = vhdef$meta1
                ) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Cell info to colour cells",
                    content = c(
                      "- Select cell info to colour cells",
                      "- Categorical covariates have a fixed colour palette",
                      "- Continuous covariates are coloured in a Blue-Yellow-Red colour scheme, which can be changed in the plot controls"
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("vh_civge_tog1", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.vh_civge_tog1 == true",
                  radioButtons("vh_civge_col1", "Colour (Continuous data):",
                    choices = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple"),
                    selected = "Blue-Yellow-Red"
                  ),
                  radioButtons("vh_civge_ord1", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Original", inline = TRUE
                  ),
                  checkboxInput("vh_civge_lab1", "Show cell info labels", value = TRUE)
                )
              )
            )
          ),
          # row 3 col 1 row 2
          fluidRow(
            class = "tab-section",
            column(
              12,
              uiOutput("vh_civge_oup1.ui")
            )
          ),
          # row 3 col 1 row 3
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("vh_civge_oup1.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("vh_civge_oup1.png", "Download PNG", class = "btn-sm"),
                downloadButton("vh_civge_oup1.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("vh_civge_oup1.svg", "Download SVG", class = "btn-sm"),
                checkboxInput("vh_civge_tog9", "Show cell numbers / statistics")
              )
            )
          ),
          conditionalPanel(
            condition = "input.vh_civge_tog9 == true",
            h4("Cell numbers / statistics"),
            radioButtons("vh_civge_splt", "Split continuous cell info into:",
              choices = c("Quartile", "Decile"),
              selected = "Decile", inline = TRUE
            ),
            dataTableOutput("vh_civge_.dt")
          )
        ), # row 3 col 1
        # row 3 col 2
        column(
          6,
          h4("Gene expression"),
          fluidRow(
            class = "tab-section",
            column(
              6,
              div(
                class = "input-panel",
                selectInput("vh_civge_inp2", "Gene name:", choices = NULL) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Gene expression to colour cells",
                    content = c(
                      "- Select gene to colour cells by gene expression",
                      "- Type in gene names for unlisted genes",
                      "- Gene expression are coloured in a White-Red colour scheme which can be changed in the plot controls"
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("vh_civge_tog2", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.vh_civge_tog2 == true",
                  radioButtons("vh_civge_col2", "Colour:",
                    choices = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple"),
                    selected = "White-Red"
                  ),
                  radioButtons("vh_civge_ord2", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Max", inline = TRUE
                  )
                )
              )
            )
          ),
          fluidRow(
            class = "tab-section",
            column(
              12,
              uiOutput("vh_civge_oup2.ui")
            )
          ),
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("vh_civge_oup2.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("vh_civge_oup2.png", "Download PNG", class = "btn-sm"),
                downloadButton("vh_civge_oup2.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("vh_civge_oup2.svg", "Download SVG", class = "btn-sm")
              )
            )
          )
        ) # row 3 col 2
      ), # row 3
      hr()
    )
  )
) # End of tab civge

,
# tab civci ----
tabPanel(
  "CellInfo vs CellInfo",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Cell info vs cell info"),
          p("Two cell infos side-by-side on low-dimensional represention.")
        ) # row 1 col 1
      ), # row 1
      # row 2 ----
      fluidRow(
        class = "tab-section",
        column(
          4,
          div(
            class = "input-panel",
            h4("Dimension Reduction"),
            selectInput("vh_civci_drX", "X-axis:",
              choices = vhconf[dimred == TRUE]$UI,
              selected = vhdef$dimred[1]
            ),
            selectInput("vh_civci_drY", "Y-axis:",
              choices = vhconf[dimred == TRUE]$UI,
              selected = vhdef$dimred[2]
            )
          )
        ), # row 2 col 2
        # row 2 col 2
        column(4,
          div(
            class = "input-panel",
            checkboxInput("vh_civci_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.vh_civci_togL == true",
              selectInput("vh_civci_sub1", "Cell info to subset:",
                choices = vhconf[grp == TRUE]$UI,
                selected = vhdef$grp1
              ),
              uiOutput("vh_civci_sub1.ui"),
              actionButton("vh_civci_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("vh_civci_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            )
          )
        ), # row 2 col 2
        # row 2 col 3
        column(
          4,
          div(
            class = "input-panel",
            checkboxInput("vh_civci_tog0", "Adjust graphics"),
            conditionalPanel(
              condition = "input.vh_civci_tog0 == true",
              sliderInput("vh_civci_siz", "Point size:",
                min = 0, max = 4, value = 1.25, step = 0.25
              ),
              radioButtons("vh_civci_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("vh_civci_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              radioButtons("vh_civci_asp", "Aspect ratio:",
                choices = c("Square", "Fixed", "Free"),
                selected = "Square", inline = TRUE
              ),
              checkboxInput("vh_civci_txt", "Show axis text", value = FALSE)
            )
          )
        ) # row 2 col 3
      ), # row 2
      # row 3 ----
      fluidRow(
        class = "tab-section",
        # row 3 col 1
        column(
          6,
          style = "border-right: 2px solid #f3f6f4",
          h4("Cell info 1"),
          # row 3 col 1 row 1
          fluidRow(
            class = "tab-section",
            column(
              6,
              div(
                class = "input-panel",
                selectInput("vh_civci_inp1", "Cell info:",
                  choices = vhconf$UI,
                  selected = vhdef$meta1
                ) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Cell info to colour cells",
                    content = c(
                      "- Select cell info to colour cells",
                      "- Categorical covariates have a fixed colour palette",
                      "- Continuous covariates are coloured in a Blue-Yellow-Red colour scheme, which can be changed in the plot controls"
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("vh_civci_tog1", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.vh_civci_tog1 == true",
                  radioButtons("vh_civci_col1", "Colour (Continuous data):",
                    choices = c(
                      "White-Red", "Blue-Yellow-Red",
                      "Yellow-Green-Purple"
                    ),
                    selected = "Blue-Yellow-Red"
                  ),
                  radioButtons("vh_civci_ord1", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Original", inline = TRUE
                  ),
                  checkboxInput("vh_civci_lab1", "Show cell info labels", value = TRUE)
                )
              )
            )
          ),
          # row 3 col 1 row 2
          fluidRow(column(12, uiOutput("vh_civci_oup1.ui"))),
          # row 3 col 1 row 3
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("vh_civci_oup1.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("vh_civci_oup1.png", "Download PNG", class = "btn-sm"),
                downloadButton("vh_civci_oup1.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("vh_civci_oup1.svg", "Download svg", class = "btn-sm")
              )
            )
          )
        ), # row 3 col 1
        # row 3 col 2
        column(
          6, h4("Cell info 2"),
          fluidRow(
            class = "tab-section",
            column(
              6,
              div(
                class = "input-panel",
                selectInput("vh_civci_inp2", "Cell info:",
                  choices = vhconf$UI,
                  selected = vhdef$meta2
                ) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Cell info to colour cells",
                    content = c(
                      "- Select cell info to colour cells",
                      "- Categorical covariates have a fixed colour palette",
                      "- Continuous covariates are coloured in a Blue-Yellow-Red colour scheme, which can be changed in the plot controls"
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("vh_civci_tog2", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.vh_civci_tog2 == true",
                  radioButtons("vh_civci_col2", "Colour (Continuous data):",
                    choices = c(
                      "White-Red", "Blue-Yellow-Red",
                      "Yellow-Green-Purple"
                    ),
                    selected = "Blue-Yellow-Red"
                  ),
                  radioButtons("vh_civci_ord2", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Original", inline = TRUE
                  ),
                  checkboxInput("vh_civci_lab2", "Show cell info labels", value = TRUE)
                )
              )
            )
          ),
          fluidRow(column(12, uiOutput("vh_civci_oup2.ui"))),
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("vh_civci_oup2.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("vh_civci_oup2.png", "Download PNG", class = "btn-sm"),
                downloadButton("vh_civci_oup2.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("vh_civci_oup2.svg", "Download SVG", class = "btn-sm")
              )
            )
          )
        ) # row 3 col 2
      ), # row 3
      hr()
    )
  )
) # End of tab civci

,
# tab gevge ----
tabPanel(
  "GeneExpr vs GeneExpr",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Gene expression vs Gene expression"),
          p("Visualise two gene expressions side-by-side on low-dimensional representions.")
        )
      ),
      # row 2 ----
      fluidRow(
        class = "tab-section",
        # row 2 col 1
        column(
          4,
          div(
            class = "input-panel",
            h4("Dimension Reduction"),
            selectInput("vh_gevge_drX", "X-axis:",
              choices = vhconf[dimred == TRUE]$UI,
              selected = vhdef$dimred[1]
            ),
            selectInput("vh_gevge_drY", "Y-axis:",
              choices = vhconf[dimred == TRUE]$UI,
              selected = vhdef$dimred[2]
            )
          )
        ), # row 1 col 1
        # row 2 col 2
        column(4,
          div(
            class = "input-panel",
            checkboxInput("vh_gevge_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.vh_gevge_togL == true",
              selectInput("vh_gevge_sub1", "Cell info to subset:",
                choices = vhconf[grp == TRUE]$UI,
                selected = vhdef$grp1
              ),
              uiOutput("vh_gevge_sub1.ui"),
              actionButton("vh_gevge_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("vh_gevge_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            )
          )
        ), # End of column
        # row 2 col 3
        column(
          4,
          div(
            class = "input-panel",
            checkboxInput("vh_gevge_tog0", "Adjust graphics"),
            conditionalPanel(
              condition = "input.vh_gevge_tog0 == true",
              sliderInput("vh_gevge_siz", "Point size:",
                min = 0, max = 4, value = 1.25, step = 0.25
              ),
              radioButtons("vh_gevge_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("vh_gevge_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              radioButtons("vh_gevge_asp", "Aspect ratio:",
                choices = c("Square", "Fixed", "Free"),
                selected = "Square", inline = TRUE
              ),
              checkboxInput("vh_gevge_txt", "Show axis text", value = FALSE)
            )
          )
        ) # row 2 col 3
      ), # row 2
      # row 3 ----
      fluidRow(
        class = "tab-section",
        column(
          6,
          style = "border-right: 2px solid #f3f6f4",
          h4("Gene expression 1"),
          fluidRow(
            class = "tab-section",
            column(
              6,
              div(
                class = "input-panel",
                selectInput("vh_gevge_inp1", "Gene name:", choices = NULL) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Gene expression to colour cells",
                    content = c(
                      "- Select gene to colour cells by gene expression",
                      "- Type in gene names for unlisted genes",
                      "- Gene expression are coloured in a White-Red colour scheme which can be changed in the plot controls"
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("vh_gevge_tog1", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.vh_gevge_tog1 == true",
                  radioButtons("vh_gevge_col1", "Colour:",
                    choices = c(
                      "White-Red", "Blue-Yellow-Red",
                      "Yellow-Green-Purple"
                    ),
                    selected = "White-Red"
                  ),
                  radioButtons("vh_gevge_ord1", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Max", inline = TRUE
                  )
                )
              )
            )
          ),
          # row 3 col 1 row 2
          fluidRow(
            class = "tab-section",
            column(12, uiOutput("vh_gevge_oup1.ui"))
          ),
          # row 3 col 1 row 3
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("vh_gevge_oup1.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("vh_gevge_oup1.png", "Download PNG", class = "btn-sm"),
                downloadButton("vh_gevge_oup1.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("vh_gevge_oup1.svg", "Download SVG", class = "btn-sm")
              )
            )
          )
        ), # row 3 col 1
        # row 3 col 2
        column(
          6, h4("Gene expression 2"),
          fluidRow(
            class = "tab-section",
            column(
              6,
              div(
                class = "input-panel",
                selectInput("vh_gevge_inp2", "Gene name:", choices = NULL) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Gene expression to colour cells",
                    content = c(
                      "- Select gene to colour cells by gene expression",
                      "- Type in gene names for unlisted genes",
                      "- Gene expression are coloured in a White-Red colour scheme which can be changed in the plot controls"
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("vh_gevge_tog2", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.vh_gevge_tog2 == true",
                  radioButtons("vh_gevge_col2", "Colour:",
                    choices = c(
                      "White-Red", "Blue-Yellow-Red",
                      "Yellow-Green-Purple"
                    ),
                    selected = "White-Red"
                  ),
                  radioButtons("vh_gevge_ord2", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Max", inline = TRUE
                  )
                )
              )
            )
          ),
          fluidRow(
            class = "tab-section",
            column(12, uiOutput("vh_gevge_oup2.ui"))
          ),
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                  numericInput("vh_gevge_oup2.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                  downloadButton("vh_gevge_oup2.png", "Download PNG", class = "btn-sm"),
                  downloadButton("vh_gevge_oup2.pdf", "Download PDF", class = "btn-sm"),
                  downloadButton("vh_gevge_oup2.svg", "Download SVG", class = "btn-sm")
              )
            )
          )
        ) # row 3 col 2
      ), # row 3
      hr()
    )
  )
) # End of tab gevge

,
# tab gem ----
tabPanel(
  "Expression",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Gene expression"),
          p("Explore gene expression on low-dimensional represention.")
        ) # row 1 col 1
      ), # row 1
      # row 2 ----
      fluidRow(
        class = "tab-section",
        column(4,
               fluidRow(
        column(
          12,
          div(
            class = "input-panel input-panel-section",
            selectInput("vh_gem_inp", "Genes:", choices = NULL, multiple = TRUE) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "List of genes to plot on bubbleplot / heatmap",
                content = c(
                  "- Input genes to plot",
                  "- Type in gene names for unlisted genes"
                )
              ),
            selectInput("vh_gem_drX", "X-axis:",
                        choices = vhconf[dimred == TRUE]$UI,
                        selected = vhdef$dimred[1]
            ),
            selectInput("vh_gem_drY", "Y-axis:",
                        choices = vhconf[dimred == TRUE]$UI,
                        selected = vhdef$dimred[2]
            ),
            checkboxInput("vh_gem_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.vh_gem_togL == true",
              selectInput("vh_gem_sub1", "Cell info to subset:",
                          choices = vhconf[grp == TRUE]$UI,
                          selected = vhdef$grp1
              ),
              uiOutput("vh_gem_sub1.ui"),
              actionButton("vh_gem_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("vh_gem_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            ),
            checkboxInput("vh_gem_tog0", "Adjust graphics"),
            conditionalPanel(
              condition = "input.vh_gem_tog0 == true",
              sliderInput("vh_gem_siz", "Point size:",
                          min = 0, max = 3, value = 0.5, step = 0.1
              ),
              radioButtons("vh_gem_psz", "Plot size:",
                           choices = c("Small", "Medium", "Large"),
                           selected = "Medium", inline = TRUE
              ),
              radioButtons("vh_gem_fsz", "Font size:",
                           choices = c("Smaller", "Small", "Medium", "Large"),
                           selected = "Small", inline = TRUE
              ),
              radioButtons("vh_gem_asp", "Aspect ratio:",
                           choices = c("Square", "Fixed", "Free"),
                           selected = "Square", inline = TRUE
              ),
              checkboxInput("vh_gem_txt", "Show axis text", value = FALSE),
              radioButtons("vh_gem_col", "Colour (Continuous data):",
                           choices = c(
                             "White-Red", "Blue-Yellow-Red",
                             "Yellow-Green-Purple"
                           ),
                           selected = "Blue-Yellow-Red"
              ),
              radioButtons("vh_gem_ord", "Plot order:",
                           choices = c("Max", "Min", "Original", "Random"),
                           selected = "Max", inline = TRUE
              ),
              numericInput("vh_gem_ncol", "Number of columns", value = 0, min = 0, step = 1)
            )
          ),
          div(
            class = "input-panel",
            fluidRow(
            column(4,
              numericInput("vh_gem_oup1.height", "Height:", min = 1, max = 50, value = 25, step = 2)
            ),
            column(4,
              numericInput("vh_gem_oup1.width", "Width:", min = 1, max = 50, value = 25, step = 2)
            ),
            column(4,
              numericInput("vh_gem_oup1.res", "Res:", min = 72, max = 600, value = 150, step = 5)
            )
            ),
            downloadButton("vh_gem_oup1.png", "Download PNG", class = "btn-sm"),
            downloadButton("vh_gem_oup1.pdf", "Download PDF", class = "btn-sm"),
            downloadButton("vh_gem_oup1.svg", "Download SVG", class = "btn-sm")
          )
        )
               )
      ),
      column(8,
             uiOutput("vh_gem_oup1.ui")
      )
      ),
      hr()
    )
  )
) # End of tab gem

,
# tab gec ----
tabPanel(
  "Gene coexpression",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Coexpression of two genes on reduced dimensions"),
          p("Visualise the coexpression of two genes on low-dimensional representions.")
        )
      ),
      # row 2 ----
      fluidRow(
        class = "tab-section",
        # row 2 col 1
        column(4,
               column(
                 12,
                 div(
                   class = "input-panel input-panel-section",
                   h4("Dimension Reduction"),
                   selectInput("vh_gec_drX", "X-axis:",
                               choices = vhconf[dimred == TRUE]$UI,
                               selected = vhdef$dimred[1]
                   ),
                   selectInput("vh_gec_drY", "Y-axis:",
                               choices = vhconf[dimred == TRUE]$UI,
                               selected = vhdef$dimred[2]
                   ),
                   selectInput("vh_gec_inp1", "Gene 1:", choices = NULL) %>%
                     helper(
                       type = "inline", size = "m", fade = TRUE,
                       title = "Gene expression to colour cells",
                       content = c(
                         "- Select gene to colour cells by gene expression",
                         "- Type in gene names for unlisted genes",
                         "- Gene expression are coloured in a White-Red colour scheme which can be changed in the plot controls"
                       )
                     ),
                   selectInput("vh_gec_inp2", "Gene 2:", choices = NULL) %>%
                     helper(
                       type = "inline", size = "m", fade = TRUE,
                       title = "Gene expression to colour cells",
                       content = c(
                         "- Select gene to colour cells by gene expression",
                         "- Type in gene names for unlisted genes",
                         "- Gene expression are coloured in a White-Blue colour scheme which can be changed in the plot controls"
                       )
                     ),
                   checkboxInput("vh_gec_togL", "Subset cells"),
                   conditionalPanel(
                    condition = "input.vh_gec_togL == true",
                    selectInput("vh_gec_sub1", "Cell info to subset:", choices = vhconf[grp == TRUE]$UI, selected = vhdef$grp1),
                     uiOutput("vh_gec_sub1.ui"),
                     actionButton("vh_gec_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
                     actionButton("vh_gec_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
                   ),
                   checkboxInput("vh_gec_tog0", "Adjust graphics"),
                   conditionalPanel(
                     condition = "input.vh_gec_tog0 == true",
                     radioButtons("vh_gec_col1", "Colour:",
                                  choices = c(
                                    "Red (Gene1); Blue (Gene2)",
                                    "Orange (Gene1); Blue (Gene2)",
                                    "Red (Gene1); Green (Gene2)",
                                    "Green (Gene1); Blue (Gene2)"
                                  ),
                                  selected = "Red (Gene1); Blue (Gene2)"
                     ),
                     radioButtons("vh_gec_ord1", "Plot order:",
                                  choices = c("Max", "Min", "Original", "Random"),
                                  selected = "Max", inline = TRUE
                     ),
                     sliderInput("vh_gec_siz", "Point size:",
                                 min = 0, max = 4, value = 1.25, step = 0.25
                     ),
                     radioButtons("vh_gec_psz", "Plot size:",
                                  choices = c("Small", "Medium", "Large"),
                                  selected = "Medium", inline = TRUE
                     ),
                     radioButtons("vh_gec_fsz", "Font size:",
                                  choices = c("Small", "Medium", "Large"),
                                  selected = "Small", inline = TRUE
                     ),
                     radioButtons("vh_gec_asp", "Aspect ratio:",
                                  choices = c("Square", "Fixed", "Free"),
                                  selected = "Square", inline = TRUE
                     ),
                     checkboxInput("vh_gec_txt", "Show axis text", value = FALSE)
                   )
                 ),
                 div(class="input-panel input-panel-section",
                     numericInput("vh_gec_oup1.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                     downloadButton("vh_gec_oup1.png", "Download PNG", class = "btn-sm"),
                     downloadButton("vh_gec_oup1.pdf", "Download PDF", class = "btn-sm"),
                     downloadButton("vh_gec_oup1.svg", "Download SVG", class = "btn-sm"),
                 ),
                 div(class="input-panel-section",
                     h4("Cell numbers"),
                     dataTableOutput("vh_gec_.dt")
                 )
               )
        ), # row 2 col 1
        # row 2 col 2
        column(
          8,
          uiOutput("vh_gec_oup1.ui"),
        )
      ), # end of row 2
      hr()
    ) # col
  ) # row
) # End of tab gec

,
# tab vio ----
tabPanel(
  "Violinplot / Boxplot",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Cell information / gene expression violin plot / box plot"),
          p("Visualise the gene expression or continuous cell information (e.g. Number of UMIs / module score) across groups of cells (e.g. libary / clusters).")
        )
      ),
      # row 2 ----
      fluidRow(
        class = "tab-section",
        # row 2 col 1
        column(
          3,
          div(
            class = "input-panel",
            style = "border-right: 2px solid #f3f6f4",
            selectInput("vh_vio_inp1", "Cell info (X-axis):",
              choices = vhconf[grp == TRUE]$UI,
              selected = vhdef$grp1
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell info to group cells",
                content = c(
                  "- Select categorical cell info to group cells by",
                  "- Single cells are grouped by this categorical covariate",
                  "- Plotted as the X-axis of the violin plot / box plot"
                )
              ),
            selectInput("vh_vio_inp2", "Cell info / Gene (Y-axis):", choices = NULL) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell Info / Gene to plot",
                content = c(
                  "- Select cell info / gene to plot on Y-axis",
                  "- Can be continuous cell info (e.g. nUMIs / scores)",
                  "- Can also be gene expression",
                  "- Type in gene names for unlisted genes"
                )
              ),
            radioButtons("vh_vio_typ", "Plot type:",
              choices = c("violin", "boxplot", "lineplot"),
              selected = "violin", inline = TRUE
            ),
            checkboxInput("vh_vio_pts", "Show data points", value = FALSE),
            checkboxInput("vh_vio_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.vh_vio_togL == true",
              selectInput("vh_vio_sub1", "Cell info to subset:",
                choices = vhconf[grp == TRUE]$UI,
                selected = vhdef$grp1
              ),
              uiOutput("vh_vio_sub1.ui"),
              actionButton("vh_vio_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("vh_vio_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            ),
            checkboxInput("vh_vio_tog", "Adjust graphics"),
            conditionalPanel(
              condition = "input.vh_vio_tog == true",
              sliderInput("vh_vio_siz", "Data point size:",
                min = 0, max = 4, value = 1.25, step = 0.25
              ),
              radioButtons("vh_vio_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              radioButtons("vh_vio_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              conditionalPanel(
              condition = "input.vh_vio_typ == 'lineplot'",
              sliderInput("vh_vio_barsz", "Line size", min = 0.05, max = 0.5, step = 0.01, value = 0.3)
              )
            ),
            selectInput("vh_vio_datatype", "Data type", choices = c("normalised", "raw"), selected = "normalised"),
          ),
          div(
            class = "input-panel",
            numericInput("vh_vio_oup.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
            downloadButton("vh_vio_oup.png", "Download PNG", class = "btn-sm"),
            downloadButton("vh_vio_oup.pdf", "Download PDF", class = "btn-sm"),
            downloadButton("vh_vio_oup.svg", "Download SVG", class = "btn-sm")
          )
        ), # row 2 col 1
        # row 2 col 2
        column(
          9, uiOutput("vh_vio_oup.ui")
        ) # row 2 col 2
      ), # row 2
      hr()
    )
  )
) # End of tab vio

,
# tab pro ----
tabPanel(
  "Proportion plot",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Proportion / cell numbers across different cell information"),
          p("Visualise the composition of single cells based on one discrete cell information across another discrete cell information.")
        )
      ),
      # row 2 ----
      fluidRow(
        class = "tab-section",
        # row 2 col 1
        column(
          3,
          div(
            class = "input-panel",
            selectInput("vh_pro_inp1", "Cell info to plot (X-axis):",
              choices = vhconf[grp == TRUE]$UI,
              selected = vhdef$grp2
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell info to group cells",
                content = c(
                  "- Select categorical cell info to group cells",
                  "- Plotted as the X-axis of the proportion plot"
                )
              ),
            selectInput("vh_pro_inp2", "Cell info to group / colour by:",
              choices = vhconf[grp == TRUE]$UI,
              selected = vhdef$grp1
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell info to group / colour cells",
                content = c(
                  "- Select categorical cell info to group / colour cells",
                  "- Proportion / cell numbers are shown in different colours"
                )
              ),
            radioButtons("vh_pro_typ", "Plot value:",
              choices = c("Proportion", "CellNumbers"),
              selected = "Proportion", inline = TRUE
            ),
            checkboxInput("vh_pro_flp", "Flip X/Y", value = FALSE),
            checkboxInput("vh_pro_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.vh_pro_togL == true",
              selectInput("vh_pro_sub1", "Cell info to subset:",
                choices = vhconf[grp == TRUE]$UI,
                selected = vhdef$grp1
              ),
              uiOutput("vh_pro_sub1.ui"),
              actionButton("vh_pro_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("vh_pro_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            ),
            checkboxInput("vh_pro_tog", "Adjust graphics"),
            conditionalPanel(
              condition = "input.vh_pro_tog == true",
              radioButtons("vh_pro_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("vh_pro_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              )
            )
          ),
          div(
            class = "input-panel",
            numericInput("vh_pro_oup.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
            downloadButton("vh_pro_oup.png", "Download PNG", class = "btn-sm"),
            downloadButton("vh_pro_oup.pdf", "Download PDF", class = "btn-sm"),
            downloadButton("vh_pro_oup.svg", "Download SVG", class = "btn-sm")
          )
        ), # row 2 col 1
        # row 2 col 2
        column(
          9, uiOutput("vh_pro_oup.ui")
        ) # row 2 col 2
      ), # row 2
      hr()
    )
  )
) # End of tab pro

,
# tab hea ----
tabPanel(
  "Bubbleplot / Heatmap",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Gene expression bubbleplot / heatmap"),
          p("Visualise the gene expression patterns of multiple genes grouped by categorical cell information (e.g. library / cluster). The normalised expression are averaged, log-transformed and then plotted.")
        )
      ),
      # row 2 ----
      fluidRow(
        class = "tab-section",
        column(
          4,
          div(
            class = "input-panel",
            selectInput("vh_hea_inp", "Genes:", choices = NULL, multiple = TRUE) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "List of genes to plot on bubbleplot / heatmap",
                content = c(
                  "- Input genes to plot",
                  "- Type in gene names for unlisted genes"
                )
              ),
            selectInput("vh_hea_grp", "Group by:",
              choices = vhconf[grp == TRUE]$UI,
              selected = vhconf[grp == TRUE]$UI[1]
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell info to group cells",
                content = c(
                  "- Select categorical cell info to group cells",
                  "- Single cells are grouped by this categorical covariate",
                  "- Plotted as the X-axis of the bubbleplot / heatmap"
                )
              ),
            radioButtons("vh_hea_plt", "Plot type:",
              choices = c("Bubbleplot", "Heatmap"),
              selected = "Bubbleplot", inline = TRUE
            ),
            checkboxInput("vh_hea_scl", "Scale gene expression", value = TRUE),
            checkboxInput("vh_hea_row", "Cluster rows (genes)", value = TRUE),
            checkboxInput("vh_hea_col", "Cluster columns (samples)", value = FALSE),
            checkboxInput("vh_hea_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.vh_hea_togL == true",
              selectInput("vh_hea_sub1", "Cell info to subset:",
                choices = vhconf[grp == TRUE]$UI,
                selected = vhdef$grp1
              ),
              uiOutput("vh_hea_sub1.ui"),
              actionButton("vh_hea_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("vh_hea_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            ),
            checkboxInput("vh_hea_tog", "Adjust graphics"),
            conditionalPanel(
              condition = "input.vh_hea_tog == true",
              radioButtons("vh_hea_cols", "Colour scheme:",
                choices = c(
                  "White-Red", "Blue-Yellow-Red",
                  "Yellow-Green-Purple"
                ),
                selected = "Blue-Yellow-Red"
              ),
              radioButtons("vh_hea_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("vh_hea_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              )
            )
          ),
          div(
            class = "input-panel",
            fluidRow(
              column(4,
                     numericInput("vh_hea_oup.height", "Height:", min = 5, max = 100, value = 18, step = 2)
              ),
              column(4,
                     numericInput("vh_hea_oup.width", "Width:", min = 5, max = 100, value = 18, step = 2)
              ),
              column(4,
                     numericInput("vh_hea_oup.res", "Res:", min = 72, max = 600, value = 150, step = 5)
              )
            ),
            downloadButton("vh_hea_oup.png", "Download PNG", class = "btn-sm"),
            downloadButton("vh_hea_oup.pdf", "Download PDF", class = "btn-sm"),
            downloadButton("vh_hea_oup.svg", "Download SVG", class = "btn-sm")
          )
        ), # row 2 col 1
        # row 2 col 2
        column(
          8, h4(htmlOutput("vh_hea_oupTxt")),
          uiOutput("vh_hea_oup.ui")
        ) # row 2 col 2
      ), # row 2
      hr()
    )
  )
) # End of tab hea

,
# tab mar ----
tabPanel(
  "Markers",
  fluidRow(
    class = "container page",
    column(
      12,
      # row 1 ----
      fluidRow(
        class = "tab-section",
        column(
          12,
          h3("Markers"),
          p("Explore markers for different clustering.")
        ) # row 1 col 1
      ), # row 1
      # row 2 ----
      fluidRow(
        class = "tab-section",
        column(3,
               fluidRow(
                 column(
                   12,
                   div(
                     class = "input-panel input-panel-section",
                     selectInput("vh_mar_cls","Select clustering:", choices = names(vhmar),selected = 1)
                   )
                 )
               )
        )
      ), # end of row 2
      # row 3 ----
      fluidRow(
      column(12,
        DTOutput("vh_mar_table")
      )
      ), # end of row 3
      hr()
    )
  )
) # End of tab mar

)
,navbarMenu("About"
,
# about ----
tabPanel(
  "About",
  fluidRow(
    class = "container page",
    column(
      12,
      includeMarkdown("about.md")
    )
  )
)
))
)
)
