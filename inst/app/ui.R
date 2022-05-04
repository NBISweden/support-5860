library(shiny)
library(shinyhelper)
library(shinythemes)
library(showtext)
library(data.table)
library(Matrix)
library(DT)
library(magrittr)
if(!exists("sc1conf")) sc1conf = readRDS("sc1conf.rds")
if(!exists("sc1def")) sc1def  = readRDS("sc1def.rds")
if(!exists("sc1gene")) sc1gene = readRDS("sc1gene.rds")
if(!exists("sc1meta")) sc1meta = readRDS("sc1meta.rds")
if(!exists("sc1mar")) sc1mar = readRDS("sc1mar.rds")
if(!exists("sc2conf")) sc2conf = readRDS("sc2conf.rds")
if(!exists("sc2def")) sc2def  = readRDS("sc2def.rds")
if(!exists("sc2gene")) sc2gene = readRDS("sc2gene.rds")
if(!exists("sc2meta")) sc2meta = readRDS("sc2meta.rds")
if(!exists("sc2mar")) sc2mar = readRDS("sc2mar.rds")
if(!exists("sc3conf")) sc3conf = readRDS("sc3conf.rds")
if(!exists("sc3def")) sc3def  = readRDS("sc3def.rds")
if(!exists("sc3gene")) sc3gene = readRDS("sc3gene.rds")
if(!exists("sc3meta")) sc3meta = readRDS("sc3meta.rds")
if(!exists("sc3mar")) sc3mar = readRDS("sc3mar.rds")
if(!exists("sc4conf")) sc4conf = readRDS("sc4conf.rds")
if(!exists("sc4def")) sc4def  = readRDS("sc4def.rds")
if(!exists("sc4gene")) sc4gene = readRDS("sc4gene.rds")
if(!exists("sc4meta")) sc4meta = readRDS("sc4meta.rds")
if(!exists("sc4mar")) sc4mar = readRDS("sc4mar.rds")

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
            selectInput("sc1_civge_drX", "X-axis:",
              choices = sc1conf[dimred == TRUE]$UI,
              selected = sc1def$dimred[1]
            ),
            selectInput("sc1_civge_drY", "Y-axis:",
              choices = sc1conf[dimred == TRUE]$UI,
              selected = sc1def$dimred[2]
            )
          )
        ), # row 1 col 1
        # row 2 col 2
        column(
          4,
          div(
            class = "input-panel",
            checkboxInput("sc1_civge_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.sc1_civge_togL == true",
              selectInput("sc1_civge_sub1", "Cell info to subset:",
                choices = sc1conf[grp == TRUE]$UI,
                selected = sc1def$grp1
              ),
              uiOutput("sc1_civge_sub1.ui"),
              actionButton("sc1_civge_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("sc1_civge_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            )
          )
        ), # End of column
        # row 2 col 3
        column(
          4,
          div(
            class = "input-panel",
            checkboxInput("sc1_civge_tog0", "Adjust graphics"),
            conditionalPanel(
              condition = "input.sc1_civge_tog0 == true",
              sliderInput("sc1_civge_siz", "Point size:",
                min = 0, max = 4, value = 1.25, step = 0.25
              ),
              radioButtons("sc1_civge_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("sc1_civge_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              radioButtons("sc1_civge_asp", "Aspect ratio:",
                choices = c("Square", "Fixed", "Free"),
                selected = "Square", inline = TRUE
              ),
              checkboxInput("sc1_civge_txt", "Show axis text", value = FALSE)
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
                selectInput("sc1_civge_inp1", "Cell info:",
                  choices = sc1conf$UI,
                  selected = sc1def$meta1
                ) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Cell info to colour cells by",
                    content = c(
                      "Select cell info to colour cells",
                      "- Categorical covariates have a fixed colour palette",
                      paste0(
                        "- Continuous covariates are coloured in a ",
                        "Blue-Yellow-Red colour scheme, which can be ",
                        "changed in the plot controls"
                      )
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("sc1_civge_tog1", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.sc1_civge_tog1 == true",
                  radioButtons("sc1_civge_col1", "Colour (Continuous data):",
                    choices = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple"),
                    selected = "Blue-Yellow-Red"
                  ),
                  radioButtons("sc1_civge_ord1", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Original", inline = TRUE
                  ),
                  checkboxInput("sc1_civge_lab1", "Show cell info labels", value = TRUE)
                )
              )
            )
          ),
          # row 3 col 1 row 2
          fluidRow(
            class = "tab-section",
            column(
              12,
              uiOutput("sc1_civge_oup1.ui")
            )
          ),
          # row 3 col 1 row 3
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("sc1_civge_oup1.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("sc1_civge_oup1.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("sc1_civge_oup1.png", "Download PNG", class = "btn-sm"),
                checkboxInput("sc1_civge_tog9", "Show cell numbers / statistics")
              )
            )
          ),
          conditionalPanel(
            condition = "input.sc1_civge_tog9 == true",
            h4("Cell numbers / statistics"),
            radioButtons("sc1_civge_splt", "Split continuous cell info into:",
              choices = c("Quartile", "Decile"),
              selected = "Decile", inline = TRUE
            ),
            dataTableOutput("sc1_civge_.dt")
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
                selectInput("sc1_civge_inp2", "Gene name:", choices = NULL) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Gene expression to colour cells by",
                    content = c(
                      "Select gene to colour cells by gene expression",
                      paste0(
                        "- Gene expression are coloured in a ",
                        "White-Red colour scheme which can be ",
                        "changed in the plot controls"
                      )
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("sc1_civge_tog2", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.sc1_civge_tog2 == true",
                  radioButtons("sc1_civge_col2", "Colour:",
                    choices = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple"),
                    selected = "White-Red"
                  ),
                  radioButtons("sc1_civge_ord2", "Plot order:",
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
              uiOutput("sc1_civge_oup2.ui")
            )
          ),
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("sc1_civge_oup2.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("sc1_civge_oup2.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("sc1_civge_oup2.png", "Download PNG", class = "btn-sm")
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
            selectInput("sc1_civci_drX", "X-axis:",
              choices = sc1conf[dimred == TRUE]$UI,
              selected = sc1def$dimred[1]
            ),
            selectInput("sc1_civci_drY", "Y-axis:",
              choices = sc1conf[dimred == TRUE]$UI,
              selected = sc1def$dimred[2]
            )
          )
        ), # row 2 col 2
        # row 2 col 2
        column(4,
          div(
            class = "input-panel",
            checkboxInput("sc1_civci_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.sc1_civci_togL == true",
              selectInput("sc1_civci_sub1", "Cell info to subset:",
                choices = sc1conf[grp == TRUE]$UI,
                selected = sc1def$grp1
              ),
              uiOutput("sc1_civci_sub1.ui"),
              actionButton("sc1_civci_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("sc1_civci_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            )
          )
        ), # row 2 col 2
        # row 2 col 3
        column(
          4,
          div(
            class = "input-panel",
            checkboxInput("sc1_civci_tog0", "Adjust graphics"),
            conditionalPanel(
              condition = "input.sc1_civci_tog0 == true",
              sliderInput("sc1_civci_siz", "Point size:",
                min = 0, max = 4, value = 1.25, step = 0.25
              ),
              radioButtons("sc1_civci_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("sc1_civci_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              radioButtons("sc1_civci_asp", "Aspect ratio:",
                choices = c("Square", "Fixed", "Free"),
                selected = "Square", inline = TRUE
              ),
              checkboxInput("sc1_civci_txt", "Show axis text", value = FALSE)
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
                selectInput("sc1_civci_inp1", "Cell info:",
                  choices = sc1conf$UI,
                  selected = sc1def$meta1
                ) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Cell info to colour cells by",
                    content = c(
                      "Select cell info to colour cells",
                      "- Categorical covariates have a fixed colour palette",
                      paste0(
                        "- Continuous covariates are coloured in a ",
                        "Blue-Yellow-Red colour scheme, which can be ",
                        "changed in the plot controls"
                      )
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("sc1_civci_tog1", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.sc1_civci_tog1 == true",
                  radioButtons("sc1_civci_col1", "Colour (Continuous data):",
                    choices = c(
                      "White-Red", "Blue-Yellow-Red",
                      "Yellow-Green-Purple"
                    ),
                    selected = "Blue-Yellow-Red"
                  ),
                  radioButtons("sc1_civci_ord1", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Original", inline = TRUE
                  ),
                  checkboxInput("sc1_civci_lab1", "Show cell info labels", value = TRUE)
                )
              )
            )
          ),
          # row 3 col 1 row 2
          fluidRow(column(12, uiOutput("sc1_civci_oup1.ui"))),
          # row 3 col 1 row 3
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("sc1_civci_oup1.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("sc1_civci_oup1.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("sc1_civci_oup1.png", "Download PNG", class = "btn-sm")
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
                selectInput("sc1_civci_inp2", "Cell info:",
                  choices = sc1conf$UI,
                  selected = sc1def$meta2
                ) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Cell info to colour cells by",
                    content = c(
                      "Select cell info to colour cells",
                      "- Categorical covariates have a fixed colour palette",
                      paste0(
                        "- Continuous covariates are coloured in a ",
                        "Blue-Yellow-Red colour scheme, which can be ",
                        "changed in the plot controls"
                      )
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("sc1_civci_tog2", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.sc1_civci_tog2 == true",
                  radioButtons("sc1_civci_col2", "Colour (Continuous data):",
                    choices = c(
                      "White-Red", "Blue-Yellow-Red",
                      "Yellow-Green-Purple"
                    ),
                    selected = "Blue-Yellow-Red"
                  ),
                  radioButtons("sc1_civci_ord2", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Original", inline = TRUE
                  ),
                  checkboxInput("sc1_civci_lab2", "Show cell info labels", value = TRUE)
                )
              )
            )
          ),
          fluidRow(column(12, uiOutput("sc1_civci_oup2.ui"))),
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("sc1_civci_oup2.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("sc1_civci_oup2.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("sc1_civci_oup2.png", "Download PNG", class = "btn-sm")
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
            selectInput("sc1_gevge_drX", "X-axis:",
              choices = sc1conf[dimred == TRUE]$UI,
              selected = sc1def$dimred[1]
            ),
            selectInput("sc1_gevge_drY", "Y-axis:",
              choices = sc1conf[dimred == TRUE]$UI,
              selected = sc1def$dimred[2]
            )
          )
        ), # row 1 col 1
        # row 2 col 2
        column(4,
          div(
            class = "input-panel",
            checkboxInput("sc1_gevge_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.sc1_gevge_togL == true",
              selectInput("sc1_gevge_sub1", "Cell info to subset:",
                choices = sc1conf[grp == TRUE]$UI,
                selected = sc1def$grp1
              ),
              uiOutput("sc1_gevge_sub1.ui"),
              actionButton("sc1_gevge_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("sc1_gevge_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            )
          )
        ), # End of column
        # row 2 col 3
        column(
          4,
          div(
            class = "input-panel",
            checkboxInput("sc1_gevge_tog0", "Adjust graphics"),
            conditionalPanel(
              condition = "input.sc1_gevge_tog0 == true",
              sliderInput("sc1_gevge_siz", "Point size:",
                min = 0, max = 4, value = 1.25, step = 0.25
              ),
              radioButtons("sc1_gevge_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("sc1_gevge_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              radioButtons("sc1_gevge_asp", "Aspect ratio:",
                choices = c("Square", "Fixed", "Free"),
                selected = "Square", inline = TRUE
              ),
              checkboxInput("sc1_gevge_txt", "Show axis text", value = FALSE)
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
                selectInput("sc1_gevge_inp1", "Gene name:", choices = NULL) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Gene expression to colour cells by",
                    content = c(
                      "Select gene to colour cells by gene expression",
                      paste0(
                        "- Gene expression are coloured in a ",
                        "White-Red colour scheme which can be ",
                        "changed in the plot controls"
                      )
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("sc1_gevge_tog1", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.sc1_gevge_tog1 == true",
                  radioButtons("sc1_gevge_col1", "Colour:",
                    choices = c(
                      "White-Red", "Blue-Yellow-Red",
                      "Yellow-Green-Purple"
                    ),
                    selected = "White-Red"
                  ),
                  radioButtons("sc1_gevge_ord1", "Plot order:",
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
            column(12, uiOutput("sc1_gevge_oup1.ui"))
          ),
          # row 3 col 1 row 3
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("sc1_gevge_oup1.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("sc1_gevge_oup1.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("sc1_gevge_oup1.png", "Download PNG", class = "btn-sm")
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
                selectInput("sc1_gevge_inp2", "Gene name:", choices = NULL) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Gene expression to colour cells by",
                    content = c(
                      "Select gene to colour cells by gene expression",
                      paste0(
                        "- Gene expression are coloured in a ",
                        "White-Red colour scheme which can be ",
                        "changed in the plot controls"
                      )
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("sc1_gevge_tog2", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.sc1_gevge_tog2 == true",
                  radioButtons("sc1_gevge_col2", "Colour:",
                    choices = c(
                      "White-Red", "Blue-Yellow-Red",
                      "Yellow-Green-Purple"
                    ),
                    selected = "White-Red"
                  ),
                  radioButtons("sc1_gevge_ord2", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Max", inline = TRUE
                  )
                )
              )
            )
          ),
          fluidRow(
            class = "tab-section",
            column(12, uiOutput("sc1_gevge_oup2.ui"))
          ),
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                  numericInput("sc1_gevge_oup2.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                  downloadButton("sc1_gevge_oup2.pdf", "Download PDF", class = "btn-sm"),
                  downloadButton("sc1_gevge_oup2.png", "Download PNG", class = "btn-sm")
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
            textAreaInput("sc1_gem_inp", "Gene names:",
                          height = "100px",
                          value = paste0(sc1def$genes, collapse = ", ")
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "List of genes to plot on bubbleplot / heatmap",
                content = c(
                  "Input genes to plot",
                  "- Maximum 16 genes (due to ploting space limitations)",
                  "- Genes should be separated by comma, semicolon or newline"
                )
              ),
            selectInput("sc1_gem_drX", "X-axis:",
                        choices = sc1conf[dimred == TRUE]$UI,
                        selected = sc1def$dimred[1]
            ),
            selectInput("sc1_gem_drY", "Y-axis:",
                        choices = sc1conf[dimred == TRUE]$UI,
                        selected = sc1def$dimred[2]
            ),
            checkboxInput("sc1_gem_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.sc1_gem_togL == true",
              selectInput("sc1_gem_sub1", "Cell info to subset:",
                          choices = sc1conf[grp == TRUE]$UI,
                          selected = sc1def$grp1
              ),
              uiOutput("sc1_gem_sub1.ui"),
              actionButton("sc1_gem_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("sc1_gem_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            ),
            checkboxInput("sc1_gem_tog0", "Adjust graphics"),
            conditionalPanel(
              condition = "input.sc1_gem_tog0 == true",
              sliderInput("sc1_gem_siz", "Point size:",
                          min = 0, max = 3, value = 0.5, step = 0.1
              ),
              radioButtons("sc1_gem_psz", "Plot size:",
                           choices = c("Small", "Medium", "Large"),
                           selected = "Medium", inline = TRUE
              ),
              radioButtons("sc1_gem_fsz", "Font size:",
                           choices = c("Smaller", "Small", "Medium", "Large"),
                           selected = "Small", inline = TRUE
              ),
              radioButtons("sc1_gem_asp", "Aspect ratio:",
                           choices = c("Square", "Fixed", "Free"),
                           selected = "Square", inline = TRUE
              ),
              checkboxInput("sc1_gem_txt", "Show axis text", value = FALSE),
              radioButtons("sc1_gem_col", "Colour (Continuous data):",
                           choices = c(
                             "White-Red", "Blue-Yellow-Red",
                             "Yellow-Green-Purple"
                           ),
                           selected = "Blue-Yellow-Red"
              ),
              radioButtons("sc1_gem_ord", "Plot order:",
                           choices = c("Max", "Min", "Original", "Random"),
                           selected = "Max", inline = TRUE
              ),
              numericInput("sc1_gem_ncol", "Number of columns", value = 0, min = 0, step = 1)
            )
          ),
          div(
            class = "input-panel",
            fluidRow(
            column(4,
              numericInput("sc1_gem_oup1.height", "Height:", min = 1, max = 50, value = 25, step = 2)
            ),
            column(4,
              numericInput("sc1_gem_oup1.width", "Width:", min = 1, max = 50, value = 25, step = 2)
            ),
            column(4,
              numericInput("sc1_gem_oup1.res", "Res:", min = 72, max = 600, value = 150, step = 5)
            )
            ),
            downloadButton("sc1_gem_oup1.pdf", "Download PDF", class = "btn-sm"),
            downloadButton("sc1_gem_oup1.png", "Download PNG", class = "btn-sm")
          )
        )
               )
      ),
      column(8,
             uiOutput("sc1_gem_oup1.ui")
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
                   selectInput("sc1_gec_drX", "X-axis:",
                               choices = sc1conf[dimred == TRUE]$UI,
                               selected = sc1def$dimred[1]
                   ),
                   selectInput("sc1_gec_drY", "Y-axis:",
                               choices = sc1conf[dimred == TRUE]$UI,
                               selected = sc1def$dimred[2]
                   ),
                   selectInput("sc1_gec_inp1", "Gene 1:", choices = NULL) %>%
                     helper(
                       type = "inline", size = "m", fade = TRUE,
                       title = "Gene expression to colour cells by",
                       content = c(
                         "Select gene to colour cells by gene expression",
                         paste0(
                           "- Gene expression are coloured in a ",
                           "White-Red colour scheme which can be ",
                           "changed in the plot controls"
                         )
                       )
                     ),
                   selectInput("sc1_gec_inp2", "Gene 2:", choices = NULL) %>%
                     helper(
                       type = "inline", size = "m", fade = TRUE,
                       title = "Gene expression to colour cells by",
                       content = c(
                         "Select gene to colour cells by gene expression",
                         paste0(
                           "- Gene expression are coloured in a ",
                           "White-Blue colour scheme which can be ",
                           "changed in the plot controls"
                         )
                       )
                     ),
                   checkboxInput("sc1_gec_togL", "Subset cells"),
                   conditionalPanel(
                    condition = "input.sc1_gec_togL == true",
                    selectInput("sc1_gec_sub1", "Cell info to subset:", choices = sc1conf[grp == TRUE]$UI, selected = sc1def$grp1),
                     uiOutput("sc1_gec_sub1.ui"),
                     actionButton("sc1_gec_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
                     actionButton("sc1_gec_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
                   ),
                   checkboxInput("sc1_gec_tog0", "Adjust graphics"),
                   conditionalPanel(
                     condition = "input.sc1_gec_tog0 == true",
                     radioButtons("sc1_gec_col1", "Colour:",
                                  choices = c(
                                    "Red (Gene1); Blue (Gene2)",
                                    "Orange (Gene1); Blue (Gene2)",
                                    "Red (Gene1); Green (Gene2)",
                                    "Green (Gene1); Blue (Gene2)"
                                  ),
                                  selected = "Red (Gene1); Blue (Gene2)"
                     ),
                     radioButtons("sc1_gec_ord1", "Plot order:",
                                  choices = c("Max", "Min", "Original", "Random"),
                                  selected = "Max", inline = TRUE
                     ),
                     sliderInput("sc1_gec_siz", "Point size:",
                                 min = 0, max = 4, value = 1.25, step = 0.25
                     ),
                     radioButtons("sc1_gec_psz", "Plot size:",
                                  choices = c("Small", "Medium", "Large"),
                                  selected = "Medium", inline = TRUE
                     ),
                     radioButtons("sc1_gec_fsz", "Font size:",
                                  choices = c("Small", "Medium", "Large"),
                                  selected = "Small", inline = TRUE
                     ),
                     radioButtons("sc1_gec_asp", "Aspect ratio:",
                                  choices = c("Square", "Fixed", "Free"),
                                  selected = "Square", inline = TRUE
                     ),
                     checkboxInput("sc1_gec_txt", "Show axis text", value = FALSE)
                   )
                 ),
                 div(class="input-panel input-panel-section",
                     numericInput("sc1_gec_oup1.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                     downloadButton("sc1_gec_oup1.pdf", "Download PDF", class = "btn-sm"),
                     downloadButton("sc1_gec_oup1.png", "Download PNG", class = "btn-sm")
                 ),
                 div(class="input-panel-section",
                     h4("Cell numbers"),
                     dataTableOutput("sc1_gec_.dt")
                 )
               )
        ), # row 2 col 1
        # row 2 col 2
        column(
          8,
          uiOutput("sc1_gec_oup1.ui"),
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
            selectInput("sc1_vio_inp1", "Cell info (X-axis):",
              choices = sc1conf[grp == TRUE]$UI,
              selected = sc1def$grp1
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell info to group cells by",
                content = c(
                  "Select categorical cell info to group cells by",
                  "- Single cells are grouped by this categorical covariate",
                  "- Plotted as the X-axis of the violin plot / box plot"
                )
              ),
            selectInput("sc1_vio_inp2", "Cell info / Gene (Y-axis):", choices = NULL) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell Info / Gene to plot",
                content = c(
                  "Select cell info / gene to plot on Y-axis",
                  "- Can be continuous cell info (e.g. nUMIs / scores)",
                  "- Can also be gene expression"
                )
              ),
            radioButtons("sc1_vio_typ", "Plot type:",
              choices = c("violin", "boxplot", "lineplot"),
              selected = "violin", inline = TRUE
            ),
            checkboxInput("sc1_vio_pts", "Show data points", value = FALSE),
            checkboxInput("sc1_vio_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.sc1_vio_togL == true",
              selectInput("sc1_vio_sub1", "Cell info to subset:",
                choices = sc1conf[grp == TRUE]$UI,
                selected = sc1def$grp1
              ),
              uiOutput("sc1_vio_sub1.ui"),
              actionButton("sc1_vio_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("sc1_vio_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            ),
            checkboxInput("sc1_vio_tog", "Adjust graphics"),
            conditionalPanel(
              condition = "input.sc1_vio_tog == true",
              sliderInput("sc1_vio_siz", "Data point size:",
                min = 0, max = 4, value = 1.25, step = 0.25
              ),
              radioButtons("sc1_vio_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              radioButtons("sc1_vio_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              conditionalPanel(
              condition = "input.sc1_vio_typ == 'lineplot'",
              sliderInput("sc1_vio_barsz", "Line size", min = 0.05, max = 0.5, step = 0.01, value = 0.3)
              )
            ),
            selectInput("sc1_vio_datatype", "Data type", choices = c("normalised", "raw"), selected = "normalised"),
          ),
          div(
            class = "input-panel",
            numericInput("sc1_vio_oup.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
            downloadButton("sc1_vio_oup.pdf", "Download PDF", class = "btn-sm"),
            downloadButton("sc1_vio_oup.png", "Download PNG", class = "btn-sm")
          )
        ), # row 2 col 1
        # row 2 col 2
        column(
          9, uiOutput("sc1_vio_oup.ui")
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
            selectInput("sc1_pro_inp1", "Cell info to plot (X-axis):",
              choices = sc1conf[grp == TRUE]$UI,
              selected = sc1def$grp2
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell info to plot cells by",
                content = c(
                  "Select categorical cell info to plot cells by",
                  "- Plotted as the X-axis of the proportion plot"
                )
              ),
            selectInput("sc1_pro_inp2", "Cell info to group / colour by:",
              choices = sc1conf[grp == TRUE]$UI,
              selected = sc1def$grp1
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell info to group / colour cells by",
                content = c(
                  "Select categorical cell info to group / colour cells by",
                  "- Proportion / cell numbers are shown in different colours"
                )
              ),
            radioButtons("sc1_pro_typ", "Plot value:",
              choices = c("Proportion", "CellNumbers"),
              selected = "Proportion", inline = TRUE
            ),
            checkboxInput("sc1_pro_flp", "Flip X/Y", value = FALSE),
            checkboxInput("sc1_pro_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.sc1_pro_togL == true",
              selectInput("sc1_pro_sub1", "Cell info to subset:",
                choices = sc1conf[grp == TRUE]$UI,
                selected = sc1def$grp1
              ),
              uiOutput("sc1_pro_sub1.ui"),
              actionButton("sc1_pro_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("sc1_pro_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            ),
            checkboxInput("sc1_pro_tog", "Adjust graphics"),
            conditionalPanel(
              condition = "input.sc1_pro_tog == true",
              radioButtons("sc1_pro_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("sc1_pro_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              )
            )
          ),
          div(
            class = "input-panel",
            numericInput("sc1_pro_oup.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
            downloadButton("sc1_pro_oup.pdf", "Download PDF", class = "btn-sm"),
            downloadButton("sc1_pro_oup.png", "Download PNG", class = "btn-sm")
          )
        ), # row 2 col 1
        # row 2 col 2
        column(
          9, uiOutput("sc1_pro_oup.ui")
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
          3,
          div(
            class = "input-panel",
            textAreaInput("sc1_hea_inp", "Gene names",
              height = "100px",
              value = paste0(sc1def$genes, collapse = ", ")
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "List of genes to plot on bubbleplot / heatmap",
                content = c(
                  "Input genes to plot",
                  "- Maximum 50 genes (due to ploting space limitations)",
                  "- Genes should be separated by comma, semicolon or newline"
                )
              ),
            selectInput("sc1_hea_grp", "Group by:",
              choices = sc1conf[grp == TRUE]$UI,
              selected = sc1conf[grp == TRUE]$UI[1]
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell info to group cells by",
                content = c(
                  "Select categorical cell info to group cells by",
                  "- Single cells are grouped by this categorical covariate",
                  "- Plotted as the X-axis of the bubbleplot / heatmap"
                )
              ),
            radioButtons("sc1_hea_plt", "Plot type:",
              choices = c("Bubbleplot", "Heatmap"),
              selected = "Bubbleplot", inline = TRUE
            ),
            checkboxInput("sc1_hea_scl", "Scale gene expression", value = TRUE),
            checkboxInput("sc1_hea_row", "Cluster rows (genes)", value = TRUE),
            checkboxInput("sc1_hea_col", "Cluster columns (samples)", value = FALSE),
            checkboxInput("sc1_hea_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.sc1_hea_togL == true",
              selectInput("sc1_hea_sub1", "Cell info to subset:",
                choices = sc1conf[grp == TRUE]$UI,
                selected = sc1def$grp1
              ),
              uiOutput("sc1_hea_sub1.ui"),
              actionButton("sc1_hea_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("sc1_hea_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            ),
            checkboxInput("sc1_hea_tog", "Adjust graphics"),
            conditionalPanel(
              condition = "input.sc1_hea_tog == true",
              radioButtons("sc1_hea_cols", "Colour scheme:",
                choices = c(
                  "White-Red", "Blue-Yellow-Red",
                  "Yellow-Green-Purple"
                ),
                selected = "Blue-Yellow-Red"
              ),
              radioButtons("sc1_hea_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("sc1_hea_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              )
            )
          ),
          div(
            class = "input-panel",
            numericInput("sc1_hea_oup.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
            downloadButton("sc1_hea_oup.pdf", "Download PDF", class = "btn-sm"),
            downloadButton("sc1_hea_oup.png", "Download PNG", class = "btn-sm")
          )
        ), # row 2 col 1
        # row 2 col 2
        column(
          9, h4(htmlOutput("sc1_hea_oupTxt")),
          uiOutput("sc1_hea_oup.ui")
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
                     selectInput("sc1_mar_cls","Select clustering:", choices = names(sc1mar),selected = 1)
                   )
                 )
               )
        )
      ), # end of row 2
      # row 3 ----
      fluidRow(
      column(12,
        DTOutput("sc1_mar_table")
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
            selectInput("sc2_civge_drX", "X-axis:",
              choices = sc2conf[dimred == TRUE]$UI,
              selected = sc2def$dimred[1]
            ),
            selectInput("sc2_civge_drY", "Y-axis:",
              choices = sc2conf[dimred == TRUE]$UI,
              selected = sc2def$dimred[2]
            )
          )
        ), # row 1 col 1
        # row 2 col 2
        column(
          4,
          div(
            class = "input-panel",
            checkboxInput("sc2_civge_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.sc2_civge_togL == true",
              selectInput("sc2_civge_sub1", "Cell info to subset:",
                choices = sc2conf[grp == TRUE]$UI,
                selected = sc2def$grp1
              ),
              uiOutput("sc2_civge_sub1.ui"),
              actionButton("sc2_civge_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("sc2_civge_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            )
          )
        ), # End of column
        # row 2 col 3
        column(
          4,
          div(
            class = "input-panel",
            checkboxInput("sc2_civge_tog0", "Adjust graphics"),
            conditionalPanel(
              condition = "input.sc2_civge_tog0 == true",
              sliderInput("sc2_civge_siz", "Point size:",
                min = 0, max = 4, value = 1.25, step = 0.25
              ),
              radioButtons("sc2_civge_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("sc2_civge_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              radioButtons("sc2_civge_asp", "Aspect ratio:",
                choices = c("Square", "Fixed", "Free"),
                selected = "Square", inline = TRUE
              ),
              checkboxInput("sc2_civge_txt", "Show axis text", value = FALSE)
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
                selectInput("sc2_civge_inp1", "Cell info:",
                  choices = sc2conf$UI,
                  selected = sc2def$meta1
                ) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Cell info to colour cells by",
                    content = c(
                      "Select cell info to colour cells",
                      "- Categorical covariates have a fixed colour palette",
                      paste0(
                        "- Continuous covariates are coloured in a ",
                        "Blue-Yellow-Red colour scheme, which can be ",
                        "changed in the plot controls"
                      )
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("sc2_civge_tog1", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.sc2_civge_tog1 == true",
                  radioButtons("sc2_civge_col1", "Colour (Continuous data):",
                    choices = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple"),
                    selected = "Blue-Yellow-Red"
                  ),
                  radioButtons("sc2_civge_ord1", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Original", inline = TRUE
                  ),
                  checkboxInput("sc2_civge_lab1", "Show cell info labels", value = TRUE)
                )
              )
            )
          ),
          # row 3 col 1 row 2
          fluidRow(
            class = "tab-section",
            column(
              12,
              uiOutput("sc2_civge_oup1.ui")
            )
          ),
          # row 3 col 1 row 3
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("sc2_civge_oup1.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("sc2_civge_oup1.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("sc2_civge_oup1.png", "Download PNG", class = "btn-sm"),
                checkboxInput("sc2_civge_tog9", "Show cell numbers / statistics")
              )
            )
          ),
          conditionalPanel(
            condition = "input.sc2_civge_tog9 == true",
            h4("Cell numbers / statistics"),
            radioButtons("sc2_civge_splt", "Split continuous cell info into:",
              choices = c("Quartile", "Decile"),
              selected = "Decile", inline = TRUE
            ),
            dataTableOutput("sc2_civge_.dt")
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
                selectInput("sc2_civge_inp2", "Gene name:", choices = NULL) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Gene expression to colour cells by",
                    content = c(
                      "Select gene to colour cells by gene expression",
                      paste0(
                        "- Gene expression are coloured in a ",
                        "White-Red colour scheme which can be ",
                        "changed in the plot controls"
                      )
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("sc2_civge_tog2", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.sc2_civge_tog2 == true",
                  radioButtons("sc2_civge_col2", "Colour:",
                    choices = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple"),
                    selected = "White-Red"
                  ),
                  radioButtons("sc2_civge_ord2", "Plot order:",
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
              uiOutput("sc2_civge_oup2.ui")
            )
          ),
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("sc2_civge_oup2.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("sc2_civge_oup2.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("sc2_civge_oup2.png", "Download PNG", class = "btn-sm")
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
            selectInput("sc2_civci_drX", "X-axis:",
              choices = sc2conf[dimred == TRUE]$UI,
              selected = sc2def$dimred[1]
            ),
            selectInput("sc2_civci_drY", "Y-axis:",
              choices = sc2conf[dimred == TRUE]$UI,
              selected = sc2def$dimred[2]
            )
          )
        ), # row 2 col 2
        # row 2 col 2
        column(4,
          div(
            class = "input-panel",
            checkboxInput("sc2_civci_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.sc2_civci_togL == true",
              selectInput("sc2_civci_sub1", "Cell info to subset:",
                choices = sc2conf[grp == TRUE]$UI,
                selected = sc2def$grp1
              ),
              uiOutput("sc2_civci_sub1.ui"),
              actionButton("sc2_civci_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("sc2_civci_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            )
          )
        ), # row 2 col 2
        # row 2 col 3
        column(
          4,
          div(
            class = "input-panel",
            checkboxInput("sc2_civci_tog0", "Adjust graphics"),
            conditionalPanel(
              condition = "input.sc2_civci_tog0 == true",
              sliderInput("sc2_civci_siz", "Point size:",
                min = 0, max = 4, value = 1.25, step = 0.25
              ),
              radioButtons("sc2_civci_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("sc2_civci_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              radioButtons("sc2_civci_asp", "Aspect ratio:",
                choices = c("Square", "Fixed", "Free"),
                selected = "Square", inline = TRUE
              ),
              checkboxInput("sc2_civci_txt", "Show axis text", value = FALSE)
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
                selectInput("sc2_civci_inp1", "Cell info:",
                  choices = sc2conf$UI,
                  selected = sc2def$meta1
                ) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Cell info to colour cells by",
                    content = c(
                      "Select cell info to colour cells",
                      "- Categorical covariates have a fixed colour palette",
                      paste0(
                        "- Continuous covariates are coloured in a ",
                        "Blue-Yellow-Red colour scheme, which can be ",
                        "changed in the plot controls"
                      )
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("sc2_civci_tog1", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.sc2_civci_tog1 == true",
                  radioButtons("sc2_civci_col1", "Colour (Continuous data):",
                    choices = c(
                      "White-Red", "Blue-Yellow-Red",
                      "Yellow-Green-Purple"
                    ),
                    selected = "Blue-Yellow-Red"
                  ),
                  radioButtons("sc2_civci_ord1", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Original", inline = TRUE
                  ),
                  checkboxInput("sc2_civci_lab1", "Show cell info labels", value = TRUE)
                )
              )
            )
          ),
          # row 3 col 1 row 2
          fluidRow(column(12, uiOutput("sc2_civci_oup1.ui"))),
          # row 3 col 1 row 3
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("sc2_civci_oup1.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("sc2_civci_oup1.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("sc2_civci_oup1.png", "Download PNG", class = "btn-sm")
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
                selectInput("sc2_civci_inp2", "Cell info:",
                  choices = sc2conf$UI,
                  selected = sc2def$meta2
                ) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Cell info to colour cells by",
                    content = c(
                      "Select cell info to colour cells",
                      "- Categorical covariates have a fixed colour palette",
                      paste0(
                        "- Continuous covariates are coloured in a ",
                        "Blue-Yellow-Red colour scheme, which can be ",
                        "changed in the plot controls"
                      )
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("sc2_civci_tog2", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.sc2_civci_tog2 == true",
                  radioButtons("sc2_civci_col2", "Colour (Continuous data):",
                    choices = c(
                      "White-Red", "Blue-Yellow-Red",
                      "Yellow-Green-Purple"
                    ),
                    selected = "Blue-Yellow-Red"
                  ),
                  radioButtons("sc2_civci_ord2", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Original", inline = TRUE
                  ),
                  checkboxInput("sc2_civci_lab2", "Show cell info labels", value = TRUE)
                )
              )
            )
          ),
          fluidRow(column(12, uiOutput("sc2_civci_oup2.ui"))),
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("sc2_civci_oup2.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("sc2_civci_oup2.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("sc2_civci_oup2.png", "Download PNG", class = "btn-sm")
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
            selectInput("sc2_gevge_drX", "X-axis:",
              choices = sc2conf[dimred == TRUE]$UI,
              selected = sc2def$dimred[1]
            ),
            selectInput("sc2_gevge_drY", "Y-axis:",
              choices = sc2conf[dimred == TRUE]$UI,
              selected = sc2def$dimred[2]
            )
          )
        ), # row 1 col 1
        # row 2 col 2
        column(4,
          div(
            class = "input-panel",
            checkboxInput("sc2_gevge_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.sc2_gevge_togL == true",
              selectInput("sc2_gevge_sub1", "Cell info to subset:",
                choices = sc2conf[grp == TRUE]$UI,
                selected = sc2def$grp1
              ),
              uiOutput("sc2_gevge_sub1.ui"),
              actionButton("sc2_gevge_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("sc2_gevge_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            )
          )
        ), # End of column
        # row 2 col 3
        column(
          4,
          div(
            class = "input-panel",
            checkboxInput("sc2_gevge_tog0", "Adjust graphics"),
            conditionalPanel(
              condition = "input.sc2_gevge_tog0 == true",
              sliderInput("sc2_gevge_siz", "Point size:",
                min = 0, max = 4, value = 1.25, step = 0.25
              ),
              radioButtons("sc2_gevge_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("sc2_gevge_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              radioButtons("sc2_gevge_asp", "Aspect ratio:",
                choices = c("Square", "Fixed", "Free"),
                selected = "Square", inline = TRUE
              ),
              checkboxInput("sc2_gevge_txt", "Show axis text", value = FALSE)
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
                selectInput("sc2_gevge_inp1", "Gene name:", choices = NULL) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Gene expression to colour cells by",
                    content = c(
                      "Select gene to colour cells by gene expression",
                      paste0(
                        "- Gene expression are coloured in a ",
                        "White-Red colour scheme which can be ",
                        "changed in the plot controls"
                      )
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("sc2_gevge_tog1", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.sc2_gevge_tog1 == true",
                  radioButtons("sc2_gevge_col1", "Colour:",
                    choices = c(
                      "White-Red", "Blue-Yellow-Red",
                      "Yellow-Green-Purple"
                    ),
                    selected = "White-Red"
                  ),
                  radioButtons("sc2_gevge_ord1", "Plot order:",
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
            column(12, uiOutput("sc2_gevge_oup1.ui"))
          ),
          # row 3 col 1 row 3
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("sc2_gevge_oup1.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("sc2_gevge_oup1.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("sc2_gevge_oup1.png", "Download PNG", class = "btn-sm")
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
                selectInput("sc2_gevge_inp2", "Gene name:", choices = NULL) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Gene expression to colour cells by",
                    content = c(
                      "Select gene to colour cells by gene expression",
                      paste0(
                        "- Gene expression are coloured in a ",
                        "White-Red colour scheme which can be ",
                        "changed in the plot controls"
                      )
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("sc2_gevge_tog2", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.sc2_gevge_tog2 == true",
                  radioButtons("sc2_gevge_col2", "Colour:",
                    choices = c(
                      "White-Red", "Blue-Yellow-Red",
                      "Yellow-Green-Purple"
                    ),
                    selected = "White-Red"
                  ),
                  radioButtons("sc2_gevge_ord2", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Max", inline = TRUE
                  )
                )
              )
            )
          ),
          fluidRow(
            class = "tab-section",
            column(12, uiOutput("sc2_gevge_oup2.ui"))
          ),
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                  numericInput("sc2_gevge_oup2.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                  downloadButton("sc2_gevge_oup2.pdf", "Download PDF", class = "btn-sm"),
                  downloadButton("sc2_gevge_oup2.png", "Download PNG", class = "btn-sm")
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
            textAreaInput("sc2_gem_inp", "Gene names:",
                          height = "100px",
                          value = paste0(sc2def$genes, collapse = ", ")
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "List of genes to plot on bubbleplot / heatmap",
                content = c(
                  "Input genes to plot",
                  "- Maximum 16 genes (due to ploting space limitations)",
                  "- Genes should be separated by comma, semicolon or newline"
                )
              ),
            selectInput("sc2_gem_drX", "X-axis:",
                        choices = sc2conf[dimred == TRUE]$UI,
                        selected = sc2def$dimred[1]
            ),
            selectInput("sc2_gem_drY", "Y-axis:",
                        choices = sc2conf[dimred == TRUE]$UI,
                        selected = sc2def$dimred[2]
            ),
            checkboxInput("sc2_gem_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.sc2_gem_togL == true",
              selectInput("sc2_gem_sub1", "Cell info to subset:",
                          choices = sc2conf[grp == TRUE]$UI,
                          selected = sc2def$grp1
              ),
              uiOutput("sc2_gem_sub1.ui"),
              actionButton("sc2_gem_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("sc2_gem_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            ),
            checkboxInput("sc2_gem_tog0", "Adjust graphics"),
            conditionalPanel(
              condition = "input.sc2_gem_tog0 == true",
              sliderInput("sc2_gem_siz", "Point size:",
                          min = 0, max = 3, value = 0.5, step = 0.1
              ),
              radioButtons("sc2_gem_psz", "Plot size:",
                           choices = c("Small", "Medium", "Large"),
                           selected = "Medium", inline = TRUE
              ),
              radioButtons("sc2_gem_fsz", "Font size:",
                           choices = c("Smaller", "Small", "Medium", "Large"),
                           selected = "Small", inline = TRUE
              ),
              radioButtons("sc2_gem_asp", "Aspect ratio:",
                           choices = c("Square", "Fixed", "Free"),
                           selected = "Square", inline = TRUE
              ),
              checkboxInput("sc2_gem_txt", "Show axis text", value = FALSE),
              radioButtons("sc2_gem_col", "Colour (Continuous data):",
                           choices = c(
                             "White-Red", "Blue-Yellow-Red",
                             "Yellow-Green-Purple"
                           ),
                           selected = "Blue-Yellow-Red"
              ),
              radioButtons("sc2_gem_ord", "Plot order:",
                           choices = c("Max", "Min", "Original", "Random"),
                           selected = "Max", inline = TRUE
              ),
              numericInput("sc2_gem_ncol", "Number of columns", value = 0, min = 0, step = 1)
            )
          ),
          div(
            class = "input-panel",
            fluidRow(
            column(4,
              numericInput("sc2_gem_oup1.height", "Height:", min = 1, max = 50, value = 25, step = 2)
            ),
            column(4,
              numericInput("sc2_gem_oup1.width", "Width:", min = 1, max = 50, value = 25, step = 2)
            ),
            column(4,
              numericInput("sc2_gem_oup1.res", "Res:", min = 72, max = 600, value = 150, step = 5)
            )
            ),
            downloadButton("sc2_gem_oup1.pdf", "Download PDF", class = "btn-sm"),
            downloadButton("sc2_gem_oup1.png", "Download PNG", class = "btn-sm")
          )
        )
               )
      ),
      column(8,
             uiOutput("sc2_gem_oup1.ui")
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
                   selectInput("sc2_gec_drX", "X-axis:",
                               choices = sc2conf[dimred == TRUE]$UI,
                               selected = sc2def$dimred[1]
                   ),
                   selectInput("sc2_gec_drY", "Y-axis:",
                               choices = sc2conf[dimred == TRUE]$UI,
                               selected = sc2def$dimred[2]
                   ),
                   selectInput("sc2_gec_inp1", "Gene 1:", choices = NULL) %>%
                     helper(
                       type = "inline", size = "m", fade = TRUE,
                       title = "Gene expression to colour cells by",
                       content = c(
                         "Select gene to colour cells by gene expression",
                         paste0(
                           "- Gene expression are coloured in a ",
                           "White-Red colour scheme which can be ",
                           "changed in the plot controls"
                         )
                       )
                     ),
                   selectInput("sc2_gec_inp2", "Gene 2:", choices = NULL) %>%
                     helper(
                       type = "inline", size = "m", fade = TRUE,
                       title = "Gene expression to colour cells by",
                       content = c(
                         "Select gene to colour cells by gene expression",
                         paste0(
                           "- Gene expression are coloured in a ",
                           "White-Blue colour scheme which can be ",
                           "changed in the plot controls"
                         )
                       )
                     ),
                   checkboxInput("sc2_gec_togL", "Subset cells"),
                   conditionalPanel(
                    condition = "input.sc2_gec_togL == true",
                    selectInput("sc2_gec_sub1", "Cell info to subset:", choices = sc2conf[grp == TRUE]$UI, selected = sc2def$grp1),
                     uiOutput("sc2_gec_sub1.ui"),
                     actionButton("sc2_gec_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
                     actionButton("sc2_gec_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
                   ),
                   checkboxInput("sc2_gec_tog0", "Adjust graphics"),
                   conditionalPanel(
                     condition = "input.sc2_gec_tog0 == true",
                     radioButtons("sc2_gec_col1", "Colour:",
                                  choices = c(
                                    "Red (Gene1); Blue (Gene2)",
                                    "Orange (Gene1); Blue (Gene2)",
                                    "Red (Gene1); Green (Gene2)",
                                    "Green (Gene1); Blue (Gene2)"
                                  ),
                                  selected = "Red (Gene1); Blue (Gene2)"
                     ),
                     radioButtons("sc2_gec_ord1", "Plot order:",
                                  choices = c("Max", "Min", "Original", "Random"),
                                  selected = "Max", inline = TRUE
                     ),
                     sliderInput("sc2_gec_siz", "Point size:",
                                 min = 0, max = 4, value = 1.25, step = 0.25
                     ),
                     radioButtons("sc2_gec_psz", "Plot size:",
                                  choices = c("Small", "Medium", "Large"),
                                  selected = "Medium", inline = TRUE
                     ),
                     radioButtons("sc2_gec_fsz", "Font size:",
                                  choices = c("Small", "Medium", "Large"),
                                  selected = "Small", inline = TRUE
                     ),
                     radioButtons("sc2_gec_asp", "Aspect ratio:",
                                  choices = c("Square", "Fixed", "Free"),
                                  selected = "Square", inline = TRUE
                     ),
                     checkboxInput("sc2_gec_txt", "Show axis text", value = FALSE)
                   )
                 ),
                 div(class="input-panel input-panel-section",
                     numericInput("sc2_gec_oup1.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                     downloadButton("sc2_gec_oup1.pdf", "Download PDF", class = "btn-sm"),
                     downloadButton("sc2_gec_oup1.png", "Download PNG", class = "btn-sm")
                 ),
                 div(class="input-panel-section",
                     h4("Cell numbers"),
                     dataTableOutput("sc2_gec_.dt")
                 )
               )
        ), # row 2 col 1
        # row 2 col 2
        column(
          8,
          uiOutput("sc2_gec_oup1.ui"),
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
            selectInput("sc2_vio_inp1", "Cell info (X-axis):",
              choices = sc2conf[grp == TRUE]$UI,
              selected = sc2def$grp1
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell info to group cells by",
                content = c(
                  "Select categorical cell info to group cells by",
                  "- Single cells are grouped by this categorical covariate",
                  "- Plotted as the X-axis of the violin plot / box plot"
                )
              ),
            selectInput("sc2_vio_inp2", "Cell info / Gene (Y-axis):", choices = NULL) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell Info / Gene to plot",
                content = c(
                  "Select cell info / gene to plot on Y-axis",
                  "- Can be continuous cell info (e.g. nUMIs / scores)",
                  "- Can also be gene expression"
                )
              ),
            radioButtons("sc2_vio_typ", "Plot type:",
              choices = c("violin", "boxplot", "lineplot"),
              selected = "violin", inline = TRUE
            ),
            checkboxInput("sc2_vio_pts", "Show data points", value = FALSE),
            checkboxInput("sc2_vio_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.sc2_vio_togL == true",
              selectInput("sc2_vio_sub1", "Cell info to subset:",
                choices = sc2conf[grp == TRUE]$UI,
                selected = sc2def$grp1
              ),
              uiOutput("sc2_vio_sub1.ui"),
              actionButton("sc2_vio_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("sc2_vio_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            ),
            checkboxInput("sc2_vio_tog", "Adjust graphics"),
            conditionalPanel(
              condition = "input.sc2_vio_tog == true",
              sliderInput("sc2_vio_siz", "Data point size:",
                min = 0, max = 4, value = 1.25, step = 0.25
              ),
              radioButtons("sc2_vio_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              radioButtons("sc2_vio_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              conditionalPanel(
              condition = "input.sc2_vio_typ == 'lineplot'",
              sliderInput("sc2_vio_barsz", "Line size", min = 0.05, max = 0.5, step = 0.01, value = 0.3)
              )
            ),
            selectInput("sc2_vio_datatype", "Data type", choices = c("normalised", "raw"), selected = "normalised"),
          ),
          div(
            class = "input-panel",
            numericInput("sc2_vio_oup.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
            downloadButton("sc2_vio_oup.pdf", "Download PDF", class = "btn-sm"),
            downloadButton("sc2_vio_oup.png", "Download PNG", class = "btn-sm")
          )
        ), # row 2 col 1
        # row 2 col 2
        column(
          9, uiOutput("sc2_vio_oup.ui")
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
            selectInput("sc2_pro_inp1", "Cell info to plot (X-axis):",
              choices = sc2conf[grp == TRUE]$UI,
              selected = sc2def$grp2
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell info to plot cells by",
                content = c(
                  "Select categorical cell info to plot cells by",
                  "- Plotted as the X-axis of the proportion plot"
                )
              ),
            selectInput("sc2_pro_inp2", "Cell info to group / colour by:",
              choices = sc2conf[grp == TRUE]$UI,
              selected = sc2def$grp1
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell info to group / colour cells by",
                content = c(
                  "Select categorical cell info to group / colour cells by",
                  "- Proportion / cell numbers are shown in different colours"
                )
              ),
            radioButtons("sc2_pro_typ", "Plot value:",
              choices = c("Proportion", "CellNumbers"),
              selected = "Proportion", inline = TRUE
            ),
            checkboxInput("sc2_pro_flp", "Flip X/Y", value = FALSE),
            checkboxInput("sc2_pro_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.sc2_pro_togL == true",
              selectInput("sc2_pro_sub1", "Cell info to subset:",
                choices = sc2conf[grp == TRUE]$UI,
                selected = sc2def$grp1
              ),
              uiOutput("sc2_pro_sub1.ui"),
              actionButton("sc2_pro_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("sc2_pro_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            ),
            checkboxInput("sc2_pro_tog", "Adjust graphics"),
            conditionalPanel(
              condition = "input.sc2_pro_tog == true",
              radioButtons("sc2_pro_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("sc2_pro_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              )
            )
          ),
          div(
            class = "input-panel",
            numericInput("sc2_pro_oup.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
            downloadButton("sc2_pro_oup.pdf", "Download PDF", class = "btn-sm"),
            downloadButton("sc2_pro_oup.png", "Download PNG", class = "btn-sm")
          )
        ), # row 2 col 1
        # row 2 col 2
        column(
          9, uiOutput("sc2_pro_oup.ui")
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
          3,
          div(
            class = "input-panel",
            textAreaInput("sc2_hea_inp", "Gene names",
              height = "100px",
              value = paste0(sc2def$genes, collapse = ", ")
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "List of genes to plot on bubbleplot / heatmap",
                content = c(
                  "Input genes to plot",
                  "- Maximum 50 genes (due to ploting space limitations)",
                  "- Genes should be separated by comma, semicolon or newline"
                )
              ),
            selectInput("sc2_hea_grp", "Group by:",
              choices = sc2conf[grp == TRUE]$UI,
              selected = sc2conf[grp == TRUE]$UI[1]
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell info to group cells by",
                content = c(
                  "Select categorical cell info to group cells by",
                  "- Single cells are grouped by this categorical covariate",
                  "- Plotted as the X-axis of the bubbleplot / heatmap"
                )
              ),
            radioButtons("sc2_hea_plt", "Plot type:",
              choices = c("Bubbleplot", "Heatmap"),
              selected = "Bubbleplot", inline = TRUE
            ),
            checkboxInput("sc2_hea_scl", "Scale gene expression", value = TRUE),
            checkboxInput("sc2_hea_row", "Cluster rows (genes)", value = TRUE),
            checkboxInput("sc2_hea_col", "Cluster columns (samples)", value = FALSE),
            checkboxInput("sc2_hea_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.sc2_hea_togL == true",
              selectInput("sc2_hea_sub1", "Cell info to subset:",
                choices = sc2conf[grp == TRUE]$UI,
                selected = sc2def$grp1
              ),
              uiOutput("sc2_hea_sub1.ui"),
              actionButton("sc2_hea_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("sc2_hea_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            ),
            checkboxInput("sc2_hea_tog", "Adjust graphics"),
            conditionalPanel(
              condition = "input.sc2_hea_tog == true",
              radioButtons("sc2_hea_cols", "Colour scheme:",
                choices = c(
                  "White-Red", "Blue-Yellow-Red",
                  "Yellow-Green-Purple"
                ),
                selected = "Blue-Yellow-Red"
              ),
              radioButtons("sc2_hea_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("sc2_hea_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              )
            )
          ),
          div(
            class = "input-panel",
            numericInput("sc2_hea_oup.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
            downloadButton("sc2_hea_oup.pdf", "Download PDF", class = "btn-sm"),
            downloadButton("sc2_hea_oup.png", "Download PNG", class = "btn-sm")
          )
        ), # row 2 col 1
        # row 2 col 2
        column(
          9, h4(htmlOutput("sc2_hea_oupTxt")),
          uiOutput("sc2_hea_oup.ui")
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
                     selectInput("sc2_mar_cls","Select clustering:", choices = names(sc2mar),selected = 1)
                   )
                 )
               )
        )
      ), # end of row 2
      # row 3 ----
      fluidRow(
      column(12,
        DTOutput("sc2_mar_table")
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
            selectInput("sc3_civge_drX", "X-axis:",
              choices = sc3conf[dimred == TRUE]$UI,
              selected = sc3def$dimred[1]
            ),
            selectInput("sc3_civge_drY", "Y-axis:",
              choices = sc3conf[dimred == TRUE]$UI,
              selected = sc3def$dimred[2]
            )
          )
        ), # row 1 col 1
        # row 2 col 2
        column(
          4,
          div(
            class = "input-panel",
            checkboxInput("sc3_civge_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.sc3_civge_togL == true",
              selectInput("sc3_civge_sub1", "Cell info to subset:",
                choices = sc3conf[grp == TRUE]$UI,
                selected = sc3def$grp1
              ),
              uiOutput("sc3_civge_sub1.ui"),
              actionButton("sc3_civge_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("sc3_civge_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            )
          )
        ), # End of column
        # row 2 col 3
        column(
          4,
          div(
            class = "input-panel",
            checkboxInput("sc3_civge_tog0", "Adjust graphics"),
            conditionalPanel(
              condition = "input.sc3_civge_tog0 == true",
              sliderInput("sc3_civge_siz", "Point size:",
                min = 0, max = 4, value = 1.25, step = 0.25
              ),
              radioButtons("sc3_civge_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("sc3_civge_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              radioButtons("sc3_civge_asp", "Aspect ratio:",
                choices = c("Square", "Fixed", "Free"),
                selected = "Square", inline = TRUE
              ),
              checkboxInput("sc3_civge_txt", "Show axis text", value = FALSE)
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
                selectInput("sc3_civge_inp1", "Cell info:",
                  choices = sc3conf$UI,
                  selected = sc3def$meta1
                ) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Cell info to colour cells by",
                    content = c(
                      "Select cell info to colour cells",
                      "- Categorical covariates have a fixed colour palette",
                      paste0(
                        "- Continuous covariates are coloured in a ",
                        "Blue-Yellow-Red colour scheme, which can be ",
                        "changed in the plot controls"
                      )
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("sc3_civge_tog1", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.sc3_civge_tog1 == true",
                  radioButtons("sc3_civge_col1", "Colour (Continuous data):",
                    choices = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple"),
                    selected = "Blue-Yellow-Red"
                  ),
                  radioButtons("sc3_civge_ord1", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Original", inline = TRUE
                  ),
                  checkboxInput("sc3_civge_lab1", "Show cell info labels", value = TRUE)
                )
              )
            )
          ),
          # row 3 col 1 row 2
          fluidRow(
            class = "tab-section",
            column(
              12,
              uiOutput("sc3_civge_oup1.ui")
            )
          ),
          # row 3 col 1 row 3
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("sc3_civge_oup1.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("sc3_civge_oup1.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("sc3_civge_oup1.png", "Download PNG", class = "btn-sm"),
                checkboxInput("sc3_civge_tog9", "Show cell numbers / statistics")
              )
            )
          ),
          conditionalPanel(
            condition = "input.sc3_civge_tog9 == true",
            h4("Cell numbers / statistics"),
            radioButtons("sc3_civge_splt", "Split continuous cell info into:",
              choices = c("Quartile", "Decile"),
              selected = "Decile", inline = TRUE
            ),
            dataTableOutput("sc3_civge_.dt")
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
                selectInput("sc3_civge_inp2", "Gene name:", choices = NULL) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Gene expression to colour cells by",
                    content = c(
                      "Select gene to colour cells by gene expression",
                      paste0(
                        "- Gene expression are coloured in a ",
                        "White-Red colour scheme which can be ",
                        "changed in the plot controls"
                      )
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("sc3_civge_tog2", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.sc3_civge_tog2 == true",
                  radioButtons("sc3_civge_col2", "Colour:",
                    choices = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple"),
                    selected = "White-Red"
                  ),
                  radioButtons("sc3_civge_ord2", "Plot order:",
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
              uiOutput("sc3_civge_oup2.ui")
            )
          ),
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("sc3_civge_oup2.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("sc3_civge_oup2.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("sc3_civge_oup2.png", "Download PNG", class = "btn-sm")
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
            selectInput("sc3_civci_drX", "X-axis:",
              choices = sc3conf[dimred == TRUE]$UI,
              selected = sc3def$dimred[1]
            ),
            selectInput("sc3_civci_drY", "Y-axis:",
              choices = sc3conf[dimred == TRUE]$UI,
              selected = sc3def$dimred[2]
            )
          )
        ), # row 2 col 2
        # row 2 col 2
        column(4,
          div(
            class = "input-panel",
            checkboxInput("sc3_civci_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.sc3_civci_togL == true",
              selectInput("sc3_civci_sub1", "Cell info to subset:",
                choices = sc3conf[grp == TRUE]$UI,
                selected = sc3def$grp1
              ),
              uiOutput("sc3_civci_sub1.ui"),
              actionButton("sc3_civci_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("sc3_civci_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            )
          )
        ), # row 2 col 2
        # row 2 col 3
        column(
          4,
          div(
            class = "input-panel",
            checkboxInput("sc3_civci_tog0", "Adjust graphics"),
            conditionalPanel(
              condition = "input.sc3_civci_tog0 == true",
              sliderInput("sc3_civci_siz", "Point size:",
                min = 0, max = 4, value = 1.25, step = 0.25
              ),
              radioButtons("sc3_civci_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("sc3_civci_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              radioButtons("sc3_civci_asp", "Aspect ratio:",
                choices = c("Square", "Fixed", "Free"),
                selected = "Square", inline = TRUE
              ),
              checkboxInput("sc3_civci_txt", "Show axis text", value = FALSE)
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
                selectInput("sc3_civci_inp1", "Cell info:",
                  choices = sc3conf$UI,
                  selected = sc3def$meta1
                ) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Cell info to colour cells by",
                    content = c(
                      "Select cell info to colour cells",
                      "- Categorical covariates have a fixed colour palette",
                      paste0(
                        "- Continuous covariates are coloured in a ",
                        "Blue-Yellow-Red colour scheme, which can be ",
                        "changed in the plot controls"
                      )
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("sc3_civci_tog1", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.sc3_civci_tog1 == true",
                  radioButtons("sc3_civci_col1", "Colour (Continuous data):",
                    choices = c(
                      "White-Red", "Blue-Yellow-Red",
                      "Yellow-Green-Purple"
                    ),
                    selected = "Blue-Yellow-Red"
                  ),
                  radioButtons("sc3_civci_ord1", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Original", inline = TRUE
                  ),
                  checkboxInput("sc3_civci_lab1", "Show cell info labels", value = TRUE)
                )
              )
            )
          ),
          # row 3 col 1 row 2
          fluidRow(column(12, uiOutput("sc3_civci_oup1.ui"))),
          # row 3 col 1 row 3
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("sc3_civci_oup1.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("sc3_civci_oup1.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("sc3_civci_oup1.png", "Download PNG", class = "btn-sm")
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
                selectInput("sc3_civci_inp2", "Cell info:",
                  choices = sc3conf$UI,
                  selected = sc3def$meta2
                ) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Cell info to colour cells by",
                    content = c(
                      "Select cell info to colour cells",
                      "- Categorical covariates have a fixed colour palette",
                      paste0(
                        "- Continuous covariates are coloured in a ",
                        "Blue-Yellow-Red colour scheme, which can be ",
                        "changed in the plot controls"
                      )
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("sc3_civci_tog2", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.sc3_civci_tog2 == true",
                  radioButtons("sc3_civci_col2", "Colour (Continuous data):",
                    choices = c(
                      "White-Red", "Blue-Yellow-Red",
                      "Yellow-Green-Purple"
                    ),
                    selected = "Blue-Yellow-Red"
                  ),
                  radioButtons("sc3_civci_ord2", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Original", inline = TRUE
                  ),
                  checkboxInput("sc3_civci_lab2", "Show cell info labels", value = TRUE)
                )
              )
            )
          ),
          fluidRow(column(12, uiOutput("sc3_civci_oup2.ui"))),
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("sc3_civci_oup2.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("sc3_civci_oup2.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("sc3_civci_oup2.png", "Download PNG", class = "btn-sm")
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
            selectInput("sc3_gevge_drX", "X-axis:",
              choices = sc3conf[dimred == TRUE]$UI,
              selected = sc3def$dimred[1]
            ),
            selectInput("sc3_gevge_drY", "Y-axis:",
              choices = sc3conf[dimred == TRUE]$UI,
              selected = sc3def$dimred[2]
            )
          )
        ), # row 1 col 1
        # row 2 col 2
        column(4,
          div(
            class = "input-panel",
            checkboxInput("sc3_gevge_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.sc3_gevge_togL == true",
              selectInput("sc3_gevge_sub1", "Cell info to subset:",
                choices = sc3conf[grp == TRUE]$UI,
                selected = sc3def$grp1
              ),
              uiOutput("sc3_gevge_sub1.ui"),
              actionButton("sc3_gevge_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("sc3_gevge_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            )
          )
        ), # End of column
        # row 2 col 3
        column(
          4,
          div(
            class = "input-panel",
            checkboxInput("sc3_gevge_tog0", "Adjust graphics"),
            conditionalPanel(
              condition = "input.sc3_gevge_tog0 == true",
              sliderInput("sc3_gevge_siz", "Point size:",
                min = 0, max = 4, value = 1.25, step = 0.25
              ),
              radioButtons("sc3_gevge_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("sc3_gevge_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              radioButtons("sc3_gevge_asp", "Aspect ratio:",
                choices = c("Square", "Fixed", "Free"),
                selected = "Square", inline = TRUE
              ),
              checkboxInput("sc3_gevge_txt", "Show axis text", value = FALSE)
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
                selectInput("sc3_gevge_inp1", "Gene name:", choices = NULL) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Gene expression to colour cells by",
                    content = c(
                      "Select gene to colour cells by gene expression",
                      paste0(
                        "- Gene expression are coloured in a ",
                        "White-Red colour scheme which can be ",
                        "changed in the plot controls"
                      )
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("sc3_gevge_tog1", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.sc3_gevge_tog1 == true",
                  radioButtons("sc3_gevge_col1", "Colour:",
                    choices = c(
                      "White-Red", "Blue-Yellow-Red",
                      "Yellow-Green-Purple"
                    ),
                    selected = "White-Red"
                  ),
                  radioButtons("sc3_gevge_ord1", "Plot order:",
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
            column(12, uiOutput("sc3_gevge_oup1.ui"))
          ),
          # row 3 col 1 row 3
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("sc3_gevge_oup1.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("sc3_gevge_oup1.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("sc3_gevge_oup1.png", "Download PNG", class = "btn-sm")
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
                selectInput("sc3_gevge_inp2", "Gene name:", choices = NULL) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Gene expression to colour cells by",
                    content = c(
                      "Select gene to colour cells by gene expression",
                      paste0(
                        "- Gene expression are coloured in a ",
                        "White-Red colour scheme which can be ",
                        "changed in the plot controls"
                      )
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("sc3_gevge_tog2", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.sc3_gevge_tog2 == true",
                  radioButtons("sc3_gevge_col2", "Colour:",
                    choices = c(
                      "White-Red", "Blue-Yellow-Red",
                      "Yellow-Green-Purple"
                    ),
                    selected = "White-Red"
                  ),
                  radioButtons("sc3_gevge_ord2", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Max", inline = TRUE
                  )
                )
              )
            )
          ),
          fluidRow(
            class = "tab-section",
            column(12, uiOutput("sc3_gevge_oup2.ui"))
          ),
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                  numericInput("sc3_gevge_oup2.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                  downloadButton("sc3_gevge_oup2.pdf", "Download PDF", class = "btn-sm"),
                  downloadButton("sc3_gevge_oup2.png", "Download PNG", class = "btn-sm")
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
            textAreaInput("sc3_gem_inp", "Gene names:",
                          height = "100px",
                          value = paste0(sc3def$genes, collapse = ", ")
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "List of genes to plot on bubbleplot / heatmap",
                content = c(
                  "Input genes to plot",
                  "- Maximum 16 genes (due to ploting space limitations)",
                  "- Genes should be separated by comma, semicolon or newline"
                )
              ),
            selectInput("sc3_gem_drX", "X-axis:",
                        choices = sc3conf[dimred == TRUE]$UI,
                        selected = sc3def$dimred[1]
            ),
            selectInput("sc3_gem_drY", "Y-axis:",
                        choices = sc3conf[dimred == TRUE]$UI,
                        selected = sc3def$dimred[2]
            ),
            checkboxInput("sc3_gem_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.sc3_gem_togL == true",
              selectInput("sc3_gem_sub1", "Cell info to subset:",
                          choices = sc3conf[grp == TRUE]$UI,
                          selected = sc3def$grp1
              ),
              uiOutput("sc3_gem_sub1.ui"),
              actionButton("sc3_gem_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("sc3_gem_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            ),
            checkboxInput("sc3_gem_tog0", "Adjust graphics"),
            conditionalPanel(
              condition = "input.sc3_gem_tog0 == true",
              sliderInput("sc3_gem_siz", "Point size:",
                          min = 0, max = 3, value = 0.5, step = 0.1
              ),
              radioButtons("sc3_gem_psz", "Plot size:",
                           choices = c("Small", "Medium", "Large"),
                           selected = "Medium", inline = TRUE
              ),
              radioButtons("sc3_gem_fsz", "Font size:",
                           choices = c("Smaller", "Small", "Medium", "Large"),
                           selected = "Small", inline = TRUE
              ),
              radioButtons("sc3_gem_asp", "Aspect ratio:",
                           choices = c("Square", "Fixed", "Free"),
                           selected = "Square", inline = TRUE
              ),
              checkboxInput("sc3_gem_txt", "Show axis text", value = FALSE),
              radioButtons("sc3_gem_col", "Colour (Continuous data):",
                           choices = c(
                             "White-Red", "Blue-Yellow-Red",
                             "Yellow-Green-Purple"
                           ),
                           selected = "Blue-Yellow-Red"
              ),
              radioButtons("sc3_gem_ord", "Plot order:",
                           choices = c("Max", "Min", "Original", "Random"),
                           selected = "Max", inline = TRUE
              ),
              numericInput("sc3_gem_ncol", "Number of columns", value = 0, min = 0, step = 1)
            )
          ),
          div(
            class = "input-panel",
            fluidRow(
            column(4,
              numericInput("sc3_gem_oup1.height", "Height:", min = 1, max = 50, value = 25, step = 2)
            ),
            column(4,
              numericInput("sc3_gem_oup1.width", "Width:", min = 1, max = 50, value = 25, step = 2)
            ),
            column(4,
              numericInput("sc3_gem_oup1.res", "Res:", min = 72, max = 600, value = 150, step = 5)
            )
            ),
            downloadButton("sc3_gem_oup1.pdf", "Download PDF", class = "btn-sm"),
            downloadButton("sc3_gem_oup1.png", "Download PNG", class = "btn-sm")
          )
        )
               )
      ),
      column(8,
             uiOutput("sc3_gem_oup1.ui")
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
                   selectInput("sc3_gec_drX", "X-axis:",
                               choices = sc3conf[dimred == TRUE]$UI,
                               selected = sc3def$dimred[1]
                   ),
                   selectInput("sc3_gec_drY", "Y-axis:",
                               choices = sc3conf[dimred == TRUE]$UI,
                               selected = sc3def$dimred[2]
                   ),
                   selectInput("sc3_gec_inp1", "Gene 1:", choices = NULL) %>%
                     helper(
                       type = "inline", size = "m", fade = TRUE,
                       title = "Gene expression to colour cells by",
                       content = c(
                         "Select gene to colour cells by gene expression",
                         paste0(
                           "- Gene expression are coloured in a ",
                           "White-Red colour scheme which can be ",
                           "changed in the plot controls"
                         )
                       )
                     ),
                   selectInput("sc3_gec_inp2", "Gene 2:", choices = NULL) %>%
                     helper(
                       type = "inline", size = "m", fade = TRUE,
                       title = "Gene expression to colour cells by",
                       content = c(
                         "Select gene to colour cells by gene expression",
                         paste0(
                           "- Gene expression are coloured in a ",
                           "White-Blue colour scheme which can be ",
                           "changed in the plot controls"
                         )
                       )
                     ),
                   checkboxInput("sc3_gec_togL", "Subset cells"),
                   conditionalPanel(
                    condition = "input.sc3_gec_togL == true",
                    selectInput("sc3_gec_sub1", "Cell info to subset:", choices = sc3conf[grp == TRUE]$UI, selected = sc3def$grp1),
                     uiOutput("sc3_gec_sub1.ui"),
                     actionButton("sc3_gec_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
                     actionButton("sc3_gec_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
                   ),
                   checkboxInput("sc3_gec_tog0", "Adjust graphics"),
                   conditionalPanel(
                     condition = "input.sc3_gec_tog0 == true",
                     radioButtons("sc3_gec_col1", "Colour:",
                                  choices = c(
                                    "Red (Gene1); Blue (Gene2)",
                                    "Orange (Gene1); Blue (Gene2)",
                                    "Red (Gene1); Green (Gene2)",
                                    "Green (Gene1); Blue (Gene2)"
                                  ),
                                  selected = "Red (Gene1); Blue (Gene2)"
                     ),
                     radioButtons("sc3_gec_ord1", "Plot order:",
                                  choices = c("Max", "Min", "Original", "Random"),
                                  selected = "Max", inline = TRUE
                     ),
                     sliderInput("sc3_gec_siz", "Point size:",
                                 min = 0, max = 4, value = 1.25, step = 0.25
                     ),
                     radioButtons("sc3_gec_psz", "Plot size:",
                                  choices = c("Small", "Medium", "Large"),
                                  selected = "Medium", inline = TRUE
                     ),
                     radioButtons("sc3_gec_fsz", "Font size:",
                                  choices = c("Small", "Medium", "Large"),
                                  selected = "Small", inline = TRUE
                     ),
                     radioButtons("sc3_gec_asp", "Aspect ratio:",
                                  choices = c("Square", "Fixed", "Free"),
                                  selected = "Square", inline = TRUE
                     ),
                     checkboxInput("sc3_gec_txt", "Show axis text", value = FALSE)
                   )
                 ),
                 div(class="input-panel input-panel-section",
                     numericInput("sc3_gec_oup1.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                     downloadButton("sc3_gec_oup1.pdf", "Download PDF", class = "btn-sm"),
                     downloadButton("sc3_gec_oup1.png", "Download PNG", class = "btn-sm")
                 ),
                 div(class="input-panel-section",
                     h4("Cell numbers"),
                     dataTableOutput("sc3_gec_.dt")
                 )
               )
        ), # row 2 col 1
        # row 2 col 2
        column(
          8,
          uiOutput("sc3_gec_oup1.ui"),
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
            selectInput("sc3_vio_inp1", "Cell info (X-axis):",
              choices = sc3conf[grp == TRUE]$UI,
              selected = sc3def$grp1
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell info to group cells by",
                content = c(
                  "Select categorical cell info to group cells by",
                  "- Single cells are grouped by this categorical covariate",
                  "- Plotted as the X-axis of the violin plot / box plot"
                )
              ),
            selectInput("sc3_vio_inp2", "Cell info / Gene (Y-axis):", choices = NULL) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell Info / Gene to plot",
                content = c(
                  "Select cell info / gene to plot on Y-axis",
                  "- Can be continuous cell info (e.g. nUMIs / scores)",
                  "- Can also be gene expression"
                )
              ),
            radioButtons("sc3_vio_typ", "Plot type:",
              choices = c("violin", "boxplot", "lineplot"),
              selected = "violin", inline = TRUE
            ),
            checkboxInput("sc3_vio_pts", "Show data points", value = FALSE),
            checkboxInput("sc3_vio_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.sc3_vio_togL == true",
              selectInput("sc3_vio_sub1", "Cell info to subset:",
                choices = sc3conf[grp == TRUE]$UI,
                selected = sc3def$grp1
              ),
              uiOutput("sc3_vio_sub1.ui"),
              actionButton("sc3_vio_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("sc3_vio_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            ),
            checkboxInput("sc3_vio_tog", "Adjust graphics"),
            conditionalPanel(
              condition = "input.sc3_vio_tog == true",
              sliderInput("sc3_vio_siz", "Data point size:",
                min = 0, max = 4, value = 1.25, step = 0.25
              ),
              radioButtons("sc3_vio_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              radioButtons("sc3_vio_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              conditionalPanel(
              condition = "input.sc3_vio_typ == 'lineplot'",
              sliderInput("sc3_vio_barsz", "Line size", min = 0.05, max = 0.5, step = 0.01, value = 0.3)
              )
            ),
            selectInput("sc3_vio_datatype", "Data type", choices = c("normalised", "raw"), selected = "normalised"),
          ),
          div(
            class = "input-panel",
            numericInput("sc3_vio_oup.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
            downloadButton("sc3_vio_oup.pdf", "Download PDF", class = "btn-sm"),
            downloadButton("sc3_vio_oup.png", "Download PNG", class = "btn-sm")
          )
        ), # row 2 col 1
        # row 2 col 2
        column(
          9, uiOutput("sc3_vio_oup.ui")
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
            selectInput("sc3_pro_inp1", "Cell info to plot (X-axis):",
              choices = sc3conf[grp == TRUE]$UI,
              selected = sc3def$grp2
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell info to plot cells by",
                content = c(
                  "Select categorical cell info to plot cells by",
                  "- Plotted as the X-axis of the proportion plot"
                )
              ),
            selectInput("sc3_pro_inp2", "Cell info to group / colour by:",
              choices = sc3conf[grp == TRUE]$UI,
              selected = sc3def$grp1
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell info to group / colour cells by",
                content = c(
                  "Select categorical cell info to group / colour cells by",
                  "- Proportion / cell numbers are shown in different colours"
                )
              ),
            radioButtons("sc3_pro_typ", "Plot value:",
              choices = c("Proportion", "CellNumbers"),
              selected = "Proportion", inline = TRUE
            ),
            checkboxInput("sc3_pro_flp", "Flip X/Y", value = FALSE),
            checkboxInput("sc3_pro_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.sc3_pro_togL == true",
              selectInput("sc3_pro_sub1", "Cell info to subset:",
                choices = sc3conf[grp == TRUE]$UI,
                selected = sc3def$grp1
              ),
              uiOutput("sc3_pro_sub1.ui"),
              actionButton("sc3_pro_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("sc3_pro_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            ),
            checkboxInput("sc3_pro_tog", "Adjust graphics"),
            conditionalPanel(
              condition = "input.sc3_pro_tog == true",
              radioButtons("sc3_pro_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("sc3_pro_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              )
            )
          ),
          div(
            class = "input-panel",
            numericInput("sc3_pro_oup.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
            downloadButton("sc3_pro_oup.pdf", "Download PDF", class = "btn-sm"),
            downloadButton("sc3_pro_oup.png", "Download PNG", class = "btn-sm")
          )
        ), # row 2 col 1
        # row 2 col 2
        column(
          9, uiOutput("sc3_pro_oup.ui")
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
          3,
          div(
            class = "input-panel",
            textAreaInput("sc3_hea_inp", "Gene names",
              height = "100px",
              value = paste0(sc3def$genes, collapse = ", ")
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "List of genes to plot on bubbleplot / heatmap",
                content = c(
                  "Input genes to plot",
                  "- Maximum 50 genes (due to ploting space limitations)",
                  "- Genes should be separated by comma, semicolon or newline"
                )
              ),
            selectInput("sc3_hea_grp", "Group by:",
              choices = sc3conf[grp == TRUE]$UI,
              selected = sc3conf[grp == TRUE]$UI[1]
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell info to group cells by",
                content = c(
                  "Select categorical cell info to group cells by",
                  "- Single cells are grouped by this categorical covariate",
                  "- Plotted as the X-axis of the bubbleplot / heatmap"
                )
              ),
            radioButtons("sc3_hea_plt", "Plot type:",
              choices = c("Bubbleplot", "Heatmap"),
              selected = "Bubbleplot", inline = TRUE
            ),
            checkboxInput("sc3_hea_scl", "Scale gene expression", value = TRUE),
            checkboxInput("sc3_hea_row", "Cluster rows (genes)", value = TRUE),
            checkboxInput("sc3_hea_col", "Cluster columns (samples)", value = FALSE),
            checkboxInput("sc3_hea_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.sc3_hea_togL == true",
              selectInput("sc3_hea_sub1", "Cell info to subset:",
                choices = sc3conf[grp == TRUE]$UI,
                selected = sc3def$grp1
              ),
              uiOutput("sc3_hea_sub1.ui"),
              actionButton("sc3_hea_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("sc3_hea_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            ),
            checkboxInput("sc3_hea_tog", "Adjust graphics"),
            conditionalPanel(
              condition = "input.sc3_hea_tog == true",
              radioButtons("sc3_hea_cols", "Colour scheme:",
                choices = c(
                  "White-Red", "Blue-Yellow-Red",
                  "Yellow-Green-Purple"
                ),
                selected = "Blue-Yellow-Red"
              ),
              radioButtons("sc3_hea_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("sc3_hea_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              )
            )
          ),
          div(
            class = "input-panel",
            numericInput("sc3_hea_oup.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
            downloadButton("sc3_hea_oup.pdf", "Download PDF", class = "btn-sm"),
            downloadButton("sc3_hea_oup.png", "Download PNG", class = "btn-sm")
          )
        ), # row 2 col 1
        # row 2 col 2
        column(
          9, h4(htmlOutput("sc3_hea_oupTxt")),
          uiOutput("sc3_hea_oup.ui")
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
                     selectInput("sc3_mar_cls","Select clustering:", choices = names(sc3mar),selected = 1)
                   )
                 )
               )
        )
      ), # end of row 2
      # row 3 ----
      fluidRow(
      column(12,
        DTOutput("sc3_mar_table")
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
            selectInput("sc4_civge_drX", "X-axis:",
              choices = sc4conf[dimred == TRUE]$UI,
              selected = sc4def$dimred[1]
            ),
            selectInput("sc4_civge_drY", "Y-axis:",
              choices = sc4conf[dimred == TRUE]$UI,
              selected = sc4def$dimred[2]
            )
          )
        ), # row 1 col 1
        # row 2 col 2
        column(
          4,
          div(
            class = "input-panel",
            checkboxInput("sc4_civge_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.sc4_civge_togL == true",
              selectInput("sc4_civge_sub1", "Cell info to subset:",
                choices = sc4conf[grp == TRUE]$UI,
                selected = sc4def$grp1
              ),
              uiOutput("sc4_civge_sub1.ui"),
              actionButton("sc4_civge_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("sc4_civge_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            )
          )
        ), # End of column
        # row 2 col 3
        column(
          4,
          div(
            class = "input-panel",
            checkboxInput("sc4_civge_tog0", "Adjust graphics"),
            conditionalPanel(
              condition = "input.sc4_civge_tog0 == true",
              sliderInput("sc4_civge_siz", "Point size:",
                min = 0, max = 4, value = 1.25, step = 0.25
              ),
              radioButtons("sc4_civge_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("sc4_civge_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              radioButtons("sc4_civge_asp", "Aspect ratio:",
                choices = c("Square", "Fixed", "Free"),
                selected = "Square", inline = TRUE
              ),
              checkboxInput("sc4_civge_txt", "Show axis text", value = FALSE)
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
                selectInput("sc4_civge_inp1", "Cell info:",
                  choices = sc4conf$UI,
                  selected = sc4def$meta1
                ) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Cell info to colour cells by",
                    content = c(
                      "Select cell info to colour cells",
                      "- Categorical covariates have a fixed colour palette",
                      paste0(
                        "- Continuous covariates are coloured in a ",
                        "Blue-Yellow-Red colour scheme, which can be ",
                        "changed in the plot controls"
                      )
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("sc4_civge_tog1", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.sc4_civge_tog1 == true",
                  radioButtons("sc4_civge_col1", "Colour (Continuous data):",
                    choices = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple"),
                    selected = "Blue-Yellow-Red"
                  ),
                  radioButtons("sc4_civge_ord1", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Original", inline = TRUE
                  ),
                  checkboxInput("sc4_civge_lab1", "Show cell info labels", value = TRUE)
                )
              )
            )
          ),
          # row 3 col 1 row 2
          fluidRow(
            class = "tab-section",
            column(
              12,
              uiOutput("sc4_civge_oup1.ui")
            )
          ),
          # row 3 col 1 row 3
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("sc4_civge_oup1.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("sc4_civge_oup1.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("sc4_civge_oup1.png", "Download PNG", class = "btn-sm"),
                checkboxInput("sc4_civge_tog9", "Show cell numbers / statistics")
              )
            )
          ),
          conditionalPanel(
            condition = "input.sc4_civge_tog9 == true",
            h4("Cell numbers / statistics"),
            radioButtons("sc4_civge_splt", "Split continuous cell info into:",
              choices = c("Quartile", "Decile"),
              selected = "Decile", inline = TRUE
            ),
            dataTableOutput("sc4_civge_.dt")
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
                selectInput("sc4_civge_inp2", "Gene name:", choices = NULL) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Gene expression to colour cells by",
                    content = c(
                      "Select gene to colour cells by gene expression",
                      paste0(
                        "- Gene expression are coloured in a ",
                        "White-Red colour scheme which can be ",
                        "changed in the plot controls"
                      )
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("sc4_civge_tog2", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.sc4_civge_tog2 == true",
                  radioButtons("sc4_civge_col2", "Colour:",
                    choices = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple"),
                    selected = "White-Red"
                  ),
                  radioButtons("sc4_civge_ord2", "Plot order:",
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
              uiOutput("sc4_civge_oup2.ui")
            )
          ),
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("sc4_civge_oup2.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("sc4_civge_oup2.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("sc4_civge_oup2.png", "Download PNG", class = "btn-sm")
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
            selectInput("sc4_civci_drX", "X-axis:",
              choices = sc4conf[dimred == TRUE]$UI,
              selected = sc4def$dimred[1]
            ),
            selectInput("sc4_civci_drY", "Y-axis:",
              choices = sc4conf[dimred == TRUE]$UI,
              selected = sc4def$dimred[2]
            )
          )
        ), # row 2 col 2
        # row 2 col 2
        column(4,
          div(
            class = "input-panel",
            checkboxInput("sc4_civci_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.sc4_civci_togL == true",
              selectInput("sc4_civci_sub1", "Cell info to subset:",
                choices = sc4conf[grp == TRUE]$UI,
                selected = sc4def$grp1
              ),
              uiOutput("sc4_civci_sub1.ui"),
              actionButton("sc4_civci_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("sc4_civci_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            )
          )
        ), # row 2 col 2
        # row 2 col 3
        column(
          4,
          div(
            class = "input-panel",
            checkboxInput("sc4_civci_tog0", "Adjust graphics"),
            conditionalPanel(
              condition = "input.sc4_civci_tog0 == true",
              sliderInput("sc4_civci_siz", "Point size:",
                min = 0, max = 4, value = 1.25, step = 0.25
              ),
              radioButtons("sc4_civci_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("sc4_civci_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              radioButtons("sc4_civci_asp", "Aspect ratio:",
                choices = c("Square", "Fixed", "Free"),
                selected = "Square", inline = TRUE
              ),
              checkboxInput("sc4_civci_txt", "Show axis text", value = FALSE)
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
                selectInput("sc4_civci_inp1", "Cell info:",
                  choices = sc4conf$UI,
                  selected = sc4def$meta1
                ) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Cell info to colour cells by",
                    content = c(
                      "Select cell info to colour cells",
                      "- Categorical covariates have a fixed colour palette",
                      paste0(
                        "- Continuous covariates are coloured in a ",
                        "Blue-Yellow-Red colour scheme, which can be ",
                        "changed in the plot controls"
                      )
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("sc4_civci_tog1", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.sc4_civci_tog1 == true",
                  radioButtons("sc4_civci_col1", "Colour (Continuous data):",
                    choices = c(
                      "White-Red", "Blue-Yellow-Red",
                      "Yellow-Green-Purple"
                    ),
                    selected = "Blue-Yellow-Red"
                  ),
                  radioButtons("sc4_civci_ord1", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Original", inline = TRUE
                  ),
                  checkboxInput("sc4_civci_lab1", "Show cell info labels", value = TRUE)
                )
              )
            )
          ),
          # row 3 col 1 row 2
          fluidRow(column(12, uiOutput("sc4_civci_oup1.ui"))),
          # row 3 col 1 row 3
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("sc4_civci_oup1.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("sc4_civci_oup1.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("sc4_civci_oup1.png", "Download PNG", class = "btn-sm")
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
                selectInput("sc4_civci_inp2", "Cell info:",
                  choices = sc4conf$UI,
                  selected = sc4def$meta2
                ) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Cell info to colour cells by",
                    content = c(
                      "Select cell info to colour cells",
                      "- Categorical covariates have a fixed colour palette",
                      paste0(
                        "- Continuous covariates are coloured in a ",
                        "Blue-Yellow-Red colour scheme, which can be ",
                        "changed in the plot controls"
                      )
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("sc4_civci_tog2", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.sc4_civci_tog2 == true",
                  radioButtons("sc4_civci_col2", "Colour (Continuous data):",
                    choices = c(
                      "White-Red", "Blue-Yellow-Red",
                      "Yellow-Green-Purple"
                    ),
                    selected = "Blue-Yellow-Red"
                  ),
                  radioButtons("sc4_civci_ord2", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Original", inline = TRUE
                  ),
                  checkboxInput("sc4_civci_lab2", "Show cell info labels", value = TRUE)
                )
              )
            )
          ),
          fluidRow(column(12, uiOutput("sc4_civci_oup2.ui"))),
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("sc4_civci_oup2.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("sc4_civci_oup2.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("sc4_civci_oup2.png", "Download PNG", class = "btn-sm")
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
            selectInput("sc4_gevge_drX", "X-axis:",
              choices = sc4conf[dimred == TRUE]$UI,
              selected = sc4def$dimred[1]
            ),
            selectInput("sc4_gevge_drY", "Y-axis:",
              choices = sc4conf[dimred == TRUE]$UI,
              selected = sc4def$dimred[2]
            )
          )
        ), # row 1 col 1
        # row 2 col 2
        column(4,
          div(
            class = "input-panel",
            checkboxInput("sc4_gevge_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.sc4_gevge_togL == true",
              selectInput("sc4_gevge_sub1", "Cell info to subset:",
                choices = sc4conf[grp == TRUE]$UI,
                selected = sc4def$grp1
              ),
              uiOutput("sc4_gevge_sub1.ui"),
              actionButton("sc4_gevge_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("sc4_gevge_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            )
          )
        ), # End of column
        # row 2 col 3
        column(
          4,
          div(
            class = "input-panel",
            checkboxInput("sc4_gevge_tog0", "Adjust graphics"),
            conditionalPanel(
              condition = "input.sc4_gevge_tog0 == true",
              sliderInput("sc4_gevge_siz", "Point size:",
                min = 0, max = 4, value = 1.25, step = 0.25
              ),
              radioButtons("sc4_gevge_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("sc4_gevge_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              radioButtons("sc4_gevge_asp", "Aspect ratio:",
                choices = c("Square", "Fixed", "Free"),
                selected = "Square", inline = TRUE
              ),
              checkboxInput("sc4_gevge_txt", "Show axis text", value = FALSE)
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
                selectInput("sc4_gevge_inp1", "Gene name:", choices = NULL) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Gene expression to colour cells by",
                    content = c(
                      "Select gene to colour cells by gene expression",
                      paste0(
                        "- Gene expression are coloured in a ",
                        "White-Red colour scheme which can be ",
                        "changed in the plot controls"
                      )
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("sc4_gevge_tog1", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.sc4_gevge_tog1 == true",
                  radioButtons("sc4_gevge_col1", "Colour:",
                    choices = c(
                      "White-Red", "Blue-Yellow-Red",
                      "Yellow-Green-Purple"
                    ),
                    selected = "White-Red"
                  ),
                  radioButtons("sc4_gevge_ord1", "Plot order:",
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
            column(12, uiOutput("sc4_gevge_oup1.ui"))
          ),
          # row 3 col 1 row 3
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                numericInput("sc4_gevge_oup1.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                downloadButton("sc4_gevge_oup1.pdf", "Download PDF", class = "btn-sm"),
                downloadButton("sc4_gevge_oup1.png", "Download PNG", class = "btn-sm")
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
                selectInput("sc4_gevge_inp2", "Gene name:", choices = NULL) %>%
                  helper(
                    type = "inline", size = "m", fade = TRUE,
                    title = "Gene expression to colour cells by",
                    content = c(
                      "Select gene to colour cells by gene expression",
                      paste0(
                        "- Gene expression are coloured in a ",
                        "White-Red colour scheme which can be ",
                        "changed in the plot controls"
                      )
                    )
                  )
              )
            ),
            column(
              6,
              div(
                class = "input-panel",
                checkboxInput("sc4_gevge_tog2", "Adjust graphics"),
                conditionalPanel(
                  condition = "input.sc4_gevge_tog2 == true",
                  radioButtons("sc4_gevge_col2", "Colour:",
                    choices = c(
                      "White-Red", "Blue-Yellow-Red",
                      "Yellow-Green-Purple"
                    ),
                    selected = "White-Red"
                  ),
                  radioButtons("sc4_gevge_ord2", "Plot order:",
                    choices = c("Max", "Min", "Original", "Random"),
                    selected = "Max", inline = TRUE
                  )
                )
              )
            )
          ),
          fluidRow(
            class = "tab-section",
            column(12, uiOutput("sc4_gevge_oup2.ui"))
          ),
          fluidRow(
            class = "tab-section",
            column(
              12,
              div(
                class = "input-panel",
                  numericInput("sc4_gevge_oup2.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                  downloadButton("sc4_gevge_oup2.pdf", "Download PDF", class = "btn-sm"),
                  downloadButton("sc4_gevge_oup2.png", "Download PNG", class = "btn-sm")
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
            textAreaInput("sc4_gem_inp", "Gene names:",
                          height = "100px",
                          value = paste0(sc4def$genes, collapse = ", ")
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "List of genes to plot on bubbleplot / heatmap",
                content = c(
                  "Input genes to plot",
                  "- Maximum 16 genes (due to ploting space limitations)",
                  "- Genes should be separated by comma, semicolon or newline"
                )
              ),
            selectInput("sc4_gem_drX", "X-axis:",
                        choices = sc4conf[dimred == TRUE]$UI,
                        selected = sc4def$dimred[1]
            ),
            selectInput("sc4_gem_drY", "Y-axis:",
                        choices = sc4conf[dimred == TRUE]$UI,
                        selected = sc4def$dimred[2]
            ),
            checkboxInput("sc4_gem_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.sc4_gem_togL == true",
              selectInput("sc4_gem_sub1", "Cell info to subset:",
                          choices = sc4conf[grp == TRUE]$UI,
                          selected = sc4def$grp1
              ),
              uiOutput("sc4_gem_sub1.ui"),
              actionButton("sc4_gem_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("sc4_gem_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            ),
            checkboxInput("sc4_gem_tog0", "Adjust graphics"),
            conditionalPanel(
              condition = "input.sc4_gem_tog0 == true",
              sliderInput("sc4_gem_siz", "Point size:",
                          min = 0, max = 3, value = 0.5, step = 0.1
              ),
              radioButtons("sc4_gem_psz", "Plot size:",
                           choices = c("Small", "Medium", "Large"),
                           selected = "Medium", inline = TRUE
              ),
              radioButtons("sc4_gem_fsz", "Font size:",
                           choices = c("Smaller", "Small", "Medium", "Large"),
                           selected = "Small", inline = TRUE
              ),
              radioButtons("sc4_gem_asp", "Aspect ratio:",
                           choices = c("Square", "Fixed", "Free"),
                           selected = "Square", inline = TRUE
              ),
              checkboxInput("sc4_gem_txt", "Show axis text", value = FALSE),
              radioButtons("sc4_gem_col", "Colour (Continuous data):",
                           choices = c(
                             "White-Red", "Blue-Yellow-Red",
                             "Yellow-Green-Purple"
                           ),
                           selected = "Blue-Yellow-Red"
              ),
              radioButtons("sc4_gem_ord", "Plot order:",
                           choices = c("Max", "Min", "Original", "Random"),
                           selected = "Max", inline = TRUE
              ),
              numericInput("sc4_gem_ncol", "Number of columns", value = 0, min = 0, step = 1)
            )
          ),
          div(
            class = "input-panel",
            fluidRow(
            column(4,
              numericInput("sc4_gem_oup1.height", "Height:", min = 1, max = 50, value = 25, step = 2)
            ),
            column(4,
              numericInput("sc4_gem_oup1.width", "Width:", min = 1, max = 50, value = 25, step = 2)
            ),
            column(4,
              numericInput("sc4_gem_oup1.res", "Res:", min = 72, max = 600, value = 150, step = 5)
            )
            ),
            downloadButton("sc4_gem_oup1.pdf", "Download PDF", class = "btn-sm"),
            downloadButton("sc4_gem_oup1.png", "Download PNG", class = "btn-sm")
          )
        )
               )
      ),
      column(8,
             uiOutput("sc4_gem_oup1.ui")
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
                   selectInput("sc4_gec_drX", "X-axis:",
                               choices = sc4conf[dimred == TRUE]$UI,
                               selected = sc4def$dimred[1]
                   ),
                   selectInput("sc4_gec_drY", "Y-axis:",
                               choices = sc4conf[dimred == TRUE]$UI,
                               selected = sc4def$dimred[2]
                   ),
                   selectInput("sc4_gec_inp1", "Gene 1:", choices = NULL) %>%
                     helper(
                       type = "inline", size = "m", fade = TRUE,
                       title = "Gene expression to colour cells by",
                       content = c(
                         "Select gene to colour cells by gene expression",
                         paste0(
                           "- Gene expression are coloured in a ",
                           "White-Red colour scheme which can be ",
                           "changed in the plot controls"
                         )
                       )
                     ),
                   selectInput("sc4_gec_inp2", "Gene 2:", choices = NULL) %>%
                     helper(
                       type = "inline", size = "m", fade = TRUE,
                       title = "Gene expression to colour cells by",
                       content = c(
                         "Select gene to colour cells by gene expression",
                         paste0(
                           "- Gene expression are coloured in a ",
                           "White-Blue colour scheme which can be ",
                           "changed in the plot controls"
                         )
                       )
                     ),
                   checkboxInput("sc4_gec_togL", "Subset cells"),
                   conditionalPanel(
                    condition = "input.sc4_gec_togL == true",
                    selectInput("sc4_gec_sub1", "Cell info to subset:", choices = sc4conf[grp == TRUE]$UI, selected = sc4def$grp1),
                     uiOutput("sc4_gec_sub1.ui"),
                     actionButton("sc4_gec_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
                     actionButton("sc4_gec_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
                   ),
                   checkboxInput("sc4_gec_tog0", "Adjust graphics"),
                   conditionalPanel(
                     condition = "input.sc4_gec_tog0 == true",
                     radioButtons("sc4_gec_col1", "Colour:",
                                  choices = c(
                                    "Red (Gene1); Blue (Gene2)",
                                    "Orange (Gene1); Blue (Gene2)",
                                    "Red (Gene1); Green (Gene2)",
                                    "Green (Gene1); Blue (Gene2)"
                                  ),
                                  selected = "Red (Gene1); Blue (Gene2)"
                     ),
                     radioButtons("sc4_gec_ord1", "Plot order:",
                                  choices = c("Max", "Min", "Original", "Random"),
                                  selected = "Max", inline = TRUE
                     ),
                     sliderInput("sc4_gec_siz", "Point size:",
                                 min = 0, max = 4, value = 1.25, step = 0.25
                     ),
                     radioButtons("sc4_gec_psz", "Plot size:",
                                  choices = c("Small", "Medium", "Large"),
                                  selected = "Medium", inline = TRUE
                     ),
                     radioButtons("sc4_gec_fsz", "Font size:",
                                  choices = c("Small", "Medium", "Large"),
                                  selected = "Small", inline = TRUE
                     ),
                     radioButtons("sc4_gec_asp", "Aspect ratio:",
                                  choices = c("Square", "Fixed", "Free"),
                                  selected = "Square", inline = TRUE
                     ),
                     checkboxInput("sc4_gec_txt", "Show axis text", value = FALSE)
                   )
                 ),
                 div(class="input-panel input-panel-section",
                     numericInput("sc4_gec_oup1.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
                     downloadButton("sc4_gec_oup1.pdf", "Download PDF", class = "btn-sm"),
                     downloadButton("sc4_gec_oup1.png", "Download PNG", class = "btn-sm")
                 ),
                 div(class="input-panel-section",
                     h4("Cell numbers"),
                     dataTableOutput("sc4_gec_.dt")
                 )
               )
        ), # row 2 col 1
        # row 2 col 2
        column(
          8,
          uiOutput("sc4_gec_oup1.ui"),
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
            selectInput("sc4_vio_inp1", "Cell info (X-axis):",
              choices = sc4conf[grp == TRUE]$UI,
              selected = sc4def$grp1
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell info to group cells by",
                content = c(
                  "Select categorical cell info to group cells by",
                  "- Single cells are grouped by this categorical covariate",
                  "- Plotted as the X-axis of the violin plot / box plot"
                )
              ),
            selectInput("sc4_vio_inp2", "Cell info / Gene (Y-axis):", choices = NULL) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell Info / Gene to plot",
                content = c(
                  "Select cell info / gene to plot on Y-axis",
                  "- Can be continuous cell info (e.g. nUMIs / scores)",
                  "- Can also be gene expression"
                )
              ),
            radioButtons("sc4_vio_typ", "Plot type:",
              choices = c("violin", "boxplot", "lineplot"),
              selected = "violin", inline = TRUE
            ),
            checkboxInput("sc4_vio_pts", "Show data points", value = FALSE),
            checkboxInput("sc4_vio_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.sc4_vio_togL == true",
              selectInput("sc4_vio_sub1", "Cell info to subset:",
                choices = sc4conf[grp == TRUE]$UI,
                selected = sc4def$grp1
              ),
              uiOutput("sc4_vio_sub1.ui"),
              actionButton("sc4_vio_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("sc4_vio_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            ),
            checkboxInput("sc4_vio_tog", "Adjust graphics"),
            conditionalPanel(
              condition = "input.sc4_vio_tog == true",
              sliderInput("sc4_vio_siz", "Data point size:",
                min = 0, max = 4, value = 1.25, step = 0.25
              ),
              radioButtons("sc4_vio_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              radioButtons("sc4_vio_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              ),
              conditionalPanel(
              condition = "input.sc4_vio_typ == 'lineplot'",
              sliderInput("sc4_vio_barsz", "Line size", min = 0.05, max = 0.5, step = 0.01, value = 0.3)
              )
            ),
            selectInput("sc4_vio_datatype", "Data type", choices = c("normalised", "raw"), selected = "normalised"),
          ),
          div(
            class = "input-panel",
            numericInput("sc4_vio_oup.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
            downloadButton("sc4_vio_oup.pdf", "Download PDF", class = "btn-sm"),
            downloadButton("sc4_vio_oup.png", "Download PNG", class = "btn-sm")
          )
        ), # row 2 col 1
        # row 2 col 2
        column(
          9, uiOutput("sc4_vio_oup.ui")
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
            selectInput("sc4_pro_inp1", "Cell info to plot (X-axis):",
              choices = sc4conf[grp == TRUE]$UI,
              selected = sc4def$grp2
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell info to plot cells by",
                content = c(
                  "Select categorical cell info to plot cells by",
                  "- Plotted as the X-axis of the proportion plot"
                )
              ),
            selectInput("sc4_pro_inp2", "Cell info to group / colour by:",
              choices = sc4conf[grp == TRUE]$UI,
              selected = sc4def$grp1
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell info to group / colour cells by",
                content = c(
                  "Select categorical cell info to group / colour cells by",
                  "- Proportion / cell numbers are shown in different colours"
                )
              ),
            radioButtons("sc4_pro_typ", "Plot value:",
              choices = c("Proportion", "CellNumbers"),
              selected = "Proportion", inline = TRUE
            ),
            checkboxInput("sc4_pro_flp", "Flip X/Y", value = FALSE),
            checkboxInput("sc4_pro_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.sc4_pro_togL == true",
              selectInput("sc4_pro_sub1", "Cell info to subset:",
                choices = sc4conf[grp == TRUE]$UI,
                selected = sc4def$grp1
              ),
              uiOutput("sc4_pro_sub1.ui"),
              actionButton("sc4_pro_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("sc4_pro_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            ),
            checkboxInput("sc4_pro_tog", "Adjust graphics"),
            conditionalPanel(
              condition = "input.sc4_pro_tog == true",
              radioButtons("sc4_pro_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("sc4_pro_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              )
            )
          ),
          div(
            class = "input-panel",
            numericInput("sc4_pro_oup.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
            downloadButton("sc4_pro_oup.pdf", "Download PDF", class = "btn-sm"),
            downloadButton("sc4_pro_oup.png", "Download PNG", class = "btn-sm")
          )
        ), # row 2 col 1
        # row 2 col 2
        column(
          9, uiOutput("sc4_pro_oup.ui")
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
          3,
          div(
            class = "input-panel",
            textAreaInput("sc4_hea_inp", "Gene names",
              height = "100px",
              value = paste0(sc4def$genes, collapse = ", ")
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "List of genes to plot on bubbleplot / heatmap",
                content = c(
                  "Input genes to plot",
                  "- Maximum 50 genes (due to ploting space limitations)",
                  "- Genes should be separated by comma, semicolon or newline"
                )
              ),
            selectInput("sc4_hea_grp", "Group by:",
              choices = sc4conf[grp == TRUE]$UI,
              selected = sc4conf[grp == TRUE]$UI[1]
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell info to group cells by",
                content = c(
                  "Select categorical cell info to group cells by",
                  "- Single cells are grouped by this categorical covariate",
                  "- Plotted as the X-axis of the bubbleplot / heatmap"
                )
              ),
            radioButtons("sc4_hea_plt", "Plot type:",
              choices = c("Bubbleplot", "Heatmap"),
              selected = "Bubbleplot", inline = TRUE
            ),
            checkboxInput("sc4_hea_scl", "Scale gene expression", value = TRUE),
            checkboxInput("sc4_hea_row", "Cluster rows (genes)", value = TRUE),
            checkboxInput("sc4_hea_col", "Cluster columns (samples)", value = FALSE),
            checkboxInput("sc4_hea_togL", "Subset cells"),
            conditionalPanel(
              condition = "input.sc4_hea_togL == true",
              selectInput("sc4_hea_sub1", "Cell info to subset:",
                choices = sc4conf[grp == TRUE]$UI,
                selected = sc4def$grp1
              ),
              uiOutput("sc4_hea_sub1.ui"),
              actionButton("sc4_hea_sub1all", "Select all groups", class = "btn btn-primary btn-sm"),
              actionButton("sc4_hea_sub1non", "Deselect all groups", class = "btn btn-primary btn-sm")
            ),
            checkboxInput("sc4_hea_tog", "Adjust graphics"),
            conditionalPanel(
              condition = "input.sc4_hea_tog == true",
              radioButtons("sc4_hea_cols", "Colour scheme:",
                choices = c(
                  "White-Red", "Blue-Yellow-Red",
                  "Yellow-Green-Purple"
                ),
                selected = "Blue-Yellow-Red"
              ),
              radioButtons("sc4_hea_psz", "Plot size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Medium", inline = TRUE
              ),
              radioButtons("sc4_hea_fsz", "Font size:",
                choices = c("Small", "Medium", "Large"),
                selected = "Small", inline = TRUE
              )
            )
          ),
          div(
            class = "input-panel",
            numericInput("sc4_hea_oup.res", "Resolution:", min = 72, max = 600, value = 150, step = 5),
            downloadButton("sc4_hea_oup.pdf", "Download PDF", class = "btn-sm"),
            downloadButton("sc4_hea_oup.png", "Download PNG", class = "btn-sm")
          )
        ), # row 2 col 1
        # row 2 col 2
        column(
          9, h4(htmlOutput("sc4_hea_oupTxt")),
          uiOutput("sc4_hea_oup.ui")
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
                     selectInput("sc4_mar_cls","Select clustering:", choices = names(sc4mar),selected = 1)
                   )
                 )
               )
        )
      ), # end of row 2
      # row 3 ----
      fluidRow(
      column(12,
        DTOutput("sc4_mar_table")
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
