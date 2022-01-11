# UI

library(shiny)
library(shinythemes)
library(DT)
source("functions.R")

fixedPage(title="Support 5860",
          theme=shinytheme("flatly"),
          tags$head(tags$link(rel="stylesheet",type="text/css",href="styles.css")),
          
          fixedRow(
            column(12,class="box-head",
                   span(tags$img(src='nbis.png',style="height:18px;"),style="vertical-align:top;display:inline-block;"),
                   span(tags$h4("•",style="margin:0px;margin-left:6px;margin-right:6px;"),style="vertical-align:top;display:inline-block;"),
                   span(tags$h4(strong("Support 5860"),style="margin:0px;"),style="vertical-align:middle;display:inline-block;"),
                   div(style="margin-top:5px;margin-bottom:5px;",
                       span(class="help-note medium-small","All results are pre-computed.")
                   )
            ),
            column(3,class="box-left",
                   selectInput("in_celltype",label="Cell type",choices=c("LEC","BEC"),selected="LEC",multiple=FALSE),
                   selectInput("in_condition",label="Condition",choices=c("Healthy","Diseased","Combined","Subcluster"),selected="Healthy",multiple=FALSE),
                   selectizeInput("in_group",label="Group (Point colour)",choices=c("Condition","Study","Cluster"),selected=1,multiple=FALSE),
                   uiOutput("ui_cluster"),
                   sliderInput("in_width","Image width",min=300,max=1500,step=10,value=600),
                   div(style="margin-top:15px;",
                       span(class='help-note small',
                            paste0(format(Sys.time(),'%Y'),' NBIS • ',fn_version())
                       )
                   )
            ),
            column(9,
                   
                   tabsetPanel(id="tabset_results",
                               tabPanel("Reduction",
                                        htmlOutput("out_reduction")
                               ),
                               tabPanel("Expression",
                                        htmlOutput("out_expression")
                               ),
                               tabPanel("Markers",
                                        div(style="margin:10px;",DTOutput("out_markers"))
                               ),
                               tabPanel("Trajectory",
                                        img(src="2489360.jpg",alt="under-construction",width="300px")
                               ),
                               tabPanel("Version",
                                        fixedRow(style="padding:15px;",
                                                 column(3,
                                                        includeMarkdown("versions.md")
                                                 )
                                        )
                               )
                               
                   )
                   
                   
                   
            )
          )
)