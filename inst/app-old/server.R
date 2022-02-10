# SERVER

library(readxl)

shinyServer(function(input,output,session) {

  #store <- reactiveValues(epath=tempdir())
  
  output$ui_cluster <- renderUI({
    if(input$in_group=="Cluster") selectInput("in_cluster","Number of clusters",choices=c(6,10,16),selected=10,multiple=FALSE)
  })
  
  fn_data <- reactive({
    celltype <- tolower(input$in_celltype)
    condition <- tolower(input$in_condition)
    group <- tolower(input$in_group)
    
    if(group=="cluster" && (!is.null(input$in_cluster))){
      fname <- paste(c(celltype,condition,group,input$in_cluster),collapse="-")
    }else{
      fname <- paste(c(celltype,condition,group),collapse="-")
    }
    
    print(fname)
    return(fname)
  })
  
  output$out_reduction <- renderUI({
    prefix <- fn_data()
    path <- paste0("images/reduction/",tolower(input$in_celltype),"/",prefix,".png")

    if(file.exists(file.path("www",path))){
      div(class="image-display",
          img(src=path,alt=prefix,width=input$in_width)
      )
    }else{
      div(class="image-display",
        p("Image not available.")
      )
    }
  })
  
  output$out_expression <- renderUI({
    prefix <- fn_data()
    path <- paste0("images/expression/",tolower(input$in_celltype),"/",prefix,".png")
    
    if(file.exists(file.path("www",path))){
      div(class="image-display",
          img(src=path,alt=prefix,width=input$in_width)
      )
    }else{
      div(class="image-display",
        p("Image not available.")
      )
    }
  })
  
  output$out_markers <- renderDT({
    path <- paste0("images/markers/",tolower(input$in_celltype),"/",fn_data(),".xlsx")
    
    if(file.exists(file.path("www",path))){
      readxl::read_xlsx(file.path("www",path))
      #as.data.frame(lapply(m,function(x) ifelse(is.numeric(x),round(x,3),x)))
    }else{
      data.frame()
    }
  })
})
