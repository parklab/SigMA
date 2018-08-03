library('DT')

myplot <- function(x) {
  plot <- gridExtra::arrangeGrob(x)
  class(plot) <- c("myplot", class(plot))
  plot
}

print.myplot <- function(x, ...) {
  grid::grid.draw(x)
}

server <- function(input, output, session){

  # run SigMA
  output_file <- eventReactive(input$do_run,{
    withProgress(message = 'Running', value = 0, {
      genomes_matrix <- make_matrix(directory = input$directory, 
                                    file_type = input$file_type)
      incProgress(0.3)
      genomes <- conv_snv_matrix_to_df(genomes_matrix)
      genomes_file = 'example.csv'
      write.table(genomes,
              genomes_file,
              sep = ',',
              row.names = F,
              col.names = T ,
              quote = F)
      incProgress(0.1)
      ind <- which(as.character(tissue_names) == input$tumor_type)
      tumor_type = names(tissue_names)[[ind]]
    
      ind <- which(as.character(platform_names) == input$data)
      data = names(platform_names)[[ind]]
      incProgress(0.1)
    
      output_file <- run(genomes_file, 
                         tumor_type = tumor_type,
                         data = data,
                         do_mva = T,
                         do_assign = T)
      incProgress(0.5)
    })
    return(output_file)
  })

  # make summary figure
  plot <- reactive({
    myplot(plot_summary(output_file()))
  })
  
  output$plot2 <- renderPlot({
    withProgress(message = 'Plotting', value = 0, {
     print(plot())
    })
  })

  
  # make a data table 
  sorted_samples <- reactive({
    df <- read.csv(output_file())
    df <- df[order(-df$Signature_3_mva),]
    df$name <- paste0('tumor', 1:dim(df)[[1]])
    df$Signature_3_mva <- round(df$Signature_3_mva, digit = 4)
    df2 <- df[, c('name', 'tumor', 'Signature_3_mva')]
    colnames(df2) <- c('tag', 'filename', 'score')

    buttons <- character(dim(df2)[[1]])
    
    for(i in 1:dim(df2)[[1]]){
      buttons[[i]] <- as.character(actionButton(
                                     df2$tag[[i]], 
                                     label = df2$tag[[i]], 
                                     onclick = 'Shiny.onInputChange(\"select_button\", this.id)'))
                  
                                                  
    }
    df2$tagname <- df2$tag
    df2$tag <- buttons
    return(df2)
  })

  output$sorted_samples <- renderDataTable(
                             datatable(sorted_samples()[, -which(colnames(sorted_samples()) == "tagname")], 
                                       escape = F, 
                                       rownames = F,
                                       options = list(dom = 'ftr', scrollY = T)) %>%
                             formatStyle(c("tag", "filename", "score"),
                                         backgroundColor = '#eaf1fc'))
                                     
                                         

  # add buttons in the data table so that                                       

  this_sample <- reactiveValues(sample = '')

  observeEvent(input$select_button, {
    selectedRow <- as.integer(c(unlist(strsplit(input$select_button, split = "tumor")))[[2]])
    this_sample$sample <<- sorted_samples()$filename[[selectedRow]]
  })

  plotd <- reactive({
    if(this_sample$sample != '')
      myplot(plot_detailed(output_file(), this_sample$sample))
  })

  output$plotd2 <- renderPlot({
    if(this_sample$sample != '')
      withProgress(message = 'Plotting', value = 0, {
        print(plotd())
        incProgress(1)
      })
  })


  observeEvent(input$do_run, {  
   removeTab(inputId = "tabs", target = "Signature 3")

   insertTab(inputId = "tabs",
      tabPanel(
        "Signature 3", 
        fluidRow(style = "padding:15px;
                             margin-left:15px;
                             margin-right:15px;
                             margin-top:15px;
                             margin-bottom:15px;",
             
          column(5, h4("Summary") ,
            plotOutput("plot2", width = "90%", height = "750px")
          ),
          column(7, style = "background-color: #eaf1fc",
            h4("Click on the tag to see more"),
            dataTableOutput("sorted_samples")
          )
        )
      ),
      target = "Info"
    )
  })

  observeEvent(input$do_run, {
    updateTabsetPanel(session, "tabs", selected = "Signature 3")
  })

  observeEvent(input$select_button, {
    removeTab(inputId = "tabs", target = "selectedTumor")

    insertTab(inputId = "tabs",
      tabPanel(
        "selectedTumor",
        fluidRow(style = "padding:15px;
                             margin-left:15px;
                             margin-right:15px;
                             margin-top:15px;
                             margin-bottom:15px;",
          column(6, 
            fluidRow(plotOutput("plotd2", width = "90%", height = "400px"))
          )
        )
      ),
      target = "Info"
    )
  })

  observeEvent(input$select_button, {
    updateTabsetPanel(session, "tabs", 
      selected = "selectedTumor")
  })

}