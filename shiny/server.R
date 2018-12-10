library('DT')
library(shinycssloaders)

myplot <- function(x) {
  plot <- gridExtra::arrangeGrob(x)
  class(plot) <- c("myplot", class(plot))
  plot
}

print.myplot <- function(x, ...) {
  grid::grid.draw(x)
}

server <- function(input, output, session){

  file_memory <- reactiveValues(
    name = '')

    
  # run SigMA
  output_file <- eventReactive(input$do_run,{
    if(recalculate$val){
      if(!is.null(input$directory1)){
        directory <- input$directory1$datapath
      }
      else if(!is.null(input$directory2)){ # & input$directory2 != ''){
        if(sum(grepl('maf', input$directory2$datapath[[1]])) != length(input$directory2$datapath) & sum(grepl('vcf', input$directory2$datapath)) != length(input$directory2$datapath)){
          error_message <- 'input directory/file can only consist of only maf or vcf files'
          showNotification(error_message, type = 'error',
                          action = a(href = "javascript:location.reload();", "Reload page"), duration = NULL)
          return('')
        }
        else if(sum(grepl('0.vcf', input$directory2$datapath[[1]])) > 0){
          directory <- input$directory2$datapath
        }
        else if(sum(grepl('0.maf', input$directory2$datapath[[1]])) > 0){
          directory <- input$directory2$datapath
        }
        else{ 
          error_message <- 'input directory is invalid'
          showNotification(error_message, type = 'error',
                          action = a(href = "javascript:location.reload();", "Reload page"), duration = NULL)
          return('')
        }
      }
      else{
        error_message <- 'input directory/file is not valid'
        showNotification(error_message, type = 'error',
                         action = a(href = "javascript:location.reload();", "Reload page"), duration = NULL)

        return('')
      }

      # check whether file format is maf or vcf
      check_input <- function(file){
        if(input$file_type == "maf"){
          header <- readLines(file, n = 1)
          if(!(grepl('Reference_Allele', header) & grepl('Tumor_Seq_Allele1|Tumor_Seq_Allele2', header) & grepl('Start_position|Start_Position', header) & grepl('End_position|End_Position', header) & grepl('Chromosome', header))){
            error_message <- 'Not a maf file'
            showNotification(error_message, type = 'error',
                             action = a(href = "javascript:location.reload();", "Reload page"), duration = NULL)
            return(F)
          }else
            return(T)
        }
        if(input$file_type == "vcf"){
          if(sum(grepl('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO', readLines(file, n = 500))) == 0){ 
            error_message <- 'Not a vcf file'
            showNotification(error_message, type = 'error',
                             action = a(href = "javascript:location.reload();", "Reload page"), duration = NULL)
            return(F)
          }else
            return(T)
        }
      }

      for(file in directory){
        if(!check_input(file)) return('')
      }

      session$sendCustomMessage(type = 'launch-modal', "modal_inprogress") # launch the modal      
      genomes_matrix <- make_matrix(directory = directory, 
                                    file_type = input$file_type)
      genomes <- conv_snv_matrix_to_df(genomes_matrix)
      genomes_file = 'example.csv'
      write.table(genomes,
              genomes_file,
              sep = ',',
              row.names = F,
              col.names = T ,
              quote = F)
      ind <- which(as.character(tissue_names) == input$tumor_type)
      tumor_type = names(tissue_names)[[ind]]
    
      ind <- which(as.character(platform_names) == input$data)
      data = names(platform_names)[[ind]]
    
      check_msi <- F
      lite_format <- F
      if(sum(grepl('check_msi', input$other_settings)) > 0){
        check_msi <- T
      }

      output_file <- run(genomes_file, 
                         tumor_type = tumor_type,
                         data = data,
                         do_mva = T,
                         do_assign = T, 
                         check_msi = check_msi)
      print('Finished running SigMA')
      file_memory$name <<- output_file
      session$sendCustomMessage(type = 'remove-modal', "modal_inprogress") # hide the modal programmatically
      return(output_file)
    }else{
      return(file_memory$name)
    }
  })
  
  recalculate <- reactiveValues(val = T)

  output$save_file <- downloadHandler(
    filename = 'SigMA_output.csv',
    content = function(file){
      if(sum(grepl('lite_format', input$other_settings)) > 0){
        df <- read.csv(output_file())
        df_output <- lite_df(df)
      }else{
        df_output <- read.csv(output_file())
      }
      write.table(df_output, 
                  file,
                  row.names = F, 
                  sep = ',',
                  quote = F)
    }
  ) 

  # make summary figure
  plot <- reactive({
      return(myplot(plot_summary(output_file())))
  })
  
  output$plot2 <- renderPlot({
    print(plot())
  })

  
  # make a data table 
  sorted_samples <- reactive({
    df <- read.csv(output_file())
    df_lite <- lite_df(df)
    df_lite$name <- paste0('tumor', 1:dim(df_lite)[[1]])
    df_lite <- df_lite[order(-df_lite$Signature_3_mva),]
    df_lite$Signature_3_mva <- round(df_lite$Signature_3_mva, digit = 4)
    df2 <- df_lite[, c('name', 'tumor', 'Signature_3_mva', 'categ')]
    colnames(df2) <- c('tag', 'filename', 'score', 'categ')

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
                             datatable(sorted_samples()[, -grep("tagname|filename", colnames(sorted_samples()))], 
                                       escape = F, 
                                       rownames = F,
                                       options = list(dom = 'ftrp', scrollY = T)) %>%
                             formatStyle(c("tag", "score", "categ"),
                                         backgroundColor = '#eaf1fc'))
                                     
                                         

  # add buttons in the data table so that                                       

  this_sample <- reactiveValues(sample = '')

  observeEvent(input$select_button, {
    ind <- which(sorted_samples()$tagname == as.character(input$select_button))
    this_sample$sample <<- sorted_samples()$filename[[ind]]
  })

  plotd <- reactive({
    if(this_sample$sample != '')
      myplot(plot_detailed(output_file(), this_sample$sample))
  })

  output$plotd2 <- renderPlot({
    if(this_sample$sample != '')
      print(plotd())
  })


  dataset <- reactiveValues(number = 0)

  observeEvent(input$do_run, {
    # cleaning tabs

    recalculate$val <<- T
    output_file()    
    recalculate$val <<- F

    removeTab(inputId = "tabs", target = "Results")
    dataset$number <<- dataset$number + 1
    removeTab(inputId = "tabs", target = "Details")   

    insertTab(inputId = "tabs",
      tabPanel(
        "Results",
        fluidRow(style = "padding:15px;
                             margin-left:15px;
                             margin-right:15px;
                             margin-top:15px;
                             margin-bottom:15px;",
             
          column(5, style = "background-color: #eaf1fc",
            h4("Click to see sample-specific information"),
            dataTableOutput("sorted_samples") %>% withSpinner(color="#dee9fc")
          ),
          column(7, h4("Summary") ,
            downloadButton("save_file", "Save data file"),
            plotOutput("plot2", width = "90%", height = "750px") %>% withSpinner(color="#dee9fc")
          )
        )
      ),
      target = "Info"
    )
  })

  observeEvent(input$do_run, {
    updateTabsetPanel(session, "tabs", selected = "Results")
  })

  observeEvent(input$select_button, {
    removeTab(inputId = "tabs", target = "Details")

    insertTab(inputId = "tabs",
      tabPanel(
        "Details",
        fluidRow(style = "padding:15px;
                             margin-left:15px;
                             margin-right:15px;
                             margin-top:15px;
                             margin-bottom:15px;",
          column(6, 
            fluidRow(
              plotOutput("plotd2", width = "90%", height = "400px") %>% withSpinner(color="#dee9fc")
            )
          )
        )
      ),
      target = "Info"
    )
  })

  observeEvent(input$select_button, {
    updateTabsetPanel(session, "tabs", 
      selected = "Details")
  })

}