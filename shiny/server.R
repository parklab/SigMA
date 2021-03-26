library('DT')
devtools::load_all(path = "../R")

options(shiny.maxRequestSize = 200*1024^2)

myplot <- function(x) {
  plot <- gridExtra::arrangeGrob(x)
  class(plot) <- c("myplot", class(plot))
  plot
}

print.myplot <- function(x, ...) {
  grid::grid.draw(x)
}

server <- function(input, output, session){
  urlMAF <- a("test MAF file", href = "https://github.com/parklab/SigMA/blob/master/test_mutations_50sample.maf")
  urlVCF <- a("test VCF directory", href = "https://github.com/parklab/SigMA/tree/master/inst/extdata/")
  urlVCFinfo <- a("VCF", href = "https://docs.gdc.cancer.gov/Data/File_Formats/VCF_Format/")
  urlMAFinfo <- a("MAF", href = "https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/")

  output$testFilesLink <- renderUI({
    tagList("To test SigMA you can download example datasets from the GitHub respository (", urlMAF, "or ", urlVCF, ").")
  })
  output$fileFormatInfo <- renderUI({
    tagList("You can upload your data as a single", 
            urlMAFinfo, 
            "file or by selecting a directory with multiple",  
            urlMAFinfo, " or", urlVCFinfo, " files.")
  })
  
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
          if(sum(grepl('Reference_Allele', header) & grepl('Tumor_Seq_Allele1|Tumor_Seq_Allele2', header) & grepl('Start_position|Start_Position', header) & grepl('End_position|End_Position', header) & grepl('Chromosome', header) & grepl('Position', header)) == 0){
            error_message <- 'Reference_Allele, Tumor_Seq_Allele1/2, Start_Position, End_Position, Chromosome and Position columns are required'
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
                                    file_type = input$file_type,
                                    is_list = T, 
                                    ref_genome_name = input$ref_genome)
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
      do_mva <- T
      if(sum(grepl('check_msi', input$other_settings)) > 0){
        check_msi <- T
      }
      if(sum(grepl('without_mva', input$other_settings)) > 0){
        do_mva <- F
      }
     
      if(!has_model(data = data, tumor_type = tumor_type)){
        error_message <- 'No MVA model for this tumor type for targetted gene panels select Without MVA option'
        showNotification(error_message, type = 'error',
                         action = a(href = "javascript:location.reload();", "Reload page"), duration = NULL)
        return(F)
      }

      output_file <- run(genomes_file, 
                         tumor_type = tumor_type,
                         data = data,
                         do_mva = do_mva,
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
  plot_trinuc_pos <- reactive({
    do_mva <- T
    if(sum(grepl('without_mva', input$other_settings)) > 0){
      do_mva <- F
    }
    return(myplot(plot_summary(output_file(), do_mva)$trinucleotide_pos))
  })
  plot_trinuc_neg <- reactive({
    do_mva <- T
    if(sum(grepl('without_mva', input$other_settings)) > 0){
      do_mva <- F
    }
    return(myplot(plot_summary(output_file(), do_mva)$trinucleotide_neg))
  })
  plot_c <- reactive({
    do_mva <- T
    if(sum(grepl('without_mva', input$other_settings)) > 0){
      do_mva <- F
    }
    return(myplot(plot_summary(output_file(), do_mva)$plot_cos))
  })
  plot_m <- reactive({
    do_mva <- T
    if(sum(grepl('without_mva', input$other_settings)) > 0){
      do_mva <- F
    }
    return(myplot(plot_summary(output_file(), do_mva)$plot_ml))
  })
  plot_s <- reactive({
    do_mva <- T
    if(sum(grepl('without_mva', input$other_settings)) > 0){
      do_mva <- F
    }
    return(myplot(plot_summary(output_file(), do_mva)$plot_score))
  })
  

  output$plot_trinucleotide_pos <- renderPlot({
    print(plot_trinuc_pos())
  })
  output$plot_trinucleotide_neg <- renderPlot({
    print(plot_trinuc_neg())
  })
  output$plot_cos <- renderPlot({
    print(plot_c())
  })
  output$plot_ml <- renderPlot({
    print(plot_m())
  })
  output$plot_score <- renderPlot({
    print(plot_s())
  })

  
  # make a data table 
  sorted_samples <- reactive({
    df <- read.csv(output_file())
    df_lite <- lite_df(df)
    df_lite$name <- paste0('tumor', 1:dim(df_lite)[[1]])
    do_mva <- T
    if(sum(grepl('without_mva', input$other_settings)) > 0){
      do_mva <- F
    }
    if(do_mva){
      df_lite <- df_lite[order(-df_lite$Signature_3_mva),]
      df_lite$Signature_3_mva <- round(df_lite$Signature_3_mva, digit = 4)
      df2 <- df_lite[, c('name', 'tumor', 'Signature_3_mva', 'categ')]
      colnames(df2) <- c('tag', 'filename', 'score', 'categ')
    }
    else{
      df_lite <- df_lite[order(-df_lite$Signature_3_ml),]
      df_lite$Signature_3_ml <- round(df_lite$Signature_3_ml, digit = 4)
      df2 <- df_lite[, c('name', 'tumor', 'Signature_3_ml', 'categ')]
      colnames(df2) <- c('tag', 'filename', 'score', 'categ')
    }
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

  output$sorted_samples <- DT::renderDataTable(
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
      myplot(plot_detailed(output_file(), this_sample$sample)$tribase)
  })
  plotd_cos <- reactive({
    if(this_sample$sample != '')
      myplot(plot_detailed(output_file(), this_sample$sample)$plot_cos)
  })
  plotd_ml <- reactive({
    if(this_sample$sample != '')
      myplot(plot_detailed(output_file(), this_sample$sample)$plot_ml)
  })
  plotd_exp <- reactive({
    if(this_sample$sample != '')
      myplot(plot_detailed(output_file(), this_sample$sample)$plot_exp)
  })


  output$plotd2 <- renderPlot({
    if(this_sample$sample != '')
      print(plotd())
  })
  output$plotd3 <- renderPlot({
    if(this_sample$sample != '')
      print(plotd_ml())
  })
  output$plotd4 <- renderPlot({
    if(this_sample$sample != '')
      print(plotd_exp())
  })
  output$plotd5 <- renderPlot({
    if(this_sample$sample != '')
      print(plotd_cos())
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
            dataTableOutput("sorted_samples") 
          ),
          column(7, h4("Summary") ,
            downloadButton("save_file", "Save data file"),
            plotOutput("plot_trinucleotide_pos", width = "90%", height = "200px"),
            plotOutput("plot_trinucleotide_neg", width = "90%", height = "150px"),
            plotOutput("plot_score", width = "90%", height = "200px")
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
              plotOutput("plotd2", width = "90%", height = "200px"),
              plotOutput("plotd3", width = "90%", height = "200px"),
              plotOutput("plotd4", width = "90%", height = "200px"),
              plotOutput("plotd5", width = "90%", height = "200px")
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