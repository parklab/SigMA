load("../R/sysdata.rda")
 
ui <- fluidPage(
  tags$head(tags$style(
    HTML('
        body,canvas,figure { 
            font-family: "Helvetica";
            background-color: #FFF;
        }')
    )),
  tags$style(HTML("
    .tabbable > .nav > li[class=active] > a[data-value = 'Input'] {
    background-color: #eaf1fc;
    color: #000000;
  }")),
  titlePanel(
    fluidRow(
      div(style="display: inline-block;
                 vertical-align:top; 
                 width: 150px; 
                 margin-left:15px;",
          img(height = 70, width = 120, src = 'sigma.jpg'))
    )
  ),
  # main panel
  fluidRow(class= "R1", 
    tabsetPanel(id = "tabs", #type= "tabs",
      tabPanel(style = "background-color:#eaf1fc",
        "Input",       
        fluidRow(
          column(4),
          column(4, style = "background-color:#dee9fc;
                             padding:15px;margin-top:
                             15px;margin-bottom:15px;",
            h4('Input directory/Input file and format'),
            fileInput("directory1", label = "Browse and select file"),
            textInput("directory", label = "OR Enter directory path"),
            selectInput(
              inputId = "file_type",
              label = "File format",
              choices = c("vcf", "maf")
            ),
            hr(),
            h4('Settings'),
            selectInput(
              inputId = "tumor_type",
              label = "Select tumor type",
              choices = sort(as.character(tissue_names))
            ),
            selectInput(
              inputId = "data",
              label = "Select sequencing platform",
              choices = as.character(platform_names)
            ),
            actionButton("do_run", "Run")
          ),
          align = "center"
        )
      ),
      tabPanel("Info", style = "padding:15px;
                                 margin-top:15px;
                                 margin-bottom:15px;",
               div(p(HTML(paste0("This software was produced by ",
                               a(href = 'https://compbio.hms.harvard.edu/index', 'Park Lab')))),        
                   p(HTML(paste0('The code is available on ', 
                                 a(href = 'https://github.com/parklab/SigMA', 
                                 'github')))))
      )
    )
    
#    align = "center"         
  )
)
