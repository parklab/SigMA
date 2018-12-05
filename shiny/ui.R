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
            h4('Load data'),
            tags$div(class="form-group shiny-input-container", 
                     tags$div(tags$label("Select a directory")),
                     tags$div(tags$label("Browse", class="btn btn-default btn-file",
                              tags$input(id = "directory2", 
                                         webkitdirectory = TRUE, 
                                         type = "file", 
                                         style="display: none;", 
                                         onchange="pressed()"))
                              ),
                    tags$div(id = paste("directory2",  "_progress", sep = ""), 
                            class = "progress progress-striped active shiny-file-input-progress", 
                            tags$div(class = "progress-bar"))
                     ),
            fileInput("directory1", label = "OR Select a file"),
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
            checkboxGroupInput("other_settings", "Options:",
                               choiceNames = c("Check for microsattelite instability",
                                               "Lite data format"),
                               choiceValues = c("check_msi", "lite_format")),
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
  ),

  # generates a modal bux upon clicking run that informs the use
  # that calculation is in progress
  tags$script("Shiny.addCustomMessageHandler('launch-modal', function(d) {$('#' + d).modal().focus();})"),
  tags$script("Shiny.addCustomMessageHandler('remove-modal', function(d) {$('#' + d).modal('hide');})"),
  tags$div(
    id = "modal_inprogress",
    class="modal fade", tabindex="-1", `data-backdrop`="static", `data-keyboard`="false",
    tags$div(
      class="modal-dialog",
      tags$div(
        class = "modal-content",
        tags$div(class="modal-header", tags$h4(class="modal-title", "Calculation in progress, please wait.")),
        tags$div(class="modal-footer", tags$button(type="button", class="btn btn-default", `data-dismiss`="modal", "Dismiss"))
      )
    )
  )
)
