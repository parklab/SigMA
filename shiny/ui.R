load("../R/sysdata.rda")
 
ui <- fluidPage(
  title = 'SigMA',
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
        "Home",
        fluidRow(
          column(8, align = 'center', style = "background-color:#eaf1fc;
                             padding:30px;margin-top:
                             10px;margin-bottom:10px;",
                 tags$div(tags$h3('Description of the algorithm'),
                          tags$br(),
                          tags$p('SigMA (Signature Multivariate Analysis): We developed SigMA to detect mutational signatures from the SNV calls of whole-genome, exome or targeted gene panel data. A detailed description of the algorithm and its performance is provided in the related manuscript.'),
                          tags$p('In brief, SigMA consists of 5 main steps. First, mutational signatures in WGS data are discovered using NMF. Second, the tumor subtypes based on their signature composition are determined with clustering, and used as a reference for panels. Third, we simulate cancer-gene panels and exomes from the WGS data. In our simulations, the labels (whether a tumor is true Signature 3-positive or -negative) are known based on the signature analysis in WGS data. In the fourth step, the likelihood measure, cosine similarity and exposure of Signature 3 with NNLS are calculated for simulated panels, exomes and WGS data. Finally, we train Gradient Boosting Classifiers (GBCs) specific for each tumor type, and sequencing platform, using the features from step 4. The GBCs yield a final combined score. We determine the thresholds on SigMA score, which corresponds to small false positive rates, using the simulated data and the true labels from the WGS analysis. The thresholds depend on tumor type and on the platform.')),
                         tags$br(),
                         img(height = 1.5*373.973, width = 1.5*366.396 , src = 'workflow.png')
          ),
          column(4, align = "center", style = "background-color:#eaf1fc;
                             padding:30px;margin-top:
                             10px;margin-bottom:10px;",
                 h3('Instructions for use'),
                 tags$div(style = "display: inline-block;
                                   vertical-align:top;",
                         img(height = 3.5*189.513, width = 3.5*85.756, src = 'howtorun.png'),
                         tags$br(),
                         tags$br(),
                         tags$p('Go to the "Get Started" tab. Browse and select a file or directory as described on the figure to upload your data.'), 
                         uiOutput("fileFormatInfo"),
                         uiOutput("testFilesLink")
                         )
                )
        )
      ),
      tabPanel(style = "background-color:#eaf1fc",
        "Get Started",       
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
