#AminoAcidsDB
#June 20th 2023
#Authors: Pawanjit Sandhu, Ryland T. Giebelhaus, and Susan J. Murch
#PlantSMART research group at UBC Okanagan
#contact: Dr. Susan J. Murch. Email: susan.murch@ubc.ca
#Website: hormonomicsDB.com

#call shiny library in
library(shiny)
library(DT)

#####

ui <- fluidPage(

  titlePanel("AminoacidDB"),
  
  sidebarLayout(
    sidebarPanel(
      strong("AminoacidDB v1.1"),
      br(),
      p("Developed by: Pawanjit Sandhu (1), Ryland T. Giebelhaus (2), Ryan Hayward (3), Tingting Zhao (4), Alix Tucker (1), Daniel Gaudet (1), Tao Huan (4), Susan J. Murch (1)"),
      p("(1) Department of Chemistry, The University of British Columbia - Okanagan"),
      p("(2) Department of Chemistry, University of Victoria"),
      p("(3) Supra Research and Development"),
      p("(4) Department of Chemistry, The University of British Columbia - Vancouver"),
      p("Contact: rgiebelhaus@uvic.ca"),
      p("AminoacidDB is a tool developed in collaboration between UVic by the Giebelhaus Laboratory, UBC Okanagan by the plantSMART
        research team, Supra Research and Development, and UBC Vancouver by the Huan laboratory to allow users to process their untargeted metabolomics
        data to putativley identify all known amino acids."),
      p("Find our published article here:	https://doi.org/10.1039/D5AN01248A"),
      strong("Currently there are: 332,154 Amino Acids in the database."),
      br(),
      textOutput("search_count"),
      textOutput("visit_count"),
    ),
    mainPanel(
      tabsetPanel(type = "tabs",
                  
                  
                  tabPanel("Instructions",
                           br(),
                           
                           strong("Instructions"),
                           #edit this text to change the instructions
                           p("Use either the 'm/z screener' which searches against our amino acids datasets or select the
                            'Custom Database Search' to upload your own dataset use our platform to perform your
                            own custom queries of your untargeted metabolomics data. View your output results in the tab
                            next to the tool you used then download your results as a .csv file."),
                           br(),
                           strong("Database Descriptions: "),
                           p(("Lotus/HMDB: Curated dataset of 6485 compounds in classes amino acids and their derivatives from Lotus and HMDB. The dataset also contains known non-protein amino acids from plants not found in either Lotus/HMDB.")),
                           p(("PubChem: 325669 compounds from PubChem containing amino and carboxylic acid functional groups, elements C, H, N, O, and S, and number of carbons ranging from 2-11.")),
                           p("Monoisotopic: Only the monoisotopic mass of amino acids from Lotus/HMDB or PubChem datasets."),
                           p("M+H: M+H adduct of amino acids in ESI+ mode."),
                           p("M-H: M-H adduct of amino acids in ESI- mode."),
                           p("Adducts: Common adducts of amino acids in ESI+ mode."),
                           br(),

                           strong("Code Availability"),
                           br(),
                           strong("Terms and Agreements"),
                           p("AminoacidDB was developed for research use only and is not intended for use in diagnostic work.
                             Despite diligent validation and error fixing, we are not responsible for any mistakes the application
                             makes in data processing. Considering this, please inform us immediately of any issues or errors that you encounter."),
                           p("We do not save any data that is uploaded to the server, it is
                             immediately deleted with every new session that you start."),
                           br(),
                           strong("Download test data here:"),
                           br(),
                           downloadButton("downloadTestData", "Download Test Data"),
                           br(),
                           strong("Download all amino acids here as a .csv:"),
                           br(),
                           downloadButton("downloadAAs", "Download All Amino Acids"),
                           tags$hr(),),

                  tabPanel("M/Z Screener",

                           br(),

                           strong("Instructions: "),
                           p("Select which datasets to search from then select a search tolerance
                             and then upload your formatted data as a .csv and allow up to 3 minutes to perform
                             the search. After this is completed select how you want your data ordered and
                             view it in the 'Screener Output' tab."),
                           

                           #select what dataset to search against
                           checkboxGroupInput("WOS_pubChem", "Choose Dataset: ",
                                              choices = list("Lotus/HMDB" = 2,
                                                             "PubChem" = 1
                                                             ),
                                              selected = 2),
                           
                           #gives checkboxes so user can select multiple datasets to search from at once
                           checkboxGroupInput("dataset", "Choose Adducts: ",
                                              choices = list("Amino Acid Monoisotopic" = 1,
                                                             "Amino Acid M+H" = 2,
                                                             "Amino Acid M-H" = 3,
                                                             "Amino Acid Adducts (ESI + Only)" = 4
                                              ),
                                              selected = 1),
                           
                           #controls for using +/- Da or +/- ppm
                           radioButtons("tolMode", "Select mass tolerance mode (Da or PPM): ",
                                        choices = list("+/- Da" = 1,
                                                       "+/- ppm" = 2),
                                        selected = 1),
                           
                          numericInput('tol', "Mass tolerance (+/- Da or PPM)", 0.02, min = 0, max = 500, step = 0.0001), #tolerance input
                  
                  fileInput('file1', 'Choose file to upload: ', #import csv button
                            accept = c(
                              'text/csv',
                              '.csv'
                            )),
           
                  downloadButton("downloadData", "Download Results"),
                  tags$hr(),),
       
                  tabPanel("Screener Output",
                           br(),
                           DT::dataTableOutput("matches")  # or verbatimTextOutput("matches") for debug
                  ),
       
                  tabPanel("All Amino Acids",
                           br(),
                           DT::dataTableOutput("allCompsDisplay")),
                  ##custom database search thing copy pasted from hormonomicsdb
                  tabPanel("Custom Database Search",
                           br(),
                           strong("Instructions: "),
                           p("Upload your custom dataset to search from then upload your experimental data and allow the tool to run, please allow
                             up to 3 minuites to run. Then view your results and download as a .csv"),
                           fileInput('dataset_input', 'Upload custom dataset: ', #import csv button
                                     accept = c(
                                       'text/csv',
                                       '.csv'
                                     )
                           ),
                           #controls for using +/- Da or +/- ppm
                           radioButtons("tolMode2", "Select mass tolerance mode (Da or PPM): ",
                                        choices = list("+/- Da" = 1,
                                                       "+/- ppm" = 2),
                                        selected = 1),
                           numericInput('tol2', "Mass tolerance (+/- Da or ppm)", 0.02, min = 0, max = 10000, step = 0.0001), #slider bar input
                           fileInput('file2', 'Choose file to upload: ', #import csv button
                                     accept = c(
                                       'text/csv',
                                       '.csv'
                                     )
                           ),
                           downloadButton("downloadData_shell", "Download Results"),
                           tags$hr(),
                  ),
                  
                  tabPanel("Custom Database Output",
                           br(),
                           DT::dataTableOutput('contents_shell'),
                  ),
                  
       
                  )
    )
  )
)


#####



