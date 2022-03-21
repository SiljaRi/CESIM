ui <- navbarPage("CESIM",
                 selected = "Home",
                 # Home_UI ####
                 tabPanel("Home",
                          fluidRow(
                              column(5,
                                     h3("Short description"),
                                     h5("Welcome to my simulator for Longitudinal Clonal Evolution Sequencing Data .
                                    With this simulator you are able to create a Clonal Evolution in detail.
                                    Choose basic settings first to get an evolution to start with. This evolution is visualized
                                    in a fishplot and a rhandsontable shows the distinct cellfractions of the clones (cluster).
                                    Afterwards you can choose more advanced settings like adding a clone or mutations
                                    to modify the evolution. The simulator 
                                    will create Single Nucleotide Variants (SNV) with variant and reference reads.
                                    If selected, also Copy Number Variations (CNV) are created with their Cancer Cell Fraction (CCF).
                                    Both variations are considered as mutations. This data can then be downloaded for use as input
                                    into several clonal reconstruction tools."),
                                     h3("Step-by-step use"),
                                     h4("Step 1: Choose Basic Settings"),
                                     h5("Go to tab Settings and choose values for all basic settings. After that click \"Run new Simulation\" 
                                   to get the simulator started. A description of every Setting can be found at More/Help."),
                                     h4("Step 2: Review the created evolution"),
                                     h5("When the simulation finishes, a fishplot and the corresponding table is shown inside the settings tab. Of note, sometimes the fishplot
                                   does not correctly show the underlying phylogeny (e.g. clone is plotted inside wrong clone at time of appearance). Be sure to doublecheck this with the phylogeny table.
                                   You can have a look at the created mutations inside the Data tab. The big table is showing all mutations (SNV and CNV together). Below you find 
                                   a seperate table of CNVs, where you can choose advanced settings for these e.g. if they should overlap other SNVs. At the very bottom of the Data tab you find the purity table
                                   which can also be altered.
                                   Now you can start editing this evolution with advanced settings or you can again click on \"Run new Simulation\",
                                   to get a different evolution. Attention: With a click on \"Run new Simulation\" all advanced settings will be 
                                   discarded!! When you found an evoultion you aimed at, you can already download the data (see step 4) or make advanced changes
                                    on the evolution (see step 3)."),
                                     h4("Step 3: Advanced Settings (optional)"),
                                     h5("If you wish to further edit the details of the evolution, you can do so with the so-called Advanced Settings.
                                   Of note, all Basic Settings can be changed, but not directly. Only the values of Mean coverage, Minimum Detection and Minimum
                                      Distance will be read when \"Update Simulation\" is pressed. All other Basic Settings will be ignored. These can be altered
                                      though changing the phylogeny, cnv or purity table. These will be desribed in More/Help. By clicking on  \"Update Simulation\"
                                      the values of the mentioned tables are read and checked. If everything is okay, new SNVs are created. This is the only data
                                      that is automatically created and not changeable. If an error occurs while checking all parameters, an error message is shown 
                                      and the data is not updated! If you download the data now it will be inconsistent. So be sure that all errors are fixed (Update successful)
                                      before downloading it."),
                                     h4("Step 4: Download the data"),
                                     h5("If you are done editing the evolution you can download the data by clicking on \"Export Data\" inside the Data tab.
                                   This should give you a zip-folder inside you Download directory containing the following components: One directory with RDS files
                                   for restore this session in the simulator for further editing (see Load session), an Excel file with all information (incl. plot), and several tsv-files 
                                   with the important tables for further use with reconstruction tools."),
                                     h4("Load session"),
                                     h5("If you already downloaded data of an earlier session you can restore this by clicking on \"Browse...\" which prompts
                                      you to choose 6 RDS files that were downloaded togehter with the data. If you choose these 6 and click open, the evoultion and
                                      all tables will be restored and you can edit it again. Note, that not all Basic Settings are loaded correctly. Better check those manually.")
                              ),
                              column(2,
                              ),
                              column(3,
                                     br(), br(), br(), br(), br(), br(), br(),
                                     h4("Linear evolution"),
                                     imageOutput("linearpng", height = 150),
                                     h4("Parallel-dependant evolution"),
                                     imageOutput("parallelpng", height= 150),
                                     h4("Independent evolution"),
                                     imageOutput("independentpng")
                              )
                          )
                 ),
                 #Settings_UI ####
                 tabPanel("Settings", 
                          fluidRow(
                              plotOutput("fishPlot"),
                              hr(),
                              column(2,
                                     h4("Basic Settings"),
                                     sliderInput('samplevalue', 'Timepoints / Samples', 
                                                 min=1, max=12, value=2, 
                                                 step=1, width = "80%", ticks = FALSE),
                                     selectInput('evolutionvalue', 'Kind of Evolution', 
                                                 c("linear","parallel-dependent","independent","random"),width = "80%",selected = "random"),
                                     selectInput('clonecount', 'Number of Clones', 1:25, width = "80%", selected = 5),
                                     selectInput('mutationcount', 'Number of Mutations (SNV+CNV)', 
                                                 1:200, selected =10, width = "80%"),
                                     sliderInput("cnvcount", "Number of CNVs", min=0, max=10, value=2, 
                                                 step=1, width = "80%", ticks = FALSE),
                                     checkboxGroupInput("cnvcheckbox", "Explicit", choices = c("Duplication", "Deletion", "LOH"), selected = c("Duplication","Deletion"))
                                     #checkboxInput('fp', 'Include false positives', width = "100%")
                              ),
                              column(2, br(), 
                                     br(), br(),
                                     actionButton("simulatebtn","Run new Simulation"),
                                     br(), br(),
                                     fileInput("upload", label=("Load session (6 RDS files)"), accept = ".RDS", multiple = TRUE, width = "70%"),
                                     # actionButton("loadbtn", "Load Session"),
                                     br(), 
                                     
                                     numericInput('coveragevalue', 'Mean coverage', min = 1, max = 1000,
                                                  value = 300, width = "60%"),
                                     sliderInput("purityvalue", "Purity", min = 1, max = 100, value = 100,
                                                 width = "60%",ticks = FALSE),
                                     numericInput('detectionvalue', 'Minimum Detection', min = 0, max = 10,
                                                  value = 2, width = "60%"),
                                     numericInput('distancevalue', 'Minimum Distance', min = 0, max = 10,
                                                  value = 4, width = "60%")
                              ),
                              column(7,
                                     h4("Advanced Settings"),
                                     br(),
                                     rHandsontableOutput("clonetable"), 
                                     br(),
                                     actionButton("updatebtn1","Update Simulation"),
                                     br(), br(),
                                     uiOutput("updatetext1"),
                                     br(), br(),
                                     fluidRow(
                                         column(2, br(), 
                                                uiOutput("uiaddsample"),
                                                br(),
                                                uiOutput("uiaddclone"),
                                                br() 
                                         ),
                                         column(3, #offset = 1,
                                                uiOutput("uideletesamplevalue"),
                                                uiOutput("uideletesample")
                                                
                                         ),
                                         column(3,
                                                uiOutput("uideleteclonevalue"),
                                                uiOutput("uideleteclone")
                                         )
                                     )
                              )
                          )
                 ),
                 # Data_UI ####
                 tabPanel("Data",  
                          h3("Mutational Data"),
                          dataTableOutput("mutationtable"), 
                          hr(),
                          fluidRow(
                              column(9,  
                                     uiOutput("uicnvtitle"), br(), 
                                     rHandsontableOutput("cnvtable"),
                                     br(),
                                     fluidRow(
                                         column(1, br(), 
                                                actionButton("addcnvbtn", "Add CNV")
                                         ),
                                         column(2,
                                                numericInput("deletecnvvalue", "Choose CNV to delete", value = 1, min = 1, width = "100%"),
                                                actionButton("deletecnvbtn", "Delete CNV")
                                         ),
                                         column(8, offset = 1,
                                                actionButton("updatebtn2", "Update Simulation"),
                                                br(), br(),
                                                uiOutput("updatetext2")
                                         )
                                     )
                              ),
                              column(3,  br(), uiOutput("uiexport2"), br(), br()
                              )
                          ),
                          fluidRow(
                              column(4, 
                                     uiOutput("uipuritytitle"), br(), 
                                     rHandsontableOutput("puritytable"),
                                     br()
                              )
                          )
                 ),
                 # More_UI ####
                 navbarMenu("More",
                            tabPanel("Help",
                                     fluidRow(
                                         column(5, 
                                                h3("Basic Settings"),
                                                br(),
                                                h4("Timepoints/Samples"),
                                                h5("Count of smples. Since we focus on longitudinal data this means samples from the same tumor at different timepoints. Default: 2-3."),
                                                br(),
                                                h4("Kinf of evolution"),
                                                h5("Choose one of three different kinds of evolution: Linear evolution - all subclones descent linear, parallel-dependent - subclones
                                               can develop parallel but have a common ancestor clone, independent - there are several clones that do not have a common ancestor.
                                               See the pictures on Home tab for a visualization. Additionally, you can choose random if you do not aim at a certain evolution."),
                                                br(),
                                                h4("Number of Clones"),
                                                h5("Ranges from 1-25. For a greater number of clones it can sometimes be impossible to create an evolution automatically. You can try it several times,
                                               the simulation of clustermeans is random and therefore can succeed by chance. This value can later be overwritten, see Phylogeny Table. Default: 4-10."),
                                                br(),
                                                h4("Number of Mutations"),
                                                h5("This number gives the count of all mutations, Single Nucleotide Variants (SNVs) and Copy Number Variants(CNVs), together. The slider 'Number of CNVs' adapts according to this input. There will be 
                                                 |Number of Mutations| - |Number of CNVs| SNVs in total. It is possible to have only SNVs, only CNVs or both. This value can be overwritten, see Phylogeny Table. Default: 10-100"),
                                                br(),
                                                h4("Number of CNVs"),
                                                h5("This slider regulates the count of CNVs as part of all mutations. Therefore it is limited by 'Number of Mutations'. This value can be overwritten, see CNV Table. Default: 2 (1 Del + 1 Dup)"),
                                                br(),
                                                h4("Explicit"),
                                                h5("Here you can choose different Types for CNVs. These types can be easily edited, see CNV Table. "),
                                                br(),
                                                h4("Run new Simulation"),
                                                h5("By clicking this button, a new simulation is run taking all basic settings into account. Advanced settings / changes are ignored.
                                                 First a phylogeny is created. If a suitable phylogeny is found, next CNVs are created and afterwards the remaining mutations (SNVs)."),
                                                br(),
                                                h4("Load session"),
                                                h5("Possibility to load a past session. Choose all downloaded 6 RDS files and the session will be restored."),
                                                br(),
                                                h4("Mean Coverage"),
                                                h5("Select a mean coverage that is used to produce the reads for the SNVs. This field can be changed and will be taken into account when updating a simulation. Default: 100-1000"),
                                                br(),
                                                h4("Purity"),
                                                h5("Select a purity value (1-100, = tumorcontent) for all samples. This can be specified later, see Purity Table. Contamination with normal cells is 100 - purityvalue. Default: 100"),
                                                br(),
                                                h4("Minimum Detection"),
                                                h5("To address the issue of sequencing errors, a mutation or clone has to have a minimum CCF of 'Minimum Detection'. Values below this limit are not allowed and will automatically be corrected
                                                 to 0. Mutations with a lower CCF are considered artifacts or sequencing errors and would be discarded in real data. Default: 2 - 5 (CCF)"),
                                                br(),
                                                h4("Minimum Distance"),
                                                h5("This value describes the distance between every pair of clones/cluster. Each clone has a CCF distance of 'Minimum Distance' in at least one sample to each other clone. This value is for 
                                                 creating unique clusters. If set to 0, two similar clones (with identical clustermeans) can exist. These clones will mostlikely be clustered into one clone when using reconstruction tools. Default: 4"),
                                         ),
                                         column(1,
                                         ),
                                         column(5,
                                                h3("Advanced Settings"),
                                                h5("Advanced Settings in this simulator are changes in the RHandsontables which will be read and processed when clicking update."),
                                                br(),
                                                h4("Phylogeny Table"),
                                                h5("The Phylogeny Table is shown in tab 'Settings'. Every row represents one clone together with its parent clone (predecessor), how many mutations are assigned to it and its CCF in all samples.
                                                 To simplify the clone assignment to the fishplot, the same colors are used. Notice that the CCF values here will be shown in relative range from 1-100 for clarity whereas every other table 
                                                 shows the absolute range from 0-1. The Phylogeny Table is an RHandsontable, so every value inside can be altered by clicking into the cell and entering a value. In this way the count of mutations, samples, the kind of
                                                 evolution and the CCFs of the clones can be costumized.
                                                 To use the altered values for simulation, click on 'Update Simulation\' Every value will be checked and a message telling if an error occured or if the update was successful will appear below the
                                                 Update Button. Below this, you find buttons to add samples and clones. Clicking on these will create a new column or row, which values have to be filled in manually. You can also delete a sample or clone by
                                                 providing its number (clone = row, sample = samplenumber(tX)), and click on 'Delete sample'/'Delete clone'. If this is executable without errors, an update will happen automatically. To be sure that your data is consistent
                                                 in all tables, click 'Update Simulation' and wait for 'Update successful' message. Note that the basic settings might not fit to your altered evolution anymore as you can overwrite samples, clones, mutations and
                                                 kind of evolution by editing the RHandsontable."),
                                                br(),
                                                h4("Update Simulation"),
                                                h5("By clicking this button, the RHandsontables will be read and checked for valid values. Then the SNV will be computed automatically. On success, the message 'Update suceessful' will be shown.
                                              Otherwise a Warning or Error message will tell you what has gone wrong.."),
                                                br(),
                                                h4("Mutation Table"),
                                                h5("This table shows all mutations, SNVs and CNVs together. For each mutation its ID, position (Chr, Start, End), genotype, id of overlapping CNV (-1 if none), id of assigned clone, as well as
                                              reference- and variant reads and its corresponding VAF (for SNVs only), as well as its CCF for every sample is listed. SNVs are created automatically and can not be costumized manually. If you wish for different
                                                 SNV data, you can hit 'Update Simulation' again. They will be created randomized. CNVs can be costumized in the CNV Table."),
                                                br(),
                                                h4("CNV Table"),
                                                h5("Similar to the Phylogeny Table the CNV Table is also an RHandsontable so every entry can be edited manually. The CNV table shows every CNV together with its fixed ID, position (Chr, Start, End),
                                                 its type (Deletion, Duplication, LOH), if it is overlapping SNVs (Y/N), how many SNVs are overlapped, in which way they overlap the SNV (see CNV Function), which clone the CNV is assigned to and its corresponding
                                                 cellfractions. With the buttons you can add and delete a CNV just like clones/samples. CNVs will also be shown in the mutation table at the top of the page. For simplicity, CNVs are not allowed to overlap each other here. "),
                                                br(),
                                                h4("CNV Function"),
                                                h5("There are 4 different ways, how a CNV can affect a SNV. These types can be described as follows:"),
                                                h5("CNV first: The CNV occured temporally before the SNV. Therefore CCF_CNV >= CCF_SNV. The SNV has the genotype B (Deletion), AAB (Duplication), AB (LOH)."),
                                                h5("SNV first (active): The SNV occured temporally before the CNV. Therefore CCF_SNV >= CCF_CNV. Active means, that the allele holding the SNV is affected by the CNV. The SNV has the genotype AB(d) (Deletion),
                                              ABB (Duplication) and AB(l) (LOH). If CCF_CNV < 100 there might be \"normal\" cells containing the SNV (genotype AB) as well. "),
                                                h5("SNV first (passive): The SNV occured temporally before the CNV. Therefore CCF_SNV >= CCF_CNV. Passive means, that the allele holding the SNV is not affected by the CNV. The SNV has the genotype B (Deletion),
                                              AAB (Duplication) and BB (LOH). If CCF_CNV < 100 there might be \"normal\" cells containing the SNV (genotype AB) as well. "),
                                                h5("Parallel: The SNV and CNV exist parallel. Therefore CCF_SNV + CCF_CNV <= 1. The SNV has the genotype AB(p) (Deletion, Duplication and LOH)."),
                                                h5("If a CNV is chosen to overlap a SNV, the reads of the SNV will be corrected for the overlapping CNV e.g. a SNV with Genotype B has VAF = CCF instead of VAF*2 = CCF."), #TODO: Formel aus der Publikation einbringen, sobald diese veröffentlicht ist.
                                                br(),
                                                h4("Purity Table"),
                                                h5("This table just gives the purity for each sample. It is an RHandsontable so you can costumize the values (range 1-100). This will overwrite the purity value in the Basic Settings. To
                                                 update your data for the new purity values, hit 'Update Simulation'."),
                                                br(),
                                         )
                                     )
                            ),
                            tabPanel("About",
                                     fluidRow(
                                         column(5,
                                                h3("About"),
                                                h5("This simulator is a part of the Master's Thesis \"Reconstruction of Clonal Evolution: 
                                              Evaluation of Available Tools and Development of a Simulator for Generating Realistic Mutational Data\"."),
                                                h5("It was developed by Silja Richter, student of Computer Science at the Westfälische Wilhelms-Universität Münster."),
                                                h5("The thesis was written in cooperation with the Institut für Medizinische Informatik Münster."),
                                                br(),
                                                h5("For visualizing the clonal evolution the fishplot package is used:"),
                                                h5("Visualizing tumor evolution with the fishplot package for R. Miller CA, McMichael J, Dang HX, Maher CA, Ding L, Ley TJ, Mardis ER, Wilson RK. BMC Genomics. doi:10.1186/s12864-016-3195-z"),
                                                br(),
                                                h5("Contact: siljarichter@wwu.de"),
                                                br(),
                                                h5("Version 1.0.0"),
                                                br(),
                                         ),
                                         column(3,
                                                imageOutput("wwulogo", height= "10%"),
                                                imageOutput("imilogo")
                                         )
                                     )
                            )
                 ),
                 #for showing simulation message in better position
                 tags$head(
                     tags$style(
                         HTML(".shiny-notification {
                       position:fixed;
                       top: calc(40%);
                       left: calc(30%);
                       width: 30em;
                       }
                       "
                         )
                     )
                 ),
                 #for show and hide()
                 useShinyjs()
)