# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# PlotsOfPhasors - A Shiny app for plotting (frequency domain) FLIM data
# Created by Joachim Goedhart (@joachimgoedhart) & Franka van der Linden, first version 2020
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Joachim Goedhart & Franka van der Linden (C) 2020
# electronic mail address: j #dot# goedhart #at# uva #dot# nl
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# To implement:
# Filter/Remove conditions
# Averages, show table with stats



library(shiny)
library(ggplot2)
library(magrittr)
library(dplyr)
# library(ggrepel)
library(shinycssloaders)
library(readr)
library(DT)
library(patchwork)
# library(RCurl)

# Load Example data
df_example <- read.csv("tau.csv", na.strings = "")
df_example2 <- read.csv("lifetime_data.csv", na.strings = "")
df_tau <- data.frame(tphi=c(0, 0.2, 0.5, 1, 2,3,4, 5,6,8, 10, 20, 100),
                     tmod=c(0, 0.2, 0.5, 1, 2,3,4, 5,6,8, 10, 20, 100)
                    )

######## Functions ##########
#Convert lifetimes to polar coordinates
lifetime_to_GS <- function(data,MHz=40) {
  tphi <- data$tphi
  tmod <- data$tmod
  omega = 2*pi*MHz/1000
  Phi = atan(omega*tphi)
  M=sqrt(1/ (1 + (omega*tmod)^2 ) )
  G = M * cos(Phi)
  S = M * sin(Phi)
  new <- data.frame(tphi=tphi, tmod=tmod, Phi=Phi, M=M, G=G, S=S)
  return(new)
}

#define a function for a circle
make_half_circle <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(G = xx, S = yy))
}

#Function to generate a pulse profile that can be used to convolve with the decay
pulse_profile <- function(shape="none", centre=25, width=11, len=256) {
  pulse <- seq(0,0,length.out=len)
  #The R index starts at 1, so we need to add 1 here to match the index with the bins
  #This seems trivial but it's important to keep track of this
  i_centre <- centre + 1
  if (i_centre < 1) {centre <- 0}
  if (width*2 > len) {width <- len/2}
  
  if (shape=="none") {
    pulse[i_centre] <- 1
  }
  if (shape=="block") {
    pulse[ceiling(i_centre-width/2):(ceiling(i_centre-width/2)+width-1)] <- 1
  }
  if (shape=="gauss") {
    x <- c(1:len)
    #Gaussian with width defined as width at half of full maximum
    pulse <- exp(-((x-i_centre)^2/(width/4*sqrt(8*log(2)))^2)) 
    
  }
  return(pulse)
}


#Several qualitative color palettes that are colorblind friendly
#Code to generate vectors in R to use these palettes

#From Paul Tol: https://personal.sron.nl/~pault/
Tol_bright <- c('#EE6677', '#228833', '#4477AA', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB')

Tol_muted <- c('#88CCEE', '#44AA99', '#117733', '#332288', '#DDCC77', '#999933','#CC6677', '#882255', '#AA4499', '#DDDDDD')

Tol_light <- c('#BBCC33', '#AAAA00', '#77AADD', '#EE8866', '#EEDD88', '#FFAABB', '#99DDFF', '#44BB99', '#DDDDDD')

#From Color Universal Design (CUD): https://jfly.uni-koeln.de/color/
Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

# Create a reactive object here that we can share between all the sessions.
vals <- reactiveValues(count=0)

# Define UI
ui <- fluidPage(
   
   # Application title
   titlePanel("PlotsOfPhasors - Plotting fluorescence lifetime data"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(width=3,
          conditionalPanel(
            condition = "input.tabs=='Plot'",

            h4("Aesthetics"),
            
            radioButtons(inputId = "data_type", label = h4("Data display"), choices = list("Dots" = "dots", "Density" = "hexagons","Iso density lines"="iso_density"), selected = "dots"),
            
            numericInput("pointSize", "Size of lines/points", min=0, max=10, value = 4),

            sliderInput("alphaInput", "Visibility of the data", 0, 1, 0.8),

            numericInput("freq", "Frequency [MHz]:", value = 40),
            
            
            
            checkboxInput("color_data", "Use color for the data", value=FALSE),
            # checkboxInput("color_stats", "Use color for the stats", value=FALSE),
            #                  selectInput("colour_list", "Colour:", choices = ""),
            # conditionalPanel(condition = "input.color_data == true || input.color_stats == true",
                             
                             conditionalPanel(condition = "input.color_data == true",
                                              radioButtons("adjustcolors", "Color palette:", choices = 
                                            list(
                                              #"Standard" = 1,
                                              "Okabe&Ito; CUD" = 6,
                                              "Tol; bright" = 2,
                                              "Tol; muted" = 3,
                                              "Tol; light" = 4,
                                              "Viridis" = 7,
                                              "User defined"=5),
                                          selected =  6),
                             conditionalPanel(condition = "input.adjustcolors == 5",
                                              textInput("user_color_list", "List of names or hexadecimal codes", value = "turquoise2,#FF2222,lawngreen")), 
                             
                             h5("",
                                a("Click here for more info on color names",
                                  href = "http://www.endmemo.com/program/R/color.php", target="_blank"))
            ),
            
            
            


            h4("Statistics"),
            
            
            # checkboxInput(inputId = "show_table",
            #               label = "Show table with stats",
            #               value = FALSE),
            
            h4("Transformation & Scaling"),
            
            checkboxInput(inputId = "hide_hemi_labels",
                          label = "Hide reference lifetimes",
                          value = FALSE),

           checkboxInput(inputId = "change_scale",
                          label = "Change scale",
                          value = FALSE),
            conditionalPanel(
              condition = "input.change_scale == true",
              
              
              textInput("range_x", "Range x-axis (min,max)", value = "")
              
            ),
            

            conditionalPanel(
              condition = "input.change_scale == true",
              textInput("range_y", "Range y-axis (min,max)", value = "")
              
            ),
            numericInput("plot_height", "Plot height (# pixels): ", value = 500),
            numericInput("plot_width", "Plot width (# pixels):", value = 1000),


            h4("Labels"),
  
            checkboxInput(inputId = "add_title",
                          label = "Add title",
                          value = FALSE),
            conditionalPanel(
              condition = "input.add_title == true",
              textInput("title", "Title:", value = "")
            ),
            
            checkboxInput(inputId = "label_axes",
                          label = "Change axis labels",
                          value = FALSE),
            conditionalPanel(
              condition = "input.label_axes == true",
              textInput("lab_x", "X-axis:", value = ""),
              textInput("lab_y", "Y-axis:", value = "")
              
            ),
            
            checkboxInput(inputId = "adj_fnt_sz",
                          label = "Change font size",
                          value = FALSE),
            conditionalPanel(
              condition = "input.adj_fnt_sz == true",
              numericInput("fnt_sz_title", "Plot title:", value = 24),
              numericInput("fnt_sz_labs", "Axis titles:", value = 24),
              numericInput("fnt_sz_ax", "Axis labels:", value = 18),
              numericInput("fnt_sz_cand", "Labels of hits:", value = 6)
              
            ),

              checkboxInput(inputId = "add_legend",
                            label = "Add legend",
                            value = FALSE),
            
            # conditionalPanel(
            #   condition = "input.add_legend == true",
            #   textInput("legend_title", "Legend title:", value = "")
            # ),
            
              NULL),
          
             conditionalPanel(
                  condition = "input.tabs=='Data'",
              h4("Data upload"),
              
              radioButtons(
                "data_input", "",
                choices = 
                  list("Example data 1" = 1,
                        "Example data 2" = 2,
                       "Paste data" = 4,
                       "Upload TXT or CSV file" = 3,
                       "URL (csv files only)" = 5
                  )
                ,
                selected =  1),
              
              conditionalPanel(
                condition = "input.data_input=='1'",p('<Description goes here>'),hr()),              
              conditionalPanel(
                condition = "input.data_input=='2'",p('Lifetime data from an EPAC sensor acquired at low, intermediate and high cAMP concentrations: https://dx.doi.org/10.1371%2Fjournal.pone.0019170'),hr()),
              
              
              conditionalPanel(
                condition = "input.data_input=='3'",
        
                fileInput("upload", NULL, multiple = FALSE),

                  radioButtons(
                    "upload_delim", "Delimiter",
                    choices =
                      list("Comma" = ",",
                           "Tab" = "\t",
                           "Semicolon" = ";",
                           "Space" = " "),
                    selected = ","),          actionButton("submit_datafile_button",
                                                           "Submit datafile")

                ),
              
              conditionalPanel(
                condition = "input.data_input=='4'",
                h5("Paste data below:"),
                tags$textarea(id = "data_paste",
                              placeholder = "Add data here",
                              rows = 10,
                              cols = 20, ""),
                actionButton("submit_data_button", "Submit data"),
                radioButtons(
                  "text_delim", "Delimiter",
                  choices = 
                    list("Tab (from Excel)" = "\t",
                         "Space" = " ",
                         "Comma" = ",",
                         "Semicolon" = ";"),
                  selected = "\t")),
              
              ### csv via URL as input      
              conditionalPanel(
                condition = "input.data_input=='5'",
                #         textInput("URL", "URL", value = "https://zenodo.org/record/2545922/files/FRET-efficiency_mTq2.csv"), 
                 textInput("URL", "URL", value = ""), 
                NULL
              ),
              h4('Select variables for plotting'),
              checkboxInput(inputId = "GS",
                            label = "G,S coordinates",
                            value = FALSE),
              selectInput("x_var", label = "Tau-phi data column", choices = "-"),
              selectInput("y_var", label = "Tau-mod data column", choices = "-"),
              selectInput("g_var", label = "Column that defines groups", choices = "-"),
              hr(),
              h4('Optional: filtering data'),
              selectInput("filter_column", "Filter based on this parameter:", choices = ""),
              selectInput("remove_these_conditions", "Deselect these conditions:", "", multiple = TRUE),
              

              NULL
              ),
          
          conditionalPanel(
            condition = "input.tabs=='Synthetic'",
            
            #Session counter: https://gist.github.com/trestletech/9926129
            h4("Settings"),
            numericInput("tau1", "Tau 1", value=4),
            numericInput("tau2", "Tau 2", value=1),
            sliderInput("fraction", "Fraction of Tau 2", min = 0, max = 1, animate = TRUE, step = 0.1, value=0),
            # numericInput("freq2", "Frequency [MHz]", value=40),
            
            selectInput("freq2", "Frequency [MHz]:", choices=c("5"=5,
                                                                "10"=10,
                                                                  "20"=20,
                                                                  "40"=40,
                                                                  "80"=80,
                                                                  "160"=160
                                                                ), selected = 40),
            
            numericInput("noise", "Noise", value=1000),
            
            checkboxInput(inputId = "scale_log_10",
                          label = "Log scale for y-axis",
                          value = FALSE),
            
            radioButtons("irf", "Shape of the pulse:", choices=c("block"="block",
                                                               "Gaussian"="gauss",
                                                               "none"="none"),
                         selected = "block"),
            
            numericInput("pulse_centre", "Centre of pulse (bin)", value=25),
            numericInput("pulse_width", "Width of pulse (bin)", value=11),
            
            numericInput("start", "Offset for phasor calculation (bin)", value=25),    
            
            hr()
            ),
          
          
          
          conditionalPanel(
            condition = "input.tabs=='About'",
            
            #Session counter: https://gist.github.com/trestletech/9926129
            h4("About"),  "There are currently", 
            verbatimTextOutput("count"),
            "session(s) connected to this app.",
            hr(),
            h4("Find our other dataViz apps at:"),a("https://huygens.science.uva.nl/", href = "https://huygens.science.uva.nl/")
          )
                   
      ),   #Close sidebarPanel

      # Show a plot
      mainPanel(
        tabsetPanel(id="tabs",
                    tabPanel("Data", h4("Data as provided"),dataTableOutput("data_uploaded")),
                    tabPanel("Plot",h3("Phasor Plot"
                                       ),
                             # downloadButton("downloadPlotPDF", "Download pdf-file"),
                             # #                          downloadButton("downloadPlotSVG", "Download svg-file"),
                             # downloadButton("downloadPlotPNG", "Download png-file"),
                             # 
                             # actionButton("settings_copy", icon = icon("clone"),
                             #              label = "Clone current setting"),

                             plotOutput("coolplot",
                                        height = 'auto')
                              ),
                    # tabPanel("iPlot", h4("iPlot"), plotlyOutput("out_plotly")),
                    tabPanel("Synthetic", h4("Time-domain and Phasor plot"), plotOutput("decayplot")),

                    tabPanel("About", includeHTML("about.html"))
                    )
        
      )   #Close mainPanel
      

   ) #Close sidebarLayout
) #Close fluidPage

server <- function(input, output, session) {
  
  # Session variable - initialize defaults
  g_var.selected <- "-"
  x_var.selected <- "tphi"
  y_var.selected <- "tmod"

  isolate(vals$count <- vals$count + 1)
  ###### DATA INPUT ###################
  
df_upload <- reactive({
  
    
    if (input$data_input == 1) {
      x_var.selected <<- "tphi"
      y_var.selected <<- "tmod"
      g_var.selected <<- "-"
      data <- df_example

    } else if (input$data_input == 2) {
      x_var.selected <<- "tau_phi"
      y_var.selected <<- "tau_mod"
      g_var.selected <<- "id"

      data <- df_example2
    } else if (input$data_input == 3) {
      transform_var_x.selected <<- "-"
      transform_var_y.selected <<- "-"
      file_in <- input$upload
      # Avoid error message while file is not uploaded yet
      if (is.null(input$upload)) {
        return(data.frame(x = "Select your datafile"))
      } else if (input$submit_datafile_button == 0) {
        return(data.frame(x = "Press 'submit datafile' button"))
      } else {
        isolate({
           data <- read.csv(file=file_in$datapath, sep = input$upload_delim, na.strings=c("",".","NA", "NaN"))
        })
      }
      
    } else if (input$data_input == 5) {

      #Read data from a URL
      #This requires RCurl
      if(input$URL == "") {
        return(data.frame(x = "Enter a full HTML address, for example: https://raw.githubusercontent.com/JoachimGoedhart/VolcaNoseR/master/elife-45916-Cdc42QL_data.csv"))
      } else if (url.exists(input$URL) == FALSE) {
        return(data.frame(x = paste("Not a valid URL: ",input$URL)))
      } else {data <- read.csv(input$URL)}
      #Read the data from textbox
    } else if (input$data_input == 4) {
      if (input$data_paste == "") {
        data <- data.frame(x = "Copy your data into the textbox,
                           select the appropriate delimiter, and
                           press 'Submit data'")
      } else {
        if (input$submit_data_button == 0) {
          return(data.frame(x = "Press 'submit data' button"))
        } else {
          isolate({
            data <- read_delim(input$data_paste,
                               delim = input$text_delim,
                               col_names = TRUE)
          })
        }
      }
}  
  
  #Replace space and dot of header names by underscore
  data <- data %>% select_all(~gsub("\\s+|\\.", "_", .))
  
  # observe({print(head(data))})
  
    return(data)
  })
  
  
df_tau_ref <- reactive({
    
  df <- lifetime_to_GS(df_tau, MHz=input$freq)
  return(df)
    
  })
  
  #### DISPLAY UPLOADED DATA (as provided) ##################
  
df_tau_upload <- reactive({
  
  df <- df_upload()
  
  x_choice <- input$x_var
  y_choice <- input$y_var
  g_choice <- input$g_var
  


  if(!input$GS) {
    df <- df %>% select(`tphi` = !!x_choice , `tmod` = !!y_choice)
    df_out <- lifetime_to_GS(df, MHz=input$freq)
  } else if (input$GS) {
    df_out <- df %>% select(`G` = !!x_choice , `S` = !!y_choice)
  }
  
  if (g_choice != "-") {
    df_treat <- df_upload() %>% select(treatment = !!g_choice)
    df_out$treatment <- df_treat$treatment
  } else {df_out$treatment <- "NA"}
  
  # observe({print(df_out)})
  return(df_out)
  
})


output$data_uploaded <- renderDataTable(
    
    #    observe({ print(input$tidyInput) })
#    df_upload(),
    df_upload(),
    rownames = FALSE,
    options = list(pageLength = 20, autoWidth = FALSE,
                   lengthMenu = c(20, 100, 1000, 10000)),
    editable = FALSE,selection = 'none'
  )
  

  ############## Export Normalized data in tidy format ###########
  
# output$downloadTransformedData <- downloadHandler(
#     
#     filename = function() {
#       paste("VolcaNoseR_transformed", ".csv", sep = "")
#     },
#     content = function(file) {
#         write.csv(df_transformed(), file, row.names = FALSE)
#     }
#   )
#   

  ##### Get Variables from the input ##############
  
observe({
  
  #Retrieve the currently selected geom and use as default, even when y_var changes

  df <- df_upload()
  
    var_names  <- names(df)
    var_list <- c("none", var_names)

    # Get the names of columns that are factors.
    nms_fact <- names(Filter(function(x) is.factor(x) || is.integer(x) || is.logical(x) || is.character(x), df))
    nms_var <- names(Filter(function(x) is.integer(x) || is.numeric(x) || is.double(x), df))
    nms_fact <- c("-",nms_fact)


    # Pre-selection works well when example 1 or 2 is selected, but may interfere with URL-loading 
    # updateSelectInput(session, "x_var", choices = nms_var, selected = "log2_FoldChange")
    # updateSelectInput(session, "y_var", choices = nms_var, selected = "minus_log10_pvalue")
    
    
    updateSelectInput(session, "x_var", choices = nms_var, selected = x_var.selected)
    updateSelectInput(session, "y_var", choices = nms_var, selected = y_var.selected)
    
    # updateSelectInput(session, "map_size", choices = mapping_list_all)
    updateSelectInput(session, "g_var", choices = nms_fact, selected = g_var.selected)
    updateSelectInput(session, "filter_column", choices = var_list, selected="none")
   
   # updateSelectizeInput(session, "user_gene_list", selected = genelist.selected)
   # 
   # updateSelectInput(session, "transform_var_x", choices = c("-",nms_var), selected = transform_var_x.selected)
   # updateSelectInput(session, "transform_var_y", choices = c("-",nms_var), selected = transform_var_y.selected)   
   

  })
  
observeEvent(input$GS, {
  
  if (input$GS) {
    updateSelectInput(session, label="G values column:","x_var")
    updateSelectInput(session, label="S values column:","y_var")
  }
  if (!input$GS) {
    updateSelectInput(session, label="Tau-phi column:","x_var")
    updateSelectInput(session, label="Tau-mod column:","y_var")
  }
  
})

  ##### Set width and height of the plot area
  width <- reactive ({ input$plot_width })
  height <- reactive ({ input$plot_height }) 
  
output$coolplot <- renderPlot(width = width, height = height, {
  
  df_upload <- as.data.frame(df_tau_upload()) %>% mutate(treatment = as.factor(treatment))
  if (input$color_data == FALSE) {
    kleur_data <- NULL
  } else {kleur_data <- "treatment"}
  
  newColors <- NULL
  
  if (input$adjustcolors == 2) {
    newColors <- Tol_bright
  } else if (input$adjustcolors == 3) {
    newColors <- Tol_muted
  } else if (input$adjustcolors == 4) {
    newColors <- Tol_light
  } else if (input$adjustcolors == 6) {
    newColors <- Okabe_Ito
  } else if (input$adjustcolors == 5) {
    newColors <- gsub("\\s","", strsplit(input$user_color_list,",")[[1]])
  }
  
  ############## Adjust X-scaling if necessary ##########
  
  #Adjust scale if range for x (min,max) is specified
  if (input$range_x != "" &&  input$change_scale == TRUE) {
    rng_x <- as.numeric(strsplit(input$range_x,",")[[1]])
    observe({ print(rng_x) })
  } else if (input$range_x == "" ||  input$change_scale == FALSE) {
    
    rng_x <- c(NULL,NULL)
  }
  
  ############## Adjust Y-scaling if necessary ##########
  
  #Adjust scale if range for y (min,max) is specified
  if (input$range_y != "" &&  input$change_scale == TRUE) {
    rng_y <- as.numeric(strsplit(input$range_y,",")[[1]])
  } else if (input$range_y == "" ||  input$change_scale == FALSE) {
    
    rng_y <- c(NULL,NULL)
  }
  
  df_ref <- as.data.frame(df_tau_ref())
  #Convert 'Change' to a factor to keep this order, necessary for getting the colors right
  # df$Change <- factor(df$Change, levels=c("Unchanged","Increased","Decreased"))
  
  #call the function, with proper centre
  polar <- make_half_circle(c(0.5,0))
  #add 0,0 to the dataframe, order the data according to x
  polar <- rbind(polar,c(0,0))
  #plot an empty polar plot
  p <- ggplot(data=polar,aes(G,S)) + geom_path() + 
    # coord_fixed(ratio=1, xlim=c(0,1), ylim=c(0,0.5)) +
    xlab("G") +ylab("S")
  
  
  if (!input$hide_hemi_labels) {
  p <-  p+geom_point(data=df_ref, aes(x=G,y=S), alpha = 1, size = 4, shape = 21, fill='white') +
    geom_text(data = subset(df_ref, G > 0.5), aes(x=G,y=S,label=tphi), nudge_x = +.015, nudge_y = 0.01, hjust=0,vjust=0, size=6)+
    geom_text(data = subset(df_ref, G < 0.5), aes(x=G,y=S,label=tphi), nudge_x = -.015, nudge_y = 0.01, hjust=1,vjust=0, size=6)
  }
    
    # This needs to go here (before annotations)
    p <- p+theme_light(base_size = 18) +

    
    #remove gridlines (if selected
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      # theme(title = element_text(family = 'NimbusMon'))
    NULL
    
    if (input$data_type == 'dots') {
      p <- p + geom_point(data=df_upload, aes_string(x='G',y='S', color=kleur_data), alpha = input$alphaInput, size = input$pointSize, shape = 16)
    } else if (input$data_type == 'iso_density') {
      p <- p + geom_density2d(data=df_upload, aes_string(x='G',y='S', color=kleur_data), alpha = input$alphaInput, size = input$pointSize)
    } else if (input$data_type == 'hexagons') {
      p <- p + geom_hex(data=df_upload, aes_string(x='G',y='S'), alpha = input$alphaInput, size = input$pointSize, bins=100)
    }
    
  p <- p + coord_cartesian(xlim=c(rng_x[1],rng_x[2]),ylim=c(rng_y[1],rng_y[2]))

    ########## Do some formatting of the lay-out ##########

  # if title specified
  if (input$add_title == TRUE) {
    #Add line break to generate some space
    title <- paste(input$title, "\n",sep="")
    p <- p + labs(title = title)
  }
  
  # # if labels specified
  if (input$label_axes)
    p <- p + labs(x = input$lab_x, y = input$lab_y)
  
  # # if font size is adjusted
  if (input$adj_fnt_sz) {
    p <- p + theme(axis.text = element_text(size=input$fnt_sz_ax))
    p <- p + theme(axis.title = element_text(size=input$fnt_sz_labs))
    p <- p + theme(plot.title = element_text(size=input$fnt_sz_title))
  }
  
  #remove legend (if selected)
  if (input$add_legend == FALSE) {  
    p <- p + theme(legend.position="none")
  }
  
  
  if (!input$color_data) {p <- p+ scale_fill_viridis_c()}
  if (input$color_data && input$adjustcolors >1 && input$adjustcolors < 7) {
    p <- p+ scale_color_manual(values=newColors)
    p <- p+ scale_fill_manual(values=newColors)
  } else if (input$adjustcolors ==7) {
    p <- p+ scale_colour_viridis_d()
    p <- p+ scale_fill_viridis_d()      
  }
  
  p
  })
  


output$decayplot <- renderPlot(width = 600, height = 800, {
  
  freq <- (as.numeric(input$freq2))
  cycle <- 1/(freq*1e6)*1e9 #in nanoseconds
  bin_time <- cycle/256
  n_bins <- 256
  bins <- c(0:(n_bins-1))
  
  
  df_decay <- data.frame(time=bins*bin_time, bin=bins)
  
  #Add decay values
  df_decay <- df_decay %>% mutate(intensity=1000*(1-input$fraction)*exp(-(time)/input$tau1)+
                                    1000*(input$fraction)*exp(-(time)/input$tau2))


  #Define an IRF
  df_decay$pulse <- pulse_profile(shape=input$irf, centre=input$pulse_centre, width=input$pulse_width, len=n_bins)
  
  #Normalize the IRF to an area of 1. In this way, it does not affect the amplitude of the decay
  df_decay$pulse <- df_decay$pulse/sum(df_decay$pulse)
  
  #Convolve pulse with theoretical excited state decay
  df_decay$intensity <- convolve(df_decay$intensity, rev(df_decay$pulse), type = "c")
  
  #Generate vector with noise
  df_decay$noise <- rpois(length(df_decay$time),input$noise)
  
  #Add noise to decay
  df_decay <- df_decay %>% mutate(int=intensity-input$noise+noise)
  
  #Calculate a pulse for display, by normalizing the intensity to the peak of the experimental data
  df_decay$pulse_plot <- df_decay$pulse/(max(df_decay$pulse))*1000
  
  #Define plot for decay
  p1 <- ggplot(df_decay, aes(x=bin, y=int))+geom_line(color="black", size=1)+theme_bw(base_size = 14)+
    geom_vline(xintercept = (input$start), size=1, alpha=0.5, linetype="longdash") +
    geom_line(aes(x=bin, y=pulse_plot), color="blue")
  # if log-scale checked specified
  if (input$scale_log_10)
    p1 <- p1 + scale_y_log10(limits=c(10,NA)) 
  
  p1 <- p1 + scale_x_continuous(sec.axis = sec_axis(~ .*bin_time, name="Time [ns]"), breaks=seq(0,256, 64))
  p1 <- p1 + labs(y="Count", x="Bin number")
  
  # G <- sum(df_decay$int*cos(2*3.141593*(df_decay$time-(input$start*bin_time))/cycle))/sum(df_decay$int)
  # S <- sum(df_decay$int*sin(2*3.141593*(df_decay$time-(input$start*bin_time))/cycle))/sum(df_decay$int)
  
  G <- sum(df_decay$int*cos(2*3.141593*(df_decay$bin-(input$start))/n_bins))/sum(df_decay$int)
  S <- sum(df_decay$int*sin(2*3.141593*(df_decay$bin-(input$start))/n_bins))/sum(df_decay$int)
    
  G_pulse <- sum(df_decay$pulse*cos(2*3.141593*(df_decay$bin-(input$start))/n_bins))/sum(df_decay$pulse)
  S_pulse <- sum(df_decay$pulse*sin(2*3.141593*(df_decay$bin-(input$start))/n_bins))/sum(df_decay$pulse)
  
  df_GS <- data.frame(G=G, S=S, G_pulse=G_pulse, S_pulse=S_pulse)
  
  observe({print(df_GS)})
  
  #define a function for a circle
  make_half_circle <- function(center = c(0,0),diameter = 1, npoints = 100){
    r = diameter / 2
    tt <- seq(0,pi,length.out = npoints)
    xx <- center[1] + r * cos(tt)
    yy <- center[2] + r * sin(tt)
    return(data.frame(x = xx, y = yy))
  }
  #call the function, with proper centre
  polar <- make_half_circle(c(0.5,0))
  #add 0,0 to the dataframe and 1,0 to close the shape
  polar <- rbind(polar,c(0,0),c(1,0))
  #plot an empty polar plot
  empty_polar <- ggplot(polar,aes(x,y)) + geom_path() + 
    coord_fixed(ratio=1, xlim=c(-.1,1.1), ylim=c(-0.1,0.6)) + xlab("G") +ylab("S") +
    scale_x_continuous(breaks=c(0,0.5,1.0)) +
    scale_y_continuous(breaks=c(0,0.5)) +
    theme_bw(base_size = 14)
  
  p2 <- empty_polar + geom_point(data=df_GS, aes(x=G, y=S), alpha=1, color="black", size=5) +
    geom_point(data=df_GS, aes(x=G_pulse, y=S_pulse), alpha=1, color="blue", size=5)
  

  #Use patchwork to show both plots:
  p1 / p2
  
  
  
})



  
  ######### DEFINE DOWNLOAD BUTTONS FOR ORDINARY PLOT ###########
  
  output$downloadPlotPDF <- downloadHandler(
    filename <- function() {
      paste("PlotsOfPhasors_", Sys.time(), ".pdf", sep = "")
    },
    content <- function(file) {
      pdf(file, width = input$plot_width/72, height = input$plot_height/72)
      plot(plot_data())
      
      dev.off()
    },
    contentType = "application/pdf" # MIME type of the image
  )
  
  
  output$downloadPlotPNG <- downloadHandler(
    filename <- function() {
      paste("PlotsOfPhasors_", Sys.time(), ".png", sep = "")
    },
    content <- function(file) {
      png(file, width = input$plot_width*4, height = input$plot_height*4, res=300)
      plot(plot_data())

      dev.off()
    },
    contentType = "application/png" # MIME type of the image
  )  
  

########### Update count #########
# Reactively update the client.
output$count <- renderText({
  vals$count
})

# When a session ends, decrement the counter.
session$onSessionEnded(function(){
  isolate(vals$count <- vals$count - 1)
})

######## The End; close server ########################
} #Close server


# Run the application 
shinyApp(ui = ui, server = server)

