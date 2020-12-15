# PlotsOfPhasors
A Shiny App for plotting fluorescence lifetime data

### Running the App

The web-tool runs from a shiny server, and can be accessed at: [https://huygens.science.uva.nl/PlotsOfPhasors/](https://huygens.science.uva.nl/PlotsOfPhasors/)

Alternatively, the app can run from R/Rstudio.

#### Preparations
Note that the app depends on several R packages that need to be installed (shiny, ggplot2, dplyr, magrittr, DT, shinycssloaders, readr). 

Run this command in R/Rstudio to download and install all the packages (only needs to be done once):
```
install.packages("shiny", "ggplot2", "dplyr", "magrittr", "DT", "shinycssloaders", "readr")
```
o The first option is running it directly from Github. In the command line (in R or Rstudio) type:
```
shiny::runGitHub('PlotsOfPhasors', 'JoachimGoedhart')
```
o The second option is download the app and to use it offline:

-download the `app.R` and all the other files from the repository

-Run RStudio and load `app.R`

-Select 'Run All' (shortcut is command-option-R on a Mac) or click on "Run App" (upper right button on the window)

This should launch a web browser with the Shiny app.


### Credits

PlotsOfPhasors 


VolcaNoseR is created by Franka van der Linden and Joachim Goedhart and maintained by Joachim Goedhart ([@joachimgoedhart](https://twitter.com/joachimgoedhart))

### Example output

Example output generated with the data from 'lifetime_data.csv':

![alt text](https://github.com/JoachimGoedhart/PlotsOfPhasors/blob/master/PlotsOfPhasors_1.png "Output")



