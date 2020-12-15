library(ggplot2)

################# Make an empty polar plot ##############
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
#add 0,0 to the dataframe, order the data according to x
polar <- rbind(polar,c(0,0))
#plot an empty polar plot
empty_polar <- ggplot(polar,aes(x,y)) + geom_path() + 
  coord_fixed(ratio=1, xlim=c(0,1), ylim=c(0,0.5)) + xlab("G") +ylab("S")


#Convert lifetimes to polar coordinates
lifetime_to_GS <- function(data,MHz=40) {
  tphi <- data$tphi
  tmod <- data$tmod
  omega = 2*pi*MHz/1000
  Phi = atan(omega*tphi)
  M=sqrt(1/ (1 + (omega*tmod)^2 ) )
  G = M * cos(Phi)
  S = M * sin(Phi)
  
  new <- data.frame(Phi=Phi, M=M, G=G, S=S)
  
  data <- cbind(data, new)
  
}

#Synthetic mono-exp set of taus
df_tau <- data.frame(tphi=c(0,0.1*2^(0:9)), tmod=c(0,0.1*2^(0:9)))

df_tau_GS <- lifetime_to_GS(df_tau,40)

empty_polar+geom_point(data=df_tau_GS, aes(x=G,y=S)) +
  geom_text(data = subset(df_tau_GS, G > 0.5), aes(x=G,y=S,label=tphi), nudge_x = .01,hjust=0,vjust=0)+
  geom_text(data = subset(df_tau_GS, G < 0.5), aes(x=G,y=S,label=tphi), nudge_x = -.01,hjust=1,vjust=0)

################ plot lifetimes on polar plot #################
setwd("/Users/Franka/surfdrive/mTQ2 paper/calculations_clean_R_code")
#I usually work with G/S data, because this gives gives less noise for low intensity cells
# the dummy datasets includes both options
# "dummy" is a histamine stimulation and "dummy2" a ionomycin stimulation
# in general my data looks much cleaner with ionomycin: all cells respond the same, and high state is sustained longer == better measurements

# dummy <- read.csv("dummy.csv")
dummy <- read.csv("dummy2.csv")


###### in case you work with G/S coordinates
empty_polar + geom_point(data=dummy, aes(x=G, y=S, color=factor(frame)), alpha=0.5)
# outliers (out of polar plot) are an artefact of addition of the histamine. 
# this change in lifetime was too fast for our machine
# you can also spot this in the G and S plot vs time 
ggplot(dummy, aes(x=frame, y=G, group=cell)) +geom_line()
ggplot(dummy, aes(x=frame, y=S, group=cell)) +geom_line()

# if you want to do a correction on the data or you want the lifetimes: use function below
correctionGS <- function(data, MHz = 40) {
  G <- data$G
  S <- data$S
  omega = 2*pi*MHz/1000
  Phi = atan(S/G)
  M = sqrt( S^2 + G^2 )
  M_corr = M/M_factor
  Phi_corr = Phi + Phi_addition
  G_corr = M_corr * cos(Phi_corr)
  S_corr = M_corr * sin(Phi_corr)
  tphi_corr = tan(Phi_corr) / omega
  tmod_corr = sqrt( (1 / M_corr^2 ) - 1 ) / omega
  new <- data.frame(Phi=Phi, M=M, Phi_corr=Phi_corr, M_corr=M_corr,
                    G_corr=G_corr, S_corr=S_corr, tphi_corr=tphi_corr, tmod_corr=tmod_corr)
  data <- cbind(data, new)
}

# set correction factors. correction factors are used in the function
# M_factor = 1 and Phi_addition = 0 for no correction
# M_factor <- 1
# Phi_addition <- 0
M_factor <- 1.0471235457417
Phi_addition <- 0.0453607956954

#correct the data
dummy_corrected_GS <- correctionGS(dummy)

#look at it
empty_polar + geom_point(data=dummy_corrected_GS, aes(x=G_corr, y=S_corr, color=factor(frame)), alpha=0.5)


#####in case you work with lifetimes

#function calculates the G and S, also does a correction of the data
correctionLifetime <- function(data, MHz = 40) {
  tphi <- data$tphi
  tmod <- data$tmod
  omega = 2*pi*MHz/1000
  Phi = atan(omega*tphi)
  M=sqrt(1/ (1 + (omega*tmod)^2 ) )
  Phi_corr = Phi + Phi_addition
  M_corr = M/M_factor
  G_corr = M_corr * cos(Phi_corr)
  S_corr = M_corr * sin(Phi_corr)
  tphi_corr = tan(Phi_corr) / omega
  tmod_corr = sqrt( (1 / M_corr^2 ) - 1 ) / omega
  new <- data.frame(Phi=Phi, M=M, Phi_corr=Phi_corr, M_corr=M_corr,
                    G_corr=G_corr, S_corr=S_corr, tphi_corr=tphi_corr, tmod_corr=tmod_corr)
  data <- cbind(data, new)
}

# set correction factors. correction factors are used in the function
# M_factor = 1 and Phi_addition = 0 for no correction
# M_factor <- 1
# Phi_addition <- 0
M_factor <- 1.0471235457417
Phi_addition <- 0.0453607956954

#correct the data
dummy_corrected_lifetime <- correctionLifetime(dummy)

#look at it
empty_polar + geom_point(data=dummy_corrected_lifetime, aes(x=G_corr, y=S_corr, color=factor(frame)), alpha=0.5)


#when using dummy (not dummy2): now this looks suddenly much better than the GS corrected data: no points outside of the polar plot
#BUT:this is due to earlier processing of my data. tmod and tphi in the dummies are averages of cells, caculated from G and S images
# I caluculated tmod from G and S for the images with this formula:
# tmod = 1/omega * sqrt( 1/(S^2+G^2)  - 1 )
# which only works if 1/(S^2+G^2) is smaller than 1. 
# pixels with 1/(S^2+G^2) is bigger than 1 were set to NaN and not used in the averaging of the cells.
## so all 'outside of polarplot' pixels do not count at all in the tphi and tmod, but do in the G and S. Tphi and tmod in this dataset are skewed.
## in visual terms: converting the tphi and tmod to G and S pixels and plotting all pixels on the polar plot will give not a circle of dots
## while plotting G and S pixels directly will give a circle.
# note: the difference between the two corrected datasets is mainly in the M, hardly in Phi.


############### calculation from G/S to fraction ##############
# you need
#       - the G and S coordinates : see above
#       - ratio of the intensity contribution of both states
#       - the extreme values, stored in a csv

# we assume a ratio of 3.02, in vivo. 
Ratio = 3.02
# we assume the following extreme values
extremes <- read.csv("extremes.csv", stringsAsFactors = FALSE, row.names = 1)



## calculate dS and dG, using the G and S from the low state
dummy_corrected_GS$dG <- dummy_corrected_GS$G_corr - extremes["min","G"]
dummy_corrected_GS$dS <- dummy_corrected_GS$S_corr - extremes["min", "S"]

##calculate inproduct (projection on the line between low and high state).
dummy_corrected_GS$IN <- dummy_corrected_GS$dG * extremes["max","dG"] + dummy_corrected_GS$dS * extremes["max","dS"]

## calculate a line fraction (inproduct divided by total length of line {= dG max state ^ 2 + dS max state ^ 2} = IN max state)
dummy_corrected_GS$lfraction <- dummy_corrected_GS$IN/(extremes["max","dG"]^2+extremes["max","dS"]^2)

##convert to true fraction
dummy_corrected_GS$Fraction <- dummy_corrected_GS$lfraction/(Ratio*(1-dummy_corrected_GS$lfraction)+dummy_corrected_GS$lfraction)


ggplot(dummy_corrected_GS, aes(x=frame, y=Fraction)) + geom_line(aes(color=factor(cell)))
ggplot(dummy_corrected_GS, aes(x=frame, y=lfraction)) + geom_line(aes(color=factor(cell)))


## now this might look a bit weird when you used dummy instead of dummy2: namely cells with fraction >1.
# I used both dummies for calculation of the correctionvalues
# so the average top of the dummies was set to the position of the high state.
# dummy contains the outliers.

