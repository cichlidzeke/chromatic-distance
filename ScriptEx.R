################################################################################
#
#                                ScriptEx.R
#
################################################################################

#04/26/2022

#This script was written by Zeke Gonzalez @ University of Maryland, College Park
#and modified from script written by Dr. Dani Escobar-Camacho and Dr. Karen Carleton.

################################################################################

#Load the necessary functions.

#===============================================================================
#   FUNCTION: Opsin_gov 
#      Calculate visual pigments absorption using Govardovskii et al 2000 templates
#           lbeg, lend = wavelength range
#           lmax = peak wavelength
#           A1_chrom = % of A1 chromophore
#-------------------------------------------------------------------------------
Opsin_gov<-function(lbeg, lend, lmax, A1_chrom)  {
  # fit A1 peak
  a1<-0.8795+0.0459*exp(-(lmax-300)^2/11940)
  lmb1<-189+0.315*lmax
  p1<--40.5+0.195*lmax  #Govardovskii calls this b but this confuses it with the alpha peak equation
  i<-lbeg:lend
  x<-lmax/i
  Salpha_A1<-1/((exp(69.7*(a1-x))+exp(28*(0.922-x))+exp(-14.9*(1.104-x))+0.674))
  Sbeta_A1<-0.26*exp(-((i-lmb1)/p1)^2)
  Speak_A1<-Salpha_A1+Sbeta_A1
  # fit A2 peak
  a2<-0.875+0.0268*(exp((lmax-665)/40.7))  
  A2<- 62.7 + 1.834*exp((lmax-625)/54.2)
  lmb2<-216.7+0.287*lmax
  p2<-317-1.149*lmax+0.00124*lmax*lmax  #Govardovskii calls this b but this confuses it with the alpha peak equation
  Salpha_A2<-1/((exp(A2*(a2-x))+exp(20.85*(0.9101-x))+exp(-10.37*(1.1123-x))+0.5343))
  Sbeta_A2<-0.26*exp(-((i-lmb1)/p1)^2)
  Speak_A2<-Salpha_A2+Sbeta_A2
  # weight by chromophore, sum and normalize
  Speaktot<-A1_chrom*Speak_A1+(100-A1_chrom)*Speak_A2
  Snorm<-Speaktot/max(Speaktot)
  Gov_results<-data.frame(i,Snorm)
  return(Gov_results)}

#==============================================================================
#   FUNCTION:  Qcatch_calc  
#     Calculate the stimulation of the S, M and L cones in terms of quantum catch
#         light = the illuminant
#         fishlens = the lens transmission
#         Vispig = 3 visual pigments which can be coexpresing and involve A1 or A2
#         targets = the set of color targets
#     This version works on target file where there is no wavelength column
# ----------------------------------------------------------------------
Qcatch_calc<-function(light, fishlens, Vispig, targets) {
  target_labels<-colnames(targets)
  number_colors<-length(target_labels)   # number of colors in set
  #
  # Calculate the von Kries correction for accommodation to the illuminant
  vonK_illum_QCs=0
  vonK_illum_QCm=0
  vonK_illum_QCl=0
  vonK_illum_QCs<-sum(light$IRRAD*fishlens$Transmission*Vispig$VPs)
  vonK_illum_QCm<-sum(light$IRRAD*fishlens$Transmission*Vispig$VPm)
  vonK_illum_QCl<-sum(light$IRRAD*fishlens$Transmission*Vispig$VPl)
  qcatchs=NULL
  qcatchm=NULL
  qcatchl=NULL
  for (j in 1:number_colors) {
    qcatchs[j]<-sum(light$IRRAD*targets[,j]*fishlens$Transmission*Vispig$VPs)/vonK_illum_QCs
    qcatchm[j]<-sum(light$IRRAD*targets[,j]*fishlens$Transmission*Vispig$VPm)/vonK_illum_QCm
    qcatchl[j]<-sum(light$IRRAD*targets[,j]*fishlens$Transmission*Vispig$VPl)/vonK_illum_QCl
  }
  Qcatchdata_results<-data.frame(qcatchs,qcatchm,qcatchl)
  return(Qcatchdata_results)  }

#==============================================================================
#  FUNCTION : JND_calc   
#     Calculate JND values between two sets of targets
#
# ----------------------------------------------------------------------
JND_calc<-function(QCset1, colorlabel1, QCset2, colorlabel2, Webers)  {
  number_colors1<-length(colorlabel1)
  number_colors2<-length(colorlabel2)
  JND_results=data.frame(matrix(nrow=number_colors1, ncol=number_colors2))
  ws<-Webers[1]
  wm<-Webers[2]
  wl<-Webers[3]
  denom<-(ws*wm)^2+(ws*wl)^2+(wm*wl)^2
  for (i in 1:number_colors1) {
    for (j in 1:number_colors2) {
      Dfs<-log(QCset1[i,1]/QCset2[j,1])
      Dfm<-log(QCset1[i,2]/QCset2[j,2])
      Dfl<-log(QCset1[i,3]/QCset2[j,3])
      JND_results[i,j]<-sqrt((ws^2*(Dfm-Dfl)^2+wm^2*(Dfs-Dfl)^2+wl^2*(Dfs-Dfm)^2)/denom)
    } } 
  rownames(JND_results)<-colorlabel1
  colnames(JND_results)<-colorlabel2
  return(JND_results)  }
#
################################################################################
# Set the working directory.

setwd("Your Directory")

### Now you need light, lens, visual pigments, and colors (reflectance).

### FIRST
# Load your illuminant file. This should have a "Wavelength" column and an "IRRAD" column.

light <- read.csv("light.csv", header=TRUE)

### SECOND
# Time to load the lens transmission for your species. 

lens <-read.table("lens.txt", header = TRUE)

### THIRD
# Now you need to put together the visual pigments for your species. 
# In this example script, I will do so for Metriaclima benetos.

# We do this by telling the program what percentage of light the short, medium, 
# and long sensitive photoreceptors will absorb for the species, which is a 
# short palette species.

# Setting the M. benetos lambda max numbers for the short (s), medium (m) 
# and long (l) pigments. I got these numbers from the comparable M. zebra in  
# Parry et al.'s 2005 paper "Mix and Match Color Vision: Tuning Spectral Sensitivity by Differential 
# Opsin Gene Expression in Lake Malawi Cichlids" from CellPress.

Mbenlambdamaxs <- 379

Mbenlambdamaxm <- 489

Mbenlambdamaxl <- 522

# Then setting the lambda minimum (min) and maximum (max) wavelengths.

lambdamin <- 350

lambdamax <- 700

# Finally, setting the percentage of the pigments that use vitamin A1.
# Metriaclima A1 percent is currently estimated at 100%.

A1percent <- 100

#Now I will use all of these variables to calculate the percent of light 
#that the short wavelength pigment will absorb from the minimum wavelength 
#(lambdamin) to the maximum wavelength (lambdamax) for M. benetos.

MbenVPs <- Opsin_gov(lambdamin, lambdamax, Mbenlambdamaxs, A1percent)

#Make column names that I like.

colnames(MbenVPs) <- c("wavelength", "percentabsorb")

#Graph it as a line in green.

plot(MbenVPs, col = "violet", type = "l")

#Now I will use all of these variables to calculate the percent of 
#light that the medium wavelength pigment will absorb from the 
#minimum wavelength (lambdamin) to the maximum wavelength (lambdamax).

MbenVPm <- Opsin_gov(lambdamin, lambdamax, Mbenlambdamaxm, A1percent)

#Rename the columns.

colnames(MbenVPm) <- c("wavelength", "percentabsorb")

#Then add the medium sensitive wavelength pigments as lines to the 
#plot I've already made of the short wavelength-sensitive pigments.

lines(MbenVPm, col = "green", type = "l")

#Then for the long wavelength sensitive.

MbenVPl <- Opsin_gov(lambdamin, lambdamax, Mbenlambdamaxl, A1percent)

#Rename the columns.

colnames(MbenVPl) <- c("wavelength", "percentabsorb")

#Then add the long sensitive wavelength pigments as lines to the 
#plot I've already made of the short/medium wavelength-sensitive pigments.

lines(MbenVPl, col = "blue", type = "l")

# Create array of peak absorbances for each pigment (short, medium, and long).
# THESE ARE GENERAL CICHLID NUMBERS.

peakabss <- 0.009

peakabsm <- 0.015

peakabsl <- 0.015

peakabs <- c(peakabss, peakabsm, peakabsl)

# Create matrix of path length for each photoreceptor type. 
# Here we use the cichlid path lengths.

Ls <- 5.5
Lm <- 26 
Ll <- 26

Mbenlength <- c(Ls, Lm, Ll)

# Create a data set with the percentage of light absorbed by each of the 
# three visual pigment (short, medium, and long) for M. aurora!

MbenVP <- data.frame(MbenVPs$wavelength, MbenVPs$percentabsorb, MbenVPm$percentabsorb, MbenVPl$percentabsorb)

#Now rename the columns.

colnames(MbenVP) <- c("Wavelength", "VPs", "VPm", "VPl")

### The visual pigment data is in MbenVP.

### FOURTH
#Now you need your colors to compare. 
#This should be a file with ONLY the colors, no wavelength column.

colors <- read.csv("Your Colors.csv", header=TRUE)

### FIFTH
# Now I need to calculate the quantum catches (von Kries corrected)! 
# We use the light, lens, visual pigments, and colors that we have put together.

Qcatch <- Qcatch_calc(light, lens, MbenVP, colors)

### SIXTH
# I need to calculate my Weber constants for the long, medium, and short cones. 

# The Weber fraction I'm going to use is 0.16, specific to the blue-grey are of 
# the color space for M. benetos (Escober-Camacho et al. 2019).

#Enter the value for the long cone used by the literature.

WeberL <- 0.16

#Enter the value for the medium cone. In this case, because the ratio of the 
#medium to long cones are 1:1, the fraction is the same, so the value is the same.

WeberM <- 0.16

#Enter the value for the short cone. Here we use Equation (3) from 
#Olsson et al., 2018. Then create an array of our Weber fractions.

WeberS <- (sqrt(2/1)*0.16)

Weber <- c(WeberS, WeberM, WeberL)

### SEVENTH

# Now that I have my Weber constants, it's time to calculate JNDs! 
# In order to run it, I need labels for the colors using my column names.

colorlabels <- c(names(colors))

# Now I'll run all of the JNDs under the light with my Weber constants.

JND <- JND_calc(Qcatch, colorlabels, Qcatch, colorlabels, Weber)

# With these JND values, now all you need to do is get rid of your duplicate values
# (from crossing the colors against themselves).

# From here you are free to calculate your summary statistics, visualize the JNDs, 
# or otherwise process your resultant chromatic distances however you'd like.
