################################################################################
#
#                            VisualFxns.R
#
################################################################################

#11/10/2022

#These functions were written by Dr. Karen Carleton @ University of Maryland, College Park.

#================================================================================
#   FUNCTION: Opsin_gov 
#      Calculate visual pigments absorption using Govardovskii et al 2000 templates
#           lbeg, lend = wavelength range
#           lmax = peak wavelength
#           A1_chrom = % of A1 chromophore
#--------------------------------------------------------------------------------
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
