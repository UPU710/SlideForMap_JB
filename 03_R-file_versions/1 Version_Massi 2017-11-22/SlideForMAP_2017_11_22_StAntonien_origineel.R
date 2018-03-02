# Finite Slope Model (elliptical)
######################################
# Author: Schwarz M., (1.10.2015)
# Bern University of Applied Sciences
######################################

#########################################################
# Loading PACKAGES
############################################################
library(pROC)
require(pROC)

library(raster)
require(raster)

library(dynatopmodel)
require(dynatopmodel)

library(SDMTools)
require(SDMTools)


#########################################################
# Loading DEM
############################################################

################
# Extent DEM


d<-read.asc("~/Documents/Work/SHL/Publikationen/2018_1_StAntonien_SlideforMAP/SlideforMAP_2017/QGIS/ClipStAntonien_22_fill.asc")


dem<-raster(d)
#dem<-crop(dem, crop1)
dem
plot(dem)



x_max<- extent(dem)[2] # 556168  #780700  #xmax(tab_2005$Koord_WE)  # m, width of the study area,  Breite Gebiet
x_min<- extent(dem)[1]#554428  #min(tab_2005$Koord_WE)

y_max<-  extent(dem)[4]    #122454 #206300 #max(tab_2005$Koord_NS)   # m, length of the study area,   länge Gebiet
y_min<- extent(dem)[3]   #120744  #min(tab_2005$Koord_NS)

ext<-extent(x_min, x_max, y_min, y_max) 
crop1<-c(x_min, x_max, y_min, y_max)
# Import random DEM

#writeRaster(dem, "DEM_crop", format = "ascii")


# calculation of contributing area and topographic wetness index
#upslope_raster<-upslope.area(dem, atb=T)
#plot(upslope_raster, 2)
#r_tpi <- raster(xmn=extent(dem)[1], ymn=extent(dem)[3], xmx=extent(dem)[2], ymx=extent(dem)[4], res=2)

#r_tpi<-raster(upslope_raster, 2)
#extent(r_tpi)<-ext
#r_tpi<-setExtent(r_tpi, ext, keepres=T, snap=T)
#plot(r_tpi, add=T)

# Importing raster of topographic wetness index, obtained with SAGA_TWI_tool -> using plugin for Qgis
twi<-read.asc("~/Documents/Work/SHL/Publikationen/2018_1_StAntonien_SlideforMAP/SlideforMAP_2017/QGIS/TWI_Clip_22.asc")
#d
TWI<-raster(twi)
TWI[TWI<0]<-0

plot(TWI)

# Importing raster of slope in degrees, obtained with SAGA_TerrainAnalysis_Morphology
slope<-read.asc("~/Documents/Work/SHL/Publikationen/2018_1_StAntonien_SlideforMAP/SlideforMAP_2017/QGIS/Slope_clip22.asc")
#d
#slope<-slope*180/pi
r_slope<-raster(slope)  #[°]
plot(r_slope)

#########################################################
# Loading data events
############################################################


# Loading data of landslides, 2005

f<-"~/Documents/Work/SHL/Publikationen/2018_1_StAntonien_SlideforMAP/SlideforMAP_2017/QGIS/landslideStAntonien22.csv"



tab_2005<-data.frame(read.table(f,header=TRUE,sep=";",dec=".",na.strings="None"))
str(tab_2005)


data_TWI<-extract(TWI,matrix(c(tab_2005$X,tab_2005$Y), ncol=2))

#data_slope<-c(tab_2005$slope)
data_slope<-(extract(r_slope, cbind(tab_2005$X,tab_2005$Y), buffer=10, fun=mean))  


hist(data_slope)                     #Wichtig Resultate/Anhang
mean(data_slope)
sd(data_slope)

plot(density(data_slope),col="red",xlab='slope',ylab='density',lwd=3) #2005 
#lines(density(data_slope),col="green",xlab='slope',ylab='density',lwd=3)   #1997 
#lines(density(data_slope),col="blue",xlab='slope',ylab='density',lwd=3) #1944
#lines(density(data_slope),col="orange",xlab='slope',ylab='density',lwd=3) #1979
#lines(density(data_slope),col="black",xlab='slope',ylab='density',lwd=3) #2011        



####################loading all landslides

f<-"~/Documents/Work/SHL/Publikationen/2018_1_StAntonien_SlideforMAP/SlideforMAP_2017/StAntonien_landlisdes_2005_all.csv"



tab_tot<-data.frame(read.table(f,header=TRUE,sep=",",dec=".",na.strings="None"))
str(tab_tot)
no_NA<-!is.na(tab_tot$Maechtigkeit)
#no_NA<-!is.na(tab_tot$Flaeche)
data_area<-tab_tot$Flaeche[no_NA]

breaks <- seq(0,max(data_area)+10,10)

hist(data_area, breaks = breaks, col="black", xlab="Landslide area [m ]", ylab="Probability density [-]")

hist_data<-hist(data_area, breaks = breaks)

freq_data<-hist_data$counts/sum(hist_data$counts)



############################################################
#loading data soil depth from fieldmesurement
############################################################

data_depth<- tab_tot$Maechtigkeit[no_NA]



data_xslope<-c(tab_tot$Neigung.Rutschflaeche[no_NA])


hist(data_xslope)
mean(data_xslope)

data_depth_r<-data_depth*cos(data_xslope*pi/180)


hist(data_depth_r)
mean(data_depth_r)   #use these values to calibrate model for the random generation of soil depth
sd(data_depth_r)



#########################################################
#Generation of random distributed landslides coordinates
############################################################
(x_max-x_min)*(y_max-y_min)

n<-100							#16000 number of total generated landsldies, Anzahl hypothetischer Rutschungen (16000Rutschungen/km2=0.016R/m2)

ls_density<-n/((x_max-x_min)*(y_max-y_min))
ls_density

x<-  runif(n,x_min,x_max)
y<-  runif(n,y_min,y_max)


plot(dem)
points(x,y,col="red", pch=10, cex=1)


#########################################################
#Fitting the landslide area distribution function, using the formulation of Malamud et al., 2004
############################################################


A_L<-seq(0.00001, max(breaks)*0.000001, 0.00001)   #classes of landslide areas are in 10 m^2 == 0.00001 km^2
length(A_L)
max(A_L)*1000000
min(A_L)*1000000

loops<-10000
rho_gamma_rand<-	runif(loops,1.4,1.7)			#1.6
a_rand<- runif(loops,0.5*10^(-4),1.5*10^(-4)) 	#9.99*10^(-5)
s_rand<-  runif(loops,1.9*10^(-10),1.9*10^(-7))    # 5.26163e-08

SSE<-s_rand

index<-1
while(index<loops+1)
{
	
	p_A_L<-(1/(a_rand[index]*gamma(rho_gamma_rand[index])))*((a_rand[index]/(A_L-s_rand[index]))^(rho_gamma_rand[index]+1))*exp(-a_rand[index]/(A_L-s_rand[index])) 
	
	freq_model<-p_A_L/sum(p_A_L)
	errors<-freq_model-freq_data
	SSE[index]<-sum(errors^2)
		
	index<- index+1
}

plot(rho_gamma_rand, SSE, ylim=c(min(SSE)*0.9,min(SSE)*1.1))
rho_gamma<-	rho_gamma_rand[which(SSE==min(SSE[!is.na(SSE)]))]  #1.5	#1.4
rho_gamma				#Wichtig resultate!!!

plot(a_rand, SSE)
a<-   a_rand[which(SSE==min(SSE[!is.na(SSE)]))]          #1.28*10^(-3)
a

plot(s_rand, SSE, ylim=c(min(SSE)*0.9,min(SSE)*1.1))
s<-   s_rand[which(SSE==min(SSE[!is.na(SSE)]))]          #0.09*10^(-4)  #1.82*10^(-4)
s







p_A_L<-(1/(a*gamma(rho_gamma)))*((a/(A_L-s))^(rho_gamma+1))*exp(-a/(A_L-s))  #probability distribution malamut et al., 2004. Landslide area is in square kilometers.
 sum(p_A_L)
 


plot(breaks[2:length(breaks)],hist_data$counts/sum(hist_data$counts), pch=19, col="green4", ylab="Density in 10m^2 classes [-]", xlab="Landslide area classes [m^2]")
lines(1000000*A_L,p_A_L/sum(p_A_L), col="red", lwd=3)
grid()



#####################################################################################################################
# plot of the cumulative distribution of landslide using the fitted probability density function
#####################################################################################################################

cum_p_A_L<-p_A_L

i<-2
while(i<length(p_A_L)+1)
{
cum_p_A_L[i]<-p_A_L[i]+cum_p_A_L[i-1]
i<-i+1

}

plot(1000000*A_L,cum_p_A_L/sum(p_A_L))#, xlim=c(0,0.005))


n_unif_landslides<-runif(n,0,1)

A_landslide<-n_unif_landslides
i<-1
while(i<length(n_unif_landslides)+1)
{
	id<-1
		while(id<length(cum_p_A_L)+1)
		{
			if(n_unif_landslides[i]>cum_p_A_L[id]/sum(p_A_L))id<-id+1
			
			if(n_unif_landslides[i]<cum_p_A_L[id]/sum(p_A_L))break
			
		}
	
A_landslide[i]<-A_L[id]
i<-i+1

}

hist(A_landslide)
hist_A_landslide<-hist(A_landslide/0.000001, breaks = seq(0,100000,100), plot=FALSE)
plot(n_unif_landslides,A_landslide)

ls_area<-1000000*A_landslide



classes_ls_area<-hist(ls_area, breaks=c(1000000*A_L,max(1000000*A_L+(0.0001/2)))-(0.0001/2),plot=FALSE)

classes_ls_area$counts
sum(classes_ls_area$counts)
classes_ls_area$mids
plot(classes_ls_area$mids,classes_ls_area$counts/length(n), pch=19, col="orange")     # random generated distribution of landslide areas
points(breaks[2:length(breaks)],n*hist_data$counts/sum(hist_data$counts), pch=19, col="green4")
grid()

hist(ls_area, breaks = seq(1,100000,100), freq = NULL,
     include.lowest = TRUE, right = TRUE,
     density = NULL, angle = 45, col = NULL, border = NULL,
     xlim = range(ls_area), ylim = NULL,
     axes = TRUE, plot = TRUE, labels = FALSE)	#


#################################################################################
#   SOIL DENSITY
#################################################################################

dens_soil<- 1.4                   				# [t/m^3]Bodendichte [g/cm^3], Literatur, Trockendichte (sandstone 1.6-2.68, mudstone 2.4-2.7) Feldaufnahme

#################################################################################
#   Calculation LADNSLIDE GEOMETRY (LENGTH, WIDTH)
#################################################################################

c<- 2       									# Quotient zwischen Breite 1 und Länge 2
ls_width<- sqrt((ls_area*4)/(pi*c))				# Rutschbreite
ls_length<-ls_width*c   						# Rutschlänge = Rutschbreite*2


#################################################################################
#   calculation of the slope of each random generated landslide
############################################################
coord_ls<-matrix(c(x,y), ncol=2)

coord_ls[,1]<-x
coord_ls[,2]<-y

slope_raster<-r_slope
slope_ls<-(extract(slope_raster, coord_ls, buffer=ls_width, fun=mean )) #, fun=mean))   #buffer=ls_width, fun=mean : extract the mean value on a radius equal to landslide width

mean(slope_ls)


#slope_ls<-atan(extract(slope_raster, coord_ls, buffer=ls_width, fun=mean))*180/pi   #buffer=ls_width, fun=mean : extract the mean value on a radius equal to landslide width

#################################################################################
#   calculation of the topografic index of each random generated landslide
############################################################

tpi_ls<-extract(TWI,coord_ls, buffer=ls_width, fun=mean)    #[,2]  define the layer of tpi, [,1]  is the contributing area

#writeRaster(tpi_ls, "~/Documents/Work/SHL/Praesentationen/2015_BAFU/SlideforMAP/tpi_zoom.asc", format="ascii")

#################################################################################
#   Estimation parameters: slope, soil depth, soil cohesion, friction angle, 
############################################################


slope<- slope_ls  						# Hangneigung [�] 
hist(slope, breaks = seq(0,90,1))
slope_rad<-2*pi*slope/360         				# [rad]


soil_depth<-rnorm(n,1.14,0.39) * (1-pnorm(slope,1.35 *mean(data_xslope),0.75*sd(data_xslope)))		# Feldaufnahme mean(data_slope), sd(data_slope) Rutschkörpermächtigkeit/-tiefe [m], pnorm serve ad evitare che ci siano frane in zone ripide!
which(soil_depth<0)
soil_depth[which(soil_depth<0)]<-0.01
plot(slope,soil_depth)                                            #Wichtig Resultate/Anhang
points(data_xslope,data_depth_r, pch=19,col="orange")								# Mittelwert für St.A.: 1.2 m
#hist(soil_depth, breaks = seq(0,5,0.1))

soil_cohesion_unsaturated<- rnorm(n,0.5,0.05)              	# Bodenkohäsion [kPa] SC:5.3(11.6),5.8(6.7)
soil_cohesion_unsaturated[which(soil_cohesion_unsaturated<0)]<-0.01							# SC:Clayey sand, SC-SM:Clayey sand-silty sand
#plot(slope,soil_cohesion_unsaturated)
hist(soil_cohesion_unsaturated, breaks = seq(0,15,0.2))		# 

phi_soil<-  rnorm(n,35,3) # slope_ls*rnorm(n,1,0.1)                  	# Reibungswinkel [°] SC:34.8(4.7), 35.3(4.7)
hist(slope_ls, breaks = seq(0,90,1))
hist(phi_soil, breaks = seq(0,90,1))
phi_soil_rad<- 2*pi*phi_soil/360  				# [rad]


#################################################################################
#   ESTIMATION LADNSLIDE VOLUME & MASS
#################################################################################

ls_volume<- ls_area*soil_depth*cos(slope_rad)					# Rutschvolumen = Rutschfläche*-Tiefe

ls_circ<- pi*(3*((ls_length/2)+(ls_width/2))-sqrt((3*(ls_length/2)+(ls_width/2))*((ls_length/2)+(3*(ls_width/2)))))  			#2*pi*ls_width         #

ls_area_circ_sup<- soil_depth*ls_circ/2			#  [m^2]

mass<- ls_volume*dens_soil						# Masse [t]

#################################################################################
#   CALIBRATION SOIL COHESION
#################################################################################

################### calibration of unsaturated cohesion
#mm_h<-0  #[mm/h],rainfall intensity
#T<-1

#R<-  0.001*mm_h/3600 #[m/s],rainfall intensity
#R

#m<- (R/T)*tpi_ls       #a steady-state shallow subsurface flow model, derived from TOPMODEL assumptions (Beven and Kirkby, 1979) , Park et al..(2013)
#hist(m)


#m[which(m>1)]<- 1


#hist(soil_cohesion_unsaturated,col="red4", breaks = seq(0,50,0.1))

#soil_cohesion <-(1-m)*soil_cohesion_unsaturated

#water_pressure<-  soil_depth*m*10 #rnorm(n,0.5,0.1) 	# [kPa], 1kPa <=> 10 cm watertable	 # Porenwasserdruck

#force_par<- (mass*9.81*sin(slope_rad))
#force_per<-(mass*9.81*cos(slope_rad))
#force_per_eff<-force_per-(ls_area*water_pressure)

#shear_res_bas<- (ls_area*soil_cohesion)+(force_per_eff*(tan(phi_soil_rad)))
#shear_res<- shear_res_bas

#SF<-shear_res/force_par

#cohesion_min<- (force_par-force_per_eff*(tan(phi_soil_rad)))/ls_area

#soil_cohesion_unsaturated <- ifelse(cohesion_min>soil_cohesion_unsaturated, cohesion_min, soil_cohesion_unsaturated) #condition that correct negative values of cohesion



#hist(cohesion_min,col="green4", breaks = seq(0,10,0.1), add=T)
#hist(soil_cohesion_unsaturated,col="green1", breaks = seq(0,50,0.1), add=T)


#min(soil_cohesion_unsaturated)



#hist(SF[SF<2],breaks = seq(0,2,0.1))


#soil_cohesion <-(1-m)*soil_cohesion_unsaturated

#shear_res_bas<- (ls_area*soil_cohesion)+(force_per_eff*(tan(phi_soil_rad)))
#shear_res<- shear_res_bas

#SF<-shear_res/force_par
#hist(SF[SF<2],breaks = seq(0,2,0.1), add=T,col="green")


#################################################################################
#   Calculations PORE WATER PRESSURE
#################################################################################

#=====Empirical calibration of Trasmissivity given the probaility distribution of calculated m
#
mm_h_calibration<-70 #[mm/h],rainfall intensity # 1h intensi�t der ereigniss, oder statistische Wert der 100jaehrliches Ereigniss zurueck berechnen

R_calibration<-  0.001*mm_h_calibration/3600 #[m/s],rainfall intensity
R_calibration
T_serie<-seq(0.00001,0.01,0.000001) #seq(0.0001,0.1,0.000001)

tpi_ls[is.na(tpi_ls)]<-0

m_sum<-T_serie
i<-1
while(i<length(T_serie)+1)
{
	m_calibration<- ifelse((R_calibration/T_serie[i])*tpi_ls<1,(R_calibration/T_serie[i])*tpi_ls,1) #(R_calibration/T_serie[i])*tpi_ls
	m_sum[i]<- sum(m_calibration) 
	i<-1+i
}

plot(T_serie,m_sum)


treshehold<- length(tpi_ls)*0.15   # this factor represents the assumed percent of fully saturated soil profiles in a catchment during an extrem rainfall, for example 0.6 == 60%.
treshehold

T_serie[max(which(m_sum>treshehold))]

abline(v=T_serie[max(which(m_sum>treshehold))],col="red" )
abline(h=treshehold,col="red" )


T<- T_serie[max(which(m_sum>treshehold))]   #[m/s],estimated transmissivity 

m<- (R_calibration/T)*tpi_ls       #a steady-state shallow subsurface flow model, derived from TOPMODEL assumptions (Beven and Kirkby, 1979) , Park et al..(2013)

hist(m)



################### calulated pore water pressure
mm_h<-35  # 35[mm/h],rainfall intensity

R<-  0.001*mm_h/3600 #[m/s],rainfall intensity
R
T<- T_serie[max(which(m_sum>treshehold))]   #[m/s],transmissivity 

m<- (R/T)*tpi_ls       #a steady-state shallow subsurface flow model, derived from TOPMODEL assumptions (Beven and Kirkby, 1979) , Park et al..(2013)


m[which(m>1)]<- 1

soil_cohesion <-(1-m)*soil_cohesion_unsaturated


water_pressure<-  soil_depth*m*10 #rnorm(n,0.5,0.1) 	# [kPa], 1kPa <=> 10 cm watertable					# Porenwasserdruck

hist(water_pressure*1000)
mean(water_pressure)*1000  #Pa
sd(water_pressure)*1000    #Pa



#=====Back calculated pore waterpressure [Pa], if event analysis dataset is available.
#c_data<-0.5                                        #kN, cohesion estimated based on soil type
#data_mass<-data_volume*dens_soil
#F_par<-data_mass*9.81*sin(data_slope*pi/180)     #kN
#F_per<-data_mass*9.81*cos(data_slope*pi/180)       #kN
#data_water_pressure<- -((F_par-(c_data*data_area))/tan(data_slope*pi/180)-F_per)/data_area    #kPa, assuming  friction angle similar to slope angles!!!
#hist(data_water_pressure)
#mean(data_water_pressure)
#sd(data_water_pressure)
#min(data_water_pressure)

#plot(data_depth,data_water_pressure )

#points(data_water_pressure, data_depth*data_tpi*(R/T), pch=19, col="red")
#plot(data_water_pressure, data_depth*data_tpi*(R/T), xlim=c(0-max(data_water_pressure[!is.na(data_water_pressure)]), max(data_water_pressure[!is.na(data_water_pressure)])))

#random<-1000
#T_random<-runif(random, 0.0001,0.001)
#SSE_pressure<-T_random

#i<-1
#while(i<random+1)
#{
#	se<-sqrt((data_water_pressure-(data_depth*data_tpi*(R/T_random[i])))^2)
#	SSE_pressure[i]<-sum(se, na.rm=T)
	
#	i<-i+1
#}

#plot(T_random, SSE_pressure)

#T_random[which(SSE_pressure==min(SSE_pressure[!is.na(SSE_pressure)]))]

#################################################################################
#   Calculations ROOT REINFORCEMENT
#################################################################################
#~/Documents/Work/EcorisQ/Support_Projekte/Morgins_VS/SlideforMAP/treefile.txt
#trees<- read.table("~/Documents/Work/EcorisQ/Support_Projekte/Morgins_VS/SlideforMAP/Ind_trees.csv", sep=";")
#str(trees)
#~/Documents/Work/SHL/Projects/2017_19_Schwemmholz_BAFU/StudySites/Diemtigtal/SlideforMAP/Rothbad_trees.csv
#tree_raster_dimension<- 10   # [m]

#r <- raster(xmn=extent(dem)[1], ymn=extent(dem)[3], xmx=extent(dem)[2], ymx=extent(dem)[4], res=tree_raster_dimension)
#r[] <- 0


#x_t <- trees[,1]
#y_t <- trees[,2]
#xy <- cbind(x_t, y_t)
#r_density_trees <- rasterize(xy, r, fun=function(x,...)length(x))

#plot(r_density_trees)

#dbh<-trees[,4]
#r_dbh <- rasterize(xy, r, field=dbh, fun=function(dbh,...)mean(dbh))
#r_dist<- sqrt(tree_raster_dimension*tree_raster_dimension/(pi*r_density_trees*3))
#plot(r_dist)


#root_reinf_values<-  2*(50*getValues(r_dbh)*dgamma(getValues(r_dist)/(0.01*getValues(r_dbh)*18.5), 5,15))   #[N], dbh in [cm], distance in [m]. Minimum root reinforcement in the cell
#r_root_reinf<-r_dist
#r_root_reinf<-setValues(r_root_reinf, root_reinf_values)


#extent(r_root_reinf)<-ext
#r_root_reinf<-setExtent(r_root_reinf, ext, keepres=T, snap=T)

#plot(r_root_reinf)

#root_reinf_ls<-(extract.data(coord_ls, r_root_reinf))

#root_reinf_ls[is.na(root_reinf_ls)]<- 0
#hist(root_reinf_ls)


#################################################################################
#   Calculations VEGETATION WEIGHT
#################################################################################



w_veg<- 0.0                      			# Gewicht Vegetation [t/m^2] Tree volume calculation -->werte in literatur suchen


#########################################################

root_lateral<-     1 #  root_reinf_ls/1000       					# Laterale Wurzelverstärkung [kN/m], 2-15
root_basal<- dgamma(soil_depth,3.1,12.57)*0.15*root_lateral  #exp((soil_depth)*-0.057191)*root_lateral                 		 		# Basale Wurzelverstärkung [kPa], 0-2, muss noch angepasst werden

plot(root_basal, soil_depth)

force_par<- (mass*9.81*sin(slope_rad))+(ls_area*w_veg*9.81*sin(slope_rad))
force_per<-(mass*9.81*cos(slope_rad))+(ls_area*w_veg*9.81*cos(slope_rad))
force_per_eff<-force_per-(ls_area*water_pressure)

shear_res_bas<- (ls_area*soil_cohesion)+(ls_area*root_basal)+(force_per_eff*(tan(phi_soil_rad)))
shear_res_lateral<- (ls_circ/2)* root_lateral    #ls_area_circ_sup*root_lateral
shear_res<- shear_res_bas+shear_res_lateral


SF<-shear_res/force_par

#points(ls_volume, SF, lwd=2, col="orange4")

hist(SF[SF<2],breaks = seq(0,2,0.1), add=T, col="red")               #Test fuer Kohaesion, im trockenen sollte es keine Rutschungen geben


################################
#plot(c(0,1000),c(0,3), pch=" ", xlab=expression("Landslide volume [m^3]"), ylab="Safety Factor [-]")
#abline(1.3,0, col="red", lwd=2)
#abline(1,0, col="red", lwd=2,lty=3)

#plot(c(0,max(n)),c(0,max(SF)), pch=" ", xlab=expression("Landslide number [N°]"), ylab="Safety Factor [-]")

########################################################


n.veg<-matrix(c(seq(1,n,1),x,y,SF,ls_area,soil_depth, slope),ncol=7)
n.veg_filtered<-n.veg[n.veg[,4]<1, ]
n.veg_filtered
sum(n.veg_filtered[,5]*n.veg_filtered[,6])


n.veg.area<-as.vector(n.veg_filtered[,5], mode="numeric")
#hist(n.veg.area, breaks = seq(1,10000,100), add=T, col="brown3")

#hist(n.veg.area, breaks = seq(1,3000,100))
hist(n.veg_filtered[,5], breaks = breaks, col="black", xlim=c(0,1500), xlab="Landslide area [m ]", ylab="Probability density [-]")
grid()


hist(n.veg_filtered[,5], add=T, col="red", breaks = breaks)
hist(n.veg_filtered[,5], add=T, col="blue", breaks = breaks)
hist(n.veg_filtered[,5], add=T, col="green3", breaks = breaks)
hist(n.veg_filtered[,5], add=T, col="orange3", breaks = breaks)

legend("topright", c("0 kN/m", "2 kN/m", "5 kN/m", "10 kN/m", "15 kN/m"), col=c("black", "red", "blue", "green3", "orange3"), pch=15 , bty='n')
##################################################********

#punkt<-ifelse(SF<1,19," ")
colo<-ifelse(SF<1,"red","NA")
colo<-ifelse(SF<1,"blue","NA")
#colo<-ifelse(SF<1,"orange","NA")

cex_1<-ifelse(SF<1,0.0004,0.00001)


plot(dem)
plot(TWI, add=T)
#plot(slope_raster)

points(x,y,pch=19, col=colo, cex=2*ls_area*cex_1 )

#SF

points(x,y,pch=19, col="gold", cex=0.1 )

#Data events

points(tab_2005$X,tab_2005$Y, pch=19, cex=2*tab_2005$area/800, col="red")

#Comparison of F-M distribution of measured and modelled landslides

hist(rep(n.veg_filtered[,5],2), breaks = breaks, col="", ylim=c(0, 80), xlim=c(0,1500), xlab="Landslide area [m ]", ylab="Probability density [-]")
grid()
hist(data_area, add=T, col="red", breaks = breaks)
hist(n.veg_filtered[,5], add=T, col="green", breaks = breaks)


hist(rep(n.veg_filtered[,7],2), breaks = seq(0,80,1), col="", ylim=c(0, 80), xlim=c(0,60), xlab="slope [�]", ylab="Probability density [-]")
grid()
hist(data_slope, add=T, col="red", breaks = seq(0,80,1))
hist(n.veg_filtered[,7], add=T, col="black", breaks = seq(0,80,1))

col_sf<-ifelse(n.veg[,4]<1,"red","green")
plot((ls_area*soil_cohesion)+(ls_area*root_basal),n.veg[,6], pch=19, col=col_sf, xlab="Cohesion [Pa]", ylab="Soil depth [m]")
points(data_xslope,data_depth, pch=19,col="orange")

#################################################################################
#   Generation of raster for landslide probability
#################################################################################



raster_resolution<-10
#r_ls <- raster(xmn=extent(dem)[1], ymn=extent(dem)[3], xmx=extent(dem)[2], ymx=extent(dem)[4], res=100)
r_ls <- raster(xmn=min(x), ymn=min(y), xmx=max(x), ymx=max(y), res=raster_resolution)

r_ls[] <- 0


x_ls <- n.veg_filtered[,2]
y_ls <- n.veg_filtered[,3]
xy_ls <- cbind(x_ls, y_ls)

r_number_ls <- rasterize(xy_ls, r_ls, fun=function(x,...)length(x))
writeRaster(r_number_ls, "Slides_without_100", overwrite=TRUE, format = "ascii")

r_number_ls_with <- rasterize(xy_ls, r_ls, fun=function(x,...)length(x))
writeRaster(r_number_ls_with, "Slides_with_forest_100", format = "ascii")

r_number_ls_optimum <- rasterize(xy_ls, r_ls, fun=function(x,...)length(x))
writeRaster(r_number_ls_optimum, "Slides_optimum_100", format = "ascii")

##############
# distribution of landslide probability without root reinforcement
##############
plot(slope_raster,legend=F, col=gray.colors(10, start = 0.3, end = 0.9, gamma = 2.2, alpha = NULL))
plot(100*(r_number_ls)/(6*ls_density*(raster_resolution^2)), add=T, col=c("orange1","orange2","orange3","orange4","green1","green2","green3","green4"))
#Data events

points(tab_2005$X,tab_2005$Y, pch=19, cex=2*tab_2005$area/800, col="red")






 ##############
# distribution of landslide probability with root reinforcement
##############
plot(slope_raster,legend=F, col=gray.colors(10, start = 0.3, end = 0.9, gamma = 2.2, alpha = NULL))
plot(100*(r_number_ls_with)/(6*ls_density*(raster_resolution^2)), add=T, col=c("orange1","orange2","orange3","orange4","green1","green2","green3","green4"))
#Data events
plot(r_number_ls_with, add=T, col=c("orange1","orange2","orange3","orange4","green1","green2","green3","green4"))
#Data events

points(tab_2005$X,tab_2005$Y, pch=19, cex=2, col="red")

##############
# distribution of stabilisation efficienty due to root reinforcement
##############
plot(slope_raster,legend=F)
#plot(slope_raster,legend=F, col=gray.colors(10, start = 0.3, end = 0.9, gamma = 2.2, alpha = NULL))
values<-getValues(r_number_ls)
values<-values[!is.na(values)]
plot(100*(r_number_ls-r_number_ls_with)/r_number_ls, add=T, col=c("lightblue4","blue2","blue3","blue4"))


#################################################################################
#ROC Curves
#################################################################################

#landslides model

mod<- getValues(r_number_ls)

mod[is.na(mod)]<-0
mod



#landslide inventory

data_x<-c(tab_2005$X)
str(data_x)

data_y<-c(tab_2005$Y)
str(data_y)

n.inv<-matrix(c(data_x, data_y),ncol=2)
n.inv

raster_resolution<-10
r_i <- raster(xmn=min(x), ymn=min(y), xmx=max(x), ymx=max(y), res=raster_resolution)

r_i[] <- 0


x_i <- n.inv[,1]
y_i <- n.inv[,2]
xy_i <- cbind(x_i, y_i)
xy_i

r_number_i <- rasterize(xy_i, r_i, fun=function(x,...)length(x))
r_number_i

inv<- getValues(r_number_i)

inv[is.na(inv)]<-0
inv



category <- inv #Inventory
prediction <- mod #model


plot(roc(category, prediction))
roc_obj <- roc(category, prediction)
str(roc_obj)
roc_obj$auc

roc_df <- data.frame(
  TPR=rev(roc_obj$sensitivities), 
  FPR=rev(1 - roc_obj$specificities))  #,labels=roc_obj$response, scores=roc_obj$predictor)
  
  
rectangle <- function(x, y, width, height, density=12, angle=-45, ...) 
  polygon(c(x,x,x+width,x+width), c(y,y+height,y+height,y), 
          density=density, angle=angle, ...)


roc_df <- transform(roc_df, 
  dFPR = c(diff(FPR), 0),
  dTPR = c(diff(TPR), 0))
  
  plot(0:10/10, 0:10/10, type='n', xlab="FPR", ylab="TPR")
abline(h=0:10/10, col="lightblue")
abline(v=0:10/10, col="lightblue")

with(roc_df, {
  mapply(rectangle, x=FPR, y=0,   
         width=dFPR, height=TPR, col="green", lwd=2)
  mapply(rectangle, x=FPR, y=TPR, 
         width=dFPR, height=dTPR, col="blue", lwd=2)

  lines(FPR, TPR, type='b', lwd=3, col="red")
})



TPR<-rev(roc_obj$sensitivities)
TPR
FPR<-rev(1 - roc_obj$specificities)
FPR

plot(FPR,TPR,xlim=c(0,1),ylim=c(0,1), yaxs="i", xaxs="i",xlab="False-positive rate (FPR)", ylab="True-positive rate( TPR)")
lines(FPR,TPR)
abline(0, 1,lty="dashed")

r<-sqrt((1-TPR)^2+FPR^2)
r
min(r)

##############################################################################################

hist(values(r_number_ls-r_number_ls_with))

plot(100*(r_number_ls-r_number_ls_optimum)/6, add=T, col=c("orange1","orange2","orange3","orange4","green1","green2","green3","green4"))


plot(slope_raster,legend=F)
plot(100*(r_number_ls/50), add=T, col=c("lightblue4","blue2","blue3","blue4"))

plot(r_number_ls, add=T, col=c("lightblue4","blue2","blue3","blue4"))

plot(100*(r_number_ls/((raster_resolution^2)*length(x_ls)/((x_max-x_min)*(y_max-y_min)))), add=T, col=c("lightblue4","blue2","blue3","blue4"))


plot(slope_raster,legend=F)
plot(r_number_ls, add=T, col=c("lightblue4","blue2","blue3","blue4"))

plot(r_number_ls_with, add=T, col=c("green1","green2","green3","green4"))

plot(r_number_ls_optimum, add=T, col=c("red1","red2","red3","red4"))

points(x_t,y_t, pch=19, col="green4", cex=dbh*0.005)






breaks_root <- seq(0,max(root_reinf_ls)+100,100)

hist(root_reinf_ls, breaks=breaks_root)





#Trees

trees<- read.table("~/Documents/Work/SHL/Praesentationen/2015_BAFU/SlideforMAP/Ind_trees.csv", sep=";")



trees_sub<-trees

c<-1
p<-1 while(c<length(trees[,1]+1))
{
	
	if(trees[c,1]>x_min && trees[c,1]<x_max && trees[c,2]>y_min  && trees[c,2]<y_max)
	{
		trees_sub[p,]<-trees_sub[c,]
		p<-p+1
	}
	
	
	c<- c+1
}
p
trees_sub<-trees_sub[1:p,]


plot(dem_zoom)
plot(dem)
plot(upslope_raster,2)
points(trees_sub[,1], trees_sub[,2], pch=19, col="green4",cex=trees[,3]/20)



root_ls<-rep(0,length(x))

i<-1  #landslide ID

while(i<n+1)
{
	
	dist<- sqrt((x[i]-trees_sub[,1])^2  + (y[i]-trees_sub[,2])^2)
	
	
	t<-1
	while(t<length(trees_sub[,1]+1))
	{
		if(dist[t]< 10)
		{
		root_ls[i]<-root_ls[i]+(50*trees_sub[t,4]*dgamma(dist[t]/(0.01*trees_sub[t,4]*18.5), 5,15))
		}	
		
		
		t<-t+1
	}
	
	
	i<-i+1
}

hist(root_ls)
plot(y,root_ls)

distanza<-seq(1,10,0.1)
r<-(50*30*dgamma(distanza/(0.3*18.5), 5,15))
plot(distanza,r)
lines(distanza,r)

(50*40*dgamma(distanza/(0.4*18.5), 5,15))


points(x,y,col="brown", cex=root_ls/10000, pch=19)







#Calculationof minimum root reinforcement based on mean Diameter and mean distnace


root_reinf_raster<-(50*trees_sub[t,4]*dgamma(dist[t]/(0.01*trees_sub[t,4]*18.5), 5,15)) 



#  calculation of the number of trees per cell (25x25 m)

library(raster)
#extent(dem)[2]

r <- raster(xmn=extent(dem)[1], ymn=extent(dem)[3], xmx=extent(dem)[2], ymx=extent(dem)[4], res=25)
r[] <- 0


x <- trees[,1]
y <- trees[,2]
xy <- cbind(x, y)
r1 <- rasterize(xy, r, fun=function(x,...)length(x))
plot(r1)
points(xy)

dbh<-trees[,4]

r2 <- rasterize(xy, r, field=dbh, fun=function(dbh,...)mean(dbh))

plot(r2)


