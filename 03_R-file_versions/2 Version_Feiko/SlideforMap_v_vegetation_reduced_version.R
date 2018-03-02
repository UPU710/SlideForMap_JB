## SlideForMap model with Feiko adaptations, Author: Schwarz M., (1.10.2015), Bern HAFL, Finite Slope Model (elliptical)

progress_display <- TRUE # decide if you want to display the model steps
#do_once <- function(){}replicate(20, do_once()) # code provided for a loop of multiple runs

#### 01: LOADING PACKAGES ####
if(progress_display){
  print("01: LOADING PACKAGES")
}
library(pROC)
library(raster)
library(dynatopmodel)
library(SDMTools)
library(velox)
library(sp)
library(rgeos) 
library(rgdal)
library(gstat)
#### 02: DEFINE INPUT PATHS AND NAMES ####
if(progress_display){
  print("02: DEFINE INPUT PATHS AND NAMES")
}
  
filePath <-  "C:/Users/vof1/Documents/SlideForCe-Feiko/study_area_3/"                   # the input path
outPath <- paste0(filePath,"R_output/")                           # the output path

# names of the input files
DEMFile <-  "DEM_clipped.asc" # the name of the digital elevation model
TWIFile <- "TWI_clipped.asc" # name of the topographic wetness index file
SlopeFile <- "Slope_clipped.asc" # name of the file with the slope values
landslidesFile <- "Hangmuren_Daten_Final.csv" # data base of existing landslides
treeFile <- "Ind_trees.csv" # file of all existing trees in the study area

# define the names of the landslide variables from the inputfile #latitude, longitude, slide depth, slide surface area and slide slope
slidevariablenames <- c("latitude", "longitude", "mmacht", "anrissflaeche", "rutschneigung") 
#### 03: DEFINE RUN SETTINGS####
if(progress_display){
  print("03: DEFINE RUN SETTINGS")
}

# decide to run without Nan's or not
noNANs <- FALSE

# setting for the landslide inventory
slides_subsetting <- TRUE

# setting for the soil depth
random_soil_depth <- FALSE
kriging_soil_depth <- TRUE

# settings for vegetation 
vegetation_included <- FALSE

# settings for the root lateral reinforcement method
global_uniform_rootlat <- TRUE
forest_uniform_rootlat <- FALSE
variable_rootlat <- FALSE
root_lateral_strength <- 5.15              # lateral root reinforcement (Kpa), only used for uniform rootlat

# settings for the vegetation weight method
global_uniform_vegweg <- TRUE
forest_uniform_vegweg <- FALSE
variable_vegweg <- FALSE
uniform_vegetation_weight <- 0.5        # uniform vegetation weight (tons/m^2), only used for uniform vegetation weights

# settings for validation 
validation <- TRUE
auc_validation <- TRUE
area_validation <- TRUE
validation_size <- 500 # if you want validation size to be equal to landslide invertory put: validation_size <- 0
validation_slope_threshold <- 20 # the slope threshold for areas included in the validation
#### 04: DEFINE INPUT PARAMETERS ####
if(progress_display){
  print("04: DEFINE INPUT PARAMETERS")
}
resolution <- 5            # the desired resolution for the total, if 0, no resampling
TWI_multiplier <- 1.236         # a multipication of the TWI
ls_density <- 0.5                       # the density of the randomly generated landslides; landslides/m^2
soil_density <- 1.702                         # the density of the soil material; g/(cm^3)
mean_soil_depth <- 1.346                     # the mean depth of the soil: m
sd_soil_depth <- 0.476                       # the standard deviation of the depth of the soil; m
mean_soil_cohesion <- 2.181                 # the mean unsaturated cohesion of the soil; Kpa
sd_soil_cohesion <- 0.239                    # the standard deviation of the unsaturated cohesion of the soil; Kpa
mean_phi <- 28.262                            # the mean angle of internal friction of the soil material; -
sd_phi <- 4.39                                # the standard deviation of the angle of internal friction of the soil material; -
hundred_year_precipitation <- 119.35           # The precipitation intensity occuring once every 100 years, to calibrate transmissivity; mm/h
catchment_saturated_fraction <- 0.216        # The fraction of the catchment assumed to be fully saturated at the 100 year event, for transmissivity calibration; -
instability_threshold_precipitation <- 34.915   # the precipitation intensity assumed to be the threshold of instability; mm/h
tree_uniformity_resolution <- 15.5            # the resolution over which the root system of the tree stretches; m
coniferfraction <- 0.5                   # the fraction of the trees assumed to be conifers; -
test_rainfall_intensity <- 1.448*hundred_year_precipitation              # the actual event precipitation that is being tested; mm/h
cohesion_minimum_value <- 0.512*mean_soil_cohesion             # the minimum residual cohesion of the soil, c`
probability_resolution <- resolution               # resolution for the output probability maps
validation_probability_threshold <- 20
validation_resolution <- 20
#### 05: LOADING DEM, TWI AND SLOPE ####
if(progress_display){
  print("05: LOADING DEM, TWI AND SLOPE")
}
dem <- raster(read.asc(paste(filePath, DEMFile, sep=""))) # turn the DEM into an R raster, DEM created using SAGA
TWI <- raster(read.asc(paste(filePath, TWIFile, sep=""))) # turn the TWI into an R raster, TWI created using SAGA
Slope <- raster(read.asc(paste(filePath, SlopeFile, sep=""))) # turn the slope into an R raster, slope created using SAGA

x_max <- extent(dem)[2]   
x_min <- extent(dem)[1]   
y_max <-  extent(dem)[4]  
y_min <- extent(dem)[3]   
ext <- extent(x_min, x_max, y_min, y_max)   # create the total extent of the area

if(resolution != 0){ # resample the input grids for the given input resolution, if 0, no resampling will be performed
  r <- raster(xmn = x_min, ymn = y_min, xmx = x_max, ymx = y_max, res=resolution)
  dem <- resample(dem,r,method="bilinear")
  TWI <- resample(TWI,r,method="bilinear")
  Slope <- resample(Slope,r,method="bilinear")
}
  
TWI[TWI<0] <- 0
Slope[Slope<0] <- 0
TWI <- TWI * TWI_multiplier # add a multipication to the TWI

if(noNANs){ # creates a square no Nan research area
  TWI[is.na(TWI)] <- 0
  dem[is.na(dem)] <- 0
  Slope[is.na(Slope)] <- 0
}
#### 06: PROCESSING LANDSLIDE INVENTORY #### 
if(progress_display){
  print("06: PROCESSING LANDSLIDE INVENTORY")
}
tab_all <- data.frame(read.csv(paste(filePath, landslidesFile, sep=""),header=TRUE,sep=";",dec=".",na.strings="None")) # load all the landslides

if(slides_subsetting){ # create a subset over the study area
  tab_all <- subset(tab_all, longitude < x_max & longitude > x_min & latitude < y_max & latitude > y_min)
}

tab_tot <- tab_all # select only the coloumn we need
tab_tot <- tab_tot[slidevariablenames]
colnames(tab_tot) <- c("Latitude", "Longitude", "Maechtigkeit", "Flaeche", "Neigung.Rutschflaeche") # rename those columns

tab_tot<- tab_tot[!is.na(tab_tot[,2]),]
tab_tot<- tab_tot[!is.na(tab_tot[,3]),]
tab_tot<- tab_tot[!is.na(tab_tot[,4]),]
tab_tot<- tab_tot[!is.na(tab_tot[,5]),]
tab_tot<- data.frame(tab_tot) # make sure only the inputs with real numbers are accepted and format back to data frame

data_area <- as.character(tab_tot$Flaeche)
data_area <- as.numeric(gsub("'", "", data_area))
data_area[is.na(data_area)] <- 10
tab_tot$Flaeche <- data_area
breaks <- seq(0,max(data_area)+10,10)
hist_data<-hist(data_area, breaks = breaks, plot=FALSE)
freq_data<-hist_data$counts/sum(hist_data$counts) # eventually create a distribution of the area of landslide
#### 07: GENERATION OF RANDOM DISTRIBUTED LANDSLIDE COORDINATES ####
if(progress_display){
  print("07: GENERATION OF RANDOM DISTRIBUTED LANDSLIDE COORDINATES")
}
n<- ceiling((x_max-x_min)*(y_max-y_min)*ls_density)							#number of total generated landsldies, Anzahl hypothetischer Rutschungen (16000Rutschungen/km2=0.016R/m2)
x<-  runif(n,x_min,x_max) # x coordinates of the random landslides
y<-  runif(n,y_min,y_max) # y coordinates of the random landslides
#### 08: FITTING THE LANDSLIDE AREA DISTRIBUTION FUNCTION ####
if(progress_display){
  print("08: FITTING THE LANDSLIDE AREA DISTRIBUTION FUNCTION")
}
A_L<-seq(0.00001, max(breaks)*0.000001, 0.00001)   #classes of landslide areas are in 10 m^2 == 0.00001 km^2
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

rho_gamma<-	rho_gamma_rand[which(SSE==min(SSE[!is.na(SSE)]))]  #1.5	#1.4

a<-   a_rand[which(SSE==min(SSE[!is.na(SSE)]))]          #1.28*10^(-3)
s<-   s_rand[which(SSE==min(SSE[!is.na(SSE)]))]          #0.09*10^(-4)  #1.82*10^(-4)
p_A_L<-(1/(a*gamma(rho_gamma)))*((a/(A_L-s))^(rho_gamma+1))*exp(-a/(A_L-s))  #probability distribution malamut et al., 2004. Landslide area is in square kilometers.
#### 09: CUMULATIVE DISTRIBUTION OF RANDOM LANDSLIDES USING FITTED PDF ####
if(progress_display){
  print("09: CUMULATIVE DISTRIBUTION OF RANDOM LANDSLIDES USING FITTED PDF")
}
cum_p_A_L<-p_A_L

i<-2
while(i<length(p_A_L)+1)
{
  cum_p_A_L[i]<-p_A_L[i]+cum_p_A_L[i-1]
  i<-i+1
}

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

ls_area<-1000000*A_landslide #eventually create a vector of all the landslide areas
#### 10: RANDOM LANDSLIDE GEOMETRY ####
if(progress_display){
  print("10: RANDOM LANDSLIDE GEOMETRY")
}
c<- 2       									# Quotient zwischen Breite 1 und Länge 2
ls_width<- sqrt((ls_area*4)/(pi*c))				# Rutschbreite
ls_length<-ls_width*c   						# Rutschlänge = Rutschbreite*2
#### 11: RANDOM LANDSLIDE SLOPE AND TWI  ####
if(progress_display){
  print("11: RANDOM LANDSLIDE SLOPE AND TWI")
}
coord_ls<-matrix(c(x,y), ncol=2) # combine the random x and y into a coordinate vector
coord_ls[,1]<-x
coord_ls[,2]<-y

# extraction of the slope and TWI of each generated landslide
spol <- gBuffer(SpatialPoints(coord_ls), width=ls_width, byid=TRUE)
spdf <- SpatialPolygonsDataFrame(spol, data.frame(id=1:length(spol)), FALSE)
velox_slope <- velox(Slope, extent=ext, res=c(resolution,resolution))
velox_TWI <- velox(TWI, extent=ext, res=c(resolution,resolution))
slope_ls <- as.vector(velox_slope$extract(spdf, fun=mean))
tpi_ls <- as.vector(velox_TWI$extract(spdf, fun=mean))
#### 12: DEFINE PARAMETERS: SOIL DENSITY, SOIL DEPTH, SOIL COHESION AND PHI ####
if(progress_display){
  print("12: DEFINE PARAMETERS: SOIL DENSITY, SOIL DEPTH, SOIL COHESION AND PHI")
}
dens_soil<- soil_density                   				# [t/m^3]Bodendichte [g/cm^3], Literatur, Trockendichte (sandstone 1.6-2.68, mudstone 2.4-2.7) Feldaufnahme
slope<- slope_ls  						# Hangneigung [�] 
slope_rad<-2*pi*slope/360         				# [rad]

soil_cohesion_unsaturated<- rnorm(n,mean_soil_cohesion,sd_soil_cohesion)              	# Bodenkohäsion [kPa] SC:5.3(11.6),5.8(6.7)
soil_cohesion_unsaturated[which(soil_cohesion_unsaturated<0)]<-0.01							# SC:Clayey sand, SC-SM:Clayey sand-silty sand

phi_soil<-  rnorm(n,mean_phi,sd_phi)          # slope_ls*rnorm(n,1,0.1)                  	# Reibungswinkel [°] SC:34.8(4.7), 35.3(4.7)
phi_soil_rad<- 2*pi*phi_soil/360  				# [rad]

if(random_soil_depth){ # create a random vector for all the soil depths, dependent on the slope
  soil_depth<-rnorm(n,mean_soil_depth,sd_soil_depth) * (1-pnorm(slope,1.35 *mean(data_xslope),0.75*sd(data_xslope)))		# Feldaufnahme mean(data_slope), sd(data_slope) Rutschkörpermächtigkeit/-tiefe [m], pnorm serve ad evitare che ci siano frane in zone ripide!
  soil_depth[which(soil_depth<0.01)]<-0.01
}

if(kriging_soil_depth){ # use kriging to get soil depth values from an interpolated raster
  slope_grid <-  read.asciigrid(paste(filePath, SlopeFile, sep=""), as.image=FALSE, plot.image=TRUE) # load the slope as a asciigrid
  names(slope_grid) <- "slope_clipped.asc"

  coordinates(tab_tot) =~ Longitude + Latitude # define the coordinates of the landslide inventory
  proj4string(tab_tot) <- proj4string(slope_grid) <- CRS("+init=epsg:21782") # put a geo referential on it
  
  coord_slides <- matrix(c(tab_tot$Longitude,tab_tot$Latitude), ncol=2) # create a vector of inventory landslide locations
  slope_slides <- extract(Slope, coord_slides, buffer=data_area, fun=mean) # extract the slope of these slides
  tab_tot$slope_clipped.asc = slope_slides

  lm.depth <- lm(Maechtigkeit~slope_clipped.asc, tab_tot) # create a linear fit through slope and soil depth from the inventory
  #plot(Maechtigkeit~slope_clipped.asc, as.data.frame(tab_tot))
  #abline(lm(Maechtigkeit~slope_clipped.asc, as.data.frame(tab_tot)))

  null.vgm <- vgm(var(tab_tot$Maechtigkeit), "Sph", sqrt(areaSpatialGrid(slope_grid))/4, nugget=0) # initial parameters
  vgm_depth_r <- fit.variogram(variogram(Maechtigkeit~slope_clipped.asc, tab_tot), model=null.vgm) # fit a variogram between the soil depth and slope from the inventory
  #plot(variogram(Maechtigkeit~slope_clipped.asc,tab_tot), vgm_depth_r, main="fitted by gstat")

  depth_uk <- krige(Maechtigkeit~slope_clipped.asc, locations= tab_tot, newdata=slope_grid, model=vgm_depth_r,debug.level=0) # use kriging to get a soil depth map
  soil_depth <- raster(depth_uk) # convert these soil depths to a real raster
  soil_depth[soil_depth<0.01] <- 0.01 # do a minimum cutoff of the values
  soil_depth[soil_depth>2.5] <- 2.5 # do a  maximum cutoff of the values

  spol <- gBuffer(SpatialPoints(coord_ls), width=ls_width, byid=TRUE)
  spdf <- SpatialPolygonsDataFrame(spol, data.frame(id=1:length(spol)), FALSE)
  velox_soil_depth <- velox(soil_depth, extent=ext, res=c(res(soil_depth)[1],res(soil_depth)[1]))
  depth_random_slides <- as.vector(velox_soil_depth$extract(spdf, fun=mean)) # extract random landslides soil depths from this raster

  linear_correction <- (90-slope_ls)/180 + 0.5 # define a small correction of the depth values
  randomize <- rnorm(length(depth_random_slides),1,sd(tab_tot$Maechtigkeit)^2)*linear_correction # add this correction through a randomizer, to create little but some randomness
  depth_ls <- depth_random_slides*randomize
  
  depth_ls[depth_ls>2.5] <- 2.5 # do a  maximum cutoff of the values
  depth_ls[depth_ls<0.01] <- 0.01 # do a minimum cutoff of the values
  soil_depth <- depth_ls
  
  #plot(slope_ls,soil_depth,ylim=c(0, 3))
  #points(tab_tot$slope_clipped.asc,tab_tot$Maechtigkeit, pch=19,col="orange")	
}
#### 13: ESTIMATION LANDSLIDE VOLUME & MASS  ####
if(progress_display){
  print("13: ESTIMATION LANDSLIDE VOLUME & MASS")
}
ls_volume<- ls_area*soil_depth*cos(slope_rad)					# Rutschvolumen = Rutschfläche*-Tiefe
ls_circ<- pi*(3*((ls_length/2)+(ls_width/2))-sqrt((3*(ls_length/2)+(ls_width/2))*((ls_length/2)+(3*(ls_width/2)))))  			#2*pi*ls_width         #
ls_area_circ_sup<- soil_depth*ls_circ/2			#  [m^2]
mass<- ls_volume*dens_soil						# Masse [t]
#### 14: EMPERICAL CALIBRATION OF TRANSMISSIVITY ####
if(progress_display){
  print("14: EMPERICAL CALIBRATION OF TRANSMISSIVITY")
}
mm_h_calibration<-hundred_year_precipitation #[mm/h],rainfall intensity # 1h intensi�t der ereigniss, oder statistische Wert der 100jaehrliches Ereigniss zurueck berechnen
R_calibration<-  0.001*mm_h_calibration/3600 #[m/s],rainfall intensity
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

treshehold<- length(tpi_ls)*catchment_saturated_fraction   # this factor represents the assumed percent of fully saturated soil profiles in a catchment during an extrem rainfall, for example 0.6 == 60%.
T<- T_serie[max(which(m_sum>treshehold))]   #[m/s],estimated transmissivity
#### 15: EMPERICAL CALIBRATION SOIL COHESION ####
if(progress_display){
  print("15: EMPIRICAL CALIBRATION SOIL COHESION")
}
mm_h<-instability_threshold_precipitation  #[mm/h],rainfall intensity, under which no instability is assumed
R<-  0.001*mm_h/3600 #[m/s],rainfall intensity
m<- (R/T)*tpi_ls
soil_cohesion <-(1-m)*soil_cohesion_unsaturated
water_pressure <- soil_depth*m*(100/9.81) #rnorm(n,0.5,0.1) 	# [kPa], 1kPa <=> 10 cm watertable	 # Porenwasserdruck
force_par<- (mass*9.81*sin(slope_rad))
force_per<-(mass*9.81*cos(slope_rad))
force_per_eff<-force_per-(ls_area*water_pressure)
shear_res_bas<- (ls_area*soil_cohesion)+(force_per_eff*(tan(phi_soil_rad)))
shear_res<- shear_res_bas
SF_ini<-shear_res/force_par

cohesion_min_reset<- (force_par-force_per_eff*(tan(phi_soil_rad)))/ls_area # calculate a minimum stability for the unstable slides, set them all good at certain
soil_cohesion_unsaturated <- ifelse(cohesion_min_reset>soil_cohesion_unsaturated, cohesion_min_reset/(1-m)*rnorm(length(cohesion_min_reset),1.3,0.15),soil_cohesion_unsaturated) #condition that corrects negative values of cohesion, with some minor randomness
#### 16: ROOT REINFORCEMENT  ####
if(vegetation_included){
  
  if(progress_display){
    print("16: ROOT REINFORCEMENT")
  }
  
  trees<- read.csv(paste(filePath, "Ind_trees.csv", sep=""), sep=";",header = TRUE) # load the treefile
  trees <- subset(trees, trees[,1] < x_max & trees[,1] > x_min & trees[,2] < y_max & trees[,2] > y_min) # subset only for the research area
  tree_raster_dimension <- tree_uniformity_resolution # define a resolution over which the tree is assumed to give one surcharge and one root reinforcement
  r <- raster(xmn=extent(dem)[1], ymn=extent(dem)[3], xmx=extent(dem)[2], ymx=extent(dem)[4], res=tree_raster_dimension) #define a raster with the new resolution
  r[] <- 0
  x_t <- trees[,1] # x coordinates of the trees
  y_t <- trees[,2] # y coordinates of the trees
  xy <- cbind(x_t, y_t) # vectorize the tree locations
  r_density_trees <- rasterize(xy, r, fun=function(x,...)length(x)) # rasterize the number of trees per cell
  dbh<-trees[,4]
  r_dbh <- rasterize(xy, r, field=dbh, fun=function(dbh,...)mean(dbh)) # rasterize the average dbh per cell

  
  # the possible selection of different scenarios
  if(variable_rootlat){
    r_dist<- sqrt(tree_raster_dimension*tree_raster_dimension/(pi*r_density_trees*3))
    root_reinf_values<-  2*(50*getValues(r_dbh)*dgamma(getValues(r_dist)/(0.01*getValues(r_dbh)*18.5), 5,15))   #[N], dbh in [cm], distance in [m]. Minimum root reinforcement in the cell
    r_root_reinf<-r_dist
    r_root_reinf<-setValues(r_root_reinf, root_reinf_values)
    extent(r_root_reinf)<-ext
    r_root_reinf<-setExtent(r_root_reinf, ext, keepres=T, snap=T)
    root_reinf_ls<-(extract.data(coord_ls, r_root_reinf))
    root_reinf_ls[is.na(root_reinf_ls)]<- 0
    root_lateral<- root_reinf_ls/1000       					# Laterale Wurzelverstärkung [kN/m], 2-15
  }

  if(forest_uniform_rootlat){
    root_lateral_raster <- r_density_trees
    root_lateral_raster[root_lateral_raster  > 0] = root_lateral_strength # create a boolean map with or without forest
    root_lateral<-(extract.data(coord_ls, root_lateral_raster)) # extract the uniform value of the root lateral only for the forests
    root_lateral[is.na(root_lateral)]<- 0
  }
  
  if(global_uniform_rootlat){
    root_lateral <- root_lateral_strength  # if a single value of the root lateral strength, this is defined
  } 

}else{
  root_lateral <- 0 # if no vegetation, the root lateral strength is 0
}
root_basal<- dgamma(soil_depth,3.1,12.57)*0.15*root_lateral  # determine the root basal strength based on the root lateral strength
#### 17: VEGETATION WEIGHT ####
if(vegetation_included){
  
  if(progress_display){
    print("17: VEGETATION WEIGHT")
  }

  if(global_uniform_vegweg){
    w_veg <- uniform_vegetation_weight # define one global vegetation weight
  } 

  if(forest_uniform_vegweg){
    vegweg_raster <- r_density_trees
    vegweg_raster[vegweg_raster  > 0] = uniform_vegetation_weight # create a boolean yes or no forest raster
    w_veg<-(extract.data(coord_ls, w_veg_raster)) # extract vegetation weight based on boolean raster
    w_veg[is.na(w_veg)]<- 0
  }
  
  if(variable_vegweg){
    velox_dem <- velox(dem, extent=ext, res=c(res(dem)[1],res(dem)[1])) # start to assign the elevation to each tree
    spol <- gBuffer(SpatialPoints(xy), width=res(dem)[1], byid=TRUE)
    spdf <- SpatialPolygonsDataFrame(spol, data.frame(id=1:length(spol)), FALSE)
    tree_elevations <- as.vector(velox_dem$extract(spdf, fun=mean)) # extract tree elevations
    height_boundary <- quantile(tree_elevations, probs = 1-coniferfraction,na.rm=TRUE) # define a boundary between conifer and broadleaf trees
    tree_elevations[which(tree_elevations<height_boundary)] <- 1 # the broadleaf fraction will be given 1
    tree_elevations[which(tree_elevations>height_boundary)] <- 2 # the conifer fraction will be given 2
    species_index <- tree_elevations # extract the elevation of each tree and tag lower growing as broadleaf and higher growing as coniferous
  
    # attempt to get the tree masses and weights. Resample the weights back to the general resolution
    heights <- trees[,3] #broadleaf, conifer
    NH <- c(25,20 + 1.5*dbh/10)   #broadleaf, conifer  # denzin species specific correction factor # start calculating tree weights according to Denzin
    korprozent <- c(0.03,0.04) #broadleaf, conifer  # denzin species specific correction procent
    volumes <- (dbh^2/1000) #+ (dbh^2/1000)*(heights-NH[species_index])*korprozent[species_index] # return the tree volume in m^3
    tree_densities_dry <- c(720,460)
    tree_densities <- tree_densities_dry*1.5 # very gambly!!! assumed half of the tree is water. fix this problem later
    masses <- volumes*tree_densities[species_index] # get the masses in kg
    
    r_masses <- rasterize(xy, r, field=masses, fun=function(masses,...)sum(masses))
    r_weights <- (r_masses/(tree_raster_dimension^2))/1000  # weights in tons/m2)
    r_weights <- resample(r_weights,dem,method="bilinear")
    r_weights[is.na(r_weights)] <- 0
    r_weights[r_weights<0] <- 0
    
    spol <- gBuffer(SpatialPoints(coord_ls), width=ls_width, byid=TRUE) # get the vegetation weight per random landslide, using the velox function once again
    spdf <- SpatialPolygonsDataFrame(spol, data.frame(id=1:length(spol)), FALSE)
    velox_weights <- velox(r_weights, extent=ext, res=c(res(dem)[1],res(dem)[1]))
    w_veg <- as.vector(velox_weights$extract(spdf, fun=mean)) # Gewicht Vegetation [t/m^2] Tree volume calculation -->werte in literatur suchen
    w_veg[is.na(w_veg)] <- 0
  }
}else{
  w_veg <- 0 # if no vegetation, no vegetation weight present
}
#### 18: CALCULATE PORE WATER PRESSURE  ####
if(progress_display){
  print("18: CALCULATE PORE WATER PRESSURE")
}
mm_h <- test_rainfall_intensity  # [mm/h],rainfall intensity
R<-  0.001*mm_h/3600 #[m/s],rainfall intensity
m<- (R/T)*tpi_ls       #a steady-state shallow subsurface flow model, derived from TOPMODEL assumptions (Beven and Kirkby, 1979) , Park et al..(2013)
m[which(m>1)]<- 1

soil_cohesion <-(1-m)*soil_cohesion_unsaturated
soil_cohesion[soil_cohesion < cohesion_minimum_value] <- cohesion_minimum_value

water_pressure<-  soil_depth*m*(100/9.81) #rnorm(n,0.5,0.1) 	# [kPa], 1kPa <=> 10 cm watertable					# Porenwasserdruck
#### 19: CALCULATE RANDOM LANDSLIDE STABILITY ####
if(progress_display){
  print("19: CALCULATE RANDOM LANDSLIDE STABILITY")
}
force_par<- (mass*9.81*sin(slope_rad))+(ls_area*w_veg*9.81*sin(slope_rad)) # define the driving force parallel to the slope
force_per<-(mass*9.81*cos(slope_rad))+(ls_area*w_veg*9.81*cos(slope_rad))
force_per_eff<-force_per-(ls_area*water_pressure)

shear_res_bas<- (ls_area*soil_cohesion)+(ls_area*root_basal)+(force_per_eff*(tan(phi_soil_rad)))
shear_res_lateral<- (ls_circ/2)* root_lateral    #ls_area_circ_sup*root_lateral
shear_res<- shear_res_bas+shear_res_lateral # define the resisting normal force
SF<-shear_res/force_par
SF <- as.vector(SF) # create a final safety factor
#### 20: PLOT UNSTABLE LANDSLIDES ####
if(progress_display){
  print("20: PLOT UNSTABLE LANDSLIDES")
}
n.veg<-matrix(c(seq(1,n,1),x,y,SF,ls_area,soil_depth, slope),ncol=7) # create the total landslides matrix
n.veg_filtered<-n.veg[n.veg[,4]<1, ] # filter the total matrix to unstable slides
n.veg_filtered<-n.veg_filtered[!is.na(n.veg_filtered[,4]),] # filter this of Nans
#### 21: RASTER FOR LANDSLIDE OCCURENCE #### 
if(progress_display){
  print("21: RASTER FOR LANDSLIDE OCCURENCE")
}
r_ls <- raster(xmn=min(x), ymn=min(y), xmx=max(x), ymx=max(y), res=probability_resolution) # create a base raster for probability
r_ls[] <- 0

x_ls <- n.veg_filtered[,2]      # ID,x,y,SF,ls_area,soil_depth,slope
y_ls <- n.veg_filtered[,3]
xy_ls <- cbind(x_ls, y_ls) #create a vector for unstable slides

if(length(xy_ls) < 1){
  xy_ls <- cbind(c(NA),c(NA))
}

r_number_ls <- rasterize(xy_ls, r_ls, fun=function(x,...)length(x)) # plot the number of unstable landslides per cell
#### 22: RASTER FOR LANDSLIDE PROBABILITY ####
if(progress_display){
  print("22: RASTER FOR LANDSLIDE PROBABILITY")
}

pdf(paste(outPath,"ls_probability_and_ls.pdf",sep="")) # save the landslide probability as pdf
plot(Slope,legend=F, col=gray.colors(10, start = 0.3, end = 0.9, gamma = 2.2, alpha = NULL))
plot(100*(r_number_ls)/(6*ls_density*(probability_resolution^2)), add=T, col=c("orange1","orange2","orange3","orange4","green1","green2","green3","green4"))

if(slides_subsetting){
  points(tab_tot$Longitude,tab_tot$Latitude, pch=19, cex=2*data_area/800, col='#FF000060', lwd=4) # plot the actual landslides
}
dev.off()

probability <- 100*(r_number_ls)/(6*ls_density*(probability_resolution^2)) # the function for landslide probability?
probability[is.na(probability[])] <- 0 # define minimum boundary of this probability function
probability[probability>100] <- 100  # define maximum boundary of this probability function
probabilities <- getValues(probability)
writeRaster(probability, paste(outPath,"ls_probability"), overwrite=TRUE, format = "ascii")
#### 23: VALIDATION ####
if(validation){
  
  if(progress_display){
    print("23: VALIDATION")
  }
  
  mod <- ifelse(probabilities < validation_probability_threshold, 0, 1) # create a raster of all unstable defined areas based on threshold
  mod[is.na(mod)]<-0
  
  data_x<-c(tab_tot$Longitude) #landslide inventory
  data_y<-c(tab_tot$Latitude)
  n.inv<-matrix(c(data_x, data_y),ncol=2)
  
  r_i <- raster(xmn=min(x), ymn=min(y), xmx=max(x), ymx=max(y), res=probability_resolution) # create a base raster for validation
  r_i[] <- 0
  
  x_i <- n.inv[,1]
  y_i <- n.inv[,2]
  xy_i <- cbind(x_i, y_i)

  if(auc_validation){ # use the auc validation
    r_number_i <- rasterize(xy_i, r_i, fun=function(x,...)length(x))
    inv<- getValues(r_number_i)
    inv[is.na(inv)]<-0
    category <- inv #Inventory
    prediction <- mod #model
    roc_obj <- roc(category, prediction)
    auc <- roc_obj$auc #area under the curve
    print(auc)
  }
  
  if(area_validation){ # a validation method using the difference in average probability over the existing landslides and the whole map
    slide_radius <- sqrt(data_area/pi)
    velox_probability <- velox(probability, extent=ext, res=c(validation_resolution,validation_resolution)) # rasteriye the probability to the valiation resolution
    spol <- gBuffer(SpatialPoints(n.inv), width=slide_radius*2, byid=TRUE)
    spdf <- SpatialPolygonsDataFrame(spol, data.frame(id=1:length(spol)), FALSE)
    inv_stab <- as.vector(velox_probability$extract(spdf, fun=mean)) # extract probabilities over the landslide areas
    inv_stab <- inv_stab[!is.na(inv_stab)]
    
    if(validation_size == 0){
      validation_size <- length(data_x) # if no size of the validation plots is given (=0), then the size of the inventory will be used
    }
  
    for (a in 1:10000){ # from 1 to a lot, to create a list of random landslide sized plof on relatively steep slopes
      xs <- runif(a, x_min, x_max)
      ys <- runif(a, y_min, y_max)
      ran_areas <- data.frame(xs,ys)
      Slopes_areas <- extract(Slope,ran_areas)
      ran_areas["Slopes"] <- Slopes_areas
      ran_areas <- ran_areas[!is.na(ran_areas[,3]),]
      
      if(nrow(ran_areas[ran_areas[,3]>validation_slope_threshold, ]) >= validation_size)
        break
    }
    
    ran_areas <- ran_areas[ran_areas[,3] > validation_slope_threshold,]
    ran_areas <- head(ran_areas,validation_size)
    areas <- cbind(ran_areas$xs,ran_areas$ys)
  
    radiusses <- rep_len(slide_radius, length.out=validation_size)
    spol <- gBuffer(SpatialPoints(areas), width=radiusses*2, byid=TRUE)
    spdf <- SpatialPolygonsDataFrame(spol, data.frame(id=1:length(spol)), FALSE)
    ran_stab <- as.vector(velox_probability$extract(spdf, fun=mean))
    ran_stab <- ran_stab[!is.na(ran_stab)]
  
    difference <- mean(inv_stab) - mean(ran_stab)
    print(difference)
  }
}


 

