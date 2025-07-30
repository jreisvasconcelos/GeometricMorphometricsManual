library(readxl)
library(StereoMorph)
library(geomorph)
library(Morpho)
library(ggplot2)


setwd("your/path/to/GeometricManual")

# Load biological dataset
BDTrac <- read_excel("Data_Trachurus_trachurus.xlsx")

# Preview structure
str(BDTrac)

# Format variables appropriately
BDTrac[1:3] <- lapply(BDTrac[1:3], as.factor) #if you want to change more than one at once
BDTrac[4:6] <- lapply(BDTrac[4:6], as.numeric)


# Digitalized landmarks
digitizeImage(image.file='Fish', shapes.file='Shapes',
              landmarks.ref='landmarks.txt')


#Read all shape files from the output folder
read_shapes <- readShapes('Shapes')


# Write 2D landmark pixel coordinates from multiple shape file to TPS
writeLMToTPS('Shapes', 'Shapes.tps', in.pixels = TRUE, 
             flip.y = TRUE, flip.x = FALSE)

#read tps file (geomorph R package)
myData<-readland.tps("Shapes.tps", specID = "ID", readcurves = FALSE, warnmsg = TRUE) #The readcurves argument is set to FALSE when we do not have any semilandmarks

## Generalized Procrustes Analysis—aligning the landmark coordinates so that data are comparable across specimens.
myGPA<-gpagen(myData, PrinAxes = FALSE, ProcD = TRUE)
summary(myGPA) #WE now have new coordinates which are the new standardized positions of each of the landmarks that you placed on our photos.

#see your new landmark coordinates for each specimen, after translating, scaling, and rotating
myGPA$coords

#We now create a data frame combining shape, centroid size, and biological metadata (e.g. sampling site and total length):
Mdf <- geomorph.data.frame(shape = myGPA$coords, cs=myGPA$Csize,
                           site=BDTrac$Capture_Location, Csize=log(myGPA$Csize), size= BDTrac$TL_cm, Size=log(BDTrac$TL_cm))


#To visualize how different specimens are from each other

#If you want to compare each specimens’ shapes to the average (i.e. consensus) shape, use mean=TRUE
plotAllSpecimens(Mdf$shape,mean=TRUE) 

#or use mean=FALSE if you don’t want to see the consensus shape
plotAllSpecimens(Mdf$shape,mean=FALSE) 


#Histogram of raw centroid size
hist(myGPA$Csize,xlab="Centroid Size",ylab="Number of fish",main="")

#Histogram of log-transformed centroid size.
hist(log(myGPA$Csize),xlab="Log transformed centroid size",
     ylab="Number of fish",main="")


# ANOVA to detect if there is an effect of the centroid size on proscrustes coordinates
factorsmodelM<-procD.lm(Mdf$shape~Mdf$Csize,iter=999) #Note here that we log transformed the centroid size. It is usually necessary to do so because as specimens get larger, the possibility for fluctuations in absolute body size increases (as how weight of elephants varies more than does the weight of mice). By log transforming the data, all data values are reduced as a percentage of themselves, and so large data values are reduced more than are small data values.
anova(modelM)

#Optional: Removing Allometric EffectsmodelM1Allometry<-procD.lm(Mdf$shape~Mdf$Csize, logsz = TRUE, iter=499, print.progress = FALSE)
summary(modelM1Allometry)
plot(modelM1Allometry)
#New ANOVA
modelM1ANOVA<-procD.lm(shape~Csize, data=Mdf, RRPP= TRUE, print.progress = FALSE, SS.type = "I")
modelM1NULL<- procD.lm(shape~1, data=Mdf, RPP= TRUE, print.progress = FALSE, SS.type = "I")
summary(modelM1ANOVA)
summary(modelM1NULL)

shape.resid1 <- arrayspecs(modelM1ANOVA$residuals,
                           p=dim(myGPAM$coords)[1], k=dim(myGPAM$coords)[2]) # size-adjusted 
shape.mean1 <- arrayspecs(modelM1NULL$fitted, p=dim(myGPAM$coords)[1], k=dim(myGPAM$coords)[2])

adj.shape1 <- shape.resid1 + array(myGPAM$consensus, dim(shape.resid1))
adj.shape1 <- shape.mean1 + shape.resid1

Mdf2 <- geomorph.data.frame(shape = adj.shape1, cs=myGPAM$Csize,
                            site=BDTrac$Capture_Location, Csize=log(myGPAM$Csize), size=BDTrac$TL_cm, Size=log(BDTrac$TL_cm))


###### As we don't need to remove 
###### ANOVA to see if there are differences in shape between sites ########
modelInteraction<-procD.lm(Mdf$shape~Mdf$site,iter=999)
anova(modelInteraction)


## Perform PCA on GPA-aligned shape data
PCA<-gm.prcomp(Mdf$shape) #residuals
summary(PCA)

# Merge scores with biological data
pc_scores<- PCA$x
pc_scores1<-cbind(BDTrac, pc_scores)
str(pc_scores1)

# Define a named vector mapping Capture_Location levels to colors
location_colors <- c("Coruna" = "#B57BA2","Gulf of Cadiz" = "#78A8CE")

# Calculate group centroids
centroids <- pc_scores1 %>%
  group_by(Capture_Location) %>%
  summarize(Comp1 = mean(Comp1), Comp2 = mean(Comp2))
centroids<-as.data.frame(centroids)

# Base R plot
plot(
  PCA,
  pch = 21,
  col = location_colors[as.character(BDTrac$Capture_Location)],
  cex = 2,
  cex.lab = 1.5,
  cex.axis = 1.5,
  alpha = 0.5
)
points(
  centroids$Comp1,
  centroids$Comp2,
  pch = 16,
  col = location_colors[names(table(BDTrac$Capture_Location))],
  cex = 3
)
legend(
  "topright",
  legend = unique(BDTrac$Capture_Location),
  col = unique(location_colors),
  pch = 16,
  cex = 1.5,
  title = NULL
)

# ggplot2 version
ggplot(pc_scores1, aes(x=Comp1, y=Comp2)) + 
  geom_point(aes(color=Capture_Location, shape=Capture_Location), size=7,shape = 16, alpha = 0.5) +
  geom_point(data = centroids, aes(x = Comp1, y = Comp2, color = as.character(Capture_Location)), shape = 16, size = 8, stroke = 3) +  # Add centroids
  xlab("PC1 (50.18%)") + ylab("PC2 (15.73%)") + 
  theme_bw() +
  scale_shape_manual(values = c(0, 1, 2)) + 
  scale_color_manual(values = c("#B57BA2", "#78A8CE")) + 
  geom_hline(yintercept=0, lwd= 0.25)+ geom_vline(xintercept=0, lwd=0.25) + 
  theme(axis.line= element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        axis.title =element_text(family="Arial", face="bold", size=18), 
        axis.text = element_text(family="Times New Roman", size=18),
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") 


# CANONICAL VARIATE ANALYSIS (to classify individuals to their original sample) 
library(candisc)

mod<-lm(cbind(Comp1, Comp2, Comp3, Comp4, Comp5, Comp6, Comp7) ~ Capture_Location, data=pc_scores1)
anova(mod) # Check degrees of freedom. Use only the first few PCs that explain most of the shape variation (e.g., 90–95% of total variance)                                    

output.candis<-candisc(mod, term = "Capture_Location")
output.candis
summary(output.candis)
plot(output.candis)

# Change color per locality
plot(output.candis,
     cex.axis = 1.7,
     cex.lab = 1.7,
     col = location_colors)


###Cluster and classify #######

### CVA con Morpho package
PCAdata<-data.frame(pc_scores1)
groups<-as.factor(PCAdata$Capture_Location)


cv<-CVA(PCAdata[7:20], PCAdata$Capture_Location, cv=TRUE, prior = c(0.3,0.7)) # cv=TRUE Leave-One-Out Cross-Validation (LOOCV); if not TRUE then a Jackknife Cross-Validation is applied
cv


# Thin-Plate Spline Visualisation
library(Morpho)

RW<-relWarps(Mdf$shape, scale = TRUE, CSinit = TRUE, alpha = 0,
             orp = TRUE, noalign = TRUE)
RW


######### To compare shape extremes or group means ########

# Example without deformation grid
deformGrid2d(
  RW$mshape,
  PCA$shapes$shapes.comp2$min,
  ngrid = 0,
  pch = 19,
  cex1 = 2,
  cex2 = 2,
  col1 = "red",
  col2 = 4,
  gridcol = "black"
)
lineplot(
  RW$mshape,
  point = c(1, 10:17, 2:7, 1, 7:10),
  col = "red",
  lwd = 2
)
lineplot(
  PCA$shapes$shapes.comp2$min,
  point = c(1, 10:17, 2:7, 1, 7:10),
  col = 4,
  lwd = 2
)

#Save image
filename <- "FishMean.png"
png(filename, width = 1000, height = 746, bg = "transparent")
par(cex = 1.4) 
deformGrid2d(RW$mshape,PCA$shapes$shapes.comp2$min,ngrid=0,pch=19, cex1=2, cex2=2, col1 = "red", col2 = 4,  gridcol = "black")
lineplot(RW$mshape,point=c(1,10:17, 2:7,1,7:10),col="red", lwd = 2)
lineplot(PCA$shapes$shapes.comp2$min,point=c(1,10:17, 2:7,1,7:10),col=4, lwd = 2)
dev.off()

# Example with deformation grid
deformGrid2d(
  RW$mshape,
  PCA$shapes$shapes.comp2$min,
  ngrid = 10,
  pch = 19,
  cex1 = 2,
  cex2 = 2,
  col1 = "red",
  col2 = "#78A8CE",
  gridcol = "black"
)
lineplot(
  RW$mshape,
  point = c(1, 10:17, 2:7, 1, 7:10),
  col = "red",
  lwd = 2
)
lineplot(
  PCA$shapes$shapes.comp2$min,
  point = c(1, 10:17, 2:7, 1, 7:10),
  col = "#78A8CE",
  lwd = 2
)
