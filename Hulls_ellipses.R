######################################################################################################################
#ISOTOPIC NICHE METRICS - CONVEX HULLS AND STANDARD ELLIPSES
######################################################################################################################

#This script is intended as material for the practical part of the course "Etude des isotopes stables et applications 
#au milieu marin", taught by Gilles Lepoint & Loïc Michel at University of Liège.

#Most of this script is based on Andrew Jackson's example scripts for his package SIBER. 
#More info at http://www.tcd.ie/Zoology/research/groups/jackson/
#If you are looking for more example scripts and guidance to use the package, 
#please visit https://github.com/AndrewLJackson/SIBER and https://cran.r-project.org/web/packages/SIBER/vignettes/.

#Prerequisites: Install the current version of R: https://cran.r-project.org/, 
#of RStudio: https://www.rstudio.com/products/rstudio/download/, 
#and of JAGS: http://mcmc-jags.sourceforge.net/.
#Create a Rstudio project associated with your working directory, and make sure your work will be saved automatically.

#Now let's get started!


######################################################################################################################
#PART 1 - CONVEX HULL METRICS ON COMMUNITIES
######################################################################################################################

#Let's install the SIBER package from the CRAN repository, and load it into our session's environment.
install.packages("SIBER", dependencies = TRUE)
library(SIBER)

#Let's import the data, that are stored in the provided CSV file.
fulldata <- read.csv("SIBER_full.csv", header=T)

#Create a key for the community and group codes: it will be useful for quick reference later
key <- read.csv("Key.csv", header=F)

#Now, let's create a SIBER object that can be manipulated by the package.
siber.full <- createSiberObject(fulldata)

#Use SIBER to take a look at the data (colours = groups, symbols = communities)
par(mar=c(5.1, 5.1, 4.1, 2.1))
plotSiberObject(siber.full,
                  ax.pad = 2, 
                  hulls = F, 
                  ellipses = F,
                  group.hulls = F,
                  bty = "o",
                  iso.order = c(1,2),
                  xlab = expression({delta}^13*C~'(\u2030)'),
                  ylab = expression({delta}^15*N~'(\u2030)')
                  )

#Add community hulls, computed using the means of every consumer group
plotCommunityHulls(siber.full, plot.args = list(col = "black", lty = 2))

#Those graphs are good to take a look at your data. It's probably not your best option to produce nice figures.
#To do that, we could use R's base functions.

#First, let's calculate the means for each species in each community. To do that, we have to split our data frame in
#two, one for each community.
SplitComm <- split(fulldata, fulldata$community)
Comm1 <- SplitComm$`1`
Comm2 <- SplitComm$`2`

#Now we can calculate the meansby group for each isotopic ratio and in each community.
means.x.comm1 <- aggregate(Comm1$iso1,list(Comm1$group),mean)$x
means.y.comm1 <- aggregate(Comm1$iso2,list(Comm1$group),mean)$x
means.x.comm2 <- aggregate(Comm2$iso1,list(Comm2$group),mean)$x
means.y.comm2 <- aggregate(Comm2$iso2,list(Comm2$group),mean)$x

#Store those means in data frames
Comm1Means <- cbind(means.x.comm1,means.y.comm1)
Comm2Means <- cbind(means.x.comm2,means.y.comm2)

#Compute the parameters of the convex hull of community 1
Hull1 <- chull(Comm1Means)
Hull1 <- c(Hull1, Hull1[1])

#Compute the parameters of the convex hull of community 2
Hull2 <- chull(Comm2Means)
Hull2 <- c(Hull2, Hull2[1])

#Plot the species means for community 1. You can easily customise this plot, see ?plot for details.
plot(Comm1Means, 
     type = "p", pch = 16, col = "black", 
     xlim = c(-26,-12), ylim = c(0,15),
     xlab = expression({delta}^13*C~'(\u2030)'),
     ylab = expression({delta}^15*N~'(\u2030)'),
     main=("Communities convex hulls")
)

#Plot the species means for community 2.
points(Comm2Means, 
       type = "p", pch = 16, col = "red")

#Plot the convex hull for community 1
lines(Comm1Means[Hull1, ], col="black", lty="solid", lwd=1.5)

#Plot the convex hull for community 2
lines(Comm2Means[Hull2, ], col="red", lty="solid", lwd=1.5)

#Add a legend
legend("bottomleft", legend=c("Anse du Lion", "Cap des Elephants"), col=c("black","red"), 
       pch=16, lty="solid", merge=TRUE, bty="n")

#Now you have a nice-looking graph that you can export to illustrate documents. Time to calculate some metrics.

#Calculate the various Layman metrics on each of the communities.
Layman.full.calc <- communityMetricsML(siber.full) 

#See those metrics in the R console
print(Layman.full.calc)

#As explained during the theoretical course, we can also use a Bayesian model to estimate these values.
#To do that, we first need to set the options for running JAGS. Details are outside the scope of this course.
parms <- list()
parms$n.iter <- 2 * 10^4   # Number of model iterations. Here, we will keep it low to limit computation time.
parms$n.burnin <- 1 * 10^3 # Number of initial discarded values
parms$n.thin <- 10     # Interval to thin the posterior
parms$n.chains <- 2        # NUmber of chains to run

# We also need to define the priors. Once again, details are outside the scope of this course.
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

#Now we will run the model, using the parameters and priors we just specified.
#This might take a while, and the run length will depend on the computer you use.
ellipses.post.full <- siberMVN(siber.full, parms, priors)

#Now we can extract the posterior means from our model, and calculate the corresponding distribution of layman 
#metrics. Once again, this takes a while.
mu.post <- extractPosteriorMeans(siber.full, ellipses.post.full)
layman.full.bayes <- bayesianLayman(mu.post)

#Let's produce density plots of the metrics to compare them across communities. Given the ranges that can be very
#different, it's more efficient to do it separately for each metric. Let's try it with TA (total area of the 
#convex hull).
par(mgp=c(3,0.6,0))
siberDensityPlot(cbind(layman.full.bayes[[1]][,"TA"], layman.full.bayes[[2]][,"TA"]),
                 xticklabels = c("Anse du Lion", "Cap des Elephants"), 
                 bty="L", ylim = c(0,60),
                 las = 1,
                 ylab = "Total area of the convex hull",
                 xlab = "Community")

#On this plot, by default, the dark, median, and light grey are the 50, 75 and 95% credibility intervals, and the
#black dot is the mode. We can also add the geometrically calculated value to compare our two methods.
points(1:2, Layman.full.calc[3,], col = "red", pch = 16)

#If we want, we can do that with other metrics. Here is an example with the standard deviation of nearest neighbour
#distance ("SDNND").
siberDensityPlot(cbind(layman.full.bayes[[1]][,"SDNND"], layman.full.bayes[[2]][,"SDNND"]),
                 xticklabels = c("Anse du Lion", "Cap des Elephants"), 
                 bty="L", ylim = c(0,1),
                 las = 1,
                 ylab = "SDNND",
                 xlab = "Community")
points(1:2, Layman.full.calc[6,], col = "red", pch = 16)

#Finally, we can use our model to compute the probability that Layman metrics' values are different between 
#communities. For TA: 
Prob.diff.TA <- sum(layman.full.bayes[[1]][,"TA"]>layman.full.bayes[[2]][,"TA"])/ NROW(layman.full.bayes[[2]][,"TA"])

#Here, this probability is inferior to 95%. It means that there are less than 95% of model runs where the TA of 
#community 1 is bigger than the one of community 2. If we do an analogy with traditional hypothesis test, it means
#that TA should be considered similar in the two communities.

#Now let's compare another metric, SDNND.
Prob.diff.SDNND <- sum(layman.full.bayes[[1]][,"SDNND"]>layman.full.bayes[[2]][,"SDNND"])/ 
                          NROW(layman.full.bayes[[2]][,"SDNND"])

#The probability is around 66%, which is much lower than for TA. Once again, values of this metric can be considered
#similar in the two communities.



######################################################################################################################
#PART 2 - STANDARD ELLIPSE METRICS ON POPULATIONS
######################################################################################################################

#For this part, and to keep it simple, we will only work with two species: the sea stars Diplasterias brucei and
#Odontaster validus. To do so, we will take only a subset of the full original data.
#First, let's split it by community
SplitComm <- split(fulldata, fulldata$community)
Comm1 <- SplitComm$`1`
Comm2 <- SplitComm$`2`

#Then, let's split those by groups.
SplitGroup1 <- split(Comm1, Comm1$group)
SplitGroup2 <- split(Comm2, Comm2$group)

#Let's extract the values for S. neumayeri (group 20) and O. validus (group 14).
DipBru1 <- SplitGroup1$`5`
DipBru2 <- SplitGroup2$`5`
OdoVal1 <- SplitGroup1$`14`
OdoVal2 <- SplitGroup2$`14`

#Collate those 4 data frames together.
DipOdo <- rbind(DipBru1,DipBru2,OdoVal1,OdoVal2)

#Now we can create our SIBER object.
siber.dipodo <- createSiberObject(DipOdo)

#And quickly visualize it.
par(mar=c(5.1, 5.1, 4.1, 2.1))
group.ellipses.args  <- list(n = 100, lty = 1, lwd = 2)
plotSiberObject(siber.dipodo,
                ax.pad = 2, 
                hulls = F, 
                ellipses = T,
                group.hulls = F,
                bty = "o",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'(\u2030)'),
                ylab = expression({delta}^15*N~'(\u2030)')
                )

#Like before, SIBER itself is not the best option to customize graphs. Instead, we can do it using base R, then
#add the ellipses.

#First, let's plot the D. brucei of community 1 as black dots.
plot(DipBru1$iso1, DipBru1$iso2, type="p", pch = 16, col = "black", 
     xlim = c(-19,-12), ylim = c(5,12),
     xlab = expression({delta}^13*C~'(\u2030)'),
     ylab = expression({delta}^15*N~'(\u2030)'),
     main=("Population ellipses")
      )

#Now, let's add the D. brucei of community 2 as red dots.
points(DipBru2$iso1, DipBru2$iso2, type="p", pch = 16, col = "red")

#Now for O. validus of communities 1 and 2, as black and red triangles.
points(OdoVal1$iso1, OdoVal1$iso2, type="p", pch = 17, col = "black")
points(OdoVal2$iso1, OdoVal2$iso2, type="p", pch = 17, col = "red")

#And let's add a legend.
legend("topleft", 
       legend=c("DB Comm1", "DB Comm 2","OV Comm 1", "OV Comm 2"), 
       col=c("black","red", "black","red"), 
       pch=c(16,16,17,17), bty="n")

#Now we can add the ellipse of D. brucei, community 1 as a solid black line.
ellipse1 <- addEllipse(siber.dipodo$ML.mu[[1]][ , , 1],
           siber.dipodo$ML.cov[[1]][ , , 1],
           m = NULL,
           n = 100,
           p.interval = NULL,
           ci.mean = FALSE,
           col = "black",
           lty = 1,
           lwd = 2)

#Add the ellipse of D. brucei, community 2 as a solid red line.
ellipse2 <- addEllipse(siber.dipodo$ML.mu[[2]][ , , 1],
           siber.dipodo$ML.cov[[2]][ , , 1],
           m = NULL,
           n = 100,
           p.interval = NULL,
           ci.mean = FALSE,
           col = "red",
           lty = 1,
           lwd = 2)

#Add the ellipse of O.validus, community 1 as a dashed black line.
ellipse3 <- addEllipse(siber.dipodo$ML.mu[[1]][ , , 2],
           siber.dipodo$ML.cov[[1]][ , , 2],
           m = NULL,
           n = 100,
           p.interval = NULL,
           ci.mean = FALSE,
           col = "black",
           lty = 2,
           lwd = 2)

#Add the ellipse of O.validus, community 2 as a dashed red line.
ellipse4 <- addEllipse(siber.dipodo$ML.mu[[2]][ , , 2],
           siber.dipodo$ML.cov[[2]][ , , 2],
           m = NULL,
           n = 100,
           p.interval = NULL,
           ci.mean = FALSE,
           col = "red",
           lty = 2,
           lwd = 2)

#Now that we have a pretty graph, let's have a look at the standard ellipse areas.
group.ML <- groupMetricsML(siber.dipodo)
print(group.ML)

#We will now use a Bayesian model to estimate SEA. To do that, we need parameters and priors like in the
#part 1 of the script. We will use the ones we defined earlier and run the model
ellipses.post.dipodo <- siberMVN(siber.dipodo, parms, priors)

#Now we can sample our posterior distribution.
SEA.B <- siberEllipses(ellipses.post.dipodo)

#Let's plot our model's output
siberDensityPlot(SEA.B, xticklabels = colnames(group.ML), 
                 xlab = c("Community | Group"),
                 ylab = expression("Standard Ellipse Area  " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "Model estimations of SEA",
                 ylims = c(0,9)
                  )

#And add red dots for the SEA computed earlier from the covariance matrix
points(1:4, group.ML[3,], col = "red", pch = 16)

#Now, let's compare SEA of D. brucei in the two communities.
Prob.diff.SEA.DB <- sum(SEA.B[,1]>SEA.B[,3])/ NROW(SEA.B[,1])

#Same thing but for O. validus. We can adapt this line to compare any pair of groups.
Prob.diff.SEA.OV <- sum(SEA.B[,2]>SEA.B[,4])/ NROW(SEA.B[,2])

#Another interesting thing to do would be to compare the overlaps between ellipses of both species in each community
#to estimate to which extent they share resources. Let's do that first for community 1:
sea.overlap.comm1 <- maxLikOverlap("1.5", "1.14", siber.dipodo, p.interval = NULL, n = 100)
sea.overlap.comm1

#We can also estimate this overlap in a relative way (i.e. as a percentage of the total niche area of the two species,
#combined)
sea.overlap.comm1.rel <- sea.overlap.comm1[3]/(sea.overlap.comm1[1]+sea.overlap.comm1[2]-sea.overlap.comm1[3])
sea.overlap.comm1.rel

#Now let's do the same thing for community 2.
sea.overlap.comm2 <- maxLikOverlap("2.5", "2.14", siber.dipodo, p.interval = NULL, n = 100)
sea.overlap.comm2
sea.overlap.comm2.rel <- sea.overlap.comm2[3]/(sea.overlap.comm2[1]+sea.overlap.comm2[2]-sea.overlap.comm2[3])
sea.overlap.comm2.rel

#Finally, we can use our model to generate Bayesian estimations of overlap, and compare the distributions to see how
#probable it is that niche overlap between the two sea stars is bigger in community 2 than in community 1.
#To do that, we first compute Bayesian overlaps for community 1.
bayes.overlap1 <- bayesianOverlap("1.5", "1.14", ellipses.post.dipodo, p.interval= NULL, draws=100, n = 100)

#Then for community 2.
bayes.overlap2 <- bayesianOverlap("2.5", "2.14", ellipses.post.dipodo, p.interval= NULL, draws=100, n = 100)

#And finally we can compare the two posterior distributions to compute our probability.
Prob.diff.overlap <- sum(bayes.overlap1[,3]>bayes.overlap2[,3])/ NROW(bayes.overlap1[,3])

#That will be it for today! Any questions?

######################################################################################################################
#END OF SCRIPT
######################################################################################################################