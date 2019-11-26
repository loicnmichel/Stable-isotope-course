######################################################################################################################
#STABLE ISOTOPE MIXING MODELS
######################################################################################################################

#This script is intended as material for the practical part of the course "Etude des isotopes stables et applications 
#au milieu marin", taught by Gilles Lepoint & Loïc Michel at University of Liège.

#Part 1 of this script is based on Andrew Parnell's example scripts for his package simmr. 
#If you are looking for more example scripts and guidance to use the package, 
#please visit https://github.com/andrewcparnell/simmr and 
#https://cran.r-project.org/web/packages/simmr/vignettes/.

#Part 2 of this script is based on Brian Stock's example scripts for his package MixSIAR. 
#If you are looking for more example scripts and guidance to use the package, 
#please visit https://github.com/brianstock/MixSIAR and 
#https://cran.r-project.org/web/packages/MixSIAR/vignettes/.

#Prerequisites: Install the current version of R: https://cran.r-project.org/, 
#of RStudio: https://www.rstudio.com/products/rstudio/download/, 
#and of JAGS: http://mcmc-jags.sourceforge.net/.
#Create a Rstudio project associated with your working directory, and make sure your work will be saved automatically.

#Now let's get started!


######################################################################################################################
#PART 1 - SIMMR
######################################################################################################################

#Install the simmr package, as well as ggplot2, a graphics package that we will use later in the script.
install.packages('simmr', dependencies = TRUE)
install.packages('ggplot2')

#Load the installed packages
library(simmr)
library(ggplot2)

#Read in the three data files: consumers, sources, and trophic enrichment factors.
consumers <- read.csv("Consumers.csv",header=TRUE)
sources <- read.csv("Sources.csv", header=TRUE)
TEFs <- read.csv("TEFs.csv",header=TRUE)

#Map your data so that simmr can use it. Here, we have 10 consumer species in two stations (hence 20 groups in total),
#2 isotope tracers, 4 sources given as means ± SDs, and we will use Mc Cutchan et al.'s TEFs for all aquatic consumers.
#We don' have any concentration dependencies (elemental contents).
simmr_input = simmr_load(mixtures=as.matrix(consumers[,6:7]),
                         group=as.character(consumers$Group),
                        source_names=as.character(sources$Source),
                        source_means=cbind(sources$Meand13C,sources$Meand15N),
                        source_sds=cbind(sources$SDd13C, sources$SDd15N),
                        correction_means=cbind(TEFs$Meand13C,TEFs$Meand15N),
                        correction_sds=cbind(TEFs$SDd13C,TEFs$SDd15N),
                        concentration_means = NULL)

#Plot your data: this is important not only to check that the data mapping is correct, but also that your model makes
#sense, i.e. that all data are in the mixing polygon.
plot(simmr_input, group=1:20, xlab=expression(paste(delta^13, "C (\u2030)",sep="")), 
     ylab=expression(paste(delta^15, "N (\u2030)",sep=""))) 

#Then you can run the model. This line performs a short and quick model run.
simmr_output_short = simmr_mcmc(simmr_input,  mcmc_control = list(iter = 10000, burn = 1000, thin = 10, n.chain = 2))

#Once the model run is done (i.e. when the progress bar is completed), you can save your model output as a .rdata
#object. This way, you will be able to access in a later session, or on a different computer. This is a good solution
#both for ease of use and reproducibility.
save(simmr_output_short, file="simmr_short.rdata")

#Should you wish to load the stored rdata file in your working environment, you can do so using
#load(file = "simmr_short.rdata")

#Here, we will use a model based on longer MCMC runs, to ensure that convergence will be reached and that model
#performance will be satisfactory.

#You could run this model and save its results using
#simmr_output = simmr_mcmc(simmr_input,  mcmc_control = list(iter = 100000, burn = 10000, thin = 100, n.chain = 4))
#save(simmr_output, file="simmr_long.rdata")

#However, the model takes a while to run (about 15 minutes on my computer). To avoid waiting for it to complete, you
#can use the .rdata file that is provided alongside this script.
load(file = "simmr_long.rdata")

#Before to explore our model's results, let's check that it performed adequately. First, let's look at convergence
#diagnostics for all consumer groups, i.e. how well the multiple iterative chains we ran converge.
#Gelman-Rubins's diagnostics values should be as close to 1 as possible. Values over 1.1 suggest that posterior
#distributions returned by each chain diverge. In this case, try longer MCMC runs.
summary(simmr_output,type="diagnostics", group=c(1:20))

#We can also plot the posterior predictive distribution, and compare them to the actual data points. This allows us to
#check how the produced model fits to your original data. Ideally, data points should lie within the predictive
#returned by your model (by default, 50%). In simmr, this has to be done separately for all consumer groups. 
#For example, for group 6 (Polychaetes Harmothoe sp. from "Anse du Lion" station):
posterior_predictive(simmr_output, group = 6, prob=0.5)

#And for group 12 (Gastropods Margarella sp. from "Anse du Lion" station):
posterior_predictive(simmr_output, group = 12, prob=0.5)

#Another useful info to evaluate model performance is to look at correlation between posterior probabilites. This can
#be done quickly for multiple groups using matrix plots. Large (close to 1) negative correlations suggest that your
#model has trouble separating two of the sources (their composition could be too close, for example), so uses either
#one or the other to produce a suitable mixture. Large positive correlations occur when your model commonly associates 
#the two sources together, which can happen with complex polygons where competing resource associations can produce
#similar mixture values. Either way, it suggests that your model struggles to determine the sources contributions' 
#efficiently. Generally speaking, lower correlations are better.
plot(simmr_output,type='matrix', group=c(1:20))

#Now that we're confident that our model performed well, let's explore its results. First, we can do summary statistics
#to look at either means and SDs or quantiles of our model estimates for the groups we're interested in (here, all).
summary(simmr_output,type="statistics", group=c(1:20))
summary(simmr_output,type="quantiles", group=c(1:20))

#Among simmr's default options, besides the matrix plots that we already saw, we can do density plots that represent
#the full distribution of our model's posterior probability density function.
plot(simmr_output,type='density', group=c(1:20),ggarg=xlim(0,1))

#Or use boxplots that offer a more synthetic view.
plot(simmr_output,type='boxplot', group=c(1:20), ggarg=ylim(0,1))

#We can also use simmr's built-in functions to compare the contribution of two sources to a consumer's diet. Sources
#are called nominatively and consumers by their group number. By default, you get a boxplot graph and the probability
#of difference is returned in the console. Although there is no global consensus on that, many people use a threshold
#of 95% of probability to consider that a trend is meaningful (analogy with frequentist p-values) Here, for Polycirrus 
#sp. in the "Anse du Lion" station:
compare_sources(simmr_output, source_names=c("SPOM", "Himantothallus"), group=16)
compare_sources(simmr_output, source_names=c("SPOM", "Sympagic"), group=16)

#You can also do that for more than two sources, but the probabilities can be harder to understand.
compare_sources(simmr_output, source_names=c("SPOM", "Sympagic", "Himantothallus","Biofilm"), group=16)

#Another thing we can do is comparing the contribution of a given source to the diet of two different consumer groups.
#Here, for example, we compare the diet of S. neumayeri in the two stations.
compare_groups(simmr_output, groups=19:20, source="SPOM")
compare_groups(simmr_output, groups=19:20, source="Himantothallus")
compare_groups(simmr_output, groups=19:20, source="Biofilm")
compare_groups(simmr_output, groups=19:20, source="Sympagic")

#Another example, where we compare the diet of A. colbecki and S. neumayeri at the "Cap des Elephants" station.
compare_groups(simmr_output, groups=c(1,19), source="SPOM")
compare_groups(simmr_output, groups=c(1,19), source="Himantothallus")
compare_groups(simmr_output, groups=c(1,19), source="Biofilm")
compare_groups(simmr_output, groups=c(1,19), source="Sympagic")

#Like for sources, you can compare multiple groups.
compare_groups(simmr_output, groups=c(2,12,20), source="Sympagic")

#You could even compare all groups at once: this is useless from the probability point of view, but can be an
#interesting way to see the importance of a source for all your consumer groups at a glance.
compare_groups(simmr_output, groups=1:20, source="Sympagic")

#Sometimes, you want to produce a more synthetic output, and therefore to combine the contributions of sources a
#a posteriori. Here, for examples, we could combine the sources "Himantothallus" and "Biofilm", as they both represent
#benthic food items.
simmr_output_combine = combine_sources(simmr_output,
                                    to_combine=c('Himantothallus','Biofilm'),
                                    new_source_name='Benthic')

#Once that's done, you can use the new "Benthic" source just like you would do for one of the original sources of the 
#model:
compare_groups(simmr_output_combine, groups=1:20, source="Benthic")

#If you want to customise the model's output more, the most effective is to extract the info you want to use manually.
#Here, we do that for the contribution of sympagic algae to Adamussium colbecki diet in both stations.
Adacol1_Sympagic <- simmr_output$output$`Adacol1`$BUGSoutput$sims.list$p[,'Sympagic']
Adacol2_Sympagic <- simmr_output$output$`Adacol2`$BUGSoutput$sims.list$p[,'Sympagic']

#Once you've done that, you can perform statistics on the output
summary(Adacol1_Sympagic)

#Or plot it
boxplot(cbind(Adacol1_Sympagic,Adacol2_Sympagic))

#Or calculate probabilities of difference in the model output
mean(Adacol1_Sympagic>Adacol2_Sympagic)

#With this manual extraction, you can represent and use your model's output in whichever way you see fit.
#Below are a a few examples.

#First, we extract data for all sources' contributions to Adamussium colbecki in site "Cap des Eléphants" and organise
#it neatly in a data frame.
Adacol1data <- c(simmr_output$output$`Adacol1`$BUGSoutput$sims.list$p[,'Sympagic'],
                 simmr_output$output$`Adacol1`$BUGSoutput$sims.list$p[,'Biofilm'],
                 simmr_output$output$`Adacol1`$BUGSoutput$sims.list$p[,'Himantothallus'],
                 simmr_output$output$`Adacol1`$BUGSoutput$sims.list$p[,'SPOM'])
Adacol1names <- c(rep("Sympagic", length(simmr_output$output$`Adacol1`$BUGSoutput$sims.list$p[,'Sympagic'])),
                  rep("Biofilm", length(simmr_output$output$`Adacol1`$BUGSoutput$sims.list$p[,'Biofilm'])),
                  rep("Macroalgae", length(simmr_output$output$`Adacol1`$BUGSoutput$sims.list$p[,'Himantothallus'])),
                  rep("SPOM", length(simmr_output$output$`Adacol1`$BUGSoutput$sims.list$p[,'SPOM'])))
Adacol1output <- data.frame(Adacol1data, Adacol1names)

#Then, we can use ggplot2 to produce a nice-looking, fully customisable, synthetic density plot.
ggplot(Adacol1output, aes(x = Adacol1data, colour = Adacol1names, fill = Adacol1names)) + 
  geom_density(alpha=0.7) + 
  theme_bw() +
  scale_x_continuous(name="Source contribution to diet", limits= c(0,1)) +
  scale_y_continuous(name=NULL, breaks=NULL) +
  scale_color_manual(values=c("#00A08A","#FF0000","#F98400","#5BBCD6")) +
  scale_fill_manual(values=c("#00A08A","#FF0000","#F98400","#5BBCD6")) +
  ggtitle(expression(italic("Adamussium colbecki")~"- Cap des Eléphants")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.background = element_rect(fill=NULL, size=NULL, linetype=NULL)) +
  guides(fill=guide_legend(title="Source"), colour=guide_legend(title="Source"))

#If you prefer, you can scale the density functions so that they all have the same maximum height. This can be useful
#when they have very different distributions that mask trends in the data.
ggplot(Adacol1output, aes(x = Adacol1data, y=..scaled.., colour = Adacol1names, fill = Adacol1names)) + 
  geom_density(alpha=0.7) + 
  theme_bw() +
  scale_x_continuous(name="Source contribution to diet", limits= c(0,1)) +
  scale_y_continuous(name=NULL, breaks=NULL) +
  scale_color_manual(values=c("#00A08A","#FF0000","#F98400","#5BBCD6")) +
  scale_fill_manual(values=c("#00A08A","#FF0000","#F98400","#5BBCD6")) +
  ggtitle(expression(italic("Adamussium colbecki")~"- Cap des Eléphants")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.background = element_rect(fill=NULL, size=NULL, linetype=NULL)) +
  guides(fill=guide_legend(title="Source"), colour=guide_legend(title="Source"))

#Alternatively, you can represent model output as boxplots.
ggplot(Adacol1output, aes(y = Adacol1data, x = Adacol1names, colour = Adacol1names, 
                      fill = Adacol1names, width=1)) +
  geom_boxplot(alpha=0.6, outlier.shape=NA) + 
  ggplot2::theme_bw() +
  scale_color_manual(values=c("#00A08A","#FF0000","#F98400","#5BBCD6")) +
  scale_fill_manual(values=c("#00A08A","#FF0000","#F98400","#5BBCD6")) +
  ggtitle(expression(italic("Adamussium colbecki")~"- Cap des Eléphants")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(name="Contribution to diet", limits=c(0,1)) +
  scale_x_discrete(name="Source") +
  guides(fill=FALSE) +
  guides(color=FALSE)

#In this other example, we will use the contributions of sympagic algae to Adamussium colbecki's diet in both sites,
#extracted earlier. First, we organise them in a data frame.
Elephants <- rep("Eléphants", length(Adacol1_Sympagic))
Lion <- rep("Lion", length(Adacol2_Sympagic))
AdaCol_Sympagic <- as.data.frame(rbind(cbind(Adacol1_Sympagic,Elephants),cbind(Adacol2_Sympagic,Lion)))
colnames(AdaCol_Sympagic)=c("Proportion","Site")
AdaCol_Sympagic$Proportion <- as.numeric(paste(AdaCol_Sympagic$Proportion))

#Then we can represent them as boxplots
ggplot(AdaCol_Sympagic, aes(y = Proportion, x = Site, colour = Site, fill = Site, width=1)) +
  geom_boxplot(alpha=0.6, outlier.shape=NA) + 
  ggplot2::theme_bw() +
  ggtitle(expression(italic("Adamussium colbecki")~"- Sympagic algae")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(name="Contribution to diet", limits=c(0,1)) +
  guides(fill=FALSE) +
  guides(color=FALSE)

#Or alternatively, as violin plots, which are more comparable to density plots.
ggplot(AdaCol_Sympagic, aes(y = Proportion, x = Site, colour = Site, fill = Site, width=0.9)) +
  geom_violin(alpha=0.7, scale="count") + 
  ggplot2::theme_bw() +
  ggtitle(expression(italic("Adamussium colbecki")~"- Sympagic algae")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(name="Contribution to diet", limits=c(0,1)) +
  guides(fill=FALSE) +
  guides(color=FALSE)

#That's it for part 1. Don't hesitate to experiment different things with your model output to find out which options
#work best in your case.


######################################################################################################################
#PART 2 - MIXSIAR
######################################################################################################################

#Install the mixSIAR package, as well as tidyr, that will be required later in the script.
install.packages("MixSIAR", dependencies = TRUE)
install.packages("tidyr")

#Load the installed packages.
library(MixSIAR)
library(tidyr)

#First, we will create a mixture object containing the consumer values. A key-innovation of MixSIAR is that it lets
#you specify up to two covariates (effects). Those effects can be random or fixed, nested or not, continuous or not.
#Here, we could create a model very similar to part 1 by using "Taxon" and "Site" as factor. Instead, we will do
#something different and use "Feeding guild" and "Taxon" as factors. Both are random, "Taxon" is nested inside "Feeding
#guild", and none of them is continuous.
mix <- load_mix_data(filename="Consumers.csv", 
                     iso_names=c("d13C","d15N"), 
                     factors=c("Taxon","Feeding_guild"), 
                     fac_random=c(TRUE,TRUE), 
                     fac_nested=c(TRUE,FALSE), 
                     cont_effects=NULL)

#We can then load the source data. We can use the raw data, or, like here, means and SDs. If relevant, you can load
#source data for each factor (e.g. different data for different sites or seasons), but that's not the case here. Like
#in part 1, we won't use concentration dependencies.
source <- load_source_data(filename="Sources.csv",
                           source_factors=NULL, 
                           conc_dep=FALSE, 
                           data_type="means", 
                           mix)

#Then we can load trophic enrichment factors.
discr <- load_discr_data(filename="TEFs.csv", mix)

#Like before, it's crucial to plot your data to make sure that everything is in order.
plot_data(filename="isospace_plot", plot_save_pdf=FALSE, plot_save_png=FALSE, mix,source,discr)

#Next, we will write our model. This command will take into account all the parameters you specified, and write the
#code to use in a text file that will later be passed on to the run_model() function. It is here that you specify the
#error structure, another MixSIAR key innovation. You can chose between process error (that will account for intra-group
#variability linked with subsampling biases or specialisation), residual error (that will account for factors such as
#inter-individual differences in assimilation efficiency and metabolic rates between consumers) or both. #Here (and 
#actually in most situations), we want to use both to have a model that is as biologically realist as possible.
model_filename <- "MixSIAR_model.txt"
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

#Once it's done, we can run our model. MixSIAR runs can be very long. Therefore, it can be a good idea to perform a
#"test" run first. These runs are very short and there is no way that they could reach convergence, it's just to check
#that everything is OK.
jags1 <- run_model(run="test", mix, source, discr, model_filename, 
                    alpha.prior = 1, resid_err, process_err)

#If your test run went fine, you can proceed with a longer run. Here is the code to perform a "normal" run. Please note
#that in this case, it took over 5 hours to run.
#jags2 <- run_model(run="normal", mix, source, discr, model_filename, 
#                   alpha.prior = 1, resid_err, process_err)
#save(jags2, file="MixSIAR_long.rdata")

#To avoid losing time during the course, load the model data that is provided alongside this script instead.
load(file = "MixSIAR_long.rdata")

#Note that if you want more flexibility, you can directly specify MCMC parameters like this:
#myparams <- list(chainLength=100000, burn=50000, thin=50, chains=3, calcDIC=TRUE)
#jagsx <- run_model(run=myparams, mix, source, discr, 
#                   model_filename, alpha.prior = 1, resid_err, process_err)

#Now that we have a model, we can check its diagnostics and explore its output. By default, MixSIAR produces a lot of
#graphs and saves them to your hard drive, so it might be a good idea to create a separate folder to keep things
#organised.
dir.create("mixsiar_output")
setwd("mixsiar_output")

#Then we can specify the output & diagnostics options. Note that even you have many different choices here.
output_options <- list(summary_save = TRUE,
                      summary_name = "summary_statistics",
                      sup_post = FALSE,
                      plot_post_save_pdf = TRUE,
                      plot_post_name = "density_plot",
                      sup_pairs = FALSE,
                      plot_pairs_save_pdf = TRUE,
                      plot_pairs_name = "matrix_plot",
                      sup_xy = TRUE,
                      plot_xy_save_pdf = TRUE,
                      plot_xy_name = "xy_plot",
                      gelman = TRUE,
                      heidel = FALSE,
                      geweke = FALSE,
                      diag_save = TRUE,
                      diag_name = "diagnostics",
                      indiv_effect = FALSE,
                      plot_post_save_png = FALSE,
                      plot_pairs_save_png = FALSE,
                      plot_xy_save_png = FALSE)

#Then we can create our diagnostics and output graphs.
output_JAGS (jags2, mix, source, output_options)

#That leaves you with many graphs and results to analyse and interpretate. If you want to go further with output
#customisation, you can extract the data of interest from the model object and use them for stats and plots.
#Here, for example, we look at summary statistics for the first source (they are order alphabetically, so here it will
#be biofilm) for the first level of factor 1, i.e. the first taxon, Adamussium colbecki.
summary(jags2$BUGSoutput$sims.list$p.fac1[,1,1])

#You can also test whether the contribution of biofilm (source 1) to Adamussium colbecki's (taxon 1) diet is higher than
#the one of Himantothallus (source 2).
mean(jags2$BUGSoutput$sims.list$p.fac1[,1,1]>jags2$BUGSoutput$sims.list$p.fac1[,1,2])

#Or test whether contribution of sympagic algae (source 4) is bigger in Sterechinus neumayeri (taxon 10) than in 
#Flabegraviera mundata (taxon 2).
mean(jags2$BUGSoutput$sims.list$p.fac1[,2,4]<jags2$BUGSoutput$sims.list$p.fac1[,10,4])

#Like before, you can also produce nice and fully customisable graphs using ggplot2. Here is an example for
#Adamussium colbecki.
#First, extract and organise the data.
post.acol.data <- c(jags2$BUGSoutput$sims.list$p.fac1[,1,1], jags2$BUGSoutput$sims.list$p.fac1[,1,2], 
                    jags2$BUGSoutput$sims.list$p.fac1[,1,3], jags2$BUGSoutput$sims.list$p.fac1[,1,4])
post.acol.names <- c(rep("Biofilm", length(jags2$BUGSoutput$sims.list$p.fac1[,1,1])),
                     rep("Macroalgae", length(jags2$BUGSoutput$sims.list$p.fac1[,1,2])),
                     rep("SPOM", length(jags2$BUGSoutput$sims.list$p.fac1[,1,3])),
                     rep("Sympagic", length(jags2$BUGSoutput$sims.list$p.fac1[,1,4])))
post.acol <- data.frame(post.acol.data, post.acol.names)

#You can then represent each source contribution to A. colbecki diet as density plots
ggplot(post.acol, aes(x = post.acol.data, y=..scaled.., colour = post.acol.names, fill = post.acol.names)) + 
  geom_density(alpha=0.5) + 
  theme_bw() +
  scale_x_continuous(name="Source contribution to diet", limits= c(0,1)) +
  scale_y_continuous(name=NULL, breaks=NULL) +
  scale_color_manual(values=c("#00A08A","#FF0000","#F98400","#5BBCD6")) +
  scale_fill_manual(values=c("#00A08A","#FF0000","#F98400","#5BBCD6")) +
  ggtitle(expression(italic("Adamussium colbecki"))) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.background = element_rect(fill=NULL, size=NULL, linetype=NULL)) +
  guides(fill=guide_legend(title="Source"), colour=guide_legend(title="Source"))

#Or as boxplots
ggplot(post.acol, aes(y = post.acol.data, x = post.acol.names, colour = post.acol.names, 
                      fill = post.acol.names, width=1)) +
  geom_boxplot(alpha=0.6, outlier.shape=NA) + 
  ggplot2::theme_bw() +
  scale_color_manual(values=c("#00A08A","#FF0000","#F98400","#5BBCD6")) +
  scale_fill_manual(values=c("#00A08A","#FF0000","#F98400","#5BBCD6")) +
  ggtitle(expression(italic("Adamussium colbecki"))) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(name="Contribution to diet", limits=c(0,1)) +
  scale_x_discrete(name="Source") +
  guides(fill=FALSE) +
  guides(color=FALSE)

#Or as violin plots
ggplot(post.acol, aes(y = post.acol.data, x = post.acol.names, colour = post.acol.names, 
                      fill = post.acol.names, width=1)) +
  geom_violin(alpha=0.6, scale="width") + 
  ggplot2::theme_bw() +
  scale_color_manual(values=c("#00A08A","#FF0000","#F98400","#5BBCD6")) +
  scale_fill_manual(values=c("#00A08A","#FF0000","#F98400","#5BBCD6")) +
  ggtitle(expression(italic("Adamussium colbecki"))) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(name="Contribution to diet", limits=c(0,1)) +
  scale_x_discrete(name="Source") +
  guides(fill=FALSE) +
  guides(color=FALSE)

#That's it for part 2. Any questions?


######################################################################################################################
#PART 3 - INDEPENDENT WORK
######################################################################################################################

#Now it's your turn to work! Using what we learned today, try to complete the following tasks / answer the following 
#questions.

                                                      
#1) Using simmr, create a model using the "short" mcmc parameters we tried earlier, and grouping the consumers by taxa,
#without site distinction.
#Tip: you can control how consumers are grouped with the "group" argument of the simmr_load() function. 

#2) Check convergence diagnostics for your model, and decide if they are satisfactory.

#3) Compare the posterior predictive distributions to data points for Adamussium colbecki and Sterechinus neumayeri.
#Using those graphs, what can you tell about the model fit for the two taxa? What could explain the trend you see?

#4) Using simmr's built-in functions, plot your model's output for Adamussium colbecki and Sterechinus neumayeri as
#density plots and as boxplots. Interpret those graphs.

#5) Test how likely it is that contributions of all 4 source differ between the two taxa. What can you conclude?

#6) Compare output of your model with the MixSIAR model we built earlier, for the two taxa. What can you tell for those
#graphs? What could explain this?


######################################################################################################################
#END OF SCRIPT
######################################################################################################################