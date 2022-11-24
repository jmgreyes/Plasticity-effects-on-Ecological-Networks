#---------------------------------------------------------------------------------------------------
# DESCRIPTION:
# A companion script to the MS 
# "The role of phenotypic plasticity in shaping ecological networks"
# by José M. Gómez, Adela González-Megías, Cristina Armas, 
# Eduardo Narbona, Luis Navarro & Francisco Perfectti
#---------------------------------------------------------------------------------------------------
## First version 18 Nov 2022.
## Hosted at https://github.com/jmgreyes/Plasticity-effects-on-Ecological-Networks
#---------------------------------------------------------------------------------------------------
#
# Starting from several bipartite ecological networks, this script calculates how 
# niche expansion (scenario 1), niche shift (scenario 2), and niche jump 
# (scenario 3) due to phenotypic plasticity of a focal species, impact several 
# network metrics (connectance, modularity, and nestedness). 
#
#---------------------------------------------------------------------------------------------------
## Loading libraries 
#---------------------------------------------------------------------------------------------------
library (bipartite) # network & modularity
library(tidyr) # to manipulate dataframes
library (gdata) # to organize R enviroment
#
#---------------------------------------------------------------------------------------------------
## Reading networks 
#---------------------------------------------------------------------------------------------------
	rm(list=ls())
#
# Reading all files at once
#
# Network files are stored in a different directory
	setwd("./Networks")
#
# We obtain a named list with all network files in that directory
	network.names <- list.files(pattern=".txt")
#
# We load all files as a list
	networks = lapply(network.names, read.table, fill=TRUE)
#
# Back to the main directory
	setwd("..")
#
# Remove the extension of the member names (AFTER DOWNLOADING THE FILES!)
	network.names <- stringr::str_remove(network.names, ".txt")
	names(networks)<-network.names
#
#---------------------------------------------------------------------------------------------------
## Creating the storing datasets 
#---------------------------------------------------------------------------------------------------
sink("nonplastic.txt", append=TRUE)
cat("network\t","focal\t", 
	"non.plastic.nestedness\t", 
	"nonplastic.modularity\t", 
	"non.plastic.connectance\t", 
	"N.modules\t", "N.partners\n")
sink()
#
#
sink("plastic.niche.expansion.txt", append=TRUE)
cat("network\t",
	"focal\t", 
	"number.interactions\t", 
	"j\t", "nestedness\t", 
	"modularity\t", "connectance\t", 
	"plastic.degree\t", 
	"N.partners\n")
sink()
#
#
sink("plastic.niche.shift.txt", append=TRUE)
cat("network\t",
	"focal\t", 
	"number.interactions\t", 
	"original.module\t", 
	"final.module\t", 
	"nestedness\t", 
	"modularity\t", 
	"connectance\t", 
	"plastic.degree\t", 
	"N.partners\n")
sink()
#
#
sink("plastic.niche.jump.txt", append=TRUE)
cat("network\t",
	"focal\t", 
	"number.interactions\t", 
	"extra.interactions\t", 
	"nestedness\t", "modularity\t", 
	"connectance\t", 
	"plastic.degree\t", 
	"N.partners\n")
sink()
#
#-----------------------------------------------------------------------------------------------------
## Calculating the NON-PLASTIC metrics
#-----------------------------------------------------------------------------------------------------
#
# We first calculate the non-plastic metrics
# For each of the 10 networks we calculated the original values of 
# connectance, nestedness and modularity before allowing species to be plastic
#
for (d in 1:length(networks))
	{
	network.name<-names(networks)[d]
	#
	net<-as.data.frame(networks[d])
	net<-ifelse(net>0, 1,0)
	rownames(net)<- paste0(rep(c("F"), 1), rep(1:dim(net)[1]))   #name rows with F1, F2, F3... 
	colnames(net) <- paste0(rep(c("P"), 1), rep(1:dim(net)[2]))  #name columns with P1, P2, P3...
	#
	#Rearranges the network according to nestedness (from generalits to specialists)
	packed.net<-nestedness(net)$packed.matrix  # NOTE: other versions needs $comm
	#
	# identify modules using computeModules 
	modularity.original<-computeModules(packed.net) 
	#
	# Generating matrices for focal and partner species with modules
	modules <- modularity.original@modules
	modules <- modules[-1, -c(1:2)]
	#
	modules.Focal<-modules[,1:dim(packed.net)[1]]
	colnames(modules.Focal) <- rownames(packed.net)
	#
	modules.Partner<-modules[,(dim(packed.net)[1]+1):dim(modules)[2]]
	colnames(modules.Partner) <- colnames(packed.net)
	#
	# Matrix as presence (1) or absence (0) of each species per module
	modules.Focal.discrete<-ifelse(modules.Focal >=1, 1, 0)
	rownames(modules.Focal.discrete)<- paste0(rep(c("M"), 1), rep(1:dim(modules.Focal.discrete)[1])) 
	#
	modules.Partner.discrete<-ifelse(modules.Partner>=1, 1,0)
	rownames(modules.Partner.discrete)<-paste0(rep(c("M"), 1), rep(1:dim(modules.Partner.discrete)[1])) 
	#
	# Number of Modules-----------------------------------------------------------------------------
	NModules <- dim(modules.Focal.discrete)[1]
	#
	# Modular membership of focals -----------------------------------------------------------------
	a<-reshape::melt(modules.Focal.discrete, id.vars=c(1:ncol(modules.Focal.discrete)),var='module')
	a[a==0]<-NA
	Focal.per.module <-a[complete.cases(a),]
	Focal.per.module$value<-NULL
	names(Focal.per.module)<- c("module", "plant")
	Focal.per.module$module <- stringr::str_remove(Focal.per.module$module, "M")
	rownames(Focal.per.module)<-Focal.per.module$plant; Focal.per.module$plant<-NULL
	Focal.per.module$module<-as.numeric(Focal.per.module$module)
	#
	# The metrics ----------------------------------------------------------------------------------
	nestedness.nonplastic<-nest.smdm(packed.net)$NODFmatrix # Nestedness
	modularity.nonplastic<-modularity.original@likelihood # Modularity
	connectance.nonplastic<-networklevel(packed.net, index="connectance") # Connectance
	number.of.partners<-dim(packed.net)[2]
	#
	write.table((data.frame(network.name, 
				round(nestedness.nonplastic,4), 
				round(modularity.nonplastic,4), 
				round(connectance.nonplastic,4), 
				NModules, number.of.partners)), 
				file="nonplastic.txt", append=TRUE, sep="\t", row.names=F, col.names=F)  
	#	
	#	
#-----------------------------------------------------------------------------------------------------
## Calculating the PLASTIC metrics
#-----------------------------------------------------------------------------------------------------
#
## SCENARIO 1 ----------------------------------------------------------------------------------------
#
# For each of the 10 networks we calculated the values of 
# connectance, nestedness and modularity after expanding niche breadth due to plasticity
#
# PROCEDURE: 	For each species in each network, we added extra links with non-connected species 
#            	according to a probability function obtained through a decreasing uniform
#            	distribution based on the degree distribution. In particular, we allowed any species 
#            	to be connected with those species obtaining a probability > 0.55 in that distribution.
#            	This method forces species to be connected with most connected species of the other set,
#            	a property shared by most ecological networks
#            	We repeated this operation 20 times per species (number of iterations)
#
#
	plastic.net<-packed.net
	focals <- labels(packed.net)[[1]]
	iterations<-20						# number of iterations
	#
	#
	for (i in focals)   					# loop for i species
		{        
		focal <- i
		sp_f<-plastic.net[focal,] 
		number.interactions<-sum(sp_f)        		# number of interactions					
		#
		for (j in 1:iterations)
	 		{                       		# loop for s iterations
    			prob<-sort(runif(length(sp_f)), decreasing = TRUE)
    			spfocal.plastic<-ifelse(sp_f, 1, ifelse(prob>0.55, 1,0))
    			plastic.net[focal,]<- spfocal.plastic
    	  		#
    			nestedness.plastic <-nest.smdm(plastic.net)$NODFmatrix               
    			modularity.plastic <-computeModules(plastic.net)@likelihood
    			connectance.plastic <-networklevel(plastic.net, index="connectance")
    			plastic.degree<-sum(spfocal.plastic)	
    			number.of.partners<-dim(packed.net)[2]		
			#
			write.table((data.frame(network.name, 
						i, 
						number.interactions, 
						j, 
						round(nestedness.plastic,4), 
						round(modularity.plastic,4), 
						round(connectance.plastic,4), 
						plastic.degree, 
						number.of.partners)), 
						file="plastic.niche.expansion.txt", 
							append=TRUE, sep="\t", row.names=F, col.names=F)
			#
       		}
   		}
   	}								# optional in case we wish to run scenario 1 alone
#	
## SCENARIO 2 ----------------------------------------------------------------------------------------
#
# For each of the 10 networks we calculated the values of 
# connectance, nestedness and modularity after allowing to shift among existing niches
#
# PROCEDURE:	For each species in each network, we removed the links with the partners of its module
#           	and established new links with all the partners defining any other module.
#           	We repeated this operation for each of the modules, once at the time.
#           	This method consider the modules to be a proxy of the interaction niche,
#           	a property already shown for many ecological networks.
#	     	  	In addition, this method forces each species to visit all other modules, and thereby
#	       	outcomes are not biased by the across-modules differences in generalization degree.
#
#
	plastic.net<-packed.net
	focals <- labels(packed.net)[[1]]
	#
	#
	for (r in focals) 			      # loop for focal species from the first to the last
		{                           
  		module.vector <- c(1: dim(modules.Partner.discrete)[1])   # vector with number of modules
  		original.module<-Focal.per.module[r,]
  		s<-module.vector[-original.module]
                #
  		for (m in s)						  # loop for modules (without the actual module)
  			{                             
    			plastic.net[r,]<-ifelse(modules.Partner.discrete[Focal.per.module[r,],]>0 & packed.net[r,]>0, 0, 
					 ifelse(modules.Partner.discrete[m,]>0, 1, packed.net[r,]))
			#
    			number.interactions <-sum(packed.net[r,])
    			plastic.degree <-sum(plastic.net[r,])
    			number.of.partners<-dim(packed.net)[2]
                      	#
    			nestedness.shift<-nest.smdm(plastic.net)$NODFmatrix 
    			modularity.shift<-computeModules(plastic.net)@likelihood
    			connectance.shift<-networklevel(plastic.net, index="connectance")
			#
			write.table((data.frame(network.name, 
						r, 
						number.interactions, 
						original.module, 
						m, 
						round(nestedness.shift,4), 
						round(modularity.shift,4), 
						round(connectance.shift,4), 
						plastic.degree, number.of.partners)), 
						file="plastic.niche.shift.txt", 
							append=TRUE, sep="\t", row.names=F, col.names=F)
			#    
  			}
		  }
	}  								 # optional in case we wish to run scenario 2 alone
#
## SCENARIO 3 ----------------------------------------------------------------------------------------
#
# For each of the 10 networks we calculated the values of 
# connectance, nestedness and modularity after allowing the species to jump to a new niche
#
# PROCEDURE: 	For each species in each network, we removed the links with the partners of its module
#            	and established new links with new partners added to the original network.
#            	We repeated this operation allowing the number of new partners to range from 1 to 5.
#            	This method consider the modules to be a proxy of the interaction niche,
#            	a property already shown for many ecological networks.
#	     		In addition, this method allows each species to establish a new niche that can range
#	     		from extreme specialisation to high generalisation, avoiding biases caused for changes
#	     		in generalization degree.
#
#
	Nmax<-5						# maximum number of invading species
	focals <- labels(packed.net)[[1]]
	#
	for (w in 1:Nmax)				# number of invading species
		{
		new.niche<-matrix(0,dim(packed.net)[1],w)
		colnames(new.niche) <- paste0(rep(c("N"), 1), rep(1:dim(new.niche)[2]))
		extra.net<-cbind(packed.net, new.niche)
		#		
		for (v in focals)
			{
  			plastic.net<-extra.net
  			plastic.net[v,1:dim(packed.net)[2]]<-0
  			plastic.net[v,(dim(packed.net)[2]+1):dim(extra.net)[2]]<-1
 			# 
  			number.interactions <-sum(packed.net[v,])
  			plastic.degree <-sum(plastic.net[v,])
  			number.of.partners<-dim(packed.net)[2]
 			# 
  			nestedness.jump<-nest.smdm(plastic.net)$NODFmatrix          
  			modularity.jump<-computeModules(plastic.net)@likelihood
  			connectance.jump<-networklevel(plastic.net, index="connectance")
  			#
			write.table((data.frame(network.name, 
						v, 
						number.interactions, 
						w, round(nestedness.jump,4), 
						round(modularity.jump,4), 
						round(connectance.jump,4), 
						plastic.degree, 
						number.of.partners)), 
						file="plastic.niche.jump.txt", 
							append=TRUE, sep="\t", row.names=F, col.names=F)
			#
			}
		}
	}  								
#
#
#-----------------------------------------------------------------------------------------------------
## Manipulating the data
#-----------------------------------------------------------------------------------------------------
# Reading the datasets
	nonplastic<-read.table("nonplastic.txt", header=T)
	plastic.expansion<-read.table("plastic.niche.expansion.txt", header=T)
	plastic.shift<-read.table("plastic.niche.shift.txt", header=T)
	plastic.jump<-read.table("plastic.niche.jump.txt", header=T)
#
#
# Categorizing each species as generalist or specialist
# Generalist: linked with more than 33% of the available partners
# Specialist: linked with less than 10% of the available partners
# 'nifunifa' are species neither generalist nor specialist
	plastic.expansion$generalization.degree<-
	ifelse(plastic.expansion$number.interactions > plastic.expansion$N.partners *0.33, "generalized", 
	ifelse(plastic.expansion$number.interactions < plastic.expansion$N.partners *0.10, "specialized", 
	"nifunifa"))					 
#
	plastic.shift$generalization.degree<-
	ifelse(plastic.shift$number.interactions > plastic.shift$N.partners *0.33, "generalized", 
	ifelse(plastic.shift$number.interactions < plastic.shift$N.partners *0.10, "specialized", "nifunifa"))
#
	plastic.jump$generalization.degree<-
	ifelse(plastic.jump$number.interactions > plastic.jump$N.partners *0.33, "generalized", 
	ifelse(plastic.jump$number.interactions < plastic.jump$N.partners *0.10, "specialized", "nifunifa"))
#
#
# Merging plastic-nonplastic datasets
	plastic.expansion <-merge(plastic.expansion, nonplastic, by="network")
	plastic.shift <-merge(plastic.shift, nonplastic, by="network")
	plastic.jump <-merge(plastic.jump, nonplastic, by="network")
#
#
# Calculating the differences as increase rates
	plastic.expansion$Diff.nest<-(plastic.expansion$nestedness-plastic.expansion$nonplastic.nestedness)
	/plastic.expansion$nonplastic.nestedness
	plastic.expansion$Diff.mod<-(plastic.expansion$modularity-plastic.expansion$nonplastic.modularity)
	/plastic.expansion$nonplastic.modularity
	plastic.expansion$Diff.con<-(plastic.expansion$connectance-plastic.expansion$nonplastic.connectance)
	/plastic.expansion$nonplastic.connectance
#
	plastic.shift$Diff.nest<-(plastic.shift$nestedness-plastic.shift$nonplastic.nestedness)
	/plastic.shift$nonplastic.nestedness
	plastic.shift$Diff.mod<-(plastic.shift$modularity-plastic.shift$nonplastic.modularity)
	/plastic.shift$nonplastic.nestedness
	plastic.shift$Diff.con<-(plastic.shift$connectance-plastic.shift$nonplastic.connectance)
	/plastic.shift$nonplastic.connectance
#
	plastic.jump$Diff.nest<-(plastic.jump$nestedness-plastic.jump$nonplastic.nestedness)
	/plastic.jump$nonplastic.nestedness
	plastic.jump$Diff.mod<-(plastic.jump$modularity-plastic.jump$nonplastic.modularity)
	/plastic.jump$nonplastic.modularity
	plastic.jump$Diff.con<-(plastic.jump$connectance-plastic.jump$nonplastic.connectance)
	/plastic.jump$nonplastic.connectance
#
	plastic.expansion$Diff.degree<-(plastic.expansion$plastic.degree-plastic.expansion$number.interactions)
	/plastic.expansion$number.interactions
	plastic.shift$Diff.degree<-(plastic.shift$plastic.degree-plastic.shift$number.interactions)
	/plastic.shift$number.interactions
	plastic.jump$Diff.degree<-(plastic.jump$plastic.degree-plastic.jump$number.interactions)
	/plastic.jump$number.interactions
#
#
# Removing the species that are neither specialist nor generalist (='nifunifa')
#
	plastic.expansion.final<-plastic.expansion[plastic.expansion$generalization.degree!="nifunifa",]
	plastic.shift.final<-plastic.shift[plastic.shift$generalization.degree!="nifunifa",]
	plastic.jump.final<-plastic.jump[plastic.jump$generalization.degree!="nifunifa",]
#
#-----------------------------------------------------------------------------------------------------
## Visualizing the outcomes
#-----------------------------------------------------------------------------------------------------
#
	par(mfrow=c(3,3))
#
	boxplot(plastic.expansion.final$Diff.con ~ plastic.expansion.final$generalization.degree)
	boxplot(plastic.shift.final$Diff.con ~ plastic.shift.final$generalization.degree)
	boxplot(plastic.jump.final$Diff.con ~ plastic.jump.final$generalization.degree)
#
#
	boxplot(plastic.expansion.final$Diff.nest~ plastic.expansion.final$generalization.degree)
	boxplot(plastic.shift.final$Diff.nest~ plastic.shift.final$generalization.degree)
	boxplot(plastic.jump.final$Diff.nest~ plastic.jump.final$generalization.degree)
#
#
	boxplot(plastic.expansion.final$Diff.mod ~ plastic.expansion.final$generalization.degree)
	boxplot(plastic.shift.final$Diff.mod ~ plastic.shift.final$generalization.degree)
	boxplot(plastic.jump.final$Diff.mod ~ plastic.jump.final$generalization.degree)
#
#-----------------------------------------------------------------------------------------------------
