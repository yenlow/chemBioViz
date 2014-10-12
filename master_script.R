#Chemical-biological read-across/Dual-space kNN of LD50 dataset with cytotoxicity as biological descriptors
#Inputs: 	descriptors read in as .xa files – available in Rcode.zip
#Dependencies: libraries (boot, caret, class, e1071, plotrix, ROCR, vegan)
#
#Dependencies: scripts (installnewpackage, readXAfile.R, sampling.R, 
#				multispaceNNobj.R, multispace_functions_AD_fold_scaled.R,
# 				validationstats_BIN.R, bootstrapSD_BIN.R)
#				 – available in Rcode.zip		
#05-June-12 Yen Low yenlow@gmail.com
############################################################################

rm(list=ls(all=T)) #clear memory

############### 1. Load required packages and scripts #####################
required.packages=c("boot","caret","class","e1071","plotrix","ROCR","vegan")
#check if required packages are already installed
new.packages=required.packages[!(required.packages %in% installed.packages()[,"Package"])]
#install only packages previously uninstalled
if(length(new.packages)) install.packages(new.packages)

#load libraries
require(plotrix) #for radial.plot
require(vegan)   #for jaccard distance
require(caret)   #for confusion matrix and prediction performance metrics
require(class)   #for cross validation
require(e1071)   #for cross validation
require(ROCR)    #for AUC calculations
require(boot)    #for calculating errors by bootstrapping

#load source scripts (scripts must be in same folder as master_script.R) 
source("readXAfile.R")  #for reading input .xa data files
source("sampling.R")	#sample n-fold CV for external validation
source("multispaceNNobj.R") #create dualspace kNN objects for dualspace kNN functions
source("multispace_functions_AD_fold_scaled.R") #dualspace kNN functions
source("validationstats_BIN.R")	#calculate model prediction performance


############### 2. Load required data files #####################
#read in data files (e.g. .x, .a, .xa files)
dschem=readXAfile("ld50_drg_n.xa") #chemical descriptors
dsgenes=readXAfile("ld50_atp_csp_n.xa") #biological descriptors
chemspace=dschem$x
genespace=dsgenes$x
hybridspace=cbind(chemspace,genespace) #concatenate to form hybrid space

dsact<-as.integer(as.character(dschem$y))	#ensure numeric vector
cpdnames<-dschem$cpds				#get compound names
rownames(hybridspace)=cpdnames	#cpdnames as row names for hybridspace

#Create and save 5-fold cross validation scheme
#random sampling with activity binning
#create matrix with 5 rows (row=cpd ID in each of the 5 external folds)
extID=genextID(ID=NULL,dsact,5)
save(file="shuffleID.RData",extID) 	#output file: cpd IDs in external folds


############### 3. Run dual-space kNN #####################
#Develop dualspace-kNN models with kchem values between krange
#					and   kbio values between k2range
#inputs: 	krange (vector between 1 to 5 for kchem values)
#		k2range (vector between 1 to 5 for kbio values; set to 0 (default) if not required)
#		extID (matrix of cpd IDs as external set; number of rows should correspond to number of external folds)
#		classlab (vector of classification labels)
#		chemspace (chemical space matrix as the first space)
#		genespace (biological space matrix as the second space; set to NULL (default) if single space)
#		cpdnames (vector of compound names)
#		outpred (name of output file containing predicted labels for each external compound)
#		outvalidstat (name of output file showing model performance)
#		similaritycutoff (minimum similarity that external cpds are to training cpds; default=0.3)
#		nfoldint (number of folds for internal cross validation for tuning kchem, kbio values; recommended nfoldint=10 (default))
#		calcErr (logical. Set to TRUE (default) if bootstrapping errors are desired. Set to FALSE to speed up)
#		ntrials (number of bootstrapping samples (1000 trials by default). Set to 100 to speed up)
#output: list object of 1. predicted values of external compounds (predtable)
#				2. prediction performance on external compounds (perf)
#				3. model object with external compounds' neighbors, neighbors' distances and classification labels (extobj)
#
#Pls be patient! Takes a few minutes on dual-core processors.
#Loops through   5-fold external cross validation
#within each is 10-fold internal cross validation (gridsearch for best kchem and kbio) - Set nfoldint=5 instead of 10 to speed up

#chemical kNN (kchemNN)
cat("\n","Single space chemical kNN model (kchemNN)","\n")
chemmod=dualspacekNN(	krange=1:5,
				chemspace=chemspace,
				outpred="chem.pred",outvalidstat="validationstats_chem.txt",
				cpdnames=cpdnames,similaritycutoff=0.3,
				extID=extID,classlab=dsact,nfoldint=10,calcErr=T,ntrials=1000)

#biological kNN (kbioNN)
cat("\n","Single space biological kNN model (kbioNN)","\n")
genemod=dualspacekNN(	krange=1:5,
				chemspace=genespace,
				outpred="gene.pred",outvalidstat="validationstats_gene.txt",
				cpdnames=cpdnames,similaritycutoff=0.3,
				extID=extID,classlab=dsact,nfoldint=10,calcErr=T,ntrials=1000)

#hybrid kNN (khybridNN)
cat("\n","Single space hybrid kNN model (khybridNN)","\n")
hybridmod=dualspacekNN(	krange=1:5,
				chemspace=hybridspace,
				outpred="hybrid.pred",outvalidstat="validationstats_hybrid.txt",
				cpdnames=cpdnames,similaritycutoff=0.3,
				extID=extID,classlab=dsact,nfoldint=10,calcErr=T,ntrials=1000)

#dual-space kNN (kchemNN + kbioNN)
cat("\n","Dual-space kNN model on chemical and biological spaces","\n")
dualmod=dualspacekNN(	krange=1:5,k2range=1:5,
				chemspace=chemspace,genespace=genespace,
				outpred="dual.pred",outvalidstat="validationstats_dual.txt",
				cpdnames=cpdnames,similaritycutoff=0.3,
				extID=extID,classlab=dsact,nfoldint=10,calcErr=T,ntrials=1000)


############### 4. Draw radial plots showing nearest neighbors of external compound #####################
#prepare kNNobj input object for radial plots
#generate kNNobj with kchem=5 and kbio=5 (of training compounds)
#calculate neighbours (default: jaccard distance)
chemNN=NNmat(chemspace,cpdnames,dsact,disttype="jaccard")
geneNN=NNmat(genespace,cpdnames,dsact)
kNNobj_k5=kNNobj(5,chemNN,geneNN) #kchem=5 and kbio=5

#output: single radial plot of cpd #20 in 1x1 grid
pdf("singlecpd.pdf",width=4,height=4)
suppressWarnings(starplotgrid(20,1,1,kNNobj_k5,label=T,byrow=T))
#warnings about axes are suppressed
dev.off()

#output: multiple radial plots of cpds #2, #20, #13, #16 in 2x3 grid
pdf("6cpds_2by3grid.pdf",width=11,height=8.5)
suppressWarnings(starplotgrid(c(2,20,13,16,18,19),2,3,kNNobj_k5,label=T,byrow=T))
#warnings about axes are suppressed
dev.off()

#generate kNNobj with kchem=5 and kbio=5 (of external compounds)
#convert modobj (or extobj) to kNNobj
#input: classlab_fold (reorder vector of external compounds' classification labels by external folds)
classlab_fold=as.numeric(as.character(dualmod$predtable$yobs))
names(classlab_fold)=rownames(dualmod$predtable)
kNNobj_dual=suppressWarnings(extobj2kNNobj(dualmod$extobj,classlab_fold))

#output: multiple radial plots of cpds #2, #20, #13, #16 in 2x2 grid
pdf("4cpds_2by2grid.pdf",width=11,height=8.5)
suppressWarnings(starplotgrid(c(2,20,13,16),2,2,kNNobj_dual[[4]],label=F,byrow=T))
#warnings about axes are suppressed
dev.off()


############## 5. save important objects into dualspace.RData ###############################################
save(	kNNobj_dual,
	dualmod,chemmod,genemod,hybridmod,
	file="dualspace.RData")

cat("\n","You have completed running master_script.R. Check output in the same folder","\n")

