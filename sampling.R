#Generates IDs for external n-fold cross validation (integer nfold>=2)
#Outputs a matrix of nfold rows with each row representing the IDs of external cpds
#
#Also generates IDs for single external set of nfold proportion (0<nfold<1)
#Outputs a list of extID and modID
#
#Stratified sampling maintains active/inactive ratio
#Inputs:	1. ID=NULL (vector of possble IDs. Optional; Default to 1:length(actvec))
#			2. actvec (vector of activity values)
#			3. nfold (number of folds in cross-validation scheme. If btwn 0 and 1, samples a ext set of that fraction. Defaults to 5)
#			4. srnd=sample.int(1000,1) (randomizer seed. Optional)
#			5. foldoverminclass=NULL (foldoverminclass=majority_class/minority_class. Optional. Defaults to original ratio)
#			6. replace=T (undersample with replacement? only valid when undersampling and n-fold CV)
#09-Aug-12 Yen Low
#
###############################################################################
##Example:
##Create an activity vector of size 100 with nminor samples
#nminor=20
#actvec=c(rep(0,nminor),rep(1,100-nminor))
##Call genextID to generate IDs for external sets and modeling sets
##Set foldoverminclass=1 to balance data set
#temp=genextID(ID=NULL,actvec,foldoverminclass=1,replace=F)
#table(unlist(temp$modID))
#Save generated IDs to text file.
#extID2txt(temp,outputroot="C:/Users/yenlow/Desktop/")


genextID<-function(ID=NULL,actvec,nfold=5,srnd=sample.int(1000,1),foldoverminclass=NULL,replace=TRUE){
    if(!is.vector(ID) && !is.null(ID)) stop("ID must be a vector of integer values representing row number of data. Otherwise, ID defaults to 1:length(actvec)")
	if(!is.vector(actvec)) stop("actvec must be a vector of activity values.")
	if(!is.null(foldoverminclass) && foldoverminclass<0) stop("foldoverminclass=majority_class/minority_class must be a number.\n Enter a positive number to downsample the majority class.")
	if(!is.logical(replace)) stop("Set replace to TRUE to sample majority class with replacement or FALSE to sample without replacement.\n Only valid for n-fold cross validation. Defaults to TRUE.\n Set to FALSE (without replacement) when minority class is 15-25% to maximize sampling of majority samples assuming 5-fold cross validation.")
		if(is.numeric(nfold)){
			if(is.null(ID))	ID=1:length(actvec)
			#stratified sampling by actvec
			set.seed(srnd)
			IDsbyact=tapply(ID,actvec,sample)
			
			#sample proportion=nfold as external set
			if(nfold>0 && nfold<1){ 
				fmod=1-nfold
				IDsbyactavail=list()
				for(k in 1:length(IDsbyact)) IDsbyactavail[[k]]=sort(sample(IDsbyact[[k]],size=floor(fmod*length(IDsbyact[[k]]))))
				#no under-/over- sampling
				modID=sort(as.numeric(unlist(IDsbyactavail)))
				extID=setdiff(ID,modID)
				
				if(!is.null(foldoverminclass)){
					IDsbyact=list()
					classsize=unlist(lapply(IDsbyactavail,length))
					minclass=which.min(classsize)
					otherclass=setdiff(1:length(classsize),minclass)
					sampsize=c()
					#undersample majority class such that majority_class/minority_class=foldoverminclass
					if(foldoverminclass>=0){
						sampsize[minclass]=classsize[minclass]
						sampsize[otherclass]=floor(foldoverminclass*sampsize[minclass])
						for(i in 1:length(classsize)) IDsbyact[[i]]=sample(IDsbyactavail[[i]],size=sampsize[i])			
					}
					modID=sort(as.numeric(unlist(IDsbyact)))
				}
				
				return(list(modID=modID,extID=extID))
				
				#nfold-cross validation
			}else if(nfold==as.integer(nfold) && nfold>=2){
				#create matrix of nfold rows with extIDs in each row
				shuffled_stratID=as.numeric(unlist(IDsbyact))		
				#Calculate number of compounds per fold and set dimensions of extIDmat
				remainder=length(ID)%%nfold
				if(remainder==0){
					extIDvec=shuffled_stratID
				}else{
					extIDvec=c(shuffled_stratID,rep(NA,(nfold-remainder)))
				}
				temp=matrix(extIDvec,byrow=F,nrow=nfold)
				extIDmat=t(apply(temp,1,sort,na.last=T))
				
				#no under-/over- sampling
				if(is.null(foldoverminclass)){
					return(extID=extIDmat)
				}else{
					modIDlist=extIDlist=list()
					for(j in 1:nfold){
						modIDavail=sort(na.omit(as.vector(extIDmat[-j,])))
						IDsbyactavail=tapply(modIDavail,actvec[modIDavail],sample)
						classsize=unlist(lapply(IDsbyactavail,length))
						minclass=which.min(classsize)
						otherclass=setdiff(1:length(classsize),minclass)
						sampsize=c()
						#undersample majority class such that majority_class/minority_class=foldoverminclass
						if(foldoverminclass>=0){
							sampsize[minclass]=classsize[minclass]
							sampsize[otherclass]=floor(foldoverminclass*sampsize[minclass])
							for(i in 1:length(classsize)) IDsbyact[[i]]=sample(IDsbyactavail[[i]],size=sampsize[i])			
							
							if(!replace && j>1){
								IDsbyactavail_noreplace=setdiff(IDsbyactavail[[otherclass]],unlist(modIDlist))
								if(sampsize[otherclass]<=length(IDsbyactavail_noreplace)){
									IDsbyact[[otherclass]]=sample(IDsbyactavail_noreplace,size=sampsize[otherclass])
								}else{
									IDsbyact[[otherclass]]=c(sample(IDsbyactavail_noreplace),
															 sample(IDsbyactavail[[otherclass]],size=sampsize[otherclass]-length(IDsbyactavail_noreplace)))
								}
							}
						}
						modIDlist[[j]]=sort(as.numeric(unlist(IDsbyact)))
						extIDlist[[j]]=as.vector(na.omit(extIDmat[j,]))
					}
					return(list(modID=modIDlist,extID=extIDlist))
				}
			}else{
				stop("nfold must be an integer above 1 for nfold-cross validation, or between 0 and 1 for a single external set of proportion=nfold" )
			}
		}else{
			stop("nfold must be an integer above 1 for nfold-cross validation, or between 0 and 1 for a single external set of proportion=nfold" )
	}
}

extID2txt<-function(extIDobj,outputroot){
	outputext=paste(outputroot,"extID.txt",sep="")
	outputmod=paste(outputroot,"modID.txt",sep="")
	cat("remove files?",outputext,outputmod,sep="\n")
	print(file.remove(outputext,outputmod))
	for(k in 1:length(temp$extID)){
		cat(extIDobj$extID[[k]],"\n",file=outputext,sep="\t",append=T)
		cat(extIDobj$modID[[k]],"\n",file=outputmod,sep="\t",append=T)
	}
}
