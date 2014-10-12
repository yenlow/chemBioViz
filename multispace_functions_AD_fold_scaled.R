# functions for Multi-space kNN read across
# requires multispaceNNobj.R, ReadXAfiles.R, sampling.R
# Changed scaling factor from theoretical max Ts to actual sum of Ts
# Added varimpt to calculate variable importance
# 
# 10-Mar-13 Yen Low
###############################################################################
if (!is.loaded("caret")) library(caret)
if (!is.loaded("vegan")) library(vegan)
if (!is.loaded("ROCR")) library(ROCR)

source("readXAfile.R")
source("sampling.R")
source("multispaceNNobj.R")

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
#		outCM (name of output file containing confusion matrices showing model performance)
#		outvalidstat (name of output file showing model performance)
#		similaritycutoff (minimum similarity that external cpds are to training cpds; default=0.3)
#		nfoldint (number of folds for internal cross validation for tuning kchem, kbio values; recommended nfoldint=10 (default))
#		calcErr (logical. Set to TRUE (default) if bootstrapping errors are desired. Set to FALSE to speed up)
#		ntrials (number of bootstrapping samples (1000 trials by default). Set to 100 to speed up)
#output: list object of 1. predicted values of external compounds (predtable)
#						2. prediction performance on external compounds (perf)
#						3. model object with external compounds' neighbors, neighbors' distances and classification labels (extobj)
#
#Pls be patient! Takes a few minutes on dual-core processors.
#Loops through   5-fold external cross validation
#within each is 10-fold internal cross validation (gridsearch for best kchem and kbio) - Set nfoldint=5 instead of 10 to speed up
dualspacekNN<-function(krange,k2range=0,extID,classlab,chemspace,genespace=NULL,
		cpdnames=cpdnames,similaritycutoff=0.3,nfoldint=10,
		optcriterion1="CCR",optcriterion2="tst",
		outpred="out.pred",outCM=NA,
		outvalidstat="validationstats.txt",calcErr=T,ntrials=ntrials){
	
	#Build dualspace-kNN models usinc all kchem values in krange and all kbio values in k2range
	modobj=modworkflow(krange=krange,k2range=k2range,extID=extID,dsact=classlab,
			chemspace=chemspace,genespace=genespace,nfoldint=nfoldint,
			cpdnames=cpdnames,similaritycutoff=similaritycutoff)
	
	#get best model with optimal kchem, kbio values yielding best test set CCR
	optkstats=getopttable(modobj$statarray,optcriterion1,optcriterion2)
	
	#apply best model (optimal k values) to external set (loop over 5 folds)
	extobj=getextvalobj(optkstats,extID,modobj$chemNNext,modobj$geneNNext,cpdnames=cpdnames,similaritycutoff=similaritycutoff)
	
	#get predictions for each external compound
	predtable=createpredtable(extobj,classlab,cpdnames,outfile=outpred)
	
	#get prediction performance statistics on external set (per fold and cumulative)
	stats=validationstats_BIN(predtable,calcErr=calcErr,ntrials=ntrials,outputfile=outCM,validationstatsfile=NA)
	statswithk=concatkvalidstat(optkstats,stats,outputfile=outvalidstat)
	print(statswithk)
	list(predtable=predtable,perf=statswithk,extobj=extobj)
}


#create NAsum function for wtedact because sum with na.rm=T is retarded
NAsum<-function(x){
    if(all(is.na(x))) NA else return(sum(x,na.rm=T))
}

#calculate similarity weighted average of activites
#returns prodmat, sums, OUTofAD
wtedact<-function(kNNdist,kNNact,similaritycutoff=0.6){
    #initialize empty matrices
    NAid=c()
    temp=kNNact[,-1]
    mode(temp)<-"numeric"
    kNNactmat=temp
    #set nontoxic to "-1" instead
    kNNactmat[temp==0]<--1
    
    kNNsimi=1-kNNdist
    prodmat=as.matrix(kNNsimi*kNNactmat)
    OUTofAD=rep(0,nrow(as.matrix(kNNdist)))

    #apply AD (similarity > similarity cutoff)
    if(is.na(similaritycutoff)){ #no similarity cutoff or AD will be applied
        sums=apply(prodmat,1,sum,na.rm=T)   
    }else if(similaritycutoff>=0 & similaritycutoff<=1){
        #set dissimilar neighbours above similairtycutoff to NA ("out of AD")
        kNNsimi_AD=kNNsimi
        kNNsimi_AD[kNNsimi<similaritycutoff]=NA
        prodmat_AD=as.matrix(kNNsimi_AD*kNNactmat)
        sums=apply(prodmat_AD,1,NAsum)
		scaling_denom=apply(abs(prodmat_AD),1,NAsum)
        
        #replace NAs in sums with sums as if without AD
        NAid=which(is.na(sums))
        #need to ensure prodmat is a matrix even with 1 col or 1 row
        sums[NAid]=apply(matrix(prodmat[NAid,],byrow=F,nrow=length(NAid),ncol=ncol(prodmat)),1,sum,na.rm=T)
		scaling_denom[NAid]=apply(matrix(abs(prodmat[NAid,]),byrow=F,nrow=length(NAid),ncol=ncol(prodmat)),1,sum,na.rm=T)
		OUTofAD[NAid]=1
		sum_scaled=sums/scaling_denom
		
    }else{
        print("Error: Enter Tanimoto similaritycutoff. Possible values: NA, between 0 and 1 inclusive, default to 0.6")
    }

    list(prodmat=prodmat,sums=sum_scaled,OUTofAD=OUTofAD)
}

#creates kNNdist,kNN,kNNact but do not return them
#calls wtedact
optk<-function(krange,k2range=0,ID,chemNN,geneNN=NULL,cpdnames,mode="range",NNlist=NULL,tunecutoff=FALSE,actcutoff=0.5,similaritycutoff=0.6){
    #initialize variables
    sens=c()
    spec=c()
    CCR=c()
    acc=c()
    auc=c()
    kval=c()
    k2val=c()
    cutofflist=c()
    row_names=cpdnames
    counter=1
    
    #creates kNNdist, kNN, kNNact to feed into wtedact
    for(k2 in k2range){
        for(k in krange){
           if(is.null(NNlist)){ #if inputs are chemNN with or without geneNN
                kNNdist=cbind(chemNN$NNdist[,2:(k+1)],geneNN$NNdist[,2:(k2+1)])
                kNN=cbind(chemNN$NN[,2:(k+1)],geneNN$NN[,2:(k2+1)])
                kNNact=cbind(chemNN$NNact[,1:(k+1)],geneNN$NNact[,2:(k2+1)])
            }else { #if input is a NNlist
                kNNdist=NNlist[[1]]$NNdist[,2:(k+1)]
                kNN=NNlist[[1]]$NN[,2:(k+1)]
                kNNact=NNlist[[1]]$NNact[,1:(k+1)]
                for (j in 2:length(NNlist)){
                    kNNdist=cbind(kNNdist,NNlist[[j]]$NNdist[,2:(k+1)])
                    kNN=cbind(kNN,NNlist[[j]]$NN[,2:(k+1)])
                    kNNact=cbind(kNNact,NNlist[[j]]$NNact[,2:(k+1)])
                }
            }
            rownames(kNNdist)=row_names
            rownames(kNN)=row_names
            rownames(kNNact)=row_names
            yact=kNNact[ID,1]
            ans=wtedact(kNNdist[ID,],kNNact[ID,],similaritycutoff=similaritycutoff)
            pred=0.5+ans$sums/2
            predobj=prediction(pred,factor(yact,levels=0:1))
            #perfobj=performance(predobj,"tpr","fpr")
            #plot(perfobj)
            if(tunecutoff==TRUE){
                accobj=performance(predobj,"acc")
                acutoff=round(accobj@x.values[[1]][which.max(accobj@y.values[[1]])],1)
            }else{
                acutoff=actcutoff
            }
            ypred=pred
            ypred[pred<acutoff]=0
            ypred[pred>=acutoff]=1
            cutofflist[counter]=acutoff
            cM=confusionMatrix(factor(ypred,levels=0:1),factor(yact,levels=0:1),"1")
            sens[counter]<-cM$byClass["Sensitivity"]
            spec[counter]<-cM$byClass["Specificity"]
            CCR[counter]<-(sens[counter]+spec[counter])/2
            acc[counter]=cM$overall["Accuracy"]
			tryobj=try(performance(predobj,"auc"))
			if (class(tryobj) != "try-error"){
				auc[counter]=tryobj@y.values[[1]]
			}else{
				auc[counter]=NA
				print("Insufficient predictions for AUC calculation")
			}
            kval[counter]=k
            k2val[counter]=k2
            counter=counter+1
        } #out of krange loop
    } #out of k2range loop
       
    stat=cbind(kval,k2val,spec,sens,CCR,acc,auc,cutofflist)
    colnames(stat)[1:2]=c("k","k2")
    
    if(mode=="range") { return(stat)
    } else if(mode=="pred") { list(stat=stat,mod=ans,pred=pred,OUTofAD=ans$OUTofAD)
    } else {print('Error: Enter mode. Possible values:"range" (for selecting optimal k),"pred" (for generating predictions from point k values)')
    }
}

#Run multi-space kNN using optk based on extID matrix
#Calls optk(mode="range")
modworkflow<-function(krange,k2range=0,extID,varID_c=NULL,varID_g=NULL,chemspace,genespace=NULL,cpdnames,
		tunecutoff=FALSE,actcutoff=0.5,nfoldint=10,similaritycutoff=0.3,dsact=dsact){

	#assign varID_c to full matrix if NULL
	if(is.null(varID_c)){
		varID_c=matrix(colnames(chemspace),ncol=nrow(extID),nrow=ncol(chemspace))
	}else if(!is.matrix(varID_c)){
		print("varID_c must be a matrix of variable names with nfold columns and nvariable rows")
	}
	if(is.null(varID_g) & !is.null(genespace)){
		varID_g=matrix(colnames(genespace),ncol=nrow(extID),nrow=ncol(genespace))
	}else if(!is.null(varID_g) & !is.matrix(varID_g)){
		print("varID_g must be a matrix of variable names with nfold columns and nvariable rows")
	}
	
	if(!is.integer(krange) & !(isTRUE(all.equal(krange,as.integer(krange))))) print("krange must be integer(s) e.g. 0 or 1:5")
	if(!is.integer(k2range) & !(isTRUE(all.equal(k2range,as.integer(k2range))))) print("k2range must be integer(s) e.g. 0 or 1:5")
	if(!is.matrix(extID)) print("extID must be a matrix with each row representing an external fold. Generate matrix with genextID")
	if(!is.numeric(chemspace)) print("chemspace must be a numeric matrix with ncpd rows and ndescriptor columns")
	if(!is.vector(cpdnames)) print("cpdnames must be a vector with ncpd names")
	if(!is.vector(dsact)) print("dsact must be a vector with ncpd activity values")
	if(!is.logical(tunecutoff)) print("tunecutoff must be boolean (TRUE/FALSE)")
	if(!(isTRUE(all.equal(nfoldint,as.integer(nfoldint))))) print("nfoldint must be an integer determining number of folds for internal cross-validation; dafaults to 10")
	if(!is.numeric(similaritycutoff) | similaritycutoff>1 | similaritycutoff<0) print("similaritycutoff must be a number between 0 and 1 inclusive")
	if(!is.numeric(actcutoff) | actcutoff>1 | actcutoff<0) print("actcutoff must be a number between 0 and 1 inclusive")
	
	print("Process may take a few minutes. Please wait...")
	flush.console()	
	
	chemNNextall=list()
	geneNNextall=list()
	
	stattrn=list()
	stattst=list()
	statext=list()
	foldname=c()

	#loop over modeling folds (ext CV loop)
	for(i in 1:nrow(extID)){
		cat(nrow(extID),"-fold external cross validation: Fold ",i,"...\n",sep="")
		flush.console()
		extIDperfold=na.omit(extID[i,])
		modID=as.vector(na.omit(as.vector(t(extID[-i,]))))
		gentestID=genextID(ID=modID,actvec=dsact[modID],nfold=nfoldint) #default 10-fold
		
		inttrn=list()
		inttst=list()
		
		#loop over training/test sets of modeling fold (int CV)
		for(j in 1:nfoldint){
			testID=na.omit(gentestID[j,])
			trainID=na.omit(as.vector(gentestID[-j,]))
			
			if(nfoldint==1){
				testID=NULL
				trainID=modID		
			}
			
			#calculate NN obj of only training data
			dschem_trn=chemspace[trainID,as.character(varID_c[,i])]
			dsgene_trn=genespace[trainID,as.character(varID_g[,i])]
			chemNNtrn=NNmat(dschem_trn,cpdnames[trainID],dsact[trainID])
			geneNNtrn=NNmat(dsgene_trn,cpdnames[trainID],dsact[trainID])
		
			#loop over test cpds 
			#(add each test cpd, one at a time, to training data for NN analysis)
			chemNN_tstcpd=list()
			geneNN_tstcpd=list()
			chemNNtst=list()
			geneNNtst=list()
			
			if(!is.null(testID)){
			for(l in 1:length(testID)){
				#calculate NN obj of training data and lth test cpd
				testcpd=testID[l]
				dschem_tst=chemspace[c(trainID,testcpd),as.character(varID_c[,i])]
				chemNN_trntst=NNmat(dschem_tst,cpdnames[c(trainID,testcpd)],dsact[c(trainID,testcpd)])
								
				#extract last row from every list element. The rows become vectors (columns)
				chemNN_tstcpd[[l]]=lapply(chemNN_trntst,`[`,length(trainID)+1,)
				
				if(!is.null(genespace)){
					dsgene_tst=genespace[c(trainID,testcpd),as.character(varID_g[,i])]
					geneNN_trntst=NNmat(dsgene_tst,cpdnames[c(trainID,testcpd)],dsact[c(trainID,testcpd)])
					geneNN_tstcpd[[l]]=lapply(geneNN_trntst,`[`,length(trainID)+1,)
				}
			}
		
			
			#reassemble NN by testIDs (NNID is re-ordered, ignore NNID)
			for(m in 2:length(chemNN_tstcpd[[1]])){
				chemNNtst[[m]]=matrix(	unlist(lapply(chemNN_tstcpd,`[`,m)),
										nrow=length(testID),byrow=T)
				if(!is.null(genespace)){
					geneNNtst[[m]]=matrix(	unlist(lapply(geneNN_tstcpd,`[`,m)),
											nrow=length(testID),byrow=T)
				}
			}
			names(chemNNtst)=names(chemNNtrn)
			if(!is.null(genespace))	names(geneNNtst)=names(geneNNtrn)

			inttst[[j]]=optk(krange=krange,k2range=k2range,ID=1:length(testID),chemNN=chemNNtst,geneNN=geneNNtst,cpdnames=cpdnames[testID],similaritycutoff=similaritycutoff)
			
		}
			#Gridsearch call opt k 
			inttrn[[j]]=optk(krange=krange,k2range=k2range,ID=1:length(trainID),chemNN=chemNNtrn,geneNN=geneNNtrn,cpdnames=cpdnames[trainID],similaritycutoff=similaritycutoff)
			if(is.null(testID))	inttst[[j]]=inttrn[[j]]
		}	#out of internal cross validation loop
		
		
		#put int val stattrn and stattst into a cvarray (averaged across internal folds)
		cvarray=array(unlist(c(inttrn,inttst)),dim=c(dim(inttst[[1]]),nfoldint,2),
				dimnames=list(1:dim(inttst[[1]])[1],
						c("k","k2","spec","sens","CCR","acc","auc","cutoff"),
						1:nfoldint,c("trn","tst")))
		#if nfoldint=1, then 4D statarray collapses to a 3D array
		if(dim(cvarray)[3]==1){
			stattrn[[i]]=cvarray[,,,"trn"]
			stattst[[i]]=cvarray[,,,"tst"]
		}else if(dim(cvarray)[1]==1){ 
			#if k, k2 are a single point, then 4D statarray collapses to a 3D array
			stattrn[[i]]=apply(cvarray[,,,"trn"],1,mean)
			stattst[[i]]=apply(cvarray[,,,"tst"],1,mean)
		}else{
			stattrn[[i]]=apply(cvarray[,,,"trn"],c(1,2),mean)
			stattst[[i]]=apply(cvarray[,,,"tst"],c(1,2),mean)
		}
		
		#calculate NN obj of only modeling set data
		dschem_mod=chemspace[modID,as.character(varID_c[,i])]
		dsgene_mod=genespace[modID,as.character(varID_g[,i])]
		chemNNmod=NNmat(dschem_mod,cpdnames[modID],dsact[modID])
		geneNNmod=NNmat(dsgene_mod,cpdnames[modID],dsact[modID])
		
		#loop over ext cpds 
		#(add each ext cpd, one at a time, to modeling data for NN analysis)
		chemNN_extcpd=list()
		geneNN_extcpd=list()
		chemNNext=list()
		geneNNext=list()
		for(p in 1:length(extIDperfold)){
			#calculate NN obj of modeling data and pth ext cpd
			extcpd=extIDperfold[p]
			dschem_ext=chemspace[c(modID,extcpd),as.character(varID_c[,i])]
			chemNN_modext=NNmat(dschem_ext,cpdnames[c(modID,extcpd)],dsact[c(modID,extcpd)])
			#extract last row from every list element. The rows become vectors (columns)
			chemNN_extcpd[[p]]=lapply(chemNN_modext,`[`,length(modID)+1,)
			
			if(!is.null(genespace)){
				dsgene_ext=genespace[c(modID,extcpd),as.character(varID_g[,i])]
				geneNN_modext=NNmat(dsgene_ext,cpdnames[c(modID,extcpd)],dsact[c(modID,extcpd)])
				geneNN_extcpd[[p]]=lapply(geneNN_modext,`[`,length(modID)+1,)
			}
		}
		
		#reassemble NN by extIDs (NNID is re-ordered, ignore NNID)
		for(q in 2:length(chemNN_extcpd[[1]])){
			chemNNext[[q]]=matrix(	unlist(lapply(chemNN_extcpd,`[`,q)),
									nrow=length(extIDperfold),byrow=T)
			if(!is.null(genespace)){
				geneNNext[[q]]=matrix(	unlist(lapply(geneNN_extcpd,`[`,q)),
										nrow=length(extIDperfold),byrow=T)
			}
		}
		names(chemNNext)=names(chemNNtrn)
		if(!is.null(genespace)) names(geneNNext)=names(geneNNtrn)
		
		#save chemNNext (one fold) into super list for all folds
		chemNNextall[[i]]=chemNNext
		geneNNextall[[i]]=geneNNext
		
		#run on external sets
		statext[[i]]=optk(krange=krange,k2range=k2range,ID=1:length(extIDperfold),
				chemNN=chemNNext,geneNN=geneNNext,cpdnames=cpdnames[extIDperfold],
				similaritycutoff=similaritycutoff)
		foldname[i]=paste("fold",i-1,sep="")
		
	}	#out of external cross validation loop
	
	#convert stats list to arrays for easy calling
	#k,stat,fold,trn/tst/ext
	#dim(stattrn[[1]])
	statarray=array(unlist(c(stattrn,stattst,statext)),
					dim=c(max(1,dim(stattrn[[1]])[1]),8,nrow(extID),3),
					dimnames=list(1:max(1,dim(stattrn[[1]])[1]),c("k","k2","spec","sens","CCR","acc","auc","cutoff"),foldname,
								c("trn","tst","ext")))
	
	list(statarray=statarray,chemNNext=chemNNextall,geneNNext=geneNNextall)
}

#get optimal k, CCRtrn, CCRtst, CCRext
getopttable<-function(statarray,ptype="CCR",optset="tst"){
    #if k, k2 are a single point, then 4D statarray collapses to a 3D array
    if(dim(statarray)[1]==1){ 
        kID=rep(1,dim(statarray)[3])
    }else{
        kID=apply(statarray[,ptype,,optset],2,which.max)
    }
    optrow=list()
    for (i in 1:dim(statarray)[3]) optrow[[i]]=statarray[kID[i],ptype,i,]
    temp=matrix(unlist(optrow),nrow=dim(statarray)[3],ncol=3,byrow=T,dimnames=dimnames(statarray)[3:4])
    opttable=cbind(statarray[kID,1:2,1,1],temp)
    rownames(opttable)=rownames(temp)
    return(opttable)
}

#get ext set predictions
#calls optk(mode="pred")
#returns stat, mod, pred, OUTofAD from optk (nested within lists)
getextvalobj<-function(optkstats,extID,chemNN,geneNN=NULL,cpdnames,similaritycutoff=0.6){
    optkval=optkstats[,1:2]
    statext=list()
    for(i in 1:nrow(extID)){
		extIDperfold=na.omit(extID[i,])
		statext[[i]]=list(optk(optkval[i,1],optkval[i,2],ID=1:length(extIDperfold),
								chemNN=chemNN[[i]],geneNN=geneNN[[i]],
								cpdnames=cpdnames[extIDperfold],mode="pred",
								tunecutoff=F,similaritycutoff=similaritycutoff))
    }
    return(statext)
}

createpredtable <- function(extobj,dsact,cpdnames,outfile=NA) {
    pred=list()
    fold=list()
    nfold=length(extobj)
    OUTofAD=list()
    
    #extract information from each fold
    #predtable must have fold (must start with 0),yobs,ypred
    #cpdnames as rownames
    for(i in 1:nfold){
        pred[[i]]=extobj[[i]][[1]]$pred
        fold[[i]]=rep(i-1,length(pred[[i]]))
        OUTofAD[[i]]=extobj[[i]][[1]]$OUTofAD
    }
    ypred=round(unlist(pred))
    foldypred=as.data.frame(cbind(unlist(fold),ypred,unlist(pred),unlist(OUTofAD)))
    colnames(foldypred)=c("fold","ypred","rawpred","OUTofAD")
    
    yobs=as.data.frame(dsact,row.names=cpdnames)
    colnames(yobs)="yobs"
    temp=merge(foldypred,yobs,by=0,all.y=T,sort=F)  
    predtable=temp[,c("fold","ypred","yobs","OUTofAD","rawpred")]
    rownames(predtable)=temp[,"Row.names"]
    
	if(!is.na(outfile)) write.table(predtable,file=outfile,sep="\t",row.names=T,col.names=NA,quote=F)
    return(predtable)
}

concatkvalidstat<-function(optkstats,validationstats,outputfile="validationstats.txt"){
    kwithvalidstat=merge(optkstats[,c("k","k2")],validationstats,by=0,all.y=T,sort=F)
    colnames(kwithvalidstat)[1]=" "
    write.table(kwithvalidstat,file=outputfile,row.names=F,,col.names=T,quote=F,sep="\t")
    return(kwithvalidstat)
}

varimpt=function(nperm,extobj,optkstats,extID,chemspace,genespace=NULL,rdmvar=NULL,varID_c=NULL,varID_g=NULL,dsact,cpdnames){
	names(dsact)=cpdnames
	if(is.null(rdmvar)) rdmvar=colnames(chemspace)
	if(is.null(genespace)) k2range=0
	else k2range=5
	imp=matrix(NA,nperm,length(rdmvar))
	localimpSim=localimp=array(NA,dim=c(nperm,length(rdmvar),nrow(chemspace)))
	for(u in 1:nperm){	#repeat number of permutations
		
		for(q in 1:length(rdmvar)){
			set.seed(sample(.Random.seed[-(1:2)],1))
			chemspace_rdm=chemspace
			chemspace_rdm[,rdmvar[q]]=sample(chemspace[,rdmvar[q]])
			
			statsall_rdm=modworkflow(5,k2range=k2range,extID=extID,varID_c=varID_c,varID_g=varID_g,
					chemspace=chemspace_rdm,genespace=genespace,cpdnames=cpdnames,tunecutoff=F,
					similaritycutoff=0.3,nfoldint=1,dsact=dsact)
			extobj_rdm=getextvalobj(optkstats,extID,statsall_rdm$chemNNext,statsall_rdm$geneNNext,cpdnames=cpdnames,similaritycutoff=0.3)
			
			DecInCCR=c()
			localIncInError=list()
			localDecInSim=list()
			for(s in 1:nrow(extID)){
				DecInCCR[s]=extobj[[s]][[1]]$stat[,"CCR"]-extobj_rdm[[s]][[1]]$stat[,"CCR"]
					yobs=as.numeric(as.character(dsact[names(extobj[[s]][[1]]$pred)]))
					localIncInError[[s]]=abs(extobj_rdm[[s]][[1]]$pred-yobs)-
										 abs(extobj[[s]][[1]]$pred-yobs)
					localDecInSim[[s]]=rowMeans(abs(extobj[[s]][[1]]$mod$prodmat)-abs(extobj_rdm[[s]][[1]]$mod$prodmat))
			}
			imp[u,q]=mean(DecInCCR,na.rm=T)
			localimp[u,q,]=unlist(localIncInError)
			localimpSim[u,q,]=unlist(localDecInSim)
		}	
	}
	colnames(imp)=rdmvar
	dimnames(localimp)[2:3]=list(rdmvar,cpdnames[na.omit(as.vector(t(extID)))])
	list(imp=imp,localimp=localimp,localimpSim=localimpSim)
}

plotImp<-function(mat,labels,topn=NULL,title=NULL){
	if(is.null(topn)) topn=ncol(mat)
	sortedmeans=sort(colMeans(mat),dec=T)
	matchID=match(names(sortedmeans),colnames(mat)) 
	stdev=apply(mat[,matchID],2,sd)
	upplimit=qt(0.975,nrow(mat))*stdev/sqrt(nrow(mat))
	if(!is.null(labels)) xlabels=labels[matchID]
	else xlabels=colnames(mat)[matchID]
	plotCI(sortedmeans,uiw=upplimit,pch=20,xlim=c(1,topn),
			ylab="importance",xaxt="n",xlab="",main=title)
	axis(1,at=1:topn,lab=xlabels[1:topn],las=2,cex.axis=0.7)
	abline(h=0,lty="dotted")
}



