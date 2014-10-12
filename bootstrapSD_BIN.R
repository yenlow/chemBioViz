#Calculate validation stats for classification models
#sensitivity, specificity, CCR, coverage (per fold and cumulative) 
#22-Apr-12 Yen Low
#Inputs:	1. predtable (must have columns "fold","ypred","yobs","OUTofAD")
#		2. nfolds (number of folds)
#Output: 1. confusionMatrices.txt
#		2. validationstats (dataframe)
###########################################################

if (!is.loaded("caret")) library(caret)
if (!is.loaded("ROCR")) library(ROCR)
if (!is.loaded("boot")) library(boot)

#write function to calculate CCR where d is a dynamic ID used by boot
ccrfun<-function(predtable,d,yobscolname="yobs",ypredcolname="ypred"){
	obs=predtable[d,yobscolname]
	pred=predtable[d,ypredcolname]
    ct=table(factor(pred,levels=0:1),factor(obs,levels=0:1))
	spec=specificity(ct,"0")
	sens=sensitivity(ct,"1")
	CCR=(spec+sens)/2
	return(CCR)
}

sensfun<-function(predtable,d,yobscolname="yobs",ypredcolname="ypred"){
    obs=predtable[d,yobscolname]
    pred=predtable[d,ypredcolname]
    ct=table(factor(pred,levels=0:1),factor(obs,levels=0:1))
    sens=sensitivity(ct,"1")
    return(sens)
}

specfun<-function(predtable,d,yobscolname="yobs",ypredcolname="ypred"){
    obs=predtable[d,yobscolname]
    pred=predtable[d,ypredcolname]
    ct=table(factor(pred,levels=0:1),factor(obs,levels=0:1))
    spec=specificity(ct,"0")
    return(spec)
}

accfun<-function(predtable,d,yobscolname="yobs",ypredcolname="ypred"){
    obs=predtable[d,yobscolname]
    pred=predtable[d,ypredcolname]
    cM=confusionMatrix(table(factor(pred,levels=0:1),factor(obs,levels=0:1)))
    acc=cM$overall["Accuracy"]
    return(acc)
}

#calculate AUC if rawpred values exists
#else AUC is meaningless (equivalent to CCR, assuming cutoff=0.5)
aucfun<-function(predtable,d,yobscolname="yobs",rawpredcolname="rawpred"){
    if(rawpredcolname %in% colnames(predtable)){
        obs=predtable[d,yobscolname]
        rawpred=predtable[d,rawpredcolname]
        predobj=try(prediction(rawpred,obs))
        if(class(predobj)!="try-error"){
            auc=performance(predobj,"auc")@y.values[[1]]
        }else{
            auc=NA
			print("Error in bootstrapSD_BIN.R line 63")
        }
        return(auc)
    }else{
        print("Error: When Ypred is not a continuous value, AUC is equivalent to CCR")
        return(NA)
    }
}
    
#Calculate original CCR
#ccrfun(predtable,1:nrow(predtable))
#sensfun(predtable,1:nrow(predtable))
#specfun(predtable,1:nrow(predtable))
#accfun(predtable,1:nrow(predtable))
#aucfun(predtable_c[predtable_c$OUTofAD==0,],1:nrow(predtable_c))

#Bootstrap CCR
getbootstrapErr<-function(predtable,fun,ntrials=1000,verbose=FALSE){
	b=try(boot(predtable,fun,R=ntrials))
	if(class(b)!="try-error"){
		if(verbose) print(b)
		flush.console()
		err=sd(as.vector(b$t))
	}else{
		err=NA
		print("Error in bootstrapSD_BIN.R line 81. No values for input")
	}
    return(err)
}


