#Calculate validation stats for classification models
#sensitivity, specificity, CCR, coverage (per fold and cumulative) 
#10-Mar-13 Yen Low
#Inputs:1. predtable (must have columns "fold","ypred","yobs") optional:"rawpred","OUTofAD"
#		2. calcErr=TRUE (if bootstrap errors are desired; set to F to speed up)
#		3. ntrials=1000 (number of bootstrap samples for calculate errors; set to 100 to speed up)
#Output:returns validationstats
#		1. confusionMatrices.txt  (suppressed with outputfile=NA)
#		2. validationstats.txt (suppressed with validationstatsfile=NA)
#######################################################################

if (!is.loaded("caret")) library(caret)
if (!is.loaded("ROCR")) library(ROCR)

source("bootstrapSD_BIN.R")

validationstats_BIN<-function(predtable,calcErr=TRUE,ntrials=1000,
								outputfile="confusionMatrices.txt",validationstatsfile="validationstats.txt"){
#initialize output variables
    sens=c()
    spec=c()
    CCR=c()
    acc=c()
    auc=c()
    
    sens_withAD=c()
    spec_withAD=c()
    CCR_withAD=c()
    acc_withAD=c()
    auc_withAD=c()
    coverage=c()
    
    foldnames=c()
    
    nfolds=length(table(predtable$fold))
    for(i in 1:nfolds){
        
        foldnames[i]=paste("fold",i-1,sep="")
        
        predtable_byfold=subset(predtable,fold==(i-1))
        predlevels=factor(predtable_byfold$ypred,levels=0:1)
        obslevels=factor(predtable_byfold$yobs,levels=0:1)        
        cM<-try(confusionMatrix(predlevels,obslevels,"1"))
		if (class(cM) != "try-error"){
			sens[i]=cM$byClass["Sensitivity"]
			spec[i]=cM$byClass["Specificity"]
			CCR[i]=(sens[i]+spec[i])/2
			acc[i]=cM$overall["Accuracy"]
		}else{
			sens[i]=NA
			spec[i]=NA
			CCR[i]=NA
			acc[i]=NA
		}
        #calculate AUC if rawpred values exists
        #else AUC is meaningless (equivalent to CCR, assuming cutoff=0.5)
        if(!is.null(predtable$rawpred)){
            predobj=try(prediction(predtable_byfold$rawpred,predtable_byfold$yobs))
			if(class(predobj) != "try-error"){
				temp2=try(performance(predobj,"auc"))
				if(class(temp2) != "try-error"){
					auc[i]=performance(predobj,"auc")@y.values[[1]]
				}else{
					auc[i]=NA
					print("Error in line 67: AUC will be NA")
				}
			}else{
				auc[i]=NA
				print("Error in line 71: AUC will be NA")
			}
        }
        
        if(!is.na(outputfile)){
			if(file.exists(outputfile)) file.remove(outputfile)
			sink(outputfile,append=T)
				cat("\n","fold ",i,": without AD","\n")
				print(cM)
				print(CCR[i])
				print(auc[i])
			sink()
		}
        
        #for AD calculations per fold
        #be careful when fewer than 2 classes after excluding out of AD
        if(!is.null(predtable$OUTofAD)){
            coverage[i]=1-sum(predtable_byfold$OUTofAD)/length(predtable_byfold$OUTofAD)
            
            predtable_byfold_withAD=subset(predtable_byfold,OUTofAD==0)
            predlevels_withAD=factor(predtable_byfold_withAD$ypred,levels=0:1)
            obslevels_withAD=factor(predtable_byfold_withAD$yobs,levels=0:1)        
            cM_withAD<-try(confusionMatrix(predlevels_withAD,obslevels_withAD,"1"))
            if (class(cM_withAD) != "try-error"){
                sens_withAD[i]<-cM_withAD$byClass["Sensitivity"]
                spec_withAD[i]<-cM_withAD$byClass["Specificity"]
                CCR_withAD[i]<-(sens_withAD[i]+spec_withAD[i])/2
                acc_withAD[i]=cM_withAD$overall["Accuracy"]
            }else{
                sens_withAD[i]=NA
                spec_withAD[i]=NA
                CCR_withAD[i]=NA
                acc_withAD[i]=NA
                print("Error in line 101: Fewer than 2 classes after excluding out of AD compounds. Performance metrics set to NA")
            }
            #calculate AUC if rawpred values exists
            #else AUC is meaningless (equivalent to CCR, assuming cutoff=0.5)
            if(!is.null(predtable$rawpred)){
                predobj_withAD=try(prediction(predtable_byfold_withAD$rawpred,predtable_byfold_withAD$yobs))
                #in case only 1 class after too many exclusions (gives error)
                if (class(predobj_withAD) != "try-error"){
                    temp3=try(performance(predobj_withAD,"auc"))
					if(class(temp3) != "try-error"){
						auc_withAD[i]=performance(predobj_withAD,"auc")@y.values[[1]]
					}else{
						auc_withAD[i]=NA
						print("Error in line 114: AUC will be NA")						
					}
                }else{
                    auc_withAD[i]=NA
                    print("Error in line 118: AUC will be NA")
                }

            }
            
            if(!is.na(outputfile)){
				sink(outputfile,append=T)
					cat("fold ",i,": with AD","\n")
					print(cM_withAD)
					print(CCR_withAD[i])
					print(auc_withAD[i])
				sink() 
			}
        }	#out of if block for AD calculations
    } 	#out of nfolds loop for per fold calculations
    
    
#calculate cumulative stats
    cMall<-try(confusionMatrix(factor(predtable$ypred,levels=0:1),
                                factor(predtable$yobs,levels=0:1),"1"))
	if(class(cMall) != "try-error"){
		sens[nfolds+1]=cMall$byClass["Sensitivity"]
		spec[nfolds+1]=cMall$byClass["Specificity"]
		CCR[nfolds+1]=(sens[nfolds+1]+spec[nfolds+1])/2
		acc[nfolds+1]=cMall$overall["Accuracy"]
	}else{
		sens[nfolds+1]=NA
		spec[nfolds+1]=NA
		CCR[nfolds+1]=NA
		acc[nfolds+1]=NA		
	}
	
    #calculate AUC if rawpred values exists
    #else AUC is meaningless (equivalent to CCR, assuming cutoff=0.5)
    if(!is.null(predtable$rawpred)){
        predobjall=try(prediction(predtable$rawpred,predtable$yobs))
		if(class(predobjall) != "try-error"){
			temp4=try(performance(predobjall,"auc"))
			if(class(temp4) != "try-error"){
				auc[nfolds+1]=performance(predobjall,"auc")@y.values[[1]]
			}else{
				auc[nfolds+1]=NA
				print("Error in line 159: AUC will be NA")
			}
		}else{
			auc[nfolds+1]=NA
			print("Error in line 163: AUC will be NA")
		}
    }
    
	if(!is.na(outputfile)){
		sink(outputfile,append=T)
			cat("cumulative: without AD","\n")
			print(cMall)
			print(CCR[nfolds+1])
			print(auc[nfolds+1])
		sink()
    }
	
    #for AD calculations (cumulative)
	if(!is.null(predtable$OUTofAD)){
		coverage[nfolds+1]=1-sum(predtable$OUTofAD)/length(predtable$OUTofAD)
		
		predtable_withAD=subset(predtable,OUTofAD==0)
		cMall_withAD<-try(confusionMatrix(factor(predtable_withAD$ypred,levels=0:1),
						factor(predtable_withAD$yobs,levels=0:1),"1"))
		if (class(cMall_withAD) != "try-error"){				  
			sens_withAD[nfolds+1]=cMall_withAD$byClass["Sensitivity"]
			spec_withAD[nfolds+1]=cMall_withAD$byClass["Specificity"]
			CCR_withAD[nfolds+1]=(sens_withAD[nfolds+1]+spec_withAD[nfolds+1])/2
			acc_withAD[nfolds+1]=cMall_withAD$overall["Accuracy"]
		}else{
			sens_withAD[nfolds+1]=NA
			spec_withAD[nfolds+1]=NA
			CCR_withAD[nfolds+1]=NA
			acc_withAD[nfolds+1]=NA
			print("Fewer than 2 classes after excluding out of AD compounds. Performance metrics set to NA")
		}
			
        #calculate AUC if rawpred values exists
        #else AUC is meaningless (equivalent to CCR, assuming cutoff=0.5)
        if(!is.null(predtable$rawpred)){
            #in case only 1 class after too many exclusions (gives error)
            predobjall_withAD=try(prediction(predtable_withAD$rawpred,predtable_withAD$yobs))
            if (class(predobjall_withAD) != "try-error"){
                temp=try(performance(predobjall_withAD,"auc"))
				if (class(temp) != "try-error"){
					auc_withAD[nfolds+1]=performance(predobjall_withAD,"auc")@y.values[[1]]
				}else{
					auc_withAD[nfolds+1]=NA
					print("Error in line 205: AUC will be NA")
				}
            }else{
                auc_withAD[nfolds+1]=NA
                print("Error in line 209: AUC will be NA")
            }
        }
        
		if(!is.na(outputfile)){
			sink(outputfile,append=T)
				cat("cumulative: with AD","\n")
				print(cMall_withAD)
				print(CCR_withAD[nfolds+1])
				print(auc_withAD[nfolds+1])
			sink()
		}
    }	#out of if block for AD calculations

#calculate errors using bootstrapping
#requires bootstrapSD_BIN.R

	sens.err=spec.err=ccr.err=acc.err=auc.err=NA
	sensAD.err=specAD.err=ccrAD.err=accAD.err=aucAD.err=NA
	
	if(calcErr){
		sens.err=getbootstrapErr(predtable,sensfun,ntrials=ntrials)
		spec.err=getbootstrapErr(predtable,specfun,ntrials=ntrials)
		ccr.err=getbootstrapErr(predtable,ccrfun,ntrials=ntrials)
		acc.err=getbootstrapErr(predtable,accfun,ntrials=ntrials)
		auc.err=getbootstrapErr(predtable,aucfun,ntrials=ntrials)
    
		if(!is.null(predtable$OUTofAD)){
			predtable_withAD=predtable[predtable$OUTofAD==0,]
			sensAD.err=getbootstrapErr(predtable_withAD,sensfun,ntrials=ntrials)
			specAD.err=getbootstrapErr(predtable_withAD,specfun,ntrials=ntrials)
			ccrAD.err=getbootstrapErr(predtable_withAD,ccrfun,ntrials=ntrials)
			accAD.err=getbootstrapErr(predtable_withAD,accfun,ntrials=ntrials)
			aucAD.err=getbootstrapErr(predtable_withAD,aucfun,ntrials=ntrials)
		}
	}
    coverage.err=NA
    errstats=matrix(c(	sens.err,	spec.err,	ccr.err,
                        acc.err,	auc.err,	coverage.err,
                        sensAD.err,	specAD.err,	ccrAD.err,
                        accAD.err,	aucAD.err),nrow=1)
    colnames(errstats)=c(	"spec","sens","CCR","acc","auc","coverage",
                            "spec_withAD","sens_withAD","CCR_withAD",
                            "acc_withAD","auc_withAD")
    rownames(errstats)="error"
    
    
#output stats
	temp=cbind(spec,sens,CCR,acc,auc,coverage,
                            spec_withAD,sens_withAD,CCR_withAD,
                            acc_withAD,auc_withAD)
    temp2=merge(t(temp),t(errstats),by=0,sort=F)
    validationstats=t(temp2)[-1,]
    mode(validationstats)="numeric"
    rownames(validationstats)=c(foldnames,"cumulative","error")
    colnames(validationstats)=temp2[,"Row.names"]
    if(!is.na(validationstatsfile)){
        write.table(validationstats,validationstatsfile,col.names=NA,row.names=T,sep="\t",quote=F)
        cat("validation stats are saved to",validationstatsfile,"\n")
    }
	return(validationstats)
       
}
