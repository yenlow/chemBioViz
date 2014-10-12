#Reads .x, .a, .xa, .txt (Dragon) data files
#22-Aug-11 Yen Low
##################################################

readXAfile <-
function(xafile,coefficient=TRUE,activity=FALSE){
#to read XA file
if(tail(unlist(strsplit(xafile,"\\.")),1)=="xa"){

	ds<-read.table(xafile,skip=1, fill=T, row.names=NULL, header=F,comment.char="",quote="",check.names=F)
	dsheader=ds[1,-((ncol(ds)-2):ncol(ds))]
    if(coefficient==TRUE){
    dsdata=ds[-c(1,nrow(ds)-1,nrow(ds)),-(1:3)]
    dsact=ds[-c(1,nrow(ds)-1,nrow(ds)),3]
    dsnames=ds[-c(1,nrow(ds)-1,nrow(ds)),2]
    dscoeff=ds[c(nrow(ds)-1,nrow(ds)),]
    coeff=matrix(dscoeff[!is.na(dscoeff)],nrow=2,byrow=F)
    mode(coeff)<-"numeric"  
    }else{
        dsdata=ds[-1,-(1:3)]
        dsact=ds[-1,3]
        dsnames=ds[-1,2]
        coeff=NULL
    }
	colnames(dsdata)<-as.vector(unlist(dsheader))
	rownames(dsdata)<-as.vector(unlist(dsnames))
	mat<-as.matrix(dsdata)
	mode(mat)<-"numeric"

	#check if dsact is binary or continuous
	if(length(table(dsact))==2) list(x=mat,y=as.factor(as.character(dsact)),cpds=as.vector(dsnames),coeff=coeff)
	else list(x=mat,y=as.numeric(as.character(dsact)),cpds=as.vector(dsnames),coeff=coeff)
}	

#to read X file
else if(tail(unlist(strsplit(xafile,"\\.")),1)=="x"){
	ds<-read.table(xafile,skip=1, fill=T, row.names=NULL, header=F,comment.char="",quote="",check.names=F)
	dsheader<-ds[1,-c(length(ds)-1,length(ds))]
	if(coefficient==TRUE){
		dsdata<-ds[-c(1,nrow(ds)-1,nrow(ds)),-(1:2)]
		dsnames=ds[-c(1,nrow(ds)-1,nrow(ds)),2]
	    dscoeff=ds[c(nrow(ds)-1,nrow(ds)),]
	    coeff=matrix(dscoeff[!is.na(dscoeff)],nrow=2,byrow=F)
	    mode(coeff)<-"numeric"
    }else{
        dsdata=ds[-1,-(1:2)]
        dsnames=ds[-1,2]
        coeff=NULL
    }
	colnames(dsdata)<-as.vector(unlist(dsheader))
	rownames(dsdata)<-as.vector(unlist(dsnames))
	mat<-as.matrix(dsdata)
	mode(mat)<-"numeric"

list(x=mat,cpds=as.vector(dsnames),coeff=coeff)
}

#to read A file
else if(tail(unlist(strsplit(xafile,"\\.")),1)=="a"){
	ds<-read.table(xafile,row.names=NULL,header=F,comment.char="",quote="", fill=T)
	if(length(table(ds$V2))==2) list(y=as.factor(as.character(ds$V2)),cpds=as.vector(ds$V1))
	else list(y=as.numeric(as.character(ds$V2)),cpds=as.vector(ds$V1))
}

#to read drg file
else if(tail(unlist(strsplit(xafile,"\\.")),1) %in% c("txt","drg")){
	ds<-read.table(xafile,skip=2,row.names=NULL,header=T,comment.char="",quote="", fill=T,check.names=F)
	dsnames=ds[,2]
	dsdata=as.matrix(ds[,-(1:2)])
	rownames(dsdata)=dsnames
	mode(dsdata)="numeric"
	#Dragon file does not contain activity column (default)
	if(!activity){
		dsact=NULL
		list(x=dsdata,cpds=as.vector(dsnames))
	}
	else { #if Dragon file contains activity column (3rd column)
		dsact=ds[,3]
		dsdata=dsdata[,-1]
		if(length(table(dsact))==2) list(x=dsdata,y=as.factor(as.character(dsact)),cpds=as.vector(dsnames))
		else list(x=dsdata,y=as.numeric(as.character(dsact)),cpds=as.vector(dsnames))
	}
}


else{print("Invalid file.Acceptable file formats include .a, .x, .xa, .drg, or .txt(Dragon)")}

}
