# Calculates nearest neighbors and generates starplots
# (uses full data set; doesn't require modeling)
#
# NNmat(spacematORdistmat,names,act,disttype="jaccard") generates NNmat object
# Input: 	1. descriptor matrix or distance matrix (must be symmetric)
#           2. cpdnames
#           3. activity matrix
#			4. distance type="jaccard","euclidean",etc (jaccard as default)
# Output:	1. NNID (first col is itself)
#           2. NNdist (first col is itself, i.e. zeros)
#           3. NN (first col is itself)
#           4. NNact (first col is itself)
#
# kNNobj(k,chemNN,geneNN,cpdnames,NNlist=NULL,mispredID=(1:nrow(chemNN[[1]]))) generates kNNobj for starplot
# Input: 	1. k
#           2. chemNN, geneNN or NNlist
#           3. cpdnames
#			4. ID of cpds of interest
# Output:	1. kNNdist (first col is itself, i.e. zeros)
#           2. kNN (first col is itself)
#           3. kNNact (first col is itself)
#
# starplotgrid(cpdvec,grid1,grid2,kNNobj,label) generates starplots
# of k nearest chemical neighbors and k nearest toxicogenomic neighbors
# Input: 	1. cpdvec (vector of desired compound IDs)
#			2. grid dimensions (grid 1 rows, grid 2 columns) 
#           3. kNN object (formed by kNNobj function)
#			4. label (logical. Should neighbors be labeled?)
#			5. byrow (TRUE as default. Fills grid by row  if TRUE)
# 
# extobj2kNNobj(extobj,yobs) convert extobj model object to kNNobj for starplots
# Input: 	1. extobj (model object/list containing cpdnames, optimal k values,
#						prodmat, OUTofAD per fold)
#			2. yobs (vector with true values for comparison. Must be named with cpds)
#
# 25-Sep-12 Yen Low
###############################################################################

if (!is.loaded("plotrix")) library(plotrix) #required for radial.plot (to display starplot)
if (!is.loaded("vegan")) library(vegan) #required for vegdist (to calculate jaccard distances)
#other distances also available (euclidean,manhattan, etc)


#get neighbours in a matrix
NNmat<-function(spacematORdistmat,names,act,disttype="jaccard"){
    #if input is a space matrix (distance matrix will be symmetric)
    #then calculate distance matrix from space matrix
    if(is.null(spacematORdistmat)){
		return(NULL)
	}else {
		if(!isSymmetric(spacematORdistmat)){
        distmat=as.matrix(vegdist(spacematORdistmat,method=disttype,upp=T,diag=T))
    	}else{
        distmat=spacematORdistmat
    	}
    NNID<-t(apply(distmat,1,order))  #sort by distances, NN in first col (itself)
    NNdist<-t(apply(distmat,1,sort))  #sort by distances, NN in first col (zeros)
    NN<-matrix(names[NNID],length(names),length(names))
    NNact<-matrix(act[NNID],length(act),length(act))
    list(NNID=NNID,NNdist=NNdist,NN=NN,NNact=NNact)
	}
}

#reverse rank of geneNN for star plots to fill from top
kNNobj<-function(k,chemNN,geneNN,NNlist=NULL,ID=NULL){
	if(is.null(ID)) ID=1:length(chemNN$NN[,1])
    if(is.null(NNlist)){ #if inputs are chemNN with or without geneNN
        kNNdist=cbind(chemNN$NNdist[ID,2:(k+1)],geneNN$NNdist[ID,(k+1):2])
        kNN=cbind(chemNN$NN[ID,2:(k+1)],geneNN$NN[ID,(k+1):2])
        kNNact=cbind(chemNN$NNact[ID,1:(k+1)],geneNN$NNact[ID,(k+1):2])
    } else { #if input is a NNlist
        kNNdist=NNlist[[1]]$NNdist[ID,2:(k+1)]
        kNN=NNlist[[1]]$NN[ID,2:(k+1)]
        kNNact=NNlist[[1]]$NNact[ID,1:(k+1)]
        for (j in 2:length(NNlist)){
            kNNdist=cbind(kNNdist,NNlist[[j]]$NNdist[ID,2:(k+1)])
            kNN=cbind(kNN,NNlist[[j]]$NN[ID,2:(k+1)])
            kNNact=cbind(kNNact,NNlist[[j]]$NNact[ID,2:(k+1)])
        }
    }
    rownames(kNNdist)=chemNN$NN[,1]
    rownames(kNN)=chemNN$NN[,1]
    rownames(kNNact)=chemNN$NN[,1]
    list(kNNdist=kNNdist,kNN=kNN,kNNact=kNNact)
}

starplotgrid<-function(cpdvec,grid1,grid2,kNNobj,label=TRUE,byrow=TRUE){
	#set text label size using cex.axis
    if(byrow) par(mfrow=c(grid1,grid2),cex.axis=0.8) else par(mfcol=c(grid1,grid2),cex.axis=0.8)
	maxlim=max(kNNobj$kNNdist[cpdvec,],na.rm=T)
	#maxlim=0.5
	for(i in cpdvec){
        #maxlim=max(kNNobj$kNNdist[i,])
        minlim=0.01
        maxk=5
		startpos=pi/2*(1-1/maxk)
		labpos=seq(0,2*pi,length=2*maxk+1)+0.1
        #removed the following options from first radial.plot statement
        #lwd=line_wdth = 2+1.3*((maxlim+minlim)/(kNNobj$kNNdist[i,]+minlim))^2.5
		#plot compound of interest
		#plot lines and labels (proportional to distance)
		if(label) NNlabel=substr(tolower(kNNobj$kNN[i,]),1,11) else NNlabel=""
        radial.plot(kNNobj$kNNdist[i,],radial.lim=c(0,maxlim),start=startpos,
                labels=NNlabel,lwd=5,label.pos=labpos,
                line.col=as.numeric(kNNobj$kNNact[i,-1])+1,clockwise=T,
				show.radial.grid=F,show.grid=F,radlab=F,label.prop=1.15,
				main=toupper(rownames(kNNobj$kNNdist)[i]),boxed.radial=F)
		#plot vertices/dots
		radial.plot(kNNobj$kNNdist[i,],radial.lim=c(0,maxlim),start=startpos,
                rp.type="s",add=T,point.symbols=16,cex=2,clockwise=T,
                point.col=as.numeric(kNNobj$kNNact[i,-1])+1,
                show.radial.grid=F,show.grid=F,radlab=F)
		#plot central ratio
		lines(c(0,0),maxlim*c(-1,1),lty="dashed",lwd=2,col="dark gray")
		points(0,0,pch=16,cex=5,col=as.numeric(kNNobj$kNNact[i,1])+1)
#        ntox=sum(as.numeric(kNNobj$kNNact[i,-1]))
#        ratio.label=paste(2*k-ntox,"/",ntox,sep="")
#        text(0,0,labels=ratio.label,col="white",font=2,cex=1)		

		#plot edges
#        radial.plot(kNNobj$kNNdist[i,],radial.lim=c(0,maxlim),start=startpos,
#                rp.type="p",add=T,clockwise=T,
#                show.radial.grid=F,show.grid=F,radlab=F)
    }
}

#create kNNobj from extobj (kNN not selected from external fold, i.e. LGO-kNN)
extobj2kNNobj<-function(extobj,yobs){
	kNNobj=optkval=yobs_wanted=prodmat=list()
	for(i in 1:length(extobj)){
		rm(kNNact) #clear kNN to avoid ambiguous reuse
		#get optimal k values
		optkval[[i]]=extobj[[i]][[1]]$stat[,c("k","k2")]
		#create matrix with 10 columns with empty columns filled with NA
		prodmat[[i]]=matrix(NA,ncol=10,nrow=nrow(extobj[[i]][[1]]$mod$prodmat))
		rownames(prodmat[[i]])=rownames(extobj[[i]][[1]]$mod$prodmat)
		prodmat[[i]][,1:optkval[[i]][1]]=extobj[[i]][[1]]$mod$prodmat[,1:optkval[[i]][1]]
		prodmat[[i]][,10:(10+1-optkval[[i]][2])]=extobj[[i]][[1]]$mod$prodmat[,(optkval[[i]][1]+1):ncol(extobj[[i]][[1]]$mod$prodmat)]
		#exclude outofAD
		prodmat[[i]]=prodmat[[i]][extobj[[i]][[1]]$OUTofAD==0,]
			
		#create objects in kNNobj
		kNNact=sign(prodmat[[i]])
		kNNact[sign(prodmat[[i]])==-1]=0
		yobs_wanted[[i]]=yobs[match(rownames(kNNact),names(yobs))]
		kNNact=cbind(yobs_wanted[[i]],kNNact)
		
		kNNobj[[i]]=list(kNNdist=(1-abs(prodmat[[i]])),kNN=NULL,kNNact=kNNact)
	}
	return(kNNobj)
}
