source("my_functions.R")

pbs.cal<-function(data.filepath,msigdb.filepath,pathway.name.prefix="KEGG_",output.dir,result.name.prefix){
	
#####################################
# pbs.cal calculates the pathway based distance scores between samples from gene expression data:
#
# Arguments:
#
#	data.filepath			A tab deliminated file that has the gene expression data file (gct file). The data
#							is in a matrix form. Rows are genes and columns are samples. The first
#							column represents the row names and the first row has the column names.
#	msigdb.filepath			A text file that defines the pathways downloaded from MsigDB.
#	pathway.name.prefix		the type of pathways in MsigDB used to calculate the distance score. 
#							default "KEGG_"
#	output.dir				the folder to save the results in.
#	result.name.prefix		the prefix to use as the initial string of the result files.
#
# Values:
#	
#	pbs.cal will generate the following two files under output.dir:
#	
#	[result.name.prefix]_clustering_classmatrix.txt		A tab delimited text file that contains a matrix representing the clustering results by all the pathways.
#														Rows are samples and columns are pathways. Both row names and column names are present.
#	[result.name.prefix]_distmatrix.txt					A tab delimited text file that contains a matrix representing the pathway based distance score between samples.
#														Both rows and columns are samples and the matrix is a symmetrix matrix.
#
#####################################


# load in the gene lists for all the pathways
msigdb.pathway.list<-my.msgidb.list.load(msigdb.filepath)
if(pathway.name.prefix==""){
pathway.list<-msigdb.pathway.list
}else{
pathway.list<-msigdb.pathway.list[substr(names(msigdb.pathway.list),1,nchar(pathway.name.prefix))%in%pathway.name.prefix]
}

library(mclust) #version:4.4
library(gplots) #version:2.17.0
library(RColorBrewer) #version:1.1.2

if(file.exists(output.dir)==F){
dir.create(output.dir,recursive=TRUE)
}


# load in the gene expression data
data.matrix<-as.matrix(read.table(data.filepath,header=T,row.names=1,sep="\t",check.names=F))

# centralize the data matrix (optional)
data.matrix<-my.normalize(data.matrix)

# calculate the distance between any two samples using all genes, pathway genes, and pathway based metrics
#cat("\tEuclidean distance for all genes......\n")
#euclid.dist.matrix.all<-as.matrix(dist(t(data.matrix)))
#output.filepath<-file.path(output.dir.1,paste("distmatrix_p",pathway.perc,"_rho",gene.perc,"_dataset",i,".txt",sep=""))
#write.table(euclid.dist.matrix.all,file=output.filepath,append=F,sep="\t",row.names=T,col.names=T,quote=F)

# calculate the distance between any two samples using the genes covered by pathways
#cat("\tEuclidean distance for pathway genes......\n")
#gene.list<-unique(unlist(pathway.list))
#gene.list<-gene.list[gene.list%in%rownames(data.matrix)]
#euclid.dist.matrix.kegg<-as.matrix(dist(t(data.matrix[gene.list,])))
#output.filepath<-file.path(output.dir.2,paste("distmatrix_p",pathway.perc,"_rho",gene.perc,"_dataset",i,".txt",sep=""))
#write.table(euclid.dist.matrix.kegg,file=output.filepath,append=F,sep="\t",row.names=T,col.names=T,quote=F)


# calculate the distance using pathways with given pathway.names.prefix
mclust.list<-list()
class.matrix<-matrix(-1,nrow=ncol(data.matrix),ncol=length(pathway.list))
rownames(class.matrix)<-colnames(data.matrix)
colnames(class.matrix)<-names(pathway.list)
bic.vect<-rep(NA,length(pathway.list))
names(bic.vect)<-names(pathway.list)
cat("\tMclust for the pathways......\n")
for(j in 1:length(pathway.list)){
	cat("                                                                                                                               \r")
	cat("\tj=",j,":",names(pathway.list)[j],"\r",sep="")

	if(sum(pathway.list[[j]]%in%rownames(data.matrix))<10){
		# remove those pathways with less than 10 genes inside. There will be no output for these pathways
		mclust.list[[j]]<-NA
	}else{

		temp<-pathway.list[[j]][pathway.list[[j]]%in%rownames(data.matrix)]
		temp<-data.matrix[temp,]
		mclust.result<-Mclust(t(temp))
		mclust.list[[j]]<-mclust.result	
		bic.vect[j]<-mclust.result$bic
		# output the classification result so that the distances between samples can be recalculated.
		cmd.out<-cbind(colnames(temp),mclust.list[[j]]$classification)
		cmd.out<-apply(cmd.out,1,paste,collapse="\t")
		cmd.out<-paste(cmd.out,collapse="\n",sep="")
		#cat(cmd.out,file=output.filepath.2,append=F)
		
		# put the clustering results in the class.matrix
		class.matrix[,j]<-mclust.list[[j]]$classification
	}

}
cat("\n")
names(mclust.list)<-names(pathway.list)
mclust.list<-mclust.list[!is.na(bic.vect)]
class.matrix<-class.matrix[,!is.na(bic.vect)]
bic.vect<-bic.vect[!is.na(bic.vect)]

temp<-apply(class.matrix,2,unique)
temp<-unlist(lapply(temp,length))
if(sum(temp>1)>0){
class.matrix.clean<-as.matrix(class.matrix[,temp>1])
}else{
class.matrix.clean<-numeric()
}

####################################
# output the classmatrix
output.filepath<-file.path(output.dir,paste(result.name.prefix,"_clustering_classmatrix.txt",sep=""))
write.table(class.matrix,file=output.filepath,append=F,sep="\t",row.names=T,col.names=T,quote=F)

####################################
# calculate the distance between samples based on the clustering results by mclust of each pathway
if(length(class.matrix.clean)>0){
kegg.dist.matrix<-matrix(-1,nrow=ncol(data.matrix),ncol=ncol(data.matrix))
rownames(kegg.dist.matrix)<-colnames(data.matrix)
colnames(kegg.dist.matrix)<-colnames(data.matrix)
for(m in 1:nrow(kegg.dist.matrix)){
for(n in 1:ncol(kegg.dist.matrix)){
if(m==n){
kegg.dist.matrix[m,n]<-0
}else{
#kegg.dist.matrix[i,j]<-1-sum((class.matrix[i,]==class.matrix[j,])*bic.vect)/sum(bic.vect)
#kegg.dist.matrix[i,j]<-1-sum(class.matrix[i,]==class.matrix[j,])/ncol(class.matrix)
kegg.dist.matrix[m,n]<-1-sum(class.matrix.clean[m,]==class.matrix.clean[n,])/ncol(class.matrix.clean)
}
}
}
}else{ 
kegg.dist.matrix<-matrix(0,nrow=ncol(data.matrix),ncol=ncol(data.matrix))
rownames(kegg.dist.matrix)<-colnames(data.matrix)
colnames(kegg.dist.matrix)<-colnames(data.matrix)
}

# output the distance matrix
output.filepath<-file.path(output.dir,paste(result.name.prefix,"_distmatrix.txt",sep=""))
write.table(kegg.dist.matrix,file=output.filepath,append=F,sep="\t",row.names=T,col.names=T,quote=F)

}
