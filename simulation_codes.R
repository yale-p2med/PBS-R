my.msgidb.list.load<-function(x){
temp<-readLines(x)
temp.names<-unname(sapply(temp,my.element.extract,splitchar="\t",index=1))
temp<-sapply(temp,my.element.remove,splitchar="\t",index=1)
temp<-sapply(temp,my.element.remove,splitchar="\t",index=1)
temp.list<-sapply(temp,strsplit,split="\t")
names(temp.list)<-temp.names
return(temp.list)
}


#=============================================================================================
# PART I: simulate the gene expression data when genes are independent
#=============================================================================================


simul.independent<-function(pathway.perc=0.2,gene.perc=0.2,delta=0.5,cluster.num=3,samplesize.percluster=40,msigdb.filepath,pathway.name.prefix="KEGG_",output.dir,simul.times=100,array.anno.filepath){

#####################################
# simul.independent generates the simulated data matrices for given parameter settings:
#
# Arguments:
#
#	pathway.perc			the percentage of pathways that are perturbed or differentially 
#							expressed across the given groups. default 0.2
#	gene.perc				the percentage of genes inside the chosen pathways to be perturbed or
#							differentially expressed. default 0.2
#	delta					the difference in the average gene expression levels between the given 
#							groups. Group 1 has mean of delta. Group 2 has mean of 2*delta, and so 
#							on. default 0.5
#	cluster.num				the total number of given groups. default 3
#	samplesize.percluster	the number of samples simulated for each given group. default 40
#	msigdb.filepath			the file that defines the pathways downloaded from MsigDB.
#	pathway.name.prefix		the type of pathways in MsigDB used to calculate the distance score. 
#							default "KEGG_"
#	output.dir				the folder to put the results in.
#	simul.times				the total number of simulations to do. default 100
#	array.anno.filepath		the csv file that has the annotation for all the probe sets or genes to be 
#							simulated. This file needs to have the same format as the example 
#							annotation file HuGene-1_0-st-v1.na32.hg19.transcript.csv. The first column
#							should have the IDs present in the msigdb.filepath.
#
# Values:
#	
#	simul.dependent will generate simul.times different set of files with each set presenting one 
#	simulated data. Each simulation include three files: simulated_data*_datamatrix.gct,
#	simulated_data*_perturbedpathways.txt and simulated_data*_perturbedgenes.txt. All files will be 
#	put under output.dir. 
# 	
#	simulated_data*_datamatrix.txt contains the simulated gene expression 
#	data matrix. Rows of the matrix represent all the probe sets in the array.anno.filepath. 
#	Columns represent the samples from cluster.num groups and each group has samplesize.percluster 
#	samples. 
#
#	simulated_data*_perturbedpathways.txt contains the names of the randomly chosen pathways to be
#	perturbed in each simulated data set.
#
#	simulated_data_perturbedgenes.txt contains the probe set IDs of the randomly chosen genes to be
#	perturbed in each simulated data set.
#
#####################################

# generate the mean of the gaussian distribution for the given groups.
mu.perturbed<-seq(from=1,to=cluster.num,by=1)*delta

# if the output.dir doesn't exist, create a new folder
if(file.exists(output.dir)==F){
dir.create(output.dir,recursive = TRUE)
}

#source(myfunction.filepath)
temp.data<-read.csv(array.anno.filepath,comment.char="#")
#temp.data<-read.table(array.anno.filepath,sep="\t",header=T,row.names=1)
gene.num<-nrow(temp.data)
ps.names<-as.matrix(temp.data)[,1]
sample.num<-cluster.num*samplesize.percluster

# load in the pathway info from given category
msigdb.pathway.list<-my.msgidb.list.load(msigdb.filepath)

if(pathway.name.prefix==""){
pathway.list<-msigdb.pathway.list
}else{
pathway.list<-msigdb.pathway.list[substr(names(msigdb.pathway.list),1,nchar(pathway.name.prefix))%in%pathway.name.prefix]
}

# simulate the gene expression data matrix for each cluster
for(simulate.index in 1:simul.times){

data.matrix<-numeric()
output.filepath<-file.path(output.dir,paste("simulated_data",simulate.index,"_datamatrix.gct",sep=""))

# randomly choose the pathways to perturb
pathway2perturb<-sample(1:length(pathway.list),floor(length(pathway.list)*pathway.perc),replace=F)
perturbed.pathway.list<-names(pathway.list)[pathway2perturb]
gene2perturb<-character()
for(i in pathway2perturb){
temp<-sample(pathway.list[[i]],floor(length(pathway.list[[i]])*gene.perc),replace=F)
gene2perturb<-c(gene2perturb,temp)
gene2perturb<-unique(gene2perturb)
}
gene2perturb<-gene2perturb[gene2perturb%in%ps.names]

for(cluster.index in 1:cluster.num){

# generate the gene expression data for all the genes
simulate.data.matrix<-matrix(rnorm(n=gene.num*samplesize.percluster,mean=0,sd=1),nrow=gene.num)
rownames(simulate.data.matrix)<-ps.names
colnames(simulate.data.matrix)<-paste("C",cluster.index,"_",1:ncol(simulate.data.matrix),sep="")

#assign the perturbed values to the chosen genes
if(length(gene2perturb)>0){
simulate.data.matrix[gene2perturb,]<-matrix(rnorm(n=length(gene2perturb)*ncol(simulate.data.matrix),mean=mu.perturbed[cluster.index],sd=1),nrow=length(gene2perturb))
}

# save the simulated data matrix
data.matrix<-cbind(data.matrix,simulate.data.matrix)

}

# output the simulated gene expression levels
write.table(data.matrix,file=output.filepath,append=F,sep="\t",row.names=T,col.names=T,quote=F)

# output the perturbed pathways and genes for each cluster
pathway.output.filepath<-file.path(output.dir,paste("simulated_data",simulate.index,"_perturbedpathways.txt",sep=""))
gene.output.filepath<-file.path(output.dir,paste("simulated_data",simulate.index,"_perturbedgenes.txt",sep=""))

cmd.out<-paste(perturbed.pathway.list,collapse="\n",sep="")
cat(cmd.out,file=pathway.output.filepath,append=F)
cat("\n",file=pathway.output.filepath,append=T)

cmd.out<-paste(gene2perturb,collapse="\n",sep="")
cat(cmd.out,file=gene.output.filepath,append=F)
cat("\n",file=gene.output.filepath,append=T)
}

}


#=============================================================================================
# PART II: simulate the gene expression data when genes are not independent
#=============================================================================================

simul.dependent<-function(pathway.perc=0.2,gene.perc=0.2,delta=0.5,r=0.8,cluster.num=3,samplesize.percluster=40,msigdb.filepath,pathway.name.prefix="KEGG_",output.dir,simul.times=100,array.anno.filepath){

#####################################
# simul.independent generates the simulated data matrices for given parameter settings:
#
# Arguments:
#
#	pathway.perc			the percentage of pathways that are perturbed or differentially 
#							expressed across the given groups. default 0.2
#	gene.perc				the percentage of genes inside the chosen pathways to be perturbed or
#							differentially expressed. default 0.2
#	delta					the difference in the average gene expression levels between the given 
#							groups. Group 1 has mean of delta. Group 2 has mean of 2*delta, and so 
#							on. default 0.5
#	r						the correlation coefficient between genes from the same pathways. default 0.8
#	cluster.num				the total number of given groups. default 3
#	samplesize.percluster	the number of samples simulated for each given group. default 40
#	msigdb.filepath			the file that defines the pathways downloaded from MsigDB.
#	pathway.name.prefix		the type of pathways in MsigDB used to calculate the distance score. 
#							default "KEGG_"
#	output.dir				the folder to put the results in.
#	simul.times				the total number of simulations to do. default 100
#	array.anno.filepath		the csv file that has the annotation for all the probe sets or genes to be 
#							simulated. This file needs to have the same format as the example 
#							annotation file HuGene-1_0-st-v1.na32.hg19.transcript.csv. The first column
#							should have the IDs present in the msigdb.filepath.
#
# Values:
#	
#	simul.dependent will generate simul.times different set of files with each set presenting one 
#	simulated data. Each simulation include three files: simulated_data*_datamatrix.gct,
#	simulated_data*_perturbedpathways.txt and simulated_data*_perturbedgenes.txt. All files will be 
#	put under output.dir. 
# 	
#	simulated_data*_datamatrix.txt contains the simulated gene expression 
#	data matrix. Rows of the matrix represent all the probe sets in the array.anno.filepath. 
#	Columns represent the samples from cluster.num groups and each group has samplesize.percluster 
#	samples. 
#
#	simulated_data*_perturbedpathways.txt contains the names of the randomly chosen pathways to be
#	perturbed in each simulated data set.
#
#	simulated_data_perturbedgenes.txt contains the probe set IDs of the randomly chosen genes to be
#	perturbed in each simulated data set.
#
#####################################

mu.perturbed<-c(1,2,3)*delta

if(file.exists(output.dir)==F){
dir.create(output.dir,recursive = TRUE)
}

# use the setting provided in the array.anno.filepath
temp.data<-read.csv(array.anno.filepath,comment.char="#")
#temp.data<-read.table(array.anno.filepath,sep="\t",header=T,row.names=1)
gene.num<-nrow(temp.data)
ps.names<-as.matrix(temp.data)[,1]
sample.num<-cluster.num*samplesize.percluster

# load in the pathway info from KEGG
msigdb.pathway.list<-my.msgidb.list.load(msigdb.filepath)
if(pathway.name.prefix==""){
pathway.list<-msigdb.pathway.list
}else{
pathway.list<-msigdb.pathway.list[substr(names(msigdb.pathway.list),1,nchar(pathway.name.prefix))%in%pathway.name.prefix]
}

library(MASS)
# simulate the gene expression data matrix for each cluster
for(simulate.index in 1:simul.times){

data.matrix<-numeric()
output.filepath<-file.path(output.dir,paste("simulated_data",simulate.index,"_datamatrix.gct",sep=""))

# randomly choose the pathways to perturb
pathway2perturb<-sample(1:length(pathway.list),floor(length(pathway.list)*pathway.perc),replace=F)
perturbed.pathway.list<-names(pathway.list)[pathway2perturb]
gene2perturb<-character()
for(i in pathway2perturb){
temp<-sample(pathway.list[[i]],floor(length(pathway.list[[i]])*gene.perc),replace=F)
gene2perturb<-c(gene2perturb,temp)
gene2perturb<-unique(gene2perturb)
}
gene2perturb<-gene2perturb[gene2perturb%in%ps.names]

for(cluster.index in 1:cluster.num){

# generate the gene expression data for all the genes
simulate.data.matrix<-matrix(rnorm(n=gene.num*samplesize.percluster,mean=0,sd=1),nrow=gene.num)
rownames(simulate.data.matrix)<-ps.names
colnames(simulate.data.matrix)<-paste("C",cluster.index,"_",1:ncol(simulate.data.matrix),sep="")

#assign the perturbed values to the chosen genes
if(length(gene2perturb)>0){
#simulate.data.matrix[gene2perturb,]<-matrix(rnorm(n=length(gene2perturb)*ncol(simulate.data.matrix),mean=mu.perturbed[cluster.index],sd=1),nrow=length(gene2perturb))
varcov<-diag(1-r,nrow=length(gene2perturb),ncol=length(gene2perturb))+r
gene_multi<-mvrnorm(n = samplesize.percluster, mu=rep(mu.perturbed[cluster.index],length(gene2perturb)), Sigma=varcov)
simulate.data.matrix[gene2perturb,]<-t(gene_multi)
}


# save the simulated data matrix
data.matrix<-cbind(data.matrix,simulate.data.matrix)

}

# output the simulated gene expression levels
write.table(data.matrix,file=output.filepath,append=F,sep="\t",row.names=T,col.names=T,quote=F)

# output the perturbed pathways and genes for each cluster
pathway.output.filepath<-file.path(output.dir,paste("simulated_data",simulate.index,"_perturbedpathways.txt",sep=""))
gene.output.filepath<-file.path(output.dir,paste("simulated_data",simulate.index,"_perturbedgenes.txt",sep=""))

cmd.out<-paste(perturbed.pathway.list,collapse="\n",sep="")
cat(cmd.out,file=pathway.output.filepath,append=F)
cat("\n",file=pathway.output.filepath,append=T)

cmd.out<-paste(gene2perturb,collapse="\n",sep="")
cat(cmd.out,file=gene.output.filepath,append=F)
cat("\n",file=gene.output.filepath,append=T)
}

}

