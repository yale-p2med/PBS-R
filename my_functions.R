
my.msgidb.list.load<-function(x){
temp<-readLines(x)
temp.names<-unname(sapply(temp,my.element.extract,splitchar="\t",index=1))
temp<-sapply(temp,my.element.remove,splitchar="\t",index=1)
temp<-sapply(temp,my.element.remove,splitchar="\t",index=1)
temp.list<-sapply(temp,strsplit,split="\t")
names(temp.list)<-temp.names
return(temp.list)
}


#-----------------------------------------------------------------------------------------------------
# this function normalizes the gene expression data so that each gene has a 0 mean and 1 sd
my.normalize<-function(temp1){
temp<-apply(temp1,1,mean)
temp<-matrix(rep(temp,ncol(temp1)),nrow=nrow(temp1),byrow=F)
serum.data.norm<-temp1-temp
temp<-apply(serum.data.norm,1,sd)
temp<-matrix(rep(temp,ncol(temp1)),nrow=nrow(temp1),byrow=F)
serum.data.norm<-serum.data.norm/temp
return(serum.data.norm)
}

#-----------------------------------------------------------------------------------------------------
# this function replaces the every character of orig with replace
#-----------------------------------------------------------------------------------------------------
my.char.replace<-function(x,orig="_",replace="-"){
result<-unlist(strsplit(x,split=orig))
result<-paste(result,collapse=replace)
return(unname(result))
}
#-----------------------------------------------------------------------------------------------------
my.list.element.extract<-function(x,index=1){
if(index>0){
return(x[index])
}else{
return(x[length(x)])
}
}
#-----------------------------------------------------------------------------------------------------
# this function split x by splitchar and return the parts with index=index
my.element.extract<-function(x,splitchar="\t",index=1){
temp<-unlist(strsplit(x,split=splitchar))
if(index>length(temp)){
cat("there are not enough elements!\n")
return(NA)
}
if(index<0){
result<-temp[length(temp)]
}else{
result<-temp[index]
}
return(unname(result))
}


#-----------------------------------------------------------------------------------------------------
# This function remote the element at position index in x seperated by splitchar
my.element.remove<-function(x,splitchar="",index=-1){
temp<-unlist(strsplit(x,split=splitchar))
if(index>length(temp)){
cat("there are not enough elements!\n")
return(NA)
}

if(index<0){
temp<-temp[-length(temp)]
}else{
temp<-temp[-index]
}
if(splitchar=="\\."){
return(unname(paste(temp,collapse=".")))
}else{
return(unname(paste(temp,collapse=splitchar)))
}
}



