###Combining annotations of samples
GenomeRegions<-read.table("GenomeRegions_segment200.txt",sep = ' ')
file <- list.files(pattern = "\\.txt") 
data <- read.table(file[1],header=F)
rownames(data)<-as.character(paste(data[,1],data[,2],data[,3],sep="_"))
data<-data[c(which(rownames(data) %in% GenomeRegions[,1])),]
data<-data[,4]
n = length(file)
for (f in 2:n){
  a<-read.table(file[f],header=F)
  rownames(a)<-as.character(paste(a[,1],a[,2],a[,3],sep="_"))
  a<-a[c(which(rownames(a) %in% GenomeRegions[,1])),]
  data<-cbind(data,a[,4])
}
Combined_annotation_segment200<-as.data.frame(data)
rownames(Combined_annotation_segment200)<-GenomeRegions[,1]
colnames(Combined_annotation_segment200)<-gsub(".txt","",file)
write.table(Combined_annotation_segment200, file ="Combined_annotation_segment200.txt",sep = ' ', quote = FALSE, row.names = TRUE,col.names = TRUE)

###Running tissue signatures
Tissues_65<-read.table("Tissues_65.txt",header=F,quote='',sep=" ")
setwd("./Tissue_specificRegions_65tissues")
TissueSpecific<-list.files(pattern = "\\.txt")
State<-c("E10","E11","E12","E13","E14","E15","E16","E17","E18","E1","E2","E3","E4","E5","E6","E7","E8","E9")
SampleContribution_18states<-vector("list", 18)
for (i in 1:18){
  CurrentState<-State[i]
  TissueSpecificFile<-read.table(TissueSpecific[i],sep=" ")
  a<-Combined_annotation_segment200[which(rownames(Combined_annotation_segment200)%in%rownames(TissueSpecificFile)),]
  Combined<-cbind(TissueSpecificFile,a[,1])  
  rownames(Combined)<-rownames(TissueSpecificFile)
  Combined[Combined==CurrentState] <- 1
  Combined[Combined!=1]<- 0
  rowN<-rownames(Combined)
  Combined <- as.data.frame(lapply(Combined,function(x) as.numeric(as.character(x)))) 
  rownames(Combined)<- rowN
  Combined$summary<- rowSums(Combined)
  rownames(Combined)<-rownames(TissueSpecificFile)
  Unique<- subset(Combined, Combined$summary == "2")
  Unique<-Unique[,-c((ncol(Unique)-1),ncol(Unique))]
  TissueContribution<-as.data.frame(colSums(Unique))  
  rownames(TissueContribution)<-Tissues_65[,1]  
  colnames(TissueContribution)<-colnames(Combined_annotation_segment200)[1]
  for (c in 2:ncol(Combined_annotation_segment200)){
    Combined<-cbind(TissueSpecificFile,a[,c])  
    rownames(Combined)<-rownames(TissueSpecificFile)
    Combined[Combined==CurrentState] <- 1
    Combined[Combined!=1]<- 0
    rowN<-rownames(Combined)
    Combined <- as.data.frame(lapply(Combined,function(x) as.numeric(as.character(x)))) 
    rownames(Combined)<- rowN
    Combined$summary<- rowSums(Combined)
    rownames(Combined)<-rownames(TissueSpecificFile)
    Unique<- subset(Combined, Combined$summary == "2")
    Unique<-Unique[,-c((ncol(Unique)-1),ncol(Unique))]
    TissueContribution_new<-as.data.frame(colSums(Unique))  
    rownames(TissueContribution_new)<-Tissues_65[,1]  
    colnames(TissueContribution_new)<-colnames(Combined_annotation_segment200)[c]
    TissueContribution<-cbind(TissueContribution,TissueContribution_new)
  }
  rownames(TissueContribution)<-paste(Tissues_65[,1],State[i],sep=".")  
  colnames(TissueContribution)<-colnames(Combined_annotation_segment200)
  SampleContribution_18states[[i]]<-TissueContribution
}
TissueContribution_E10<-SampleContribution_18states[[1]]
for (i in 2:18){
  TissueContribution_new<-SampleContribution_18states[[i]]
  TissueContribution_allStages<-rbind(TissueContribution,TissueContribution_new)
}
write.table(TissueContribution_allStages,file="Tissue_signatures_18states.txt")