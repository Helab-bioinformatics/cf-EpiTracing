GenomeBin<-read.table("Ordered_chromatinLocation_segment200.txt",header=F,quote='',sep=" ")
Statefile <- list.files(pattern = "\\.txt") #Combined possibility matrix of tissues and primary cells for each chromatin states
for (s in 1:18){
  Mydata<-read.table(Statefile[s],header=T,quote='',sep=" ")
  Judgement09<-as.data.frame(Mydata > 0.9)
  Judgement01<-as.data.frame(Mydata < 0.1)
  rownames(Judgement09)<-as.character(paste(GenomeBin[,1],GenomeBin[,2],GenomeBin[,3],sep="_"))
  rownames(Judgement01)<-as.character(paste(GenomeBin[,1],GenomeBin[,2],GenomeBin[,3],sep="_"))
  rownames(Mydata)<-as.character(paste(GenomeBin[,1],GenomeBin[,2],GenomeBin[,3],sep="_"))
  Judgement09[Judgement09==TRUE]<- 1
  Judgement09[Judgement09!=TRUE]<- 0
  Judgement01[Judgement01==TRUE]<- 1
  Judgement01[Judgement01!=TRUE]<- 0
  Judgement09$summary<- rowSums(Judgement09)
  rownames(Judgement09)<-rownames(State1)
  Judgement09<-subset(Judgement09, Judgement09$summary ==1)
  Judgement01$summary<- rowSums(Judgement01)
  rownames(Judgement01)<-rownames(State1)
  Judgement01<-subset(Judgement01, Judgement01$summary ==64)
  Judgement<-intersect(rownames(Judgement09),rownames(Judgement01))
  Mydata<-Mydata[c(which(rownames(Mydata)%in% Judgement)),]
  write.table(Mydata, file =paste("SpecificRegions_of_State",as.character(Statefile[s])),sep = '', quote = FALSE, row.names = TRUE)
}