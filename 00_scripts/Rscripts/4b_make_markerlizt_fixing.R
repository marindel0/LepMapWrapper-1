setwd("/home/zjons/Github_play/LepMapWrapper/LM-Run-221025@2256/")


argv <- commandArgs(T)
CONTIGFILE <- "04_filtering/contig_H19-109-gp-mina2-maxa4-mm09-mac10-NoCont_Miss-0.1.txt"
POSFILE <- argv[2]
MARKERLIST <- argv[3]
i <- argv[4]




contig<-read.table(CONTIGFILE, stringsAsFactors=F )
pos<-read.table(POSFILE, header=T, stringsAsFactors=F)
marker_list<-cbind(c(1:(dim(contig)[1]-6)), contig[7:dim(contig)[1],1], pos[7:dim(contig)[1],1], paste(contig[7:dim(contig)[1],1], pos[7:dim(contig)[1],1], sep="_"))
head(marker_list)
colnames(marker_list)<-c("marker_id", "contig", "pos", "contig_pos")
print(paste("after filtering, there are", dim(marker_list)[1], "markers in the family", i))
write.table( marker_list, file=MARKERLIST, row.names=F, quote=F)
