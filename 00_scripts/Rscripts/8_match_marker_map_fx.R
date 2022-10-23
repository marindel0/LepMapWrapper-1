#this script is called by the match_marker_and_plot function 
#it uses the marker list from 04_filtering to assign a contig & pos to the marker id given by lepmap, 
#it makes the map more easily explorable, and correspondance with the genome used to align the markers
#last step also outputs a plot with the female position which allow looking at the repartition of markers along each LG
#Rscript 01_scripts/Rscripts/match_marker_map.R "0.1" "cov6" "aa" "10" "2" "7"

suppressPackageStartupMessages(library(dplyr))

argv <- commandArgs(T)
MARKERLIST <- argv[1]               #MARKERLIST_FILE="04_filtering/data_f_"$CORENAME"_miss-"$MISS".markerlist"
NB_CHR <- as.numeric(argv[2])
REFINE_STEPS <- as.numeric(argv[3])
ORDERFILES <- argv[4]
OUTFILEBASENAME <- argv[5]
OUTDIRBASE <- argv[6]

Nmax<-200 #nb maximum of marker above which we consider it is a plateau and remove it

cat("Markerlist: ", MARKERLIST,"\n")
cat("NB_CHR: ", NB_CHR,"\n")
cat("REFINE_STEPS: ", REFINE_STEPS,"\n")
cat("Orderfiles: ", ORDERFILES,"\n")
cat("Orderfilebasename:", OUTFILEBASENAME,"\n")
cat("Output to:", OUTDIRBASE, "\n\n")


marker_list<-read.table(MARKERLIST, stringsAsFactors=F, header=T)
contig_map<-matrix(ncol=7)
contig_map_filtered<-matrix(ncol=7)
colnames(contig_map)=colnames(contig_map_filtered)=c("CHR", "marker_id", "male_pos", "female_pos", "contig", "pos", "contig_pos")

for (j in 1 : NB_CHR)
{
  cat(paste0("File:",ORDERFILES,j,".",(REFINE_STEPS + 1),".txt --> "))
  order<-read.table(paste0(ORDERFILES,j,".",(REFINE_STEPS + 1),".txt"), stringsAsFactors=F) [,1:3]
  head(order)
  order_LG<-cbind(rep(paste0("LG", j), dim(order)[1]), order)
  colnames(order_LG)<-c("CHR", "marker_id", "male_pos", "female_pos")
  head ( order_LG)

  order_LG_marker<-left_join(order_LG, marker_list)
  head(order_LG_marker)
  contig_map<-rbind(contig_map, order_LG_marker)
  #write.table(order_LG_marker, paste0("07_order_LG/contig_order_",OUTFILEBASENAME".LG",j,".",(REFINE_STEPS + 1),".txt"), row.names=F, quote=F, sep="\t")

  plateau<-which(table(factor(order_LG_marker$female_pos))>=Nmax)
  order_LG_marker_filtered<-order_LG_marker
  if (length(plateau)>=1)
	{
    for (k in 1:length(names(plateau)))
  	 {
		    order_LG_marker_filtered<-order_LG_marker_filtered[-which(order_LG_marker_filtered$female_pos==names(plateau)[k]),]
	   }
  }
  contig_map_filtered<-rbind(contig_map_filtered, order_LG_marker_filtered)
}

write.table(contig_map[2: dim(contig_map)[1],], paste0(OUTDIRBASE,"01_maps/contig_order_",OUTFILEBASENAME,".LGall.refined",(REFINE_STEPS + 1),".txt") , row.names=F,quote=F, sep="\t")
write.table(contig_map_filtered[2: dim(contig_map_filtered)[1],], paste0(OUTDIRBASE,"01_maps/contig_order_",OUTFILEBASENAME,".LGall.refined",(REFINE_STEPS + 1),".f.txt") , row.names=F,quote=F, sep="\t")

jpeg(paste0(OUTDIRBASE,"02_plots/contig_order_",OUTFILEBASENAME,".LGall.refined",(REFINE_STEPS + 1),".jpeg"))
plot(contig_map[2:dim(contig_map)[1],4], ylab="female_pos")
dev.off()
	 
jpeg(paste0(OUTDIRBASE,"02_plots/contig_order_",OUTFILEBASENAME,".LGall.refined",(REFINE_STEPS + 1),".f.jpeg"))
plot(contig_map_filtered[2:dim(contig_map_filtered)[1],4], ylab="female_pos")
dev.off()