library(pheatmap)
all<-read.table("FIG_S9.txt",header=TRUE,sep="\t",row.names=1)
pheatmap(all, cellwidth=6,cluster_col = FALSE, cluster_row= FALSE, cellheight=6,fontsize_row=6,fontsize_col=6, color = colorRampPalette(c("black", "gold2"))(50), filename="fig_s9.pdf")

