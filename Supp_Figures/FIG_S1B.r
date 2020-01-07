library(ggplot2)
df <- read.table("FIG_1E_LEFT.txt",header=T,row.names=1,sep="\t")
g3=ggplot(df)+geom_smooth(aes(x=x,y=y, colour=colour),method='auto',alpha=.3, size=3)+geom_point(data=df, aes(x=x,y=y, colour=colour), alpha=.5, size=3)+scale_colour_brewer(type = "div", palette = "Set3")
plot(g3)+theme_mininal()
