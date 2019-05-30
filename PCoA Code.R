#PCOA
#install.packages("vegan")
#install.packages("ggplot2")
library(vegan)
otu<- read.table("otu.txt", row.names=1,header=T)#read the out table
otu.hel<-decostand(otu, "hellinger")#data transformation
otu.hel=otu
group<- read.table("groups.txt", row.names=1,header=T)#group table, each sample is assigned to a given group (e.g Dry or Wet)
head(group)
otu<-t(otu)#transpose the out table
bray_curtis<-vegdist(otu)
attach(group)
pcoa = cmdscale(bray_curtis, k=3, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
points = as.data.frame(pcoa$points) # get coordinate string, format to dataframme
colnames(points) = c("x", "y", "z") 
eig = pcoa$eig
points = cbind(points, group[match(rownames(points), rownames(group)), ])
# Make the figure with ggplot2
library(ggplot2)
p = ggplot(points, aes(x=x, y=y, color=group$Season)) +
  geom_point(alpha=.7, size=4) + 
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       title="bray_curtis PCoA of CRT")
z=p+theme_bw()+scale_y_continuous(breaks=c(-0.2, -0.1, 0, 0.1, 0.2))#scaling y axis
z+scale_x_continuous(breaks=c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3))#scaling x axis
