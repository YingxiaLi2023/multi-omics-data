########################################################

# Set the working directory to the directory 'multi-omics-data' 
# of the electronic appendix (outcomment the following line
# and replace 'pathtomulti-omics-data' by the path to 'multi-omics-data'
# on your computer):

## setwd("pathtomulti-omics-data/multi-omics-data/Results")

########################################################

##### load data 
library("grid")
library("gridExtra")
library("ggplot2")
library("ggpubr")
library("patchwork")
load("./rda_files/resultsumsum.RData")

mnames <- c("rna","mirna","methy","mutation","cnv",
            
            "rna_mirna","rna_cnv","rna_methy","rna_mutation","mirna_methy", 
            "mirna_mutation", "miran_cnv", "methy_mutation", "methy_cnv", "mutation_cnv",
            
            "miran_methy_cnv", "miran_mutation_cnv", "rna_mirna_methy","rna_mirna_mutation","rna_mirna_cnv",
            "rna_methy_mutation","rna_methy_cnv", "rna_mutation_cnv",
            "mirna_methy_mutation", "methy_mutation_cnv",
            
            "rna_mirna_methy_mutation","rna_mirna_methy_cnv", "rna_mirna_mutation_cnv", 
            "rna_methy_mutation_cnv", "mirna_methy_mutation_cnv",
            
            "rna_mirna_methy_mutation_cnv")

###### prepare the data of heatmap 
a <- rep(mnames,5)

b <- c(rep("rna", 31),rep("mirna", 31),rep("methy", 31),rep("mutation", 31),rep("cnv", 31))

c <- c("1","0","0","0","0",
       "1","1","1","1","0",
       "0","0","0","0","0",
       "0","0","1","1","1",
       "1","1","1","0","0",
       "1","1","1","1","0","1",
       
       "0","1","0","0","0",
       "1","0","0","0","1",
       "1","1","0","0","0",
       "1","1","1","1","1",
       "0","0","0","1","0",
       "1","1","1","0","1","1",
       
       "0","0","1","0","0",
       "0","0","1","0","1",
       "0","0","1","1","0",
       "1","0","1","0","0",
       "1","1","0","1","1",
       "1","1","0","1","1","1",
       
       "0","0","0","1","0",
       "0","0","0","1","0",
       "1","0","1","0","1",
       "0","1","0","1","0",
       "1","0","1","1","1",
       "1","0","1","1","1","1",
       
       "0","0","0","0","1",
       "0","1","0","0","0",
       "0","1","0","1","1",
       "1","1","0","0","1",
       "0","1","1","0","1",
       "0","1","1","1","1","1"
       
)

d <- data.frame(as.factor(a),as.factor(b),as.character(c))
colnames(d) <- c("a","b","c")

###### non parameter methods #####
#### the cindex of block forest
resultscindex_bf <- cbind(resultsumsum[,c(1,2)],resultsumsum$cindex_bf)
colnames(resultscindex_bf) <- c("comb", "dat","cindex_bf")

resultswide <- reshape(resultscindex_bf , idvar = c("dat"), timevar = "comb", direction = "wide")
resultswide1 <- resultswide[,-1]

# Calculate ranks of the methods:
resultranks <- t(apply(resultswide1, 1, function(x) rank(-x)))

colnames(resultranks) <- gsub("cindex_bf.", "", colnames(resultranks))
names(resultswide1) <- gsub("cindex_bf.", "", names(resultswide1))

resultrankstemp <- data.frame(resultranks)


resultranks2 <- reshape(resultrankstemp, varying=colnames(resultranks), 
                        v.names="rank", 
                        timevar="combin", times=colnames(resultranks),
                        direction="long")
resultswide2 <- reshape(resultswide1, varying=names(resultswide1), 
                        v.names="val", 
                        timevar="combin", times=names(resultswide1),
                        direction="long")
resultranks2$combin <- factor(resultranks2$combin, levels=mnames)
resultswide2$combin <- factor(resultswide2$combin, levels=mnames)
set.seed(1234)

means <- aggregate(rank ~  combin, resultranks2, mean)
means <- means[order(means$rank),]

meansmetric <- aggregate(val ~  combin, resultswide2, mean)

p1 <- ggplot(data=resultranks2, aes(x=combin, y=rank)) + 
  theme_bw() + 
  ggtitle('b     bf')+
  geom_boxplot() + 
  geom_point(data=means, aes(x=combin, y=rank), shape=23, fill="blue", size=2)+#, stroke=1.5) + # Adding diamond shapes
  scale_x_discrete(limits=means$combin) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y = element_text(size = 13),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 17))+
  labs(x = " ", y = "Ranks")

p1_2 <- ggplot(data=resultswide2, aes(x=combin, y=val)) + 
  theme_bw() + 
  ggtitle('b     bf')+
  geom_boxplot() + 
  geom_line(aes(group=factor(id), color=factor(id)), size=0.3) + 
  geom_point(data=meansmetric, aes(x=combin, y=val), shape=23, fill="blue", size=2)+#, stroke=1.5) + # Adding diamond shapes
  scale_x_discrete(limits=means$combin) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y = element_text(size = 13),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 17),
        legend.position="none")+
  labs(x = " ", y = "cindex")

p2 <- ggplot(d, aes(x=a, y=b, fill=c)) + 
  geom_tile(color = "white",lwd = 1.5,linetype = 1) +
  scale_fill_manual(values = c("grey70","darkviolet"))+
  labs(x = " ", y = "") +
  coord_fixed() + 
  scale_y_discrete(limits=c("mutation", "methy","cnv","mirna","rna"),
                   labels= c("mut", "met","cnv","mirna","rna")) +
  scale_x_discrete(limits=means$combin) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.title.y = element_text(size = 13),
        legend.position = 'none') +
  annotate("rect", xmin = 28.5, xmax = 29.5, ymin = -Inf, ymax = Inf, color="red", fill=NA, linewidth=1.5)

vertdist <- -0.33
heightboxplot <- 0.35

## draw figure
pall_cindex_bf <- ggarrange(p1, NULL, p2,  nrow = 3, align="v",heights = c(heightboxplot, vertdist, 1) )
pall_raw_cindex_bf <- ggarrange(p1_2, NULL, p2,  nrow = 3, align="v",heights = c(heightboxplot, vertdist, 1) )


#### the ibrier of block forest
resultsibrier_bf <- cbind(resultsumsum[,c(1,2)],resultsumsum$ibrier_bf)
colnames(resultsibrier_bf) <- c("comb", "dat","ibrier_bf")

resultswide <- reshape(resultsibrier_bf , idvar = c("dat"), timevar = "comb", direction = "wide")
resultswide1 <- resultswide[,-1]

# Calculate ranks of the methods:
resultranks <- t(apply(resultswide1, 1, function(x) rank(x)))

colnames(resultranks) <- gsub("ibrier_bf.", "", colnames(resultranks))
names(resultswide1) <- gsub("ibrier_bf.", "", names(resultswide1))

resultrankstemp <- data.frame(resultranks)


resultranks2 <- reshape(resultrankstemp, varying=colnames(resultranks), 
                        v.names="rank", 
                        timevar="combin", times=colnames(resultranks),
                        direction="long")
resultswide2 <- reshape(resultswide1, varying=names(resultswide1), 
                        v.names="val", 
                        timevar="combin", times=names(resultswide1),
                        direction="long")
resultranks2$combin <- factor(resultranks2$combin, levels=mnames)
resultswide2$combin <- factor(resultswide2$combin, levels=mnames)
set.seed(1234)

means <- aggregate(rank ~  combin, resultranks2, mean)
means <- means[order(means$rank),]

meansmetric <- aggregate(val ~  combin, resultswide2, mean)

p1 <- ggplot(data=resultranks2, aes(x=combin, y=rank)) + 
  theme_bw() + 
  ggtitle('b     bf')+
  geom_boxplot() + 
  geom_point(data=means, aes(x=combin, y=rank), shape=23, fill="blue", size=2) +
  scale_x_discrete(limits=means$combin) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y = element_text(size = 13),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 17))+
  labs(x = " ", y = "Ranks")

p1_2 <- ggplot(data=resultswide2, aes(x=combin, y=val)) + 
  theme_bw() + 
  ggtitle('b     bf')+
  geom_boxplot() + 
  geom_line(aes(group=factor(id), color=factor(id)), size=0.3) + 
  geom_point(data=meansmetric, aes(x=combin, y=val), shape=23, fill="blue", size=2)+#, stroke=1.5) + # Adding diamond shapes
  scale_x_discrete(limits=means$combin) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y = element_text(size = 13),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 17),
        legend.position = "none")+
  labs(x = " ", y = "ibrier")
  

p2 <- ggplot(d, aes(x=a, y=b, fill=c)) + 
  geom_tile(color = "white",lwd = 1.5,linetype = 1) +
  scale_fill_manual(values = c("grey70","darkviolet"))+
  labs(x = " ", y = "") +
  coord_fixed() + 
  scale_y_discrete(limits=c("mutation", "methy","cnv","mirna","rna"),
                   labels= c("mut", "met","cnv","mirna","rna")) +
  scale_x_discrete(limits=means$combin) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 13),
        legend.position = 'none') +
  annotate("rect", xmin = 14.5, xmax = 15.5, ymin = -Inf, ymax = Inf, color="red", fill=NA, linewidth=1.5)

## draw figure
pall_ibrier_bf <- ggarrange(p1, NULL, p2,  nrow = 3, align="v",heights = c(heightboxplot, vertdist, 1) )
pall_raw_ibrier_bf <- ggarrange(p1_2, NULL, p2,  nrow = 3, align="v",heights = c(heightboxplot, vertdist, 1) )


#### the cindex of random forest 
resultscindex_rf <- cbind(resultsumsum[,c(1,2)],resultsumsum$cindex_rf)
colnames(resultscindex_rf) <- c("comb", "dat","cindex_rf")

resultswide <- reshape(resultscindex_rf , idvar = c("dat"), timevar = "comb", direction = "wide")
resultswide1 <- resultswide[,-1]

# Calculate ranks of the methods:
resultranks <- t(apply(resultswide1, 1, function(x) rank(-x)))

colnames(resultranks) <- gsub("cindex_rf.", "", colnames(resultranks))
names(resultswide1) <- gsub("cindex_rf.", "", names(resultswide1))

resultrankstemp <- data.frame(resultranks)


resultranks2 <- reshape(resultrankstemp, varying=colnames(resultranks), 
                        v.names="rank", 
                        timevar="combin", times=colnames(resultranks),
                        direction="long")
resultswide2 <- reshape(resultswide1, varying=names(resultswide1), 
                        v.names="val", 
                        timevar="combin", times=names(resultswide1),
                        direction="long")
resultranks2$combin <- factor(resultranks2$combin, levels=mnames)
resultswide2$combin <- factor(resultswide2$combin, levels=mnames)
set.seed(1234)

means <- aggregate(rank ~  combin, resultranks2, mean)
means <- means[order(means$rank),]

meansmetric <- aggregate(val ~  combin, resultswide2, mean)

p1 <- ggplot(data=resultranks2, aes(x=combin, y=rank)) + 
  theme_bw() + 
  geom_boxplot() + 
  geom_point(data=means, aes(x=combin, y=rank), shape=23, fill="blue", size=2) +
  ggtitle('a     rsf')+
  scale_x_discrete(limits=means$combin) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y = element_text(size = 13),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 17))+
  labs(x = " ", y = "Ranks")
  
  p1_2 <- ggplot(data=resultswide2, aes(x=combin, y=val)) + 
  theme_bw() + 
  ggtitle('a     rsf')+
  geom_boxplot() + 
  geom_line(aes(group=factor(id), color=factor(id)), size=0.3) + 
  geom_point(data=meansmetric, aes(x=combin, y=val), shape=23, fill="blue", size=2)+#, stroke=1.5) + # Adding diamond shapes
  scale_x_discrete(limits=means$combin) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y = element_text(size = 13),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 17),
        legend.position = "none")+
  labs(x = " ", y = "cindex")
  

p2 <- ggplot(d, aes(x=a, y=b, fill=c)) + 
  geom_tile(color = "white",lwd = 1.5,linetype = 1) +
  scale_fill_manual(values = c("grey70","darkviolet"))+
  labs(x = " ", y = "") +
  coord_fixed() + 
  scale_y_discrete(limits=c("mutation", "methy","cnv","mirna","rna"),
                   labels= c("mut", "met","cnv","mirna","rna")) +
  scale_x_discrete(limits=means$combin) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.title.y = element_text(size = 13),
        legend.position = 'none') +
  annotate("rect", xmin = 24.5, xmax = 25.5, ymin = -Inf, ymax = Inf, color="red", fill=NA, linewidth=1.5)

## draw figure
pall_cindex_rf <- ggarrange(p1, NULL, p2,  nrow = 3, align="v",heights = c(heightboxplot, vertdist, 1) )
pall_raw_cindex_rf <- ggarrange(p1_2, NULL, p2,  nrow = 3, align="v",heights = c(heightboxplot, vertdist, 1) )


#### the ibrier of random forest 
resultsibrier_rf <- cbind(resultsumsum[,c(1,2)],resultsumsum$ibrier_rf)
colnames(resultsibrier_rf) <- c("comb", "dat","ibrier_rf")

resultswide <- reshape(resultsibrier_rf , idvar = c("dat"), timevar = "comb", direction = "wide")
resultswide1 <- resultswide[,-1]

# Calculate ranks of the methods:
resultranks <- t(apply(resultswide1, 1, function(x) rank(x)))

colnames(resultranks) <- gsub("ibrier_rf.", "", colnames(resultranks))
names(resultswide1) <- gsub("ibrier_rf.", "", names(resultswide1))

resultrankstemp <- data.frame(resultranks)


resultranks2 <- reshape(resultrankstemp, varying=colnames(resultranks), 
                        v.names="rank", 
                        timevar="combin", times=colnames(resultranks),
                        direction="long")
resultswide2 <- reshape(resultswide1, varying=names(resultswide1), 
                        v.names="val", 
                        timevar="combin", times=names(resultswide1),
                        direction="long")
resultranks2$combin <- factor(resultranks2$combin, levels=mnames)
resultswide2$combin <- factor(resultswide2$combin, levels=mnames)
set.seed(1234)

means <- aggregate(rank ~  combin, resultranks2, mean)
means <- means[order(means$rank),]

meansmetric <- aggregate(val ~  combin, resultswide2, mean)

p1 <- ggplot(data=resultranks2, aes(x=combin, y=rank)) + 
  theme_bw() + 
  geom_boxplot() + 
  ggtitle('a     rsf')+
  geom_point(data=means, aes(x=combin, y=rank), shape=23, fill="blue", size=2) +
  scale_x_discrete(limits=means$combin) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y = element_text(size = 13),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 17))+
  labs(x = " ", y = "Ranks")

p1_2 <- ggplot(data=resultswide2, aes(x=combin, y=val)) + 
  theme_bw() + 
  ggtitle('a     rsf')+
  geom_boxplot() + 
  geom_line(aes(group=factor(id), color=factor(id)), size=0.3) + 
  geom_point(data=meansmetric, aes(x=combin, y=val), shape=23, fill="blue", size=2)+#, stroke=1.5) + # Adding diamond shapes
  scale_x_discrete(limits=means$combin) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y = element_text(size = 13),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 17),
        legend.position = "none")+
  labs(x = " ", y = "ibrier")
  
p2 <- ggplot(d, aes(x=a, y=b, fill=c)) + 
  geom_tile(color = "white",lwd = 1.5,linetype = 1) +
  scale_fill_manual(values = c("grey70","darkviolet"))+
  labs(x = " ", y = "") +
  coord_fixed() + 
  scale_y_discrete(limits=c("mutation", "methy","cnv","mirna","rna"),
                   labels= c("mut", "met","cnv","mirna","rna")) +
  scale_x_discrete(limits=means$combin) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.title.y = element_text(size = 13),
        legend.position = 'none') +
  annotate("rect", xmin = 12.5, xmax = 13.5, ymin = -Inf, ymax = Inf, color="red", fill=NA, linewidth=1.5)


## draw figure
pall_ibrier_rf <- ggarrange(p1, NULL, p2,  nrow = 3, align="v",heights = c(heightboxplot, vertdist, 1) )
pall_raw_ibrier_rf <- ggarrange(p1_2, NULL, p2,  nrow = 3, align="v",heights = c(heightboxplot, vertdist, 1) )


###### parameter methods #####
#### the cindex of lasso
resultscindex_lasso <- cbind(resultsumsum[,c(1,2)],resultsumsum$cindex_lasso)
colnames(resultscindex_lasso) <- c("comb", "dat","cindex_lasso")

resultswide <- reshape(resultscindex_lasso , idvar = c("dat"), timevar = "comb", direction = "wide")
resultswide1 <- resultswide[,-1]

# Calculate ranks of the methods:
resultranks <- t(apply(resultswide1, 1, function(x) rank(-x)))

colnames(resultranks) <- gsub("cindex_lasso.", "", colnames(resultranks))
names(resultswide1) <- gsub("cindex_lasso.", "", names(resultswide1))

resultrankstemp <- data.frame(resultranks)


resultranks2 <- reshape(resultrankstemp, varying=colnames(resultranks), 
                        v.names="rank", 
                        timevar="combin", times=colnames(resultranks),
                        direction="long")
resultswide2 <- reshape(resultswide1, varying=names(resultswide1), 
                        v.names="val", 
                        timevar="combin", times=names(resultswide1),
                        direction="long")
resultranks2$combin <- factor(resultranks2$combin, levels=mnames)
resultswide2$combin <- factor(resultswide2$combin, levels=mnames)
set.seed(1234)

means <- aggregate(rank ~  combin, resultranks2, mean)
means <- means[order(means$rank),]

meansmetric <- aggregate(val ~  combin, resultswide2, mean)

p1 <- ggplot(data=resultranks2, aes(x=combin, y=rank)) + 
  theme_bw() + 
  geom_boxplot() + 
  ggtitle('a     lasso')+
  geom_point(data=means, aes(x=combin, y=rank), shape=23, fill="blue", size=2) +
  scale_x_discrete(limits=means$combin) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y = element_text(size = 13),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 17))+
  labs(x = " ", y = "Ranks")

p1_2 <- ggplot(data=resultswide2, aes(x=combin, y=val)) + 
  theme_bw() + 
  ggtitle('a     lasso')+
  geom_boxplot() + 
  geom_line(aes(group=factor(id), color=factor(id)), size=0.3) + 
  geom_point(data=meansmetric, aes(x=combin, y=val), shape=23, fill="blue", size=2)+#, stroke=1.5) + # Adding diamond shapes
  scale_x_discrete(limits=means$combin) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y = element_text(size = 13),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 17),
        legend.position = "none")+
  labs(x = " ", y = "cindex")
  

p2 <- ggplot(d, aes(x=a, y=b, fill=c)) + 
  geom_tile(color = "white",lwd = 1.5,linetype = 1) +
  scale_fill_manual(values = c("grey70","darkviolet"))+
  labs(x = " ", y = "") +
  coord_fixed() + 
  scale_y_discrete(limits=c("mutation", "methy","cnv","mirna","rna"),
                   labels= c("mut", "met","cnv","mirna","rna")) +
  scale_x_discrete(limits=means$combin) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.ticks.x = element_blank(),
        legend.position = 'none') +
  annotate("rect", xmin = 20.5, xmax = 21.5, ymin = -Inf, ymax = Inf, color="red", fill=NA, linewidth=1.5)

## draw figure
pall_cindex_lasso <- ggarrange(p1, NULL, p2,  nrow = 3, align="v",heights = c(1, -0.32, 1) )
pall_raw_cindex_lasso <- ggarrange(p1_2, NULL, p2,  nrow = 3, align="v",heights = c(1, -0.32, 1) )


#### the ibrier of lasso
resultsibrier_lasso <- cbind(resultsumsum[,c(1,2)],resultsumsum$ibrier_lasso)
colnames(resultsibrier_lasso) <- c("comb", "dat","ibrier_lasso")

resultswide <- reshape(resultsibrier_lasso , idvar = c("dat"), timevar = "comb", direction = "wide")
resultswide1 <- resultswide[,-1]

# Calculate ranks of the methods:
resultranks <- t(apply(resultswide1, 1, function(x) rank(x)))

colnames(resultranks) <- gsub("ibrier_lasso.", "", colnames(resultranks))
names(resultswide1) <- gsub("ibrier_lasso.", "", names(resultswide1))

resultrankstemp <- data.frame(resultranks)


resultranks2 <- reshape(resultrankstemp, varying=colnames(resultranks), 
                        v.names="rank", 
                        timevar="combin", times=colnames(resultranks),
                        direction="long")
resultswide2 <- reshape(resultswide1, varying=names(resultswide1), 
                        v.names="val", 
                        timevar="combin", times=names(resultswide1),
                        direction="long")
resultranks2$combin <- factor(resultranks2$combin, levels=mnames)
resultswide2$combin <- factor(resultswide2$combin, levels=mnames)
set.seed(1234)

means <- aggregate(rank ~  combin, resultranks2, mean)
means <- means[order(means$rank),]

meansmetric <- aggregate(val ~  combin, resultswide2, mean)

p1 <- ggplot(data=resultranks2, aes(x=combin, y=rank)) + 
  theme_bw() + 
  geom_boxplot() + 
  ggtitle('a     lasso')+
  geom_point(data=means, aes(x=combin, y=rank), shape=23, fill="blue", size=2) +
  scale_x_discrete(limits=means$combin) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y = element_text(size = 13),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 17))+
  labs(x = " ", y = "Ranks")

p1_2 <- ggplot(data=resultswide2, aes(x=combin, y=val)) + 
  theme_bw() + 
  ggtitle('a     lasso')+
  geom_boxplot() + 
  geom_line(aes(group=factor(id), color=factor(id)), size=0.3) + 
  geom_point(data=meansmetric, aes(x=combin, y=val), shape=23, fill="blue", size=2)+#, stroke=1.5) + # Adding diamond shapes
  scale_x_discrete(limits=means$combin) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y = element_text(size = 13),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 17),
        legend.position = "none")+
  labs(x = " ", y = "ibrier")
  
p2 <- ggplot(d, aes(x=a, y=b, fill=c)) + 
  geom_tile(color = "white",lwd = 1.5,linetype = 1) +
  scale_fill_manual(values = c("grey70","darkviolet"))+
  labs(x = " ", y = "") +
  coord_fixed() + 
  scale_y_discrete(limits=c("mutation", "methy","cnv","mirna","rna"),
                   labels= c("mut", "met","cnv","mirna","rna")) +
  scale_x_discrete(limits=means$combin) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.ticks.x = element_blank(),
        legend.position = 'none') +
  annotate("rect", xmin = 8.5, xmax = 9.5, ymin = -Inf, ymax = Inf, color="red", fill=NA, linewidth=1.5)

## draw figure
pall_ibrier_lasso <- ggarrange(p1, NULL, p2,  nrow = 3, align="v",heights = c(1, -0.32, 1) )
pall_raw_ibrier_lasso <- ggarrange(p1_2, NULL, p2,  nrow = 3, align="v",heights = c(1, -0.32, 1) )



#### the cindex of prioritylasso
resultscindex_prioritylasso <- cbind(resultsumsum[,c(1,2)],resultsumsum$cindex_prioritylasso)
colnames(resultscindex_prioritylasso) <- c("comb", "dat","cindex_prioritylasso")

resultswide <- reshape(resultscindex_prioritylasso , idvar = c("dat"), timevar = "comb", direction = "wide")
resultswide1 <- resultswide[,-1]

# Calculate ranks of the methods:
resultranks <- t(apply(resultswide1, 1, function(x) rank(-x)))

colnames(resultranks) <- gsub("cindex_prioritylasso.", "", colnames(resultranks))
names(resultswide1) <- gsub("cindex_prioritylasso.", "", names(resultswide1))

resultrankstemp <- data.frame(resultranks)


resultranks2 <- reshape(resultrankstemp, varying=colnames(resultranks), 
                        v.names="rank", 
                        timevar="combin", times=colnames(resultranks),
                        direction="long")
resultswide2 <- reshape(resultswide1, varying=names(resultswide1), 
                        v.names="val", 
                        timevar="combin", times=names(resultswide1),
                        direction="long")
resultranks2$combin <- factor(resultranks2$combin, levels=mnames)
resultswide2$combin <- factor(resultswide2$combin, levels=mnames)
set.seed(1234)

means <- aggregate(rank ~  combin, resultranks2, mean)
means <- means[order(means$rank),]
means$combin

meansmetric <- aggregate(val ~  combin, resultswide2, mean)

p1 <- ggplot(data=resultranks2, aes(x=combin, y=rank)) + 
  theme_bw() + 
  geom_boxplot() +
  ggtitle('b     prioritylasso')+
  geom_point(data=means, aes(x=combin, y=rank), shape=23, fill="blue", size=2) +
  scale_x_discrete(limits=means$combin) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y = element_text(size = 13),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 17))+
  labs(x = " ", y = "Ranks")

p1_2 <- ggplot(data=resultswide2, aes(x=combin, y=val)) + 
  theme_bw() + 
  ggtitle('b     prioritylasso')+
  geom_boxplot() + 
  geom_line(aes(group=factor(id), color=factor(id)), size=0.3) + 
  geom_point(data=meansmetric, aes(x=combin, y=val), shape=23, fill="blue", size=2)+#, stroke=1.5) + # Adding diamond shapes
  scale_x_discrete(limits=means$combin) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y = element_text(size = 13),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 17),
        legend.position = "none")+
  labs(x = " ", y = "cindex")
  
p2 <- ggplot(d, aes(x=a, y=b, fill=c)) + 
  geom_tile(color = "white",lwd = 1.5,linetype = 1) +
  scale_fill_manual(values = c("grey70","darkviolet"))+
  labs(x = " ", y = "") +
  coord_fixed() + 
  scale_y_discrete(limits=c("mutation", "methy","cnv","mirna","rna"),
                   labels= c("mut", "met","cnv","mirna","rna")) +
  scale_x_discrete(limits=means$combin) +
  theme(panel.background = element_blank(), # Make panel background transparent
        axis.title.y = element_text(size = 13),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 17),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 13),
        legend.position = 'none') +
  annotate("rect", xmin = 30.5, xmax = 31.5, ymin = 0.4, ymax = Inf, color="red", fill=NA, linewidth=1.5)

# Adding an arrow and text
p2 <- p2 + 
  annotate("segment", x = 0.5, xend = 31.5, y = -0, yend = -0, 
           arrow = arrow(type = "closed", ends = "first", length = unit(0.2, "inches"), angle = 15), 
           color = "black", size = 1) +
  annotate("text", x = 20, y = -1, label = "Better mean performance", 
           hjust = 1, vjust = -0.5, color = "black", size = 5)

## draw figure
pall_cindex_prioritylasso <- ggarrange(p1, NULL, p2,  nrow = 3, align="v",heights = c(1, -0.22, 1) )#c(1, -0.32, 1) )
pall_raw_cindex_prioritylasso <- ggarrange(p1_2, NULL, p2,  nrow = 3, align="v",heights = c(1, -0.22, 1) )#c(1, -0.32, 1) )



#### the ibrier of prioritylasso
resultsibrier_prioritylasso <- cbind(resultsumsum[,c(1,2)],resultsumsum$ibrier_prioritylasso)
colnames(resultsibrier_prioritylasso) <- c("comb", "dat","ibrier_prioritylasso")

resultswide <- reshape(resultsibrier_prioritylasso , idvar = c("dat"), timevar = "comb", direction = "wide")
resultswide1 <- resultswide[,-1]

# Calculate ranks of the methods:
resultranks <- t(apply(resultswide1, 1, function(x) rank(x)))

colnames(resultranks) <- gsub("ibrier_prioritylasso.", "", colnames(resultranks))
names(resultswide1) <- gsub("ibrier_prioritylasso.", "", names(resultswide1))

resultrankstemp <- data.frame(resultranks)


resultranks2 <- reshape(resultrankstemp, varying=colnames(resultranks), 
                        v.names="rank", 
                        timevar="combin", times=colnames(resultranks),
                        direction="long")
resultswide2 <- reshape(resultswide1, varying=names(resultswide1), 
                        v.names="val", 
                        timevar="combin", times=names(resultswide1),
                        direction="long")
resultranks2$combin <- factor(resultranks2$combin, levels=mnames)
resultswide2$combin <- factor(resultswide2$combin, levels=mnames)
set.seed(1234)

means <- aggregate(rank ~  combin, resultranks2, mean)
means <- means[order(means$rank),]

meansmetric <- aggregate(val ~  combin, resultswide2, mean)

p1 <- ggplot(data=resultranks2, aes(x=combin, y=rank)) + 
  theme_bw() + 
  geom_boxplot() + 
  ggtitle('b     prioritylasso')+
  geom_point(data=means, aes(x=combin, y=rank), shape=23, fill="blue", size=2) +
  scale_x_discrete(limits=means$combin) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y = element_text(size = 13),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 17))+
  labs(x = " ", y = "Ranks")

p1_2 <- ggplot(data=resultswide2, aes(x=combin, y=val)) + 
  theme_bw() + 
  ggtitle('b     prioritylasso')+
  geom_boxplot() + 
  geom_line(aes(group=factor(id), color=factor(id)), size=0.3) + 
  geom_point(data=meansmetric, aes(x=combin, y=val), shape=23, fill="blue", size=2)+#, stroke=1.5) + # Adding diamond shapes
  scale_x_discrete(limits=means$combin) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y = element_text(size = 13),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 17),
        legend.position = "none")+
  labs(x = " ", y = "ibrier")
  
  
p2 <- ggplot(d, aes(x=a, y=b, fill=c)) + 
  geom_tile(color = "white",lwd = 1.5,linetype = 1) +
  scale_fill_manual(values = c("grey70","darkviolet"))+
  labs(x = " ", y = "") +
  coord_fixed() + 
  scale_y_discrete(limits=c("mutation", "methy","cnv","mirna","rna"),
                   labels= c("mut", "met","cnv","mirna","rna")) +
  scale_x_discrete(limits=means$combin) +
  theme(panel.background = element_blank(), # Make panel background transparent
        axis.title.y = element_text(size = 13),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 17),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 13),
        legend.position = 'none') +
  annotate("rect", xmin = 30.5, xmax = 31.5, ymin = 0.4, ymax = Inf, color="red", fill=NA, linewidth=1.5)

# Adding an arrow and text
p2 <- p2 + 
  annotate("segment", x = 0.5, xend = 31.5, y = -0, yend = -0, 
           arrow = arrow(type = "closed", ends = "first", length = unit(0.2, "inches"), angle = 15), 
           color = "black", size = 1) +
  annotate("text", x = 20, y = -1, label = "Better mean performance", 
           hjust = 1, vjust = -0.5, color = "black", size = 5)

## draw figure
pall_ibrier_prioritylasso <- ggarrange(p1, NULL, p2,  nrow = 3, align="v",heights = c(1, -0.22, 1) )#c(1, -0.32, 1) )
pall_raw_ibrier_prioritylasso <- ggarrange(p1_2, NULL, p2,  nrow = 3, align="v",heights = c(1, -0.22, 1) )#c(1, -0.32, 1) )




#### the cindex of ipflasso
resultscindex_ipflasso <- cbind(resultsumsum[,c(1,2)],resultsumsum$cindex_ipflasso)
colnames(resultscindex_ipflasso) <- c("comb", "dat","cindex_ipflasso")

resultswide <- reshape(resultscindex_ipflasso , idvar = c("dat"), timevar = "comb", direction = "wide")
resultswide1 <- resultswide[,-1]

# Calculate ranks of the methods:
resultranks <- t(apply(resultswide1, 1, function(x) rank(-x)))

colnames(resultranks) <- gsub("cindex_ipflasso.", "", colnames(resultranks))
names(resultswide1) <- gsub("cindex_ipflasso.", "", names(resultswide1))

resultrankstemp <- data.frame(resultranks)


resultranks2 <- reshape(resultrankstemp, varying=colnames(resultranks), 
                        v.names="rank", 
                        timevar="combin", times=colnames(resultranks),
                        direction="long")
resultswide2 <- reshape(resultswide1, varying=names(resultswide1), 
                        v.names="val", 
                        timevar="combin", times=names(resultswide1),
                        direction="long")
resultranks2$combin <- factor(resultranks2$combin, levels=mnames)
resultswide2$combin <- factor(resultswide2$combin, levels=mnames)
set.seed(1234)

means <- aggregate(rank ~  combin, resultranks2, mean)
means <- means[order(means$rank),]

meansmetric <- aggregate(val ~  combin, resultswide2, mean)

p1 <- ggplot(data=resultranks2, aes(x=combin, y=rank)) + 
  theme_bw() + 
  geom_boxplot() + 
  ggtitle('c     ipflasso')+
  geom_point(data=means, aes(x=combin, y=rank), shape=23, fill="blue", size=2) +
  scale_x_discrete(limits=means$combin) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y = element_text(size = 13),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 17))+
  labs(x = " ", y = "Ranks")

p1_2 <- ggplot(data=resultswide2, aes(x=combin, y=val)) + 
  theme_bw() + 
  ggtitle('c     ipflasso')+
  geom_boxplot() + 
  geom_line(aes(group=factor(id), color=factor(id)), size=0.3) + 
  geom_point(data=meansmetric, aes(x=combin, y=val), shape=23, fill="blue", size=2)+#, stroke=1.5) + # Adding diamond shapes
  scale_x_discrete(limits=means$combin) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y = element_text(size = 13),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 17),
        legend.position = "none")+
  labs(x = " ", y = "cindex")
  
  
p2 <- ggplot(d, aes(x=a, y=b, fill=c)) + 
  geom_tile(color = "white",lwd = 1.5,linetype = 1) +
  scale_fill_manual(values = c("grey70","darkviolet"))+
  labs(x = " ", y = "") +
  coord_fixed() + 
  scale_y_discrete(limits=c("mutation", "methy","cnv","mirna","rna"),
                   labels= c("mut", "met","cnv","mirna","rna")) +
  scale_x_discrete(limits=means$combin) +
  theme(panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 13),
        axis.title.y = element_text(size = 13),
        legend.position = 'none') +
  annotate("rect", xmin = 28.5, xmax = 29.5, ymin = 0.4, ymax = Inf, color="red", fill=NA, linewidth=1.5)

# Adding an arrow and text
p2 <- p2 + 
  annotate("segment", x = 0.5, xend = 31.5, y = -0, yend = -0, 
           arrow = arrow(type = "closed", ends = "first", length = unit(0.2, "inches"), angle = 15), 
           color = "black", size = 1) +
  annotate("text", x = 20, y = -1, label = "Better mean performance", 
           hjust = 1, vjust = -0.5, color = "black", size = 5)

## draw figure
pall_cindex_ipflasso <- ggarrange(p1, NULL, p2,  nrow = 3, align="v",heights = c(0.35, -0.13, 1) )# c(1, -0.05, 1) ) # c(1, -0.25, 1) )
pall_raw_cindex_ipflasso <- ggarrange(p1_2, NULL, p2,  nrow = 3, align="v",heights = c(0.35, -0.13, 1) )# c(1, -0.05, 1) ) # c(1, -0.25, 1) )



#### the ibrier of ipflasso
resultsibrier_ipflasso <- cbind(resultsumsum[,c(1,2)],resultsumsum$ibrier_ipflasso)
colnames(resultsibrier_ipflasso) <- c("comb", "dat","ibrier_ipflasso")

resultswide <- reshape(resultsibrier_ipflasso , idvar = c("dat"), timevar = "comb", direction = "wide")
resultswide1 <- resultswide[,-1]

# Calculate ranks of the methods:
resultranks <- t(apply(resultswide1, 1, function(x) rank(x)))

colnames(resultranks) <- gsub("ibrier_ipflasso.", "", colnames(resultranks))
names(resultswide1) <- gsub("ibrier_ipflasso.", "", names(resultswide1))

resultrankstemp <- data.frame(resultranks)


resultranks2 <- reshape(resultrankstemp, varying=colnames(resultranks), 
                        v.names="rank", 
                        timevar="combin", times=colnames(resultranks),
                        direction="long")
resultswide2 <- reshape(resultswide1, varying=names(resultswide1), 
                        v.names="val", 
                        timevar="combin", times=names(resultswide1),
                        direction="long")
resultranks2$combin <- factor(resultranks2$combin, levels=mnames)
resultswide2$combin <- factor(resultswide2$combin, levels=mnames)
set.seed(1234)

means <- aggregate(rank ~  combin, resultranks2, mean)
means <- means[order(means$rank),]
means$combin

meansmetric <- aggregate(val ~  combin, resultswide2, mean)

p1 <- ggplot(data=resultranks2, aes(x=combin, y=rank)) + 
  theme_bw() + 
  geom_boxplot() + 
  ggtitle('c     ipflasso')+
  geom_point(data=means, aes(x=combin, y=rank), shape=23, fill="blue", size=2) +
  scale_x_discrete(limits=means$combin) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y = element_text(size = 13),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 17))+
  labs(x = " ", y = "Ranks")

p1_2 <- ggplot(data=resultswide2, aes(x=combin, y=val)) + 
  theme_bw() + 
  ggtitle('c     ipflasso')+
  geom_boxplot() + 
  geom_line(aes(group=factor(id), color=factor(id)), size=0.3) + 
  geom_point(data=meansmetric, aes(x=combin, y=val), shape=23, fill="blue", size=2)+#, stroke=1.5) + # Adding diamond shapes
  scale_x_discrete(limits=means$combin) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y = element_text(size = 13),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 17),
        legend.position = "none")+
  labs(x = " ", y = "ibrier")
  
p2 <- ggplot(d, aes(x=a, y=b, fill=c)) + 
  geom_tile(color = "white",lwd = 1.5,linetype = 1) +
  scale_fill_manual(values = c("grey70","darkviolet"))+
  labs(x = " ", y = "") +
  coord_fixed() + 
  scale_y_discrete(limits=c("mutation", "methy","cnv","mirna","rna"),
                   labels= c("mut", "met","cnv","mirna","rna")) +
  scale_x_discrete(limits=means$combin) +
  theme(panel.background = element_blank(), # Make panel background transparent
    axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 13),
        legend.position = 'none') +
  annotate("rect", xmin = 29.5, xmax = 30.5, ymin = 0.4, ymax = Inf, color="red", fill=NA, linewidth=1.5)

# Adding an arrow and text
p2 <- p2 + 
  annotate("segment", x = 0.5, xend = 31.5, y = -0, yend = -0, 
           arrow = arrow(type = "closed", ends = "first", length = unit(0.2, "inches"), angle = 15), 
           color = "black", size = 1) +
  annotate("text", x = 20, y = -1, label = "Better mean performance", 
           hjust = 1, vjust = -0.5, color = "black", size = 5)

## draw figure
pall_ibrier_ipflasso <- ggarrange(p1, NULL, p2,  nrow = 3, align="v",heights = c(0.35, -0.13, 1) ) # c(0.35, -0.33, 1) )
pall_raw_ibrier_ipflasso <- ggarrange(p1_2, NULL, p2,  nrow = 3, align="v",heights = c(0.35, -0.13, 1) ) # c(0.35, -0.33, 1) )



#### combine figures ####
p3_cindex <- ggarrange(pall_cindex_rf,  NULL, pall_cindex_bf,  NULL, pall_cindex_ipflasso,
                             nrow = 5, align="v", heights = c(1, -0.35, 1,-0.35,1)) #c(1, -0.14, 1,-0.14,1))
ggsave(file="./figures/figureS6_raw.png", 
       p3_cindex, width=11*0.865, height=14*0.865)

p3_cindex <- ggarrange(pall_raw_cindex_rf,  NULL, pall_raw_cindex_bf,  NULL, pall_raw_cindex_ipflasso,
                             nrow = 5, align="v", heights = c(1, -0.35, 1,-0.35,1)) #c(1, -0.14, 1,-0.14,1))
ggsave(file="./figures/figureS3_raw.png", 
       p3_cindex, width=11*0.865, height=14*0.865)


p2_cindex <- ggarrange(pall_cindex_lasso,  NULL, pall_cindex_prioritylasso,
                       nrow = 3, align="v", heights = c(1, -0.2, 1))
ggsave(file="./figures/figureS7_raw.png", 
       p2_cindex, width=8, height=12)

p2_cindex <- ggarrange(pall_raw_cindex_lasso,  NULL, pall_raw_cindex_prioritylasso,
                       nrow = 3, align="v", heights = c(1, -0.2, 1))
ggsave(file="./figures/figureS4_raw.png", 
       p2_cindex, width=8, height=12)


p3_ibrier <- ggarrange(pall_ibrier_rf,  NULL, pall_ibrier_bf,  NULL, pall_ibrier_ipflasso,
                       nrow = 5, align="v", heights = c(1, -0.35, 1,-0.35,1))
ggsave(file="./figures/figure1_raw.png", 
       p3_ibrier, width=11*0.865, height=14*0.865)

p3_ibrier <- ggarrange(pall_raw_ibrier_rf,  NULL, pall_raw_ibrier_bf,  NULL, pall_raw_ibrier_ipflasso,
                       nrow = 5, align="v", heights = c(1, -0.35, 1,-0.35,1))
ggsave(file="./figures/figureS1_raw.png", 
       p3_ibrier, width=11*0.865, height=14*0.865)


p2_ibrier <- ggarrange(pall_ibrier_lasso,  NULL, pall_ibrier_prioritylasso,
                       nrow = 3, align="v", heights = c(1, -0.2, 1))
ggsave(file="./figures/figureS5_raw.png", 
       p2_ibrier, width=8, height=12)

p2_ibrier <- ggarrange(pall_raw_ibrier_lasso,  NULL, pall_raw_ibrier_prioritylasso,
                       nrow = 3, align="v", heights = c(1, -0.2, 1))
ggsave(file="./figures/figureS2_raw.png", 
       p2_ibrier, width=8, height=12)