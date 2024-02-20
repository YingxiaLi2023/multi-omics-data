########################################################

# Set the working directory to the directory 'multi-omics-data' 
# of the electronic appendix (outcomment the following line
# and replace 'pathtomulti-omics-data' by the path to 'multi-omics-data'
# on your computer):

## setwd("pathtomulti-omics-data/multi-omics-data/Results")
setwd("Z:/Projects/SideProjects/Yingxia/Multi-Omics-Importance/Neu_wichtig/Final_Code_on_GitHub/multi-omics-data/Results")

########################################################


##### load data 
library("grid")
library("gridExtra")
library("ggplot2")
library("ggpubr")
library("patchwork")
load("./rda_files/resultsumsum.RData")

mnames <- c("rna","mirna","cnv","methy","mutation",
            
            "rna_mirna","rna_cnv","rna_methy","rna_mutation","miran_cnv", 
            
            "mirna_methy", "mirna_mutation", "methy_cnv", "mutation_cnv", "methy_mutation",
            
            "rna_mirna_cnv", "rna_mirna_methy","rna_mirna_mutation",
            "rna_methy_cnv", "rna_mutation_cnv",
            
            "rna_methy_mutation","miran_methy_cnv", "miran_mutation_cnv", 
            "mirna_methy_mutation", "methy_mutation_cnv",
            
            "rna_mirna_methy_cnv","rna_mirna_mutation_cnv", "rna_mirna_methy_mutation", 
            "rna_methy_mutation_cnv", "mirna_methy_mutation_cnv",
            
            "rna_mirna_methy_mutation_cnv")

###### prepare the data of heatmap 
a <- rep(mnames,5)

b <- c(rep("rna", 31),rep("mirna", 31),rep("cnv", 31),rep("methy", 31),rep("mutation", 31))

c <- c("1","0","0","0","0",
       "1","1","1","1","0",
       "0","0","0","0","0",
       "1","1","1","1","1",
       "1","0","0","0","0",
       "1","1","1","1","0","1",
       
       "0","1","0","0","0",
       "1","0","0","0","1",
       "1","1","0","0","0",
       "1","1","1","0","0",
       "0","1","1","1","0",
       "1","1","1","0","1","1",
       
       "0","0","1","0","0",
       "0","1","0","0","1",
       "0","0","1","1","0",
       "1","0","0","1","1",
       "0","1","1","0","1",
       "1","1","0","1","1","1",
       
       "0","0","0","1","0",
       "0","0","1","0","0",
       "1","0","1","0","1",
       "0","1","0","1","0",
       "1","1","0","1","1",
       "1","0","1","1","1","1",
       
       "0","0","0","0","1",
       "0","0","0","1","0",
       "0","1","0","1","1",
       "0","0","1","0","1",
       "1","0","1","1","1",
       "0","1","1","1","1","1"
       
)

d <- data.frame(as.factor(a),as.factor(b),as.character(c))
colnames(d) <- c("a","b","c")

p2 <- ggplot(d, aes(x=a, y=b, fill=c)) + 
  geom_tile(color = "white",lwd = 1.5,linetype = 1) +
  scale_fill_manual(values = c("grey70","darkviolet"))+
  labs(x = " ", y = "") +
  coord_fixed() + 
  scale_y_discrete(limits=c("mutation", "methy","cnv","mirna","rna"),
                   labels= c("mut", "met","cnv","mirna","rna")) +
  scale_x_discrete(limits=mnames) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(colour = "black", size = 12),
        legend.position = 'none')

###### cindex #####
# bf
resultscindex_bf <- cbind(resultsumsum[,c(1,2)],resultsumsum$cindex_bf)
colnames(resultscindex_bf) <- c("comb", "dat","cindex_bf")

cindex_bf <- ggplot(data=resultscindex_bf, aes(x=comb, y=cindex_bf)) + 
  theme_bw() + 
  ggtitle('b')+
  geom_boxplot() + 
  scale_x_discrete(limits=mnames)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
       axis.ticks.x=element_blank())+
  labs(x = " ", y = "Cross-validated cindex ")

# rf
resultscindex_rf <- cbind(resultsumsum[,c(1,2)],resultsumsum$cindex_rf)
colnames(resultscindex_rf) <- c("comb", "dat","cindex_rf")

cindex_rf <- ggplot(data=resultscindex_rf, aes(x=comb, y=cindex_rf)) + 
  theme_bw() + 
  ggtitle('a')+
  geom_boxplot() + 
  scale_x_discrete(limits=mnames)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  labs(x = " ", y = "Cross-validated cindex ")

# lasso
resultscindex_lasso <- cbind(resultsumsum[,c(1,2)],resultsumsum$cindex_lasso)
colnames(resultscindex_lasso) <- c("comb", "dat","cindex_lasso")

cindex_lasso <- ggplot(data=resultscindex_lasso, aes(x=comb, y=cindex_lasso)) + 
  theme_bw() + 
  ggtitle('c')+
  geom_boxplot() + 
  scale_x_discrete(limits=mnames)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  labs(x = " ", y = "Cross-validated cindex ")

# ipflasso
resultscindex_ipflasso <- cbind(resultsumsum[,c(1,2)],resultsumsum$cindex_ipflasso)
colnames(resultscindex_ipflasso) <- c("comb", "dat","cindex_ipflasso")

cindex_ipflasso <- ggplot(data=resultscindex_ipflasso, aes(x=comb, y=cindex_ipflasso)) + 
  theme_bw() + 
  ggtitle('d')+
  geom_boxplot() + 
  scale_x_discrete(limits=mnames)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  labs(x = " ", y = "Cross-validated cindex ")

# prioritylasso
resultscindex_prioritylasso <- cbind(resultsumsum[,c(1,2)],resultsumsum$cindex_prioritylasso)
colnames(resultscindex_prioritylasso) <- c("comb", "dat","cindex_prioritylasso")

cindex_prioritylasso <- ggplot(data=resultscindex_prioritylasso, aes(x=comb, y=cindex_prioritylasso)) + 
  theme_bw() + 
  ggtitle('e')+
  geom_boxplot() + 
  scale_x_discrete(limits=mnames)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  labs(x = " ", y = "Cross-validated cindex ")

##heatmap
p2 <- ggplot(d, aes(x=a, y=b, fill=c)) + 
  geom_tile(color = "white",lwd = 1.5,linetype = 1) +
  scale_fill_manual(values = c("grey70","darkviolet"))+
  labs(x = " ", y = "") +
  coord_fixed() + 
  scale_y_discrete(limits=c("mutation", "methy","cnv","mirna","rna"),
                   labels= c("mut", "met","cnv","mirna","rna")) +
  scale_x_discrete(limits=mnames) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(colour = "black", size = 12),
        legend.position = 'none')



###### ibrier #####
# bf
resultsibrier_bf <- cbind(resultsumsum[,c(1,2)],resultsumsum$ibrier_bf)
colnames(resultsibrier_bf) <- c("comb", "dat","ibrier_bf")

ibrier_bf <- ggplot(data=resultsibrier_bf, aes(x=comb, y=ibrier_bf)) + 
  theme_bw() + 
  ggtitle('b')+
  geom_boxplot() + 
  scale_x_discrete(limits=mnames)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  labs(x = " ", y = "Cross-validated ibrier ")

# rf
resultsibrier_rf <- cbind(resultsumsum[,c(1,2)],resultsumsum$ibrier_rf)
colnames(resultsibrier_rf) <- c("comb", "dat","ibrier_rf")

ibrier_rf <- ggplot(data=resultsibrier_rf, aes(x=comb, y=ibrier_rf)) + 
  theme_bw() + 
  ggtitle('a')+
  geom_boxplot() + 
  scale_x_discrete(limits=mnames)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  labs(x = " ", y = "Cross-validated ibrier ")

# lasso
resultsibrier_lasso <- cbind(resultsumsum[,c(1,2)],resultsumsum$ibrier_lasso)
colnames(resultsibrier_lasso) <- c("comb", "dat","ibrier_lasso")

ibrier_lasso <- ggplot(data=resultsibrier_lasso, aes(x=comb, y=ibrier_lasso)) + 
  theme_bw() + 
  ggtitle('c')+
  geom_boxplot() + 
  scale_x_discrete(limits=mnames)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  labs(x = " ", y = "Cross-validated ibrier ")

# ipflasso
resultsibrier_ipflasso <- cbind(resultsumsum[,c(1,2)],resultsumsum$ibrier_ipflasso)
colnames(resultsibrier_ipflasso) <- c("comb", "dat","ibrier_ipflasso")

ibrier_ipflasso <- ggplot(data=resultsibrier_ipflasso, aes(x=comb, y=ibrier_ipflasso)) + 
  theme_bw() + 
  ggtitle('d')+
  geom_boxplot() + 
  scale_x_discrete(limits=mnames)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  labs(x = " ", y = "Cross-validated ibrier ")

# prioritylasso
resultsibrier_prioritylasso <- cbind(resultsumsum[,c(1,2)],resultsumsum$ibrier_prioritylasso)
colnames(resultsibrier_prioritylasso) <- c("comb", "dat","ibrier_prioritylasso")

ibrier_prioritylasso <- ggplot(data=resultsibrier_prioritylasso, aes(x=comb, y=ibrier_prioritylasso)) + 
  theme_bw() + 
  ggtitle('e')+
  geom_boxplot() + 
  scale_x_discrete(limits=mnames)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  labs(x = " ", y = "Cross-validated ibrier ")

#### combine figures ####
p5_cindex_value <- ggarrange(cindex_rf,  NULL, cindex_bf,  NULL, cindex_lasso,NULL, 
                             cindex_ipflasso,NULL, cindex_prioritylasso,NULL, p2,
                       nrow = 11, align="v", heights = c(1, 0, 1,0,1,0, 1,0, 1,-0.2, 1))
ggsave(file="./figures/figureS2.png", 
       p5_cindex_value, width=8, height=14)

p5_ibrier_value <- ggarrange(ibrier_rf,  NULL, ibrier_bf,  NULL, ibrier_lasso,NULL, 
                             ibrier_ipflasso,NULL, ibrier_prioritylasso,NULL, p2,
                             nrow = 11, align="v", heights = c(1, 0, 1,0,1,0, 1,0, 1,-0.2, 1))
ggsave(file="./figures/figureS1.png", 
       p5_ibrier_value, width=8, height=14)


