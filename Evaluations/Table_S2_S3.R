
########################################################

# Set the working directory to the directory 'multi-omics-data' 
# of the electronic appendix (outcomment the following line
# and replace 'pathtomulti-omics-data' by the path to 'multi-omics-data'
# on your computer):

 setwd("pathtomulti-omics-data/multi-omics-data/Results")
#setwd("Z:/Projects/SideProjects/Yingxia/Multi-Omics-Importance/Neu_wichtig/Final_Code_on_GitHub/multi-omics-data/Results")
#setwd("D:/IBE/paper3/LRZ_Jul_Results")
load("./resultsumsum.RData")
########################################################


##### load data 
library("grid")
library("gridExtra")
library("ggplot2")
library("ggpubr")
library("patchwork")
library(tidyr)
library(flextable)
library(officer)
library(dplyr)
#load("./rda_files/resultsumsum.RData")

####step1: heatmap #####
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

# prepare the data of heatmap 
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

# Plot
p2 <- ggplot(d, aes(x = a, y = b, fill = c)) + 
  geom_tile(color = "white", lwd = 1.5, linetype = 1) +
  scale_fill_manual(values = c("grey70", "darkviolet")) +
  labs(x = " ", y = "") +
  coord_fixed() + 
  scale_y_discrete(limits = c("mutation", "methy", "cnv", "mirna", "rna"),
                   labels = c("mut", "met", "cnv", "mirna", "rna")) +
  scale_x_discrete(limits = mnames) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = 'none') +
  annotate("text", x = 1, y = 5, label = "rna", size = 6, hjust = 0.5) +
  annotate("text", x = 1, y = 4, label = "mirna", size = 6, hjust = 0.5) +
  annotate("text", x = 1, y = 3, label = "cnv", size = 6, hjust = 0.5) +
  annotate("text", x = 1, y = 2, label = "met", size = 6, hjust = 0.5) +
  annotate("text", x = 1, y = 1, label = "mut", size = 6, hjust = 0.5)

p2


ggsave(file="./heatmap2.png", 
       p2, width=24, height=4)




######step2: the table of cindex #####
cindex <- resultsumsum[,c(1,2,4,6,8,10,12)]

cindex_bf <- cindex[,c(1,2,3)]
cindex_bf <- cindex_bf %>% pivot_wider(names_from = comb, values_from = cindex_bf)
cindex_rf <- cindex[,c(1,2,4)]
cindex_rf <- cindex_rf %>% pivot_wider(names_from = comb, values_from = cindex_rf)
cindex_lasso <- cindex[,c(1,2,5)]
cindex_lasso <- cindex_lasso %>% pivot_wider(names_from = comb, values_from = cindex_lasso)
cindex_ipflasso <- cindex[,c(1,2,6)]
cindex_ipflasso <- cindex_ipflasso %>% pivot_wider(names_from = comb, values_from = cindex_ipflasso)
cindex_prioritylasso <- cindex[,c(1,2,7)]
cindex_prioritylasso <- cindex_prioritylasso %>% pivot_wider(names_from = comb, values_from = cindex_prioritylasso)
cindex_c <- rbind(cindex_bf,cindex_rf,cindex_lasso,cindex_ipflasso,cindex_prioritylasso)

cindex_c <- cindex_c %>%
  mutate(dat = sub("\\.RData$", "", dat))

cindex_c <- cindex_c %>% arrange(dat)

cindex_c <- cindex_c %>% select(dat, one_of(mnames))

cindex_c <- cindex_c %>%
  mutate(blank = "", .before = rna)

cindex_c <- cindex_c %>%
  mutate(across(where(is.numeric), ~ round(., 2)))

# flextable
ft <- flextable(cindex_c)

ft <- set_header_labels(ft, values = rep("", ncol(cindex_c)))

ft <- border_remove(ft)
ft <- border_outer(ft, part = "all", border = fp_border(color="black", width = 1))
ft <- border_inner_v(ft, part = "all", border = fp_border(color="black", width = 1))
ft <- border_inner_h(ft, part = "all", border = fp_border(color="black", width = 1))


ft <- set_table_properties(ft, width = 1.0, layout = "autofit") 


doc <- read_docx() %>% 
  body_add_par(value = "", style = "Normal") %>%
  body_end_section_continuous() %>%
  body_add_flextable(value = ft) %>%
  body_end_section_landscape()


print(doc, target = "./table_cindex.docx")

######step3: the table of cindex #####
ibrier <- resultsumsum[,c(1,2,3,5,7,9,11)]

ibrier_bf <- ibrier[,c(1,2,3)]
ibrier_bf <- ibrier_bf %>% pivot_wider(names_from = comb, values_from = ibrier_bf)
ibrier_rf <- ibrier[,c(1,2,4)]
ibrier_rf <- ibrier_rf %>% pivot_wider(names_from = comb, values_from = ibrier_rf)
ibrier_lasso <- ibrier[,c(1,2,5)]
ibrier_lasso <- ibrier_lasso %>% pivot_wider(names_from = comb, values_from = ibrier_lasso)
ibrier_ipflasso <- ibrier[,c(1,2,6)]
ibrier_ipflasso <- ibrier_ipflasso %>% pivot_wider(names_from = comb, values_from = ibrier_ipflasso)
ibrier_prioritylasso <- ibrier[,c(1,2,7)]
ibrier_prioritylasso <- ibrier_prioritylasso %>% pivot_wider(names_from = comb, values_from = ibrier_prioritylasso)
ibrier_c <- rbind(ibrier_bf,ibrier_rf,ibrier_lasso,ibrier_ipflasso,ibrier_prioritylasso)


ibrier_c <- ibrier_c %>%
  mutate(dat = sub("\\.RData$", "", dat))

ibrier_c <- ibrier_c %>% arrange(dat)

ibrier_c <- ibrier_c %>% select(dat, one_of(mnames))

ibrier_c <- ibrier_c %>%
  mutate(blank = "", .before = rna)

ibrier_c <- ibrier_c %>%
  mutate(across(where(is.numeric), ~ round(., 2)))

# flextable
ft <- flextable(ibrier_c)

ft <- set_header_labels(ft, values = rep("", ncol(ibrier_c)))


ft <- border_remove(ft)
ft <- border_outer(ft, part = "all", border = fp_border(color="black", width = 1))
ft <- border_inner_v(ft, part = "all", border = fp_border(color="black", width = 1))
ft <- border_inner_h(ft, part = "all", border = fp_border(color="black", width = 1))


ft <- set_table_properties(ft, width = 1.0, layout = "autofit")  


doc <- read_docx() %>% 
  body_add_par(value = "", style = "Normal") %>%
  body_end_section_continuous() %>%
  body_add_flextable(value = ft) %>%
  body_end_section_landscape()


print(doc, target = "./table_ibrier.docx")







                  
