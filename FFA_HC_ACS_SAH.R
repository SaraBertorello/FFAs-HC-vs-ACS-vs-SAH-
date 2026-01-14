###################### PREPARING DATA #########################
{library(Hmisc)
  library(ggpubr)
  library(ggplot2)
  library(vegan)
  library(pairwiseAdonis)
  library(mixOmics)
  library(ecodist)
  library(reshape2)
  library(readxl)
  library(FSA)
}


options(scipen=100)
#setwd('../Acidi_messicana')
FFA <- read_excel("FFA_Simone.xlsx", sheet=4) 
colnames(FFA)
unique(FFA$Disease)

FFA[! colnames(FFA) %in% c("ID","Disease")]<-apply(FFA[!colnames(FFA)%in% c("ID","Disease")], MARGIN = 2, as.numeric)
FFA$Disease<-factor(FFA$Disease, levels = c("HC","ACS","SAH"))

colnames(FFA)
SCFA<-FFA[,1:9] # from acetic to valeric 
MCFA<-FFA[,c(1,2,11:15)] # to dodecanoic
LCFA<-FFA[,c(1,2,17:19)] # to Octadecanoic, before the raw total columns
colnames(LCFA)

SCFA_raw<-as.data.frame(SCFA)
MCFA_raw<-as.data.frame(MCFA)
LCFA_raw<- as.data.frame(LCFA)

### raw graph and k-w and dunn test files
list <- list("SCFA" = SCFA_raw, "MCFA" = MCFA_raw, "LCFA" = LCFA_raw)
names(list)

for (y in 1:length(list)) {
  temp <- list[[y]]
  results <- NULL
  posthoc_results <- NULL
  
  for (x in colnames(temp)[!colnames(temp) %in% c("ID", "Disease")]) {
    # Kruskal-Wallis Test
    test <- kruskal.test(temp[[x]] ~ temp$Disease)
    kruskal <- cbind(paste(x), test$statistic, test$p.value)
    results <- rbind.data.frame(results, kruskal)
    
    # Post-hoc Dunn Test
    dunn <- dunnTest(temp[[x]] ~ temp$Disease, method = "bh")
    dunn_res <- cbind(variable = paste(x), dunn$res)
    posthoc_results <- rbind.data.frame(posthoc_results, dunn_res)
  }
  
  colnames(results) <- c(paste(names(list)[y]), "Chi-squared", "p.value")
  results$'p.adj(BH)' <- p.adjust(results$p.value, method = "BH")
  row.names(results) <- NULL
  write.csv2(results, file = paste0("FFA/", names(list)[y], "_raw_abund_KW.csv"), row.names = F)
  write.csv2(posthoc_results, file = paste0("FFA/", names(list)[y], "_raw_abund_Dunn.csv"), row.names = F)
  table <- melt(data = temp, id.vars = c("ID", "Disease"))
  if (y == 1) {
    table_SCFA <- table 
  }
  
  # Plot
  ggplot(data = table, mapping = aes(x = variable, y = value, fill = Disease)) +
    scale_fill_manual(values = c("HC" = "#32CD32", "ACS" = "#FF4500", "SAH" = 'yellow3')) +
    theme(axis.text.x = element_text(angle = -25, size = 8, vjust = 1, hjust = 0.05)) +
    geom_boxplot() +
    labs(title = paste("Raw", names(list)[y], "quantities", sep = " "), x = "", y = "Quantity")
  
  ggsave(filename = paste0("FFA/", names(list)[y], "_raw_abundance.png"), height = 3, width = 6, dpi = 300)
  
  # Save each category
  if (y == 2) {
    table_MCFA <- table
  }
  if (y == 3) {
    table_LCFA <- table
  }
}

####################### Kruskal-wallis TESTs ####################
###### loop for each FFA
library(stringr)
### Stool
Stool<-FFA[,c(1:9,11:15,17:19)]
Stool$Disease<-trimws(Stool$Disease)
Stool$ID<-trimws(Stool$ID)
head(Stool)
Stool<-as.data.frame(Stool)
Stool<-reshape::melt(Stool)
Stool$Disease
Stool$variable<-factor(Stool$variable, levels=unique(Stool$variable))
Stool$Disease<-factor(Stool$Disease, levels = c("HC","ACS","SAH"))
Stool$value<-Stool$value*100 # percentual

for(i in unique(Stool$variable)){
  target<-Stool[Stool$variable==i , ]
  
  # statistics
  test <- dunn.test::dunn.test(target$value, target$Disease, method = "bh",list = T) # pairwise
  test <- cbind.data.frame(test$comparison, p_adj_bh=round(test$P.adjusted,3))
  #test$p_adj_bh[test$p_adj_bh==0.000]<-0.001 # if...
  specular_vector<-test$p_adj_bh # can not modify the numeric vector with characters at each step, otherwise every number will became a character itself!
  test$p_adj_bh[specular_vector <=0.001]<-"***"
  test$p_adj_bh[specular_vector >0.001 & specular_vector <=0.01]<-"**"
  test$p_adj_bh[specular_vector >0.01 & specular_vector <=0.05]<-"*"
  
  if(length(which( grepl("*",test$p_adj_bh, fixed=T) )) > 0 ) { # then if there are results (otherwise errors with empty tables!)
    test<-test[ grepl("*",test$p_adj_bh, fixed=T) , ]  # mainteining only significant rows
    comparisons<-unlist(str_split(test$`test$comparison`, pattern ="-")) # it gives a unique vector, not two columns
    comparisons<-gsub(" ","",comparisons)
    # NB: c(T,F) --> only even indexes,   c(F,T) --> only pair indexes   (see the output of unlist-str_split functions, the even indexes are the first pseudo-column elements)
    test$group1<-comparisons[c(T,F)] # 'group1' and 'group2' HAVE to be the col names according to stat_pvalue_manual function
    test$group2<-comparisons[c(F,T)]
    # test # a 'manual' check --> OK
    test<-test[ seq(length(test$`test$comparison`),1,-1) , ] # inverting the order of comparisons to enhance the plot aesthetics (the last row here is the one with most width), moreover using and inverse seq instead of c(3,2,1) to adapt each loop to the row numbers (=number of signif results)
    test$`test$comparison`<-NULL
    
    kruskal_res<-kruskal.test(target$value~target$Disease)
    kruskal_res<-round(kruskal_res$p.value, 3)
    kruskal_res[kruskal_res==0]<-0.001 # if...
    
    ggboxplot(data=target, x="Disease", y="value", fill="Disease", 
              width=0.6, size= 0.15, alpha= 0.2, outlier.shape = NA) +
      stat_pvalue_manual(test, label="p_adj_bh",
                         y.position = max(target$value),
                         step.increase = 0.065,
                         size = 2.5,
                         bracket.size = 0.15
      ) +
      scale_fill_manual(values=c("HC" = "#32CD32", "ACS" = "#FF4500", "SAH" = 'yellow3')) +
      scale_color_manual(values=c("HC" = "#32CD32", "ACS" = "#FF4500", "SAH" = 'yellow3')) +
      geom_point(position = position_jitterdodge(seed = 1994, dodge.width = 0.6, jitter.width = 0.9),
                 aes(color=Disease), size= 0.3, alpha= 0.5) +  
      theme_classic2(base_size = 7) + 
      theme(strip.text.x=element_text(size=12,colour="black"),
            axis.text.x = element_text(angle = 0, vjust=0.5, hjust=0.5, size=7), 
            axis.text.y = element_text(size=4.5),
            plot.margin = margin(5,2,0,2),
            plot.title= element_text(size=10, hjust = 0.5, vjust=1.8) ,
            plot.subtitle = element_text(size=7.4, vjust=1) ,
            legend.key.size=unit(0.8,"cm"), 
            legend.text=element_text(size=7),
            panel.grid.major.y = element_line(size=0.12, color="gray"),
            panel.grid.minor.y = element_line(size=0.03, color="gray"),
            panel.grid.major.x = element_blank()
      ) +
      guides( color="none", fill="none" ) +
      labs(y="Percentual quantity", x="",
           title = i,
           subtitle = paste("Kruskal-Wallis p-value:", kruskal_res)) + 
      scale_x_discrete(expand=c(0.2, 0))
    
    ggsave(file=paste0("FFA/Stool_",i,"_with_Dunnett_pairwise_corrected_though_Benjamini_H.png"), width = 2.5, height = 2.6, dpi=300)
  }
  
  if(length(which( grepl("*",test$p_adj_bh, fixed=T) )) == 0 ) { # otherwise, if there are not pairwise results...
    
    kruskal_res<-kruskal.test(target$value~target$Disease)
    kruskal_res<-round(kruskal_res$p.value, 3)
    kruskal_res[kruskal_res==0]<-0.001 # if...
    
    ggboxplot(data=target, x="Disease", y="value", fill="Disease", 
              width=0.6, size= 0.15, alpha= 0.2, outlier.shape = NA) +
      scale_fill_manual(values=c("HC" = "#32CD32", "ACS" = "#FF4500", "SAH" = 'yellow3')) +
      scale_color_manual(values=c("HC" = "#32CD32", "ACS" = "#FF4500", "SAH" = 'yellow3')) +
      geom_point(position = position_jitterdodge(seed = 1994, dodge.width = 0.6, jitter.width = 0.9),
                 aes(color=Disease), size= 0.3, alpha= 0.5) +  
      theme_classic2(base_size = 7) + 
      theme(strip.text.x=element_text(size=12,colour="black"),
            axis.text.x = element_text(angle = 0, vjust=0.5, hjust=0.5, size=7), 
            axis.text.y = element_text(size=5),
            plot.margin = margin(5,2,0,2),
            plot.title= element_text(size=10, hjust = 0.5, vjust=1.8) ,
            plot.subtitle = element_text(size=7.4, vjust=1) ,
            legend.key.size=unit(0.8,"cm"), 
            legend.text=element_text(size=7),
            panel.grid.major.y = element_line(size=0.12, color="gray"),
            panel.grid.minor.y = element_line(size=0.03, color="gray"),
            panel.grid.major.x = element_blank()
      ) +
      guides( color="none", fill="none" ) +
      labs(y="Percentual quantity", x="",
           title = i,
           subtitle = paste("Kruskal-Wallis p-value:", kruskal_res)) + 
      scale_x_discrete(expand=c(0.2, 0))
    
    ggsave(file=paste0("FFA/Stool_",i,"_with_Dunnett_pairwise_corrected_though_Benjamini_H.png"), width = 2.5, height = 2.6, dpi=300)
  }
  
}

#total 
colnames(FFA)
SCFA_t<-FFA[,c(1,2,10)] # from acetic to valeric 
MCFA_t<-FFA[,c(1,2,16)] # to dodecanoic
LCFA_t<-FFA[,c(1,2,20)] # to Octadecanoic, before the raw total columns
colnames(LCFA_t)

SCFA_raw_t<-as.data.frame(SCFA_t)
MCFA_raw_t<-as.data.frame(MCFA_t)
LCFA_raw_t<- as.data.frame(LCFA_t)

list <- list("SCFA" = SCFA_raw_t, "MCFA" = MCFA_raw_t, "LCFA" = LCFA_raw_t)
names(list)

for (y in 1:length(list)) {
  temp <- list[[y]]
  results <- NULL
  posthoc_results <- NULL
  
  for (x in colnames(temp)[!colnames(temp) %in% c("ID", "Disease")]) {
    # Kruskal-Wallis Test
    test <- kruskal.test(temp[[x]] ~ temp$Disease)
    kruskal <- cbind(paste(x), test$statistic, test$p.value)
    results <- rbind.data.frame(results, kruskal)
    
    # Post-hoc Dunn Test
    dunn <- dunnTest(temp[[x]] ~ temp$Disease, method = "bh")
    dunn_res <- cbind(variable = paste(x), dunn$res)
    posthoc_results <- rbind.data.frame(posthoc_results, dunn_res)
  }
  
  colnames(results) <- c(paste(names(list)[y]), "Chi-squared", "p.value")
  results$'p.adj(BH)' <- p.adjust(results$p.value, method = "BH")
  row.names(results) <- NULL
  
  # Kruskal-Wallis
  write.csv2(results, file = paste0("FFA/", names(list)[y], "_raw_abund_KW.csv"), row.names = F)
  
  # post hoc Dunn
  write.csv2(posthoc_results, file = paste0("FFA/", names(list)[y], "_raw_abund_Dunn.csv"), row.names = F)
  table <- melt(data = temp, id.vars = c("ID", "Disease"))
  if (y == 1) {
    table_SCFA <- table 
  }
  
  # Plot
  ggplot(data = table, mapping = aes(x = variable, y = value, fill = Disease)) +
    scale_fill_manual(values = c("HC" = "#32CD32", "ACS" = "#FF4500", "SAH" = 'yellow3')) +
    theme(axis.text.x = element_text(angle = -25, size = 8, vjust = 1, hjust = 0.05)) +
    geom_boxplot() +
    labs(title = paste("Raw", names(list)[y], "quantities", sep = " "), x = "", y = "Quantity")
  
  ggsave(filename = paste0("FFA/", names(list)[y], "_raw_abundance.png"), height = 3, width = 6, dpi = 300)
  
  # Salvataggio per ogni categoria
  if (y == 2) {
    table_MCFA <- table
  }
  if (y == 3) {
    table_LCFA <- table
  }
}

Stool<-FFA[,c(1,2,10,16,20)]
Stool$Disease<-trimws(Stool$Disease)
Stool$ID<-trimws(Stool$ID)
head(Stool)
Stool<-as.data.frame(Stool)
Stool<-reshape::melt(Stool)
Stool$Disease
Stool$variable<-factor(Stool$variable, levels=unique(Stool$variable))
Stool$Disease<-factor(Stool$Disease, levels = c("HC","ACS","SAH"))
Stool$value<-Stool$value*100 # percentual

for(i in unique(Stool$variable)){
  target<-Stool[Stool$variable==i , ]
  
  # statistics
  test <- dunn.test::dunn.test(target$value, target$Disease, method = "bh",list = T) # pairwise
  test <- cbind.data.frame(test$comparison, p_adj_bh=round(test$P.adjusted,3))
  #test$p_adj_bh[test$p_adj_bh==0.000]<-0.001 # if...
  specular_vector<-test$p_adj_bh # can not modify the numeric vector with characters at each step, otherwise every number will became a character itself!
  test$p_adj_bh[specular_vector <=0.001]<-"***"
  test$p_adj_bh[specular_vector >0.001 & specular_vector <=0.01]<-"**"
  test$p_adj_bh[specular_vector >0.01 & specular_vector <=0.05]<-"*"
  
  if(length(which( grepl("*",test$p_adj_bh, fixed=T) )) > 0 ) { # then if there are results (otherwise errors with empty tables!)
    test<-test[ grepl("*",test$p_adj_bh, fixed=T) , ]  # mainteining only significant rows
    comparisons<-unlist(str_split(test$`test$comparison`, pattern ="-")) # it gives a unique vector, not two columns
    comparisons<-gsub(" ","",comparisons)
    # NB: c(T,F) --> only even indexes,   c(F,T) --> only pair indexes   (see the output of unlist-str_split functions, the even indexes are the first pseudo-column elements)
    test$group1<-comparisons[c(T,F)] # 'group1' and 'group2' HAVE to be the col names according to stat_pvalue_manual function
    test$group2<-comparisons[c(F,T)]
    # test # a 'manual' check --> OK
    test<-test[ seq(length(test$`test$comparison`),1,-1) , ] # inverting the order of comparisons to enhance the plot aesthetics (the last row here is the one with most width), moreover using and inverse seq instead of c(3,2,1) to adapt each loop to the row numbers (=number of signif results)
    test$`test$comparison`<-NULL
    
    kruskal_res<-kruskal.test(target$value~target$Disease)
    kruskal_res<-round(kruskal_res$p.value, 3)
    kruskal_res[kruskal_res==0]<-0.001 # if...
    
    ggboxplot(data=target, x="Disease", y="value", fill="Disease", 
              width=0.6, size= 0.15, alpha= 0.2, outlier.shape = NA) +
      stat_pvalue_manual(test, label="p_adj_bh",
                         y.position = max(target$value),
                         step.increase = 0.065,
                         size = 2.5,
                         bracket.size = 0.15
      ) +
      scale_fill_manual(values=c("HC" = "#32CD32", "ACS" = "#FF4500", "SAH" = 'yellow3')) +
      scale_color_manual(values=c("HC" = "#32CD32", "ACS" = "#FF4500", "SAH" = 'yellow3')) +
      geom_point(position = position_jitterdodge(seed = 1994, dodge.width = 0.6, jitter.width = 0.9),
                 aes(color=Disease), size= 0.3, alpha= 0.5) +  
      theme_classic2(base_size = 7) + 
      theme(strip.text.x=element_text(size=12,colour="black"),
            axis.text.x = element_text(angle = 0, vjust=0.5, hjust=0.5, size=7), 
            axis.text.y = element_text(size=4.5),
            plot.margin = margin(5,2,0,2),
            plot.title= element_text(size=10, hjust = 0.5, vjust=1.8) ,
            plot.subtitle = element_text(size=7.4, vjust=1) ,
            legend.key.size=unit(0.8,"cm"), 
            legend.text=element_text(size=7),
            panel.grid.major.y = element_line(size=0.12, color="gray"),
            panel.grid.minor.y = element_line(size=0.03, color="gray"),
            panel.grid.major.x = element_blank()
      ) +
      guides( color="none", fill="none" ) +
      labs(y="Percentual quantity", x="",
           title = i,
           subtitle = paste("Kruskal-Wallis p-value:", kruskal_res)) + 
      scale_x_discrete(expand=c(0.2, 0))
    
    ggsave(file=paste0("FFA/Stool_",i,"_with_Dunnett_pairwise_corrected_though_Benjamini_H.png"), width = 2.5, height = 2.6, dpi=300)
  }
  
  if(length(which( grepl("*",test$p_adj_bh, fixed=T) )) == 0 ) { # otherwise, if there are not pairwise results...
    
    kruskal_res<-kruskal.test(target$value~target$Disease)
    kruskal_res<-round(kruskal_res$p.value, 3)
    kruskal_res[kruskal_res==0]<-0.001 # if...
    
    ggboxplot(data=target, x="Disease", y="value", fill="Disease", 
              width=0.6, size= 0.15, alpha= 0.2, outlier.shape = NA) +
      scale_fill_manual(values=c("HC" = "#32CD32", "ACS" = "#FF4500", "SAH" = 'yellow3')) +
      scale_color_manual(values=c("HC" = "#32CD32", "ACS" = "#FF4500", "SAH" = 'yellow3')) +
      geom_point(position = position_jitterdodge(seed = 1994, dodge.width = 0.6, jitter.width = 0.9),
                 aes(color=Disease), size= 0.3, alpha= 0.5) +  
      theme_classic2(base_size = 7) + 
      theme(strip.text.x=element_text(size=12,colour="black"),
            axis.text.x = element_text(angle = 0, vjust=0.5, hjust=0.5, size=7), 
            axis.text.y = element_text(size=5),
            plot.margin = margin(5,2,0,2),
            plot.title= element_text(size=10, hjust = 0.5, vjust=1.8) ,
            plot.subtitle = element_text(size=7.4, vjust=1) ,
            legend.key.size=unit(0.8,"cm"), 
            legend.text=element_text(size=7),
            panel.grid.major.y = element_line(size=0.12, color="gray"),
            panel.grid.minor.y = element_line(size=0.03, color="gray"),
            panel.grid.major.x = element_blank()
      ) +
      guides( color="none", fill="none" ) +
      labs(y="Percentual quantity", x="",
           title = i,
           subtitle = paste("Kruskal-Wallis p-value:", kruskal_res)) + 
      scale_x_discrete(expand=c(0.2, 0))
    
    ggsave(file=paste0("FFA/Stool_",i,"_with_Dunnett_pairwise_corrected_though_Benjamini_H.png"), width = 2.5, height = 2.6, dpi=300)
  }
  
}

#PCoA
list_PCoA <- list("SCFA" = SCFA_raw, "MCFA" = MCFA_raw, "LCFA" = LCFA_raw)

p_values <- NULL
for (x in 1:length(list_PCoA)) {
  rm(data, coord)
  data <- as.data.frame(list_PCoA[[x]])
  row.names(data) <- data$ID
  
  if (x == 3) {
    data <- data[!is.na(data$Octadecanoic),] # rimuovere NA in LCFA
  }
  
  dist <- vegdist(data[, !colnames(data) %in% c("ID", "Disease")], method = "bray") # Bray-Curtis dissimilarity
  coord <- cmdscale(dist, k = 2, eig = TRUE) # PCoA

  coord_df <- as.data.frame(coord$points)
  colnames(coord_df) <- c("PC1", "PC2")
  coord_df$ID <- row.names(coord_df)
  coord_df$Disease <- data$Disease
  
  ggplot(data = coord_df, aes(x = PC1, y = PC2, color = Disease)) + 
    geom_point(size = 3, alpha = 1) +
    theme_classic() + 
    scale_color_manual(values = c("HC" = alpha("#32CD32", 0.5), "ACS" = alpha("#FF4500", 0.5), "SAH" = alpha ('yellow3', 0.5))) +
    stat_ellipse(size = 0.2) +
    geom_text(aes(label = ID), size = 2, color = "black") +
    labs(
      x = paste("PC1:", round(coord$eig[1] / sum(coord$eig) * 100, digits = 2), "%"),
      y = paste("PC2:", round(coord$eig[2] / sum(coord$eig) * 100, digits = 2), "%"),
      color = "", 
      title = paste("PCoA computed on", names(list_PCoA)[x], "quantities\n of HC, ASC and SAH subjects\n with Bray-Curtis dissimilarity")
    )
  
  ggsave(
    file = paste("PCoA_BRAY_", names(list_PCoA)[x], ".png", sep = "_"), 
    height = 5, 
    width = 6.5, 
    dpi = 300
  )     
  # again but with no names
  ggplot(data = coord_df, aes(x = PC1, y = PC2, color = Disease)) + 
    geom_point(size = 3, alpha = 1) +
    theme_classic() + 
    scale_color_manual(values = c("HC" = alpha("#32CD32", 0.5), "ACS" = alpha("#FF4500", 0.5), "SAH" = alpha ('yellow3', 0.5))) +
    stat_ellipse(size = 0.2) +
    labs(
      x = paste("PC1:", round(coord$eig[1] / sum(coord$eig) * 100, digits = 2), "%"),
      y = paste("PC2:", round(coord$eig[2] / sum(coord$eig) * 100, digits = 2), "%"),
      color = "", 
      title = paste("PCoA computed on", names(list_PCoA)[x], "quantities\n of HC, ASC and SAH subjects\n with Bray-Curtis dissimilarity")
    )
  
  ggsave(
    file = paste("PCoA_BRAY_", names(list_PCoA)[x], "_noname.png", sep = "_"), 
    height = 5, 
    width = 6.5, 
    dpi = 300
  )                        
  
  # checking for significant dispersion
  suppressWarnings(rm(p_value_diver, p_value_disp))
  p_value_diver<-adonis(dist ~ data$Disease, data = data, permutations = 9999)
  p_value_disp<-permutest(betadisper(dist, data$Disease), permutations = 9999)
  p_value_disp<-as.data.frame(p_value_disp$tab[1,])
  colnames(p_value_disp)<-colnames(p_value_diver$aov.tab[1,])
  
  new_values<-rbind(p_value_diver$aov.tab[1,],p_value_disp)
  row.names(new_values)<-c(paste("Beta_diversity",names(list_PCoA)[x], sep = "_"),
                           paste("Beta_dispersion",names(list_PCoA)[x], sep = "_"))
  p_values<-rbind(p_values,new_values)
}  

p_values$adjusted_BH_p_values<-p.adjust(p_values$`Pr(>F)`, method = "BH")
p_values
write.csv2(p_values, file="p_values_BRAY.csv", quote = F,row.names = T)


pairwise_results <- pairwise.adonis2(dist ~ Disease, 
                                     data = SCFA_raw,
                                     permutations = 9999, method='fdr')
print(pairwise_results)
write.csv2(pairwise_results, file="p_values_BRAY_SCFA.csv", quote = F,row.names = T)


pairwise_results <- pairwise.adonis2(dist ~ Disease, 
                                     data = MCFA_raw,
                                     permutations = 9999, method='fdr')
print(pairwise_results)
write.csv2(pairwise_results, file="p_values_BRAY_MCFA.csv", quote = F,row.names = T)


pairwise_results <- pairwise.adonis2(dist ~ Disease, 
                                     data = LCFA_raw,
                                     permutations = 9999, method='fdr')
print(pairwise_results)
write.csv2(pairwise_results, file="p_values_BRAY_LCFA.csv", quote = F,row.names = T)

















