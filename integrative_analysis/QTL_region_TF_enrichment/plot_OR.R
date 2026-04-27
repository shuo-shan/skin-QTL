library(ggplot2)
library(dplyr)

##### this part is for promoter and enhancers
# Sample data
data <- data.frame(
  odds_ratio = c(3.09, 3.75, 3.7, 1.71, 2.66, 3.34, 3.04, 3.72, 3.76, 4.32, 6.79, 7.05, 2.99, 5.6, 5.66, 1.96, 4.39, 6.63),
  CI_lower = c(1.65, 2.54, 2.66, 0.76, 2.09, 2.68, 1.44, 2.53, 2.59, 1.61, 4.31, 4.57, 0.96, 3.17, 3.5, 0.49, 3.11, 4.96),
  CI_upper = c(5.78, 5.54, 5.16, 3.82, 3.38, 4.16, 6.42, 5.49, 5.47, 11.56, 10.69, 10.85, 9.3, 9.91, 9.14, 7.87, 6.2, 8.86),
  tag = c("MEL_reQTL", "MEL_PBSeQTL", "MEL_IFNeQTL", "FRB_reQTL", "FRB_PBSeQTL", "FRB_IFNeQTL", "KRT_reQTL", "KRT_PBSeQTL", "KRT_IFNeQTL", "KRT_reQTL", "KRT_PBSeQTL", "KRT_IFNeQTL", "MEL_reQTL", "MEL_PBSeQTL", "MEL_IFNeQTL", "FRB_reQTL", "FRB_PBSeQTL", "FRB_IFNeQTL"),
  region_type = c("enhancer", "enhancer", "enhancer", "enhancer", "enhancer", "enhancer", "enhancer", "enhancer", "enhancer", "promoter", "promoter", "promoter", "promoter", "promoter", "promoter", "promoter", "promoter", "promoter")
)

# Adding a midpoint for plotting
data <- data %>%
  mutate(midpoint = odds_ratio) %>%
  mutate(separate(data, tag, into = c("celltype","QTL_type"), sep="_", remove=F))

# Plot
p <- ggplot(data, aes(x = midpoint, y = tag, color = region_type)) +
  geom_pointrange(aes(xmin = CI_lower, xmax = CI_upper),
                  position = position_dodge(width = 0.3), 
                  size = 1.2) +
  geom_vline(xintercept = 1.0, linetype = "dotted", linewidth = 1) +
  scale_x_log10(breaks = c(0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10),
                minor_breaks = NULL) +
  labs(x = "Odds Ratio", y = "QTL Type") 

# Show plot
print(p)

#### this part is for TF binding sites
data <- data.table::fread("~/Downloads/nl/human/skin/eQTLs/integrative_analysis/QTL_region_TF_enrichment/ODDS_RATIO_95CI.txt")
colnames(data) <- c("tag","region","odds_ratio","CI_lower","CI_upper","QTL_in_region","QTL_notIn_region","nonQTL_in_region","nonQTL_notIn_region","fisher_test_pval")
data <- data %>%
  mutate(midpoint = odds_ratio) %>%
  mutate(separate(data, tag, into = c("celltype","QTL_type"), sep="_", remove=F)) %>%
  as.data.frame()
data$region <- gsub("_peak","",data$region)

p.reQTL.mel <- data %>% 
  dplyr::filter(QTL_type=="reQTL") %>%
  dplyr::filter(celltype=="MEL") %>%
  dplyr::filter(fisher_test_pval < 0.05) %>%
  arrange(odds_ratio) %>%
  na.omit() %>%
  mutate(region = factor(region, levels = unique(region))) %>%
  ggplot(., aes(y = midpoint, x = region)) +
  geom_pointrange(aes(ymin = CI_lower, ymax = CI_upper),
                  position = position_dodge(width = 0.3), 
                  size = 1.2) +
  geom_hline(yintercept = 1.0, linetype = "dotted", linewidth = 1) +
  scale_y_log10(breaks = c(0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10),
                minor_breaks = NULL) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5, face="bold")) +
  labs(x = "", y = "Odds Ratio with 95% CI") +
  ggtitle("MEL reQTL")


  
