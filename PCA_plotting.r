#!/usr/bin/env Rscript
#chmod +x script.r

#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# https://stackoverflow.com/questions/3452086/getting-path-of-an-r-script/35842176#35842176

common_pca <- read.table("common_PCA.eigenvec")
print(common_pca)

# https://www.biostars.org/p/487372/

# going to brute-force assigning colors, one to the first 50 (presumably from pop A) and one to the second, since deme id's not written in by the treefile

library(dplyr)
# popA <- (common_pca[2:52,])
# popB <- (common_pca[53:101,])
# new_col = (popA$popB)
# new_col


# df$deme_id <- ifelse(df[2:52,] = 1, df[52:102] =2)

# common_pca$deme_id <- NA
# common_pca
#common_pca$row_number <- common_pca %>% mutate(row_number = 1:n())
# common_pca$row_number <- common_pca %>% mutate(1:nrow(common_pca[,1]))
common_pca$row_number <- 1:nrow(common_pca)

#common_pca
common_pca$deme_id <- ifelse(common_pca$row_number >= 0 & common_pca$row_number <= 500, 1, 2) # galen's suggestion
common_pca


new_df <- data.frame(matrix(ncol = 1, nrow = 500))
new_df$deme_id <- new_df %>% mutate(ifelse(common_pca$row_number >= 0 & common_pca$row_number <= 500, 1, 2)) # galen's suggestion
#new_df
# if (common_pca$row_number <= 50){common_pca$deme_id = 1}
# else {common_pca$deme_id = 2}



#common_pca %>% mutate(if df$row_number <= 50 {common_pca$newrow =1}
#     else {common_pca$newrow = 2}
#     )

# plot(common_pca$V3, common_pca$V4, color = common_pca$deme_id)
# plot(common_pca$V3, common_pca$V4)

library(ggplot2)
plot <- ggplot(common_pca, aes(x = V3, y = V4, color = deme_id)) +
    geom_point()
plot <- plot + xlab("PC1")+ ylab("PC2") + theme(axis.text=element_text(size=20), axis.title =element_text(size=20))

ggsave("common_PCA_plot.jpg", plot = plot)

rare_pca <- read.table("rare_PCA.eigenvec")
rare_pca$row_number <- 1:nrow(rare_pca)
rare_pca
rare_pca$deme_id <- ifelse(rare_pca$row_number >= 0 & rare_pca$row_number <= 500, 1, 2) # galen's suggestion
rare_pca


new_df_2 <- data.frame(matrix(ncol = 1, nrow = 500))
new_df_2$deme_id <- new_df_2 %>% mutate(ifelse(rare_pca$row_number >= 0 & rare_pca$row_number <= 500, 1, 2)) # galen's suggestion
new_df_2

library(ggplot2)
plot2 <- ggplot(rare_pca, aes(x = V3, y = V4, color = deme_id)) +
    geom_point() 
plot2 <- plot2 + xlab("PC1")+ ylab("PC2") +  theme(axis.text=element_text(size=20), axis.title =element_text(size=20))
plot2

ggsave("rare_pca_plot.jpg", plot = plot2)

library(patchwork)

combined_plot <- (plot / plot2) + plot_annotation(title = "perpetual structure") & theme(plot.title = element_text(hjust = 0.5, size = 20))
combined_plot
ggsave("PCAs.png", combined_plot)













