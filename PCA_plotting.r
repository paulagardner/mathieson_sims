setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# https://stackoverflow.com/questions/3452086/getting-path-of-an-r-script/35842176#35842176

ret_pca <- read.table("common_PCA.eigenvec")
print(ret_pca)
plot(ret_pca$V3, ret_pca$V4)
# https://www.biostars.org/p/487372/

# going to brute-force assigning colors, one to the first 50 (presumably from pop A) and one to the second, since deme id's not written in by the treefile


library(dplyr)
# popA <- (ret_pca[2:52,])
# popB <- (ret_pca[53:101,])
# new_col = (popA$popB)
# new_col


# df$deme_id <- ifelse(df[2:52,] = 1, df[52:102] =2)

# ret_pca$deme_id <- NA
# ret_pca
# ret_pca$row_number <- ret_pca %>% mutate(row_number = 1:n())
# ret_pca$row_number <- ret_pca %>% mutate(1:nrow(ret_pca[,1]))
ret_pca$row_number <- 1:nrow(ret_pca)
ret_pca
ret_pca$deme_id <- ifelse(ret_pca$row_number >= 0 & ret_pca$row_number <= 50, 1, 2) # galen's suggestion
ret_pca


new_df <- data.frame(matrix(ncol = 1, nrow = 100))
new_df$deme_id <- new_df %>% mutate(ifelse(ret_pca$row_number >= 0 & ret_pca$row_number <= 50, 1, 2)) # galen's suggestion
new_df
# if (ret_pca$row_number <= 50){ret_pca$deme_id = 1}
# else {ret_pca$deme_id = 2}




# ret_pca %>% mutate(if df$row_number <= 50 {ret_pca$newrow =1}
#     else {ret_pca$newrow = 2}
#     )

# plot(ret_pca$V3, ret_pca$V4, color = ret_pca$deme_id)
# plot(ret_pca$V3, ret_pca$V4)

library(ggplot2)
plot <- ggplot(ret_pca, aes(x = V3, y = V4, color = deme_id)) +
    geom_point()
plot