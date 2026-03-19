pacman::p_load("vroom", "tidyverse")

setwd("~/Desktop/local_work/fomo_networks/DEGS/")
files <- list.files(path = "~/Desktop/local_work/fomo_networks/DEGS/")

X<- lapply(files, vroom::vroom)
X<- lapply(X, as.data.frame)
names(X) <- tools::file_path_sans_ext(files)

# Filtrar DEGs
DEGS <- lapply(X, function(df) {
  df %>% dplyr::filter(abs(logFC) >= 0.5 | adj.P.Val < 0.05)
})

gene_of_interest <- "ENSG00000078668"

presence <- sapply(DEGS, function(df) {
  gene_of_interest %in% df$gene
})

presence

matches <- lapply(DEGS, function(df) {
  df[df$gene == gene_of_interest, ]
})

# Quitar los vacíos
matches <- matches[sapply(matches, nrow) > 0]

matches

library(dplyr)
library(ggplot2)

gene_of_interest <- "ENSG00000078668"

df_plot <- bind_rows(DEGS, .id = "comparison") %>%
  filter(gene == gene_of_interest)

ggplot(df_plot, aes(x = comparison, y = logFC)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  theme_minimal()
df_all <- bind_rows(X, .id = "comparison")

ggplot(df_all, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(alpha = 0.4) +
  geom_point(
    data = df_all %>% filter(gene == gene_of_interest),
    size = 3
  ) +
  facet_wrap(~comparison) +
  theme_minimal()

