#!/usr/bin/env Rscript

#core_modules.R

#Packages
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}

pacman::p_load(
 "data.table", 
  "tidyverse",
  "tools",
  "optparse",
  "ggplot2",
  "igraph", 
  "cowplot", 
  "pheatmap",
 "clusterProfiler"
)

#Parser
option_list <- list(
  make_option(c("-i", "--input_dir"), type = "character", default = getwd(),
              help = "Directorio con archivos *_nodes_summary.csv", metavar = "character"),
  make_option(c("-o", "--output_dir"), type = "character", default = "./results_jaccard",
              help = "Directorio de salida [default: %default]", metavar = "character"),
  make_option(c("-p", "--pattern"), type = "character", default = "_nodes_summary\\.csv$",
              help = "Regex para seleccionar archivos [default: %default]", metavar = "character"), 
  make_option(c("-r", "--enrich_dir"), type = "character", default = "./results_comm",
              help = "Directorio de enriquecimientos [default: %default]", metavar = "character"))

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
#Inputs
opt$input_dir <- "~/Desktop/local_work/fomo_networks/results_topos_louvain"
opt$output_dir <- "~/Desktop/local_work/fomo_networks/results_core_modules"
opt$enrich_dir <- "~/Desktop/local_work/fomo_networks/results_topos_louvain/results_comm/"

input_dir <- opt$input_dir
output_dir <- opt$output_dir
pattern <- opt$pattern
enrich_dir <- opt$enrich_dir

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
setwd(output_dir)

#Functions --- ---

#Extract metadata
extract_info <- function(filename) {
  base <- basename(filename)
  if (grepl("_counts_", base)) {
    region <- sub("_counts_.*", "", base)
    phenotype <- ifelse(grepl("_AD_", base), "AD", "Control")
    return(data.frame(filename = filename, region = region, phenotype = phenotype))
  } else {
    return(NULL)
  }
}

#Load modules with unique names
load_modules <- function(file) {
  df <- vroom::vroom(file)
  if (!all(c("node", "membership") %in% names(df))) {
    stop("File ", file, " must contain columns 'node' and 'membership'")
  }
  region <- sub("_counts_.*", "", basename(file))
  region <- gsub("^(Mayo_|ROSMAP_)", "", region)
  phenotype <- ifelse(grepl("_AD_", file), "AD", "Control")
  df$unique_module_id <- paste(region, phenotype, df$membership, sep = "_")
  split(df$node, df$unique_module_id)
}

#Jaccard index
jaccard_index <- function(a, b) {
  if (length(a) == 0 && length(b) == 0) return(1)
  if (length(a) == 0 || length(b) == 0) return(0)
  length(intersect(a, b)) / length(union(a, b))
}

#Function to compute NMI for each region
compute_nmi_region <- function(region_name, meta_df, min_genes = 3) {
  submeta <- meta_df %>% dplyr::filter(region == region_name)
  
  #We need both phenos
  if (!all(c("AD", "Control") %in% submeta$phenotype)) {
    message("Skip region ", region_name, " because it does not have both phenotypes")
    return(NULL)
  }
  
  file_ad   <- submeta$filename[submeta$phenotype == "AD"]
  file_ctrl <- submeta$filename[submeta$phenotype == "Control"]
  
  #Only the necessary cols
  ad   <- data.table::fread(file_ad, select = c("node", "membership"))
  ctrl <- data.table::fread(file_ctrl, select = c("node", "membership"))
  
  #Filter genes by minimal number of genes
  if (!is.null(min_genes)) {
    keep_ad   <- names(which(table(ad$membership)   >= min_genes))
    keep_ctrl <- names(which(table(ctrl$membership) >= min_genes))
    
    ad   <- ad   %>% dplyr::filter(membership %in% keep_ad)
    ctrl <- ctrl %>% dplyr::filter(membership %in% keep_ctrl)
  }
  
  #Merge
  merged <- dplyr::inner_join(ad, ctrl, by = "node", suffix = c("_AD", "_CTRL"))
  
  if (nrow(merged) == 0) {
    message("Region ", region_name, ": 0 genes after filtering. NMI = NA.")
    return(
      tibble::tibble(
        region = region_name,
        n_genes_common = 0L,
        NMI = NA_real_
      )
    )
  }
  
  #NMI w/igraph
  nmi_val <- igraph::compare(
    merged$membership_AD,
    merged$membership_CTRL,
    method = "nmi"
  )
  
  tibble::tibble(
    Region = region_name,
    n_genes_common = nrow(merged),
    NMI = nmi_val
  )
}


#Build like an enrich result

create_enrichResult <- function(df,
                                ontology = "MF",
                                organism = "Homo sapiens",
                                pAdjustMethod = "BH",
                                pCutoff = 0.05,
                                qCutoff = 0.2) {
  
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("clusterProfiler package is required")
  }
  
  #Check required columns
  required_cols <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", 
                     "p.adjust", "qvalue", "geneID", "Count")
  
  if (!all(required_cols %in% colnames(df))) {
    stop("La tabla debe contener las columnas: ", paste(required_cols, collapse = ", "))
  }
  
  #Vector of unique analyzed genes
  gene <- unique(unlist(strsplit(df$geneID, split = "/")))
  
  #Get TERM2GENE list
  geneSets <- df %>%
    tidyr::separate_rows(geneID, sep = "/") %>%
    dplyr::group_by(ID) %>%
    dplyr::summarise(genes = list(unique(geneID)), .groups = "drop") %>%
    tibble::deframe()
  
  #Deplete unwanted columns
  df_core <- df[, ..required_cols]
  
  #Creatr enrichResult object
  enrich_obj <- methods::new(
    "enrichResult",
    result = df_core,
    pvalueCutoff = pCutoff,
    pAdjustMethod = pAdjustMethod,
    qvalueCutoff = qCutoff,
    gene = gene,
    universe = character(0),
    geneSets = geneSets,
    organism = organism,
    ontology = ontology,
    keytype = "SYMBOL",   
    readable = FALSE
  )
  
  return(enrich_obj)
}


#Read and Prepare Files
files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
meta <- do.call(rbind, lapply(files, extract_info))
module_list <- lapply(files, load_modules)
#names(module_list) <- basename(files)

#Extract unique names from regions
meta$region <- gsub("^(Mayo_|ROSMAP_)", "", meta$region)
regions <- unique(meta$region)

#Read modules
modules <- unlist(module_list, recursive = FALSE)
# Filter modules with at least N genes
min_genes <- 3
modules <- modules[sapply(modules, length) >= min_genes]
#Recalculate `module_list` by filtering for minimum size
module_list <- lapply(module_list, function(mods) {
  mods[sapply(mods, length) >= min_genes]
})

#Read enrichments
enrich_files <- list.files(enrich_dir, pattern = "_enrichment\\.csv$", full.names = TRUE)
enrich_data <- lapply(enrich_files, fread)

#Add enrichment metadata
enrich_data <- Map(function(df, filename) {
  fname <- basename(filename)
  region <- sub("^(Mayo_|ROSMAP_)", "", sub("_counts_.*", "", fname))
  phenotype <- ifelse(grepl("_AD_", fname), "AD", "Control")
  df$Region <- region
  df$Phenotype <- phenotype
  return(df)
}, enrich_data, enrich_files)

#Bind enrichment tables
enrich_all <- rbindlist(enrich_data)

enrich_all <- enrich_all %>%
  mutate(
    CommunityID = as.character(CommunityID),
    Module = paste(Region, Phenotype, CommunityID, sep = "_"))

################################## PART 0 ################################## 
#How similar is the overall modularity between controls and AD in each region??

#Using NMI per region

#Get NMI for all regions
nmi_results <- purrr::map_dfr(
  regions,
  ~compute_nmi_region(.x, meta_df = meta, min_genes = min_genes)
)

#Save results 
#vroom::vroom_write(nmi_results, file = file.path(output_dir, "NMI_AD_vs_Control_by_region.csv"))

print(nmi_results)

################################## PART 1 ################################## 

#FIRST: COMPARE NETWORKS IN EACH REGION PER PHENOTYPE TO FIND AD-exclusive and ctrl-exclusive modules in each region
#Comparison between phenotypes per region

#Count modules per network
module_counts <- tibble(
  Region = meta$region,
  Phenotype = meta$phenotype,
  N_modules = sapply(module_list, length),
  Network = paste0(Region, "_", Phenotype)
) %>%
  group_by(Region) %>%
  mutate(Total_region_modules = sum(N_modules)) %>%
  ungroup()

module_gene_counts <- tibble(
  Module = names(unlist(module_list, recursive = FALSE)),
  N_genes = sapply(unlist(module_list, recursive = FALSE), length)
)

#A module is considered unique to a phenotype (e.g., AD) if it has no similar module in the opposite phenotype in the same region.
#This means:
#If an X_AD module does not have any Y_Control modules with Jaccard_Index ≥ upp_thres → it is AD-exclusive.
#If a Y_Control module does not have any X_AD modules with Jaccard_Index ≥ upp_thres → it is Control-exclusive.

#Check if there is at least one row with Classification == ‘Similar’ → then it is NOT exclusive.
#If no row has ‘Similar’, then this Module_AD is AD-exclusive.

#Thresholds
upp_thres <- 0.7
low_thres <- 0.2

#Loop to calculate Jaccards per region
jaccards <- list()
for (current_region in regions) {
  #Obtain module names for this region
  mods_control <- modules[grep(paste0(current_region, "_Control"), names(modules))]
  mods_ad <- modules[grep(paste0(current_region, "_AD"), names(modules))]
  
  #Compare all modules Control vs AD
  for (ctrl_name in names(mods_control)) {
    for (ad_name in names(mods_ad)) {
      jaccard_val <- jaccard_index(mods_control[[ctrl_name]], mods_ad[[ad_name]])
      
      clasificacion <- ifelse(
        jaccard_val >= upp_thres, "Similar",
        ifelse(jaccard_val <= low_thres, "Dissimilar", "Intermediate")
      )
      
      jaccards[[length(jaccards) + 1]] <- data.frame(
        Region = current_region,
        Module_Control = ctrl_name,
        Module_AD = ad_name,
        Jaccard_Index = jaccard_val,
        Classification = clasificacion,
        stringsAsFactors = FALSE
      )
    }
  }
}

#Bind results
jaccards.tb <- rbindlist(jaccards)
#Factors in order
ordered_regions <- c("HCN", "PCC", "TC", "DLPFC","CRB")
jaccards.tb$Region <- factor(jaccards.tb$Region, levels = ordered_regions)
print(table(jaccards.tb$Classification))
#Save results <- RESULT #1
#vroom::vroom_write(jaccards.tb, file = file.path(output_dir, "jaccards_all_regions.csv"))
message("First Jaccard comparison completed and table saved!")

#Reorder modules to have a summary table 
classified_AD <- jaccards.tb %>%
  group_by(Region, Module_AD) %>%
  summarise(
    has_similar = any(Classification == "Similar"),
    has_intermediate = any(Classification == "Intermediate"),
    .groups = "drop"
  ) %>%
  mutate(
    Classification = case_when(
      has_similar ~ "Similar",
      has_intermediate ~ "Intermediate",
      TRUE ~ "AD_exclusive"
    )
  ) %>%
  dplyr::rename(Module = Module_AD) %>%
  mutate(Phenotype = "AD")

classified_CTRL <- jaccards.tb %>%
  group_by(Region, Module_Control) %>%
  summarise(
    has_similar = any(Classification == "Similar"),
    has_intermediate = any(Classification == "Intermediate"),
    .groups = "drop"
  ) %>%
  mutate(
    Classification = case_when(
      has_similar ~ "Similar",
      has_intermediate ~ "Intermediate",
      TRUE ~ "Control_exclusive"
    )
  ) %>%
   dplyr::rename(Module = Module_Control) %>%
  mutate(Phenotype = "Control")

#Bind tables to get a classification of all modules
classified_all <- bind_rows(classified_AD, classified_CTRL)

#Save table
#vroom::vroom_write(classified_all, file = file.path(output_dir, "module_classification.csv"))

#Make a summary per classification and region
summary_all <- classified_all %>%
  group_by(Region, Phenotype, Classification) %>%
  summarise(n_modules = n(), .groups = "drop") %>%
  mutate(Network = paste0(Region, "_", Phenotype)) %>% 
  dplyr::select(Network, n_modules, Classification) %>%
  left_join(module_counts, by = "Network") %>%
  mutate(local_proportion = n_modules / N_modules * 100) %>% 
  mutate(global_proportion = n_modules / Total_region_modules * 100) %>%
  mutate(across(c(local_proportion, global_proportion), round, 1))

#Save result
#vroom::vroom_write(summary_all, file = file.path(output_dir, "modules_classification_summary.csv"))

#Filter exclusive modules --- ---
candidates_exclusive_modules <- jaccards.tb %>%
  filter(Classification %in% c("Dissimilar", "Similar"))

#Save result
#vroom::vroom_write(candidates_exclusive_modules, file = file.path(output_dir, "candidates_exclusive_modules.csv"))

#Count by region and phenotype
candidates_exclusive_modules.c <- candidates_exclusive_modules %>%
  group_by(Region, Classification) %>%
  summarise(n = n(), .groups = "drop")

#Save result
#vroom::vroom_write(candidates_exclusive_modules.c, file = file.path(output_dir, "n_candidates_exclusive_modules.csv"))

#Extract AD and Control modules by region

exclusive_AD <- jaccards.tb %>%
  group_by(Region, Module_AD) %>%
  summarise(
    has_similar = any(Classification == "Similar"),
    .groups = "drop"
  ) %>%
  filter(!has_similar) %>%
  dplyr::select(Region, Module_AD) %>%
  mutate(Phenotype = "AD_exclusive")

exclusive_CTRL <- jaccards.tb %>%
  group_by(Region, Module_Control) %>%
  summarise(
    has_similar = any(Classification == "Similar"),
    .groups = "drop"
  ) %>%
  filter(!has_similar) %>%
  dplyr::select(Region, Module_Control) %>%
  mutate(Phenotype = "Control_exclusive")

#Bind results
names(exclusive_AD)[names(exclusive_AD) == "Module_AD"] <- "Module"
names(exclusive_CTRL)[names(exclusive_CTRL) == "Module_Control"] <- "Module"

#Bind
exclusive_modules <- rbind(exclusive_AD, exclusive_CTRL)

#Add number of genes per module
exclusive_modules <- exclusive_modules %>%
  left_join(module_gene_counts, by = "Module") %>%
  arrange(Region, Phenotype, desc(N_genes))

#Save results
#vroom::vroom_write(exclusive_modules, file = file.path(output_dir, "exclusive_modules_AD_CTRL.csv"))

#check

exclusive_summary <- summary_all %>% filter(Classification %in% c("AD_exclusive", "Control_exclusive"))

#Save results
#vroom::vroom_write(exclusive_summary, file = file.path(output_dir, "exclusive_summary.csv"))

#Exclusivity scores

#Exclusive AD modules with Jaccard == 0
exclusive_AD_0 <- jaccards.tb %>%
  group_by(Region, Module_AD) %>%
  summarise(
    Min_Jaccard = min(Jaccard_Index),
    .groups = "drop"
  ) %>%
  filter(Min_Jaccard == 0) %>%
  mutate(Phenotype = "AD_exclusive") %>%
  dplyr::rename(Module = Module_AD)

#Exclusive Control Modules with Jaccard == 0
exclusive_CTRL_0 <- jaccards.tb %>%
  group_by(Region, Module_Control) %>%
  summarise(
    Min_Jaccard = min(Jaccard_Index),
    .groups = "drop"
  ) %>%
  filter(Min_Jaccard == 0) %>%
  mutate(Phenotype = "Control_exclusive") %>%
  dplyr::rename(Module = Module_Control)

#Bind both
exclusive_0_all <- bind_rows(exclusive_AD_0, exclusive_CTRL_0)

#Bind with enrich_all
enrich_exclusive_0 <- enrich_all %>%
  filter(Module %in% exclusive_0_all$Module) %>%
  mutate(Region_Pheno = paste(Region, Phenotype, sep = "_"))

enrich_split_list <- enrich_exclusive_0 %>%
  split(.$Region_Pheno)

# Si quieres los términos más comunes o más significativos por grupo
enrich_summary <- enrich_exclusive_0 %>%
  filter(p.adjust < 0.05) %>%
  group_by(Phenotype, Region, Description) %>%
  summarise(
    n_modules = n_distinct(Module),
    min_p = min(p.adjust),
    example_genes = geneID[which.min(p.adjust)],
    .groups = "drop"
  ) %>%
  arrange(Phenotype, Region, min_p)

#Apply function
# ego <- create_enrichResult(enrich_split_list[["CRB_AD"]], ontology = "BP")
# 
# x<-clusterProfiler::simplify(ego)

# enrich_obj_list <- map(enrich_split_list, ~ create_enrichResult(.x, ontology = "BP"))
# 
# enrich_simplified_list <- map(enrich_obj_list, ~ simplify(
#   .x,
#   cutoff = 0.7,
#   by = "p.adjust",
#   select_fun = min,
#   measure = "Wang"
# ))

################################## PART 2 ################################## 

#Phenotypic disagreement
#Are there AD-exclusive modules that are present (conserved) in several regions?
#Are there Control-exclusive modules that also appear in multiple regions?
#You take the AD-exclusive modules from each region, for example:
#
# Mayo_CRB_AD_, Mayo_PCC_AD_3, etc.
# 
# And then you compare each pair between regions using, for example, the Jaccard index between their genes.
# 
# If two modules from different regions have a Jaccard index ≥ threshold (e.g. 0.8), you consider them ‘conserved between regions’.

#Goal: Identify exclusive modules (AD or Control) that appear in more than one region with high similarity (Jaccard_Index ≥ 0.8).

#Jaccard between exclusive modules in different regions (within each phenotype)

find_conserved_modules <- function(exclusive_df, modules_list, threshold = 0.6) {
  conserved_list <- list()
  #Filter per phenotype
  phenos <- unique(exclusive_df$Phenotype)
  for (pheno in phenos) {
    sub_df <- exclusive_df %>% filter(Phenotype == pheno)
    #Group per modules in diff regions
    for (i in 1:(nrow(sub_df) - 1)) {
      for (j in (i + 1):nrow(sub_df)) {
        mod1 <- sub_df$Module[i]
        mod2 <- sub_df$Module[j]
        
        reg1 <- sub_df$Region[i]
        reg2 <- sub_df$Region[j]
        
        #Skip if it is the same region
        if (reg1 == reg2) next
        
        genes1 <- modules[[mod1]]
        genes2 <- modules[[mod2]]
        
        jaccard_val <- jaccard_index(genes1, genes2)
        
        if (jaccard_val >= threshold) {
          conserved_list[[length(conserved_list) + 1]] <- data.frame(
            Module_1 = mod1,
            Region_1 = reg1,
            Module_2 = mod2,
            Region_2 = reg2,
            Phenotype = pheno,
            Jaccard_Index = jaccard_val,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  
  bind_rows(conserved_list)
}

#Run
conserved_modules_across_regions <- find_conserved_modules(
  exclusive_modules,
  modules,
  threshold = upp_thres
)

conserved_modules_across_regions

#Save
#vroom::vroom_write(conserved_modules_across_regions, file = file.path(output_dir, "conserved_exclusive_modules_across_regions.csv"))

#What biological functions have the conserved modules?

#what are those modules?
conserved_modules <- unique(c(
  conserved_modules_across_regions$Module_1,
  conserved_modules_across_regions$Module_2
))

#Filter the enrichments for conserved modules
enrich_conserved <- enrich_all %>%
  filter(Module %in% conserved_modules)  %>%
  filter(p.adjust < 0.05)

#Summary of conserved enriched biological functions
conserved_bio_summary <- enrich_conserved %>%
  group_by(Description, Phenotype) %>%
  summarise(
    N_modules = n_distinct(Module),
    Regions = paste(unique(Region), collapse = ", "),
    Min_pval = min(p.adjust),
    Example_genes = geneID[which.min(p.adjust)],
    .groups = "drop"
  ) %>%
  arrange(desc(N_modules), Min_pval)

conserved_bio_summary_f <- conserved_bio_summary %>% 
  filter(N_modules >=2)

enrich_conserved_f <- enrich_conserved %>% 
  filter(Description %in% conserved_bio_summary_f$Description)

#Concept plot like if we were cavernicolas
#Get edges
gene_term_edges <- enrich_conserved_f %>%
  separate_rows(geneID, sep = "/") %>%
  dplyr::select(Description, geneID, p.adjust)

#Bipartite graph
edges <- gene_term_edges %>%
  distinct(Description, geneID)

nodes <- unique(c(edges$Description, edges$geneID))
node_type <- ifelse(nodes %in% enrich_conserved_f$Description, "GO term", "Gene")

unique_terms <- enrich_conserved_f %>%
  group_by(Description) %>%
  summarise(p.adjust = min(p.adjust))
graph_df <- data.frame(name = nodes, type = node_type)

#Asign node type
graph_df <- graph_df %>%
  mutate(label_type = ifelse(type == "GO term", "Pathway", "Gene")) %>%
  left_join(unique_terms, by = c("name" = "Description")) %>%
  mutate(logp = ifelse(label_type == "Pathway", -log10(p.adjust), NA))

#Create graph
g <- graph_from_data_frame(edges, vertices = graph_df, directed = FALSE)

#Save graph
#igraph::write_graph(g, file = "cnetplot_conserved_modules.graphml", format = "graphml")

########################## PLOTTING ########################## 

#Number of modules per network
#1. Hoy many modules has each network?

diffs_modules <- module_counts %>%
  dplyr::select(Region, Phenotype, N_modules) %>%
  pivot_wider(
    names_from = Phenotype,
    values_from = N_modules
  ) %>%
  mutate(diff = AD - Control)

diffs_modules

#Times that
positives <- sum(diffs_modules$diff > 0)

#Nonsense
binom.test(positives, n = 5, alternative = "greater")
wilcox.test(diffs_modules$AD, diffs_modules$Control, paired = TRUE, alternative = "greater")

#Plot module counts
module_counts.p <- ggplot(module_counts, aes(x = Region, y = N_modules, fill = Phenotype)) +
  geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.3) + 
  facet_wrap(~Phenotype) +
  theme_cowplot() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("AD" = "firebrick", "Control" = "navyblue")) +
  scale_y_continuous(breaks = c(0:max(module_counts$N_modules))) +
  labs(
    title = "",
    y = "N modules",
    x = ""
  )

#Vis
module_counts.p

#2.How similar are the modules between AD and control within each region?

nmi_mat <- nmi_results %>%
  dplyr::select(Region, NMI) %>%
  arrange(NMI) %>%
  mutate(Region = factor(Region, levels = unique(Region)))

mini_nmiheat <- ggplot(nmi_mat, aes(x = "AD vs Control", 
                                   y = Region,
                                   fill = NMI)) +
 geom_tile(color = "white") +
 geom_text(aes(label = round(NMI, 3)), color = "black", size = 4) +
  scale_fill_gradient(low = "cornflowerblue", high = "firebrick") +
  labs(title = "NMI",
    x = "",
    y = "") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(), 
    legend.position="none")

#Vis
mini_nmiheat

#What about the means and medias

jaccard_summary <- jaccards.tb %>%
  dplyr::group_by(Region) %>%
  dplyr::summarise(
    mean_jaccard = mean(Jaccard_Index, na.rm = TRUE),
    median_jaccard = median(Jaccard_Index, na.rm = TRUE),
    n_pairs = dplyr::n()
  ) %>%
  arrange(mean_jaccard)%>%
  mutate(Region = factor(Region, levels = unique(Region)))

#Plot
mini_meanjacc_heat <- ggplot(jaccard_summary, aes(x = "", y = Region, fill = mean_jaccard)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(mean_jaccard, 3)),  color = "black", size = 4) +
  scale_fill_gradient(low = "cornflowerblue", high = "firebrick") +
  theme_minimal(base_size = 14) +
  labs(title = "Jaccard mean",
       x = "",
       y = "") +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(), 
    legend.position="none")

mini_meanjacc_heat

#Plot
mini_median_jacc_heat <- ggplot(jaccard_summary, aes(x = "", y = Region, fill = median_jaccard)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(median_jaccard, 3)),  color = "black", size = 4) +
  scale_fill_gradient(low = "cornflowerblue", high = "firebrick") +
  theme_minimal(base_size = 14) +
  labs(title = "Jaccard median",
       x = "",
       y = "") +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(), 
    legend.position="none")

mini_median_jacc_heat#lol

#Plot panel

miniheats <- cowplot::plot_grid(mini_nmiheat, mini_meanjacc_heat,nrow = 1)
miniheats

# #Save histograms
# ggsave(filename = file.path(output_dir, "miniheats.pdf"),
#        plot = miniheats,
#        width = 4,
#        height = 7)

#Histogram of Jaccard values

jaccard_hist.p <- ggplot(jaccards.tb, aes(x = Jaccard_Index)) +
  geom_histogram(binwidth = 0.05, fill = "navyblue", color = "black", alpha = 0.7) +
  #geom_density(aes(y = ..count.. * 0.005), adjust = 10,  na.rm = TRUE, color = "orange", size = 1) +
  geom_vline(xintercept = low_thres, linetype = "dashed", color = "red", size = 1) +
  geom_vline(xintercept = upp_thres, linetype = "dashed", color = "darkgreen", size = 1) +
  scale_y_continuous(trans = "log10", limits = c(1, NA)) +  # evita problemas con ceros
  labs(
    title = "Global",
    x = "Jaccard Index",
    y = "Frequency (log10)"
  ) +
  #facet_wrap(~Region) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    axis.text = element_text(size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

jaccard_hist.p

#Plot per region (then grid them)
jaccard_region.p <- lapply(levels(jaccards.tb$Region), function(reg) {
  ggplot(jaccards.tb %>% filter(Region == reg), aes(x = Jaccard_Index)) +
    geom_histogram(binwidth = 0.04, fill = "navyblue", color = "black", alpha = 0.7) +
    #geom_density(aes(y = ..count.. * 0.005), color = "orange", size = 1) +
    geom_vline(xintercept = low_thres, linetype = "dashed", color = "red", size = 1) +
    geom_vline(xintercept = upp_thres, linetype = "dashed", color = "darkgreen", size = 1) +
    scale_y_continuous(trans = "log10", limits = c(1, NA)) +
    labs(
      title = reg,
      x = "Jaccard Index",
      y = "Frequency (log10)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )
})

#Panel
jaccard.p <- plot_grid(jaccard_hist.p, 
                       jaccard_region.p[[1]], 
                       jaccard_region.p[[2]],
                       jaccard_region.p[[3]],
                       jaccard_region.p[[4]],
                       jaccard_region.p[[5]],
                       ncol = 3)
jaccard.p

#Save histograms
# ggsave(filename = file.path(output_dir, "jaccard_histogram.pdf"),
#        plot = jaccard.p,
#        width = 18,
#        height = 10)

#Visualize local proportion

#Get cool colors
fill_colors <- c(
  "Control_exclusive" = rgb(0.1, 0.4, 0.3, 0.6),   #dark green
  "Intermediate"      = rgb(0.4, 0.6, 0.3, 0.6),   #olive green
  "Similar"           = rgb(0.4, 0.3, 0.6, 0.6),   #soft violeta
  "AD_exclusive"      = rgb(0.6, 0.2, 0.5, 0.6)    #magenta
)

#Plot stacked barplot per local proportion

local.p <- ggplot(summary_all, aes(x = Region, y = local_proportion, fill = Classification)) +
  geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.3) +
  facet_wrap(~Phenotype) +
  scale_fill_manual(values = fill_colors) +
  labs(
    title = "Local",
    y = "%",
    x = ""
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

#Vis 
local.p

#Visualize global proportion

#Plot proportions per region (global proportion)
global.p <- ggplot(summary_all, aes(x = Region, y = global_proportion, fill = Classification)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(
    values = fill_colors,
    labels = c(
      "AD exclusive",
      "Control exclusive",
      "Intermediate",
      "Similar"
    )
  ) +
  labs(
    title = "Global",
    x = "",
    y = "%"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank()
  )

#Vis
global.p 

#Grid plot
proportions <- plot_grid(local.p, global.p, labels = c("A", "B"))
proportions

#Save histograms
# ggsave(filename = file.path(output_dir, "proportions.pdf"),
#        plot = proportions,
#        width = 13,
#        height = 7)

#Barplot
barplot.p <- ggplot(candidates_exclusive_modules.c,
                    aes(x = Region, y = n, fill = Classification)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = paste0("n: ", n)),
            position = position_dodge(width = 0.8),
            vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = fill_colors) +
  labs(title = "Exclusive Modules by Region and Phenotype",
       x = "Brain Region",
       y = "Number of Exclusive Modules",
       fill = "Classification") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title.x = element_text(size = 13, face = "bold"),
    axis.title.y = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(angle = 0, size = 11, margin = margin(t = -5)),
    axis.text.y = element_text(size = 11),
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.title = element_text(face = "bold"),
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

barplot.p

#What 

