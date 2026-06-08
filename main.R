library(readxl)
library(GEOquery)
library(limma)
library(WGCNA)
library(dplyr)
library(readr)
library(fastcluster)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(cowplot)

options(stringsAsFactors = FALSE, warn = -1)
allowWGCNAThreads()

#=====================================================================================
# Command line arguments
#=====================================================================================
args <- commandArgs(trailingOnly = TRUE)
dataset_arg <- grep("^--dataset=", args, value = TRUE)
if (length(dataset_arg) != 1) {
  stop(
    "Missing dataset argument.\n",
    "Usage: Rscript main.R --dataset=GSE159090"
  )
}

selected_dataset <- sub("^--dataset=", "", dataset_arg)

#=====================================================================================
# Project directories
# Assumption: the script is launched from the repository root
#=====================================================================================
workingDir_base <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
config_file <- file.path(workingDir_base, "dataset_config.csv")
plot_utils  <- file.path(workingDir_base, "plot_utility.R")
if (!file.exists(config_file)) {
  stop("Missing dataset_config.csv in repository root: ", workingDir_base)
}
if (!file.exists(plot_utils)) {
  stop("Missing plot_utility.R in repository root: ", workingDir_base)
}
source(plot_utils)

#===================================================================================
# Run-specific output directory
#===================================================================================
OUT_DIR_DAY <- format(Sys.Date(), "%Y_%m_%d")
config_df <- read.csv(
  config_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  colClasses = c(original_gsms = "character")
)
cfg <- config_df %>%
  filter(dataset_id == selected_dataset)
if (nrow(cfg) != 1) {
  stop("Configuration error: expected exactly one row for dataset_id = ", selected_dataset)
}
dataset_subdir <- cfg$data_subdir[1]
dataset_in     <- cfg$dataset_in[1]
trait_in       <- cfg$trait_in[1]
original_gsms  <- as.character(cfg$original_gsms[1])
min_group_size <- 30
max_block_size <- 100000
DATASET_DIR <- file.path(workingDir_base, dataset_subdir)
OUT_DATASET_DIR <- file.path(workingDir_base, OUT_DIR_DAY, dataset_subdir)
dir.create(OUT_DATASET_DIR, showWarnings = FALSE, recursive = TRUE)
dataset_in <- file.path(DATASET_DIR, dataset_in)
trait_in   <- file.path(DATASET_DIR, trait_in)
if (!file.exists(dataset_in)) {
  stop("Dataset file not found: ", dataset_in)
}
if (!file.exists(trait_in)) {
  stop("Trait file not found: ", trait_in)
}
setwd(OUT_DATASET_DIR)
message("Selected dataset: ", selected_dataset)
message("Working directory: ", workingDir_base)
message("Output directory: ", OUT_DATASET_DIR)

#=====================================================================================
# Global hyperparameters
#=====================================================================================
MAD                     <- TRUE
MAD_PERCENTAGE_DISCARD  <- 0.5
ALPHA                   <- 0.05
ALPHA_LIMMA             <- 0.05
GS_threshold_init       <- 0.9
MM_threshold_init       <- 0.9
options(stringsAsFactors = FALSE)

#=====================================================================================
# Function: run full pipeline with the first k_remove samples dropped
#=====================================================================================
run_leave_k_out_pipeline <- function(dataset_in, trait_in, gsms_sub, keep_idx, k= NULL, iter = NULL) {
  # load expression + trait tables 
  expr0   <- read.csv(dataset_in, row.names = 1, check.names = FALSE)
  traits0 <- read.csv(trait_in,   row.names = 1, check.names = FALSE)
  total   <- ncol(expr0)
  
  # prepare output folder
  ts     <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  removed_idx <- setdiff(seq_len(total), keep_idx)
  if (!is.null(k) && !is.null(iter)) {
    folder <- paste0(ts, "_k", k, "_iter", iter, "_removed_", paste(removed_idx, collapse = "_"))
  } else {
    folder <- paste0(ts, "_removed_", paste(removed_idx, collapse = "_"))
  }
  dir.create(folder)
  
  expr     <- expr0[, keep_idx, drop = FALSE]
  chars    <- strsplit(gsms_sub, "")[[1]]
  classes  <- as.numeric(ifelse(chars == "0", 1, 2)) # 1 = control, 2 = case
  
  # subset traits to match expr
  if (nrow(traits0) == total) {
    traits <- traits0[keep_idx, , drop = FALSE]
  } else if (ncol(traits0) == total) {
    traits <- traits0[, keep_idx, drop = FALSE]
  } else {
    traits <- data.frame(Group = classes,
                         row.names = colnames(expr))
  }
  
  message("=== leave-k-out: dropped ", length(setdiff(seq_len(total), keep_idx)),
          " samples → kept ", ncol(expr),
          " | GSMS = ", gsms_sub)
  
  #===================================================================================
  # Cleaning if needed
  #===================================================================================
  datExpr0 <- as.data.frame(t(expr))
  datExpr  <- t(datExpr0)[complete.cases(t(datExpr0)), ]
  datExpr  <- t(datExpr)
  dim(datExpr)
  gsg <- goodSamplesGenes(datExpr0, verbose = 3)
  if (!gsg$allOK) stop("Some samples/genes flagged by goodSamplesGenes()")
  
  #===================================================================================
  # PCA
  #===================================================================================
  pca_result <- prcomp(datExpr, scale. = TRUE) 
  pca_scores <- pca_result$x[, 1:2]
  labels <- character(length(classes)) 
  counters <- list()
  for (i in seq_along(classes)) {
    value <- classes[i]
    if (!is.null(counters[[as.character(value)]])) {
      counters[[as.character(value)]] <- counters[[as.character(value)]] + 1
    } else {
      counters[[as.character(value)]] <- 1
    }
    labels[i] <- counters[[as.character(value)]]
  }

  save_PCA(pca_scores = pca_scores,  classes = classes, labels = labels, output_file = paste0(folder, "/PCA.pdf"))
  eigenvalues_plot <- pca_result$sdev^2
  save_variance_plot(eigenvalues_plot, output_file = paste0(folder, "/PCA_variance_plot.pdf"))
  
  #===================================================================================
  # Sample Clustering
  #===================================================================================
  sampleTree = fastcluster::hclust(dist(datExpr), method = "average");
  save_hierarchical_clustering(sample_tree = sampleTree,  output_file = paste0(folder, "/sampleClustering.pdf"))
  
  #===================================================================================
  # LIMMA, differential expression analysis
  #===================================================================================
  gset = ExpressionSet(data.matrix(t(datExpr)))
  sml <- strsplit(gsms_sub, split="")[[1]]
  
  # assign samples to groups and set up design matrix
  gs <- factor(sml)
  groups <- make.names(c("chow","hfd"))
  levels(gs) <- groups
  gset$group <- gs
  design <- model.matrix(~group + 0, gset)
  colnames(design) <- levels(gs)
  gset <- gset[complete.cases(exprs(gset)), ] # skip missing values
  fit <- lmFit(gset, design)  # fit linear model
  # set up contrasts of interest and recalculate model coefficients
  cts <- c(paste(groups[1],"-",groups[2],sep=""))
  cont.matrix <- makeContrasts(contrasts=cts, levels=design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  # compute statistics and table of top significant genes
  fit2 <- eBayes(fit2, 0.01)
  tT <- topTable(fit2, adjust="fdr", sort.by="B", number=1000000)
  write.table(tT, file=paste(folder,'/apoe_out.tsv', sep = ''), sep="\t")
  
  #===================================================================================
  # parse LIMMA results and save DEGs
  #===================================================================================
  read_genes <- function(out_table, ALPHA_LIMMA, TH_FC) {
    # Filter the dataframe based on adjusted p-value and log fold change
    df_limma <- out_table %>%
      filter(adj.P.Val < ALPHA_LIMMA, abs(logFC) > TH_FC)
    genes = rownames(df_limma)
    # Create a list of unique Entrez IDs (removing NAs)
    list_limma <- unique(as.integer(genes)) # EntrezIDs are integers
    list_limma = na.omit(list_limma)
    return(list_limma)
  }
  
  FC_values <- c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5)
  all_DEGs   <- vector("list", length(FC_values))
  names(all_DEGs) <- as.character(FC_values)
  for (TH_FC in FC_values) {
    save_volcano_plot(tT, ALPHA_LIMMA, TH_FC, paste(folder, "/volcano", paste(ALPHA_LIMMA*100, TH_FC*100, sep = '_'), ".pdf", sep=''), width = 6, height = 4)
    limma_genes <- read_genes(tT, ALPHA_LIMMA, TH_FC)
    writeLines(as.character(limma_genes), con = paste(folder, '/deg', paste(ALPHA_LIMMA*100, TH_FC*100, sep = '_'), '.txt', sep = ''))
    # Store in our “all_DEGs” list
    all_DEGs[[as.character(TH_FC)]] <- limma_genes
  }
  
  #===================================================================================
  # WGCNA: selection of genes by MAD
  #===================================================================================
  if (MAD) {
    mad_values <- apply(datExpr, 2, mad)
    threshold <- quantile(mad_values, MAD_PERCENTAGE_DISCARD)
    datExpr <- datExpr[, mad_values >= threshold]
    dim(datExpr)
  }
  
  #===================================================================================
  # WGCNA: soft thresholding & plotting
  #===================================================================================
  networkType <- "unsigned" 
  powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
  sft = pickSoftThreshold(datExpr, 
                          powerVector = powers,
                          networkType = networkType,
                          verbose = 5, 
                          blockSize = max_block_size)
  
  save_sft_plots(sft = sft, powers = powers, output_file = paste0(folder, "/sft_plot", MAD, MAD*MAD_PERCENTAGE_DISCARD, ".pdf"))
  dim(datExpr)
  powerestimate <- sft$powerEstimate
  if (is.na(powerestimate) || powerestimate < 1 || powerestimate > 30) powerestimate <- 20

  #===================================================================================
  # WGCNA: network construction
  #===================================================================================
  net = blockwiseModules(datExpr, power = powerestimate,
                         maxBlockSize = max_block_size, 
                         networkType = networkType,
                         TOMType = networkType, minModuleSize = min_group_size,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE, pamRespectsDendro = TRUE,
                         saveTOMs = FALSE,
                         saveTOMFileBase = "APO-TOM", 
                         verbose = 3)
  
  table(net$colors)
  moduleColors = labels2colors(net$colors)
  mergedColors = labels2colors(net$colors)
  write.csv(table(mergedColors), paste(paste(folder, "/clusters", sep=''), MAD, MAD*MAD_PERCENTAGE_DISCARD, ".csv", sep = '_'), row.names = FALSE)
  save_dendrogram_plot(net = net, output_file = paste0(folder, "/dendrogram_plot", MAD, MAD*MAD_PERCENTAGE_DISCARD, ".pdf"))

  #===================================================================================
  # WGCNA: module eigengenes analysis
  #===================================================================================
  MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)

  #===================================================================================
  # WGCNA: correlation of modules with traits
  #===================================================================================
  nSamples <- nrow(datExpr)
  datTraits <- traits
  # if traits came in with samples as columns rather than rows, transpose
  if (ncol(datTraits) == nSamples && all(colnames(datTraits) == rownames(datExpr))) {
    datTraits <- as.data.frame(t(datTraits))
  }
  # ensure rownames align exactly
  if (!all(rownames(datTraits) == rownames(datExpr))) {
    datTraits <- datTraits[rownames(datExpr), , drop = FALSE]
  }
  stopifnot(nrow(datTraits) == nSamples)
  # compute correlations and p‐values
  moduleTraitCor    <- cor(MEs,  datTraits, use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
  # build the text matrix
  textMatrix <- paste0(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")")
  dim(textMatrix) <- dim(moduleTraitCor)
  # save the heatmap
  save_heatmap_plot(matrix=moduleTraitCor, xLabels=names(datTraits), yLabels=names(MEs), ySymbols=names(MEs), textMatrix=textMatrix, output_file=file.path(folder, "heatmap_plot.pdf")
  )
  
  #===================================================================================
  # WGCNA: gene significance & module membership
  #===================================================================================
  group = as.data.frame(datTraits$Group)
  names(group) = "group"
  modNames = substring(names(MEs), 3)
  geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
  names(geneModuleMembership) = paste("MM", modNames, sep="")
  names(MMPvalue) = paste("p.MM", modNames, sep="")
  geneTraitSignificance = as.data.frame(cor(datExpr, group, use = "p"))
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
  names(geneTraitSignificance) = paste("GS.", names(group), sep="")
  names(GSPvalue) = paste("p.GS.", names(group), sep="")
  
  #===================================================================================
  # WGCNA: hub gene selection
  #===================================================================================
  wgcna_all <- list()
  for (FRACTION in c(0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95)) {
    GS_threshold <- FRACTION
    MM_threshold <- FRACTION
    df_traitp <- data.frame(moduleTraitPvalue)
    significant_modules <- rownames(df_traitp)[df_traitp$Group < ALPHA]
    hub_genes_list <- list()
    significant_gene_counts <- data.frame(Module = character(), Square = integer(), Circle = integer(), Ellipse = integer(), stringsAsFactors = FALSE)
    
    for (MEmodule in significant_modules) {
      module = sub("^ME", "", MEmodule)
      column = match(module, modNames)
      moduleGenes = moduleColors == module
      abs_mm = abs(geneModuleMembership[[paste0("MM", module)]])
      abs_gs = abs(geneTraitSignificance[[paste0("GS.group")]])
      gene_names = colnames(datExpr)
      module_filter = moduleColors == module
      
      # ---- SELECTIONS ----
      # Square
      square_filter = gene_names[module_filter & (abs_mm >= MM_threshold) & (abs_gs >= GS_threshold)]
      # Circle
      squared_distance = (abs_mm - 1)^2 + (abs_gs - 1)^2
      r_squared = 2 * (1 - FRACTION)^2
      circle_filter = gene_names[module_filter & (squared_distance <= r_squared)]
      # Ellipse (main)
      a = 1 - MM_threshold  # semi-axis for MM (x)
      b = 1 - GS_threshold  # semi-axis for GS (y)
      ellipse_distance = ((abs_mm - 1)^2) / (a^2) + ((abs_gs - 1)^2) / (b^2)
      ellipse_filter = gene_names[module_filter & ellipse_distance <= 1]
      
      # Print comparison
      cat(sprintf("Module: %s | Square: %d | Circle: %d | Ellipse: %d genes\n", module, length(square_filter), length(circle_filter), length(ellipse_filter)))
      significant_gene_counts <- rbind(significant_gene_counts, data.frame(Module = module, Square = length(square_filter), Circle = length(circle_filter), Ellipse = length(ellipse_filter)))
      
      # ---- SAVE PLOT ----
      save_scatter_plot_with_selections(
        geneModuleMembership, geneTraitSignificance, moduleGenes, column, module,
        paste(paste(folder,"/cor_compare",sep=''), paste(GS_threshold*100,MM_threshold*100,sep = '_'), module,".pdf",sep='_'),
        MM_threshold, GS_threshold, r_squared, a, b
      )
      # ---- SELECT GENES FOR OUTPUT (Ellipse by default) ----
      hub_module = datExpr[, ellipse_filter, drop = FALSE]
      write.csv(hub_module, paste(paste(folder,"/hub_genes",sep=''), module, paste(GS_threshold*100,MM_threshold*100,sep = '_'), ".csv", sep='_'), row.names = FALSE)
      hub_genes_list <- append(hub_genes_list, colnames(hub_module))
    }
    # Save gene counts comparison
    write.csv(significant_gene_counts, paste0(folder, "/selection_counts_", FRACTION, ".csv"), row.names = FALSE)
    
    #=================================================================================
    # Save merged module‐trait + counts
    #=================================================================================
    moduleTraitCor <- data.frame(moduleTraitCor)
    moduleTraitPvalue_new <- data.frame(moduleTraitPvalue)
    moduleTraitCor$Module <- rownames(moduleTraitCor)
    moduleTraitPvalue_new$Module <- rownames(moduleTraitPvalue_new)
    pval_columns <- setdiff(colnames(moduleTraitPvalue_new), "Module")
    colnames(moduleTraitPvalue_new)[colnames(moduleTraitPvalue_new) %in% pval_columns] <- paste0("pval_", pval_columns)
    merged_df <- merge(moduleTraitCor, moduleTraitPvalue_new, by = "Module")
    rownames(merged_df) <- merged_df$Module
    merged_df$Module <- NULL
    module_counts_df <- as.data.frame(table(mergedColors))
    colnames(module_counts_df) <- c("Module", "GeneCount")
    merged_df$Module <- sub("ME", "", rownames(merged_df))
    merged_df <- merge(merged_df, module_counts_df, by = "Module", all.x = TRUE)
    merged_df <- merge(merged_df, significant_gene_counts, by = "Module", all.x = TRUE)
    write.csv(merged_df,paste(paste(folder, "/merged_df", sep=''),paste(GS_threshold*100, MM_threshold*100, sep = '_'),".csv", sep = '_'),row.names = FALSE)
    list_wgcna <- unique(as.integer(unlist(hub_genes_list)))
    write.csv(list_wgcna, paste(paste(folder, "/wgcna_genes", sep=''), paste(GS_threshold*100, MM_threshold*100, sep = '_'), ".csv", sep = '_'), row.names = FALSE)
    wgcna_all[[as.character(FRACTION * 100)]] <- list_wgcna
  }
  
  return(list(kept=ncol(expr), folder=folder, gene_lists=list(DEA_all_FC=all_DEGs, WGCNA_all_frac=wgcna_all)))
}

#=================================================================================
# Main: iterate k_remove = 0,1,…,n-1
#=================================================================================
chars        <- strsplit(original_gsms, "")[[1]]
n_controls   <- sum(chars == "0")
n_cases      <- sum(chars == "1")
max_k        <- floor(min(n_controls, n_cases)/2)
n_iter       <- 10      # Number of subsampling iterations
N            <- length(chars)
control_idx  <- which(chars == "0")
case_idx     <- which(chars == "1")
subsampling_results <- vector("list", max_k + 1) # +1 for k = 0

# --- No-removal baseline (k = 0) ---
keep_idx_all <- seq_len(N)
gsms_sub_all <- paste0(chars[keep_idx_all], collapse = "")
subsampling_results[[1]] <- run_leave_k_out_pipeline(dataset_in, trait_in, gsms_sub_all, keep_idx_all, k = 0, iter = 1)

# --- Random subsampling for k = 1 to max_k ---
set.seed(12345)
for (k in 1:max_k) {
  subsampling_results[[k + 1]] <- vector("list", n_iter)
  for (iter in 1:n_iter) {
    remove_controls <- sample(control_idx, k, replace = FALSE)
    remove_cases    <- sample(case_idx, k, replace = FALSE)
    remove_idx      <- sort(c(remove_controls, remove_cases))
    keep_idx        <- setdiff(seq_len(N), remove_idx)
    gsms_sub        <- paste0(chars[keep_idx], collapse = "")
    subsampling_results[[k + 1]][[iter]] <- run_leave_k_out_pipeline(dataset_in, trait_in, gsms_sub, keep_idx, k, iter)
  }
}