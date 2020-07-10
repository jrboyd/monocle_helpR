#should be run from monocle directory
stopifnot(basename(getwd()) == "monocle")
rm(list = ls())

library(monocle3)
library(Seurat)
library(ggplot2)
library(data.table)
library(magrittr)
if(!require(ssvRecipes)){
    devtools::install_github("jrboyd/ssvRecipes")
}
source("functions_monocle.R")
source("functions_steps.R")
source("shiny_gadgets.R")

seurat_obj_rds = "../Bcell_combined.refixed_120919.monocle_ready.rds"
stopifnot(file.exists(seurat_obj_rds))
#this monocle object will be created if needed


#if you want to start fresh, alter this path or move existing
output_root_dir = "output"
dir.create(output_root_dir, showWarnings = FALSE, recursive = TRUE)



#don't mess with this please
prefix_processed = "processed"
prefix_subset = "subset_branch"
file_loaded = "monocle_loaded.Rds"
file_processed = "monocle_processed.Rds"
file_subset = "monocle_selected.Rds"

monocle_obj_rds = file.path(output_root_dir, file_loaded)

## Step 1 parameters
#Loading
#where to load raw data from, only used once if loaded data not found in Step1
GEN = "mm10"
data_path = "~/../joeboyd/BCELL_scRNAseq/"
pattern = "monocle$"
data_dirs = dir(data_path, pattern = pattern, full.names = TRUE)
data_dirs = normalizePath(data_dirs)
data_dirs
#dir2prefix should translate directory paths to prefixes that Seurat uses for cell barcodes.
dir2prefix = basename(data_dirs); names(dir2prefix) = data_dirs
# n = 3
# dir2prefix = (stringi::stri_split_fixed(dir2prefix, "_", n = n) %>% unlist %>% matrix(ncol = n, byrow = TRUE))[,2]
dir2prefix = c("df", 'wt')
names(dir2prefix) = data_dirs

#Processing
#monocle has an alignment procedure analogous to integration in Seurat
align_apply = FALSE
#number of PC to use for reduction
#I have it at 5 here because beyond that PCs explain less that 5% of the variance
num_pca = 10
#increase to remove small branches and decrease to see more
minimal_branch_len = 8

proc_desc = paste0(ifelse(align_apply, "align", "unalign"), ".", 
                   num_pca, "_PC.",
                   minimal_branch_len, "_branchLen")

get_proc_dir = function(name, out_dir = output_root_dir, prefix = prefix_processed, suffix = proc_desc){
    file.path(out_dir,
              paste0(
                  prefix,
                  ".",
                  name,
                  ".", 
                  suffix)
    )
}

default_proc_dir = get_proc_dir("full")

dir.create(default_proc_dir, showWarnings = FALSE, recursive = TRUE)

default_monocle_processed_rds = file.path(default_proc_dir,
                                          paste0(file_processed))

cluster_identities = read.table("cluster_idents.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")
signature_genes = strsplit(cluster_identities$signature_genes, ",")

rename_clust = cluster_identities$bio_name
names(rename_clust) = cluster_identities$seurat_clusters
rename_clust = factor(rename_clust, levels = unique(rename_clust))

tmp.col = unique(cluster_identities[, c("bio_name", "color")])
color_clust = unique(tmp.col$color)
names(color_clust) = unique(tmp.col$bio_name)
stopifnot(!any(duplicated(names(color_clust))))
stopifnot(!any(duplicated(color_clust)))

color_genotype = c( "wt" = "gray30", "df4" = "red")

#' should return a named vector of cluster names
get_cluster_rename = function(){
    rename_clust
}


#' should return a named vector of valid R colors, so "red" or "#FF0000".
#' names should be all possible seurat_names in monocle object
get_clusters_colors = function(){
    color_clust
}

get_treatments_colors = function(){
    color_genotype
}

## Step 3
min_morans = .2 #Feel free to explore different values of Morans test statistic. I find values below .1 pretty noisy looking
max_q = .05