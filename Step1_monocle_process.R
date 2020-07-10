#load and process data with monocle

source("Step0_setup.R")

seu = readRDS(seurat_obj_rds)

#if seurat doesn't have cell cycle info
# seu = CellCycleScoring(seu, 
# s.features = gene_name_hs2mm(Seurat::cc.genes.updated.2019$s.genes), 
# g2m.features = gene_name_hs2mm(Seurat::cc.genes.updated.2019$g2m.genes))

message("looking for loaded monocle data at:\n", monocle_obj_rds)

if(!file.exists(monocle_obj_rds)){
    if(length(data_dirs) < 1) stop("No data directories located. Check `data_path` and `pattern`.")
    monocle_loaded = lapply(data_dirs, function(d){
        #load data
        mon = monocle3::load_cellranger_data(pipestance_path = d, genome = GEN)
        #creates cell barcodes compatible with Seurat
        prefix = dir2prefix[d]
        colnames(mon) = paste0(prefix, "_", sub("-.+", "", colnames(mon)))
        mon
    })
    total_loaded = sum(sapply(monocle_loaded, ncol))
    #combine all monocle objects
    if(length(monocle_loaded) > 1){
        mon = cbind(monocle_loaded[[1]], monocle_loaded[[2]])    
    }else{
        mon = monocle_loaded
    }
    if(length(monocle_loaded) > 2){
        for(i in seq(3, length(monocle_loaded))){
            mon = cbind(mon, monocle_loaded[[i]])
        }
    }
    #check we combined everything
    stopifnot(ncol(mon) == total_loaded)
    
    saveRDS(mon, file = monocle_obj_rds)
}else{
    mon = readRDS(monocle_obj_rds)
}

#see how seurat and monocle match
seqsetvis::ssvFeatureVenn(list(monocle = colnames(mon), seurat = colnames(seu))) +
    labs(title = "intersection between seurat and monocle barcodes")
n_before = ncol(mon)
mon = mon[, intersect(colnames(mon), colnames(seu))]
n_after = ncol(mon)

message(n_before - n_after , " cells removed based on Seurat object. ", n_after, " cells in processed Monocle object.")

message("looking for loaded monocle data at:\n", default_monocle_processed_rds)

run_process = TRUE
if(file.exists(default_monocle_processed_rds)){
    run_process = shiny_overwrite()
    message(run_process)
}
if(run_process){
    #if you want to bring over cluster or other info from Seurat, do it here
    pData(mon)$seurat_names = get_cluster_rename()[seu$seurat_clusters[rownames(pData(mon))]]
    pData(mon)$seurat_clusters = seu$seurat_clusters[rownames(pData(mon))]
    pData(mon)$genotype = seu$treatment[rownames(pData(mon))]
    pData(mon)$rep = seu$rep[rownames(pData(mon))]
    pData(mon)$Phase = seu$Phase[rownames(pData(mon))]
    pData(mon)$G2M.Score = seu$G2M.Score[rownames(pData(mon))]
    pData(mon)$S.Score = seu$S.Score[rownames(pData(mon))]
    mon = mon[, colnames(seu)]
    table(pData(mon)$seurat_clusters)
    
    
    mon <- preprocess_cds(mon, num_dim = num_pca)
    if(align_apply){
        mon <- align_cds(mon, alignment_group = "genotype")    
    }
    mon = reduce_dimension(mon, reduction_method = "UMAP", cores = 20)
    mon = cluster_cells(mon)
    #minimal_branch_len is quite important and controls the detail of branches
    mon <- learn_graph(mon, learn_graph_control = list(minimal_branch_len = minimal_branch_len), use_partition = TRUE)
    plot_cells(mon)
    saveRDS(mon, default_monocle_processed_rds)
    message('saved ', default_monocle_processed_rds)
}else{
    mon = readRDS(default_monocle_processed_rds)
    message('loaded ', default_monocle_processed_rds)
}

processed_plots(mon, default_proc_dir)

# p = seurat_cluster_choose_plot(mon, show.legend = TRUE, col_scale = get_clusters_colors(), color_by = "seurat_names", facet_by = "genotype") + labs(title = "Seurat clusters")
# xp = 2
# ggsave("seurat_genotypes.png", p+facet_grid(seurat_names~genotype), width = 4*xp, height = 16*xp)
