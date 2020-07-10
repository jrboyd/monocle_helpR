#select specific cells and also root nodes.
#should be run from monocle directory
source("Step0_setup.R")

# stopifnot(file.exists(monocle_processed_rds))
# mon = readRDS(monocle_processed_rds)

#patch to bring in cell clusters
stopifnot(file.exists(seurat_obj_rds))
seu = readRDS(seurat_obj_rds)
# if(is.null(pData(mon)$Phase)) pData(mon)$Phase = seu$Phase[rownames(pData(mon))]
# if(is.null(pData(mon)$S.Score)) pData(mon)$S.Score = seu$S.Score[rownames(pData(mon))]
# if(is.null(pData(mon)$G2M.Score)) pData(mon)$G2M.Score = seu$G2M.Score[rownames(pData(mon))]

select_dir = shiny_choose_branch(output_path = output_root_dir)
Sys.sleep(1)
select_name = sub(paste0(prefix_subset, "."), "", basename(select_dir))
stopifnot(nchar(select_dir) > 0)
dir.create(select_dir, showWarnings = FALSE)
monocle_selected_rds = file.path(select_dir, file_subset)

monocle_processed_rds = file.path(dirname(select_dir), file_processed)
mon = readRDS(monocle_processed_rds)

run_select = TRUE
#if selection rds exists, should it be overwritten?
if(file.exists(monocle_selected_rds)){
    run_select = shiny_overwrite(msg = "Selection data found, overwrite?")
    Sys.sleep(1)
    message(run_select)
}

#this is where the analysis really starts
if(run_select){
    message("Run selection")
    #select a region of interest (or assign sel_branch = mon to analyze everything)
    if(select_name != "full"){
        FUN = function(mon, not_used){
            seurat_cluster_choose_plot.alpha(mon, color_by = "seurat_names", show.legend = TRUE, col_scale = get_clusters_colors())
        }
        sel_branch = my_choose_cells(mon, plot_FUN = FUN)    
    }else{
        sel_branch = mon
    }
    
    #select root node
    sel_branch <- order_cells(sel_branch)
    #if you like your results, rename this file so you can repeat them.
    saveRDS(sel_branch, file = monocle_selected_rds)
    
}else{
    message("Load selection")
    sel_branch = readRDS(monocle_selected_rds)
}

# plot_cells(sel_branch, color_cells_by = "pseudotime")
barplot_cnt_per_seurat(mon, sel_branch) +
    facet_wrap(~genotype)
