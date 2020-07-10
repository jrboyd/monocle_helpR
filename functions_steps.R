run_step3 = function(monocle_selected_rds, max_q, min_morans, force = FALSE){
    stopifnot(file.exists(monocle_selected_rds))
    sel_branch = readRDS(monocle_selected_rds)
    #change this for every analysis
    analyze_out_dir = dirname(monocle_selected_rds)
    stopifnot(dir.exists(analyze_out_dir))
    
    res_file = function(f)file.path(analyze_out_dir, f)
    
    heatmap_gl_file = res_file("module_heatmap_geneLists.csv")
    if(file.exists(heatmap_gl_file) & !force){
        stop("previous anlaysis found, use force = TRUE to overwrite.")
    }
    
    monocle_processed_rds = file.path(monocle_selected_rds %>% dirname %>% dirname, file_processed)
    stopifnot(file.exists(monocle_processed_rds))
    
    mon = readRDS(monocle_processed_rds)
    meta_dt = as.data.table(pData(mon))
    meta_dt$id = colnames(mon)
    meta_dt[, in_selection := id %in% colnames(sel_branch)]
    
    ggplot(meta_dt[, .(.N), .(genotype, in_selection, seurat_names)], 
           aes(x = seurat_names, y = N, fill = in_selection)) +
        geom_bar(stat = "identity") +
        facet_wrap(~genotype)
    
    ggplot(meta_dt[, .(.N), .(genotype, in_selection, seurat_names)], 
           aes(x = genotype, y = N, fill = in_selection)) +
        geom_bar(stat = "identity") +
        facet_wrap(~seurat_names, scales = "free_y")
    
    
    cnt_dt = meta_dt[, .(.N), .(seurat_names, seurat_clusters, in_selection)][order(in_selection)][order(seurat_clusters)]
    cnt_dt = dcast(cnt_dt, seurat_names+seurat_clusters~in_selection, fill = 0, value.var = "N")
    if(is.null(cnt_dt[["FALSE"]])) cnt_dt[["FALSE"]] = 0
    setnames(cnt_dt, c("FALSE", "TRUE"), c("not_selected", "in_selection"))
    fwrite(cnt_dt, res_file("selection_count_per_cluster.csv"))
    p_sel_per_cluster = ggplot(meta_dt, aes(x = seurat_names, fill = in_selection)) +
        geom_bar() +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
        scale_fill_manual(values = c("FALSE" = "gray", 'TRUE' = "blue")) +
        labs(title = "Selection breakdown by Seurat cluster")
    ggsave(res_file("selection_count_per_cluster.png"), p_sel_per_cluster, width = 4.65, height = 3.5)
    ggsave(res_file("selection_count_per_cluster.pdf"), p_sel_per_cluster, width = 4.65, height = 3.5)
    
    .p = my_plot_cells(mon, return_data = TRUE)
    p_dt = as.data.table(.p$data_df)
    .p_sel = my_plot_cells(sel_branch, color_cells_by = "pseudotime", return_data = TRUE) 
    p_dt_sel = as.data.table(.p_sel$data_df)
    p_seg = as.data.table(.p$edge_df)
    p_dt_sel$pseudotime = pseudotime(sel_branch)[p_dt_sel$sample_name]
    p_pseudo = ggplot(p_dt_sel, aes(x = data_dim_1, y = data_dim_2, color = pseudotime))+
        annotate("point", x = p_dt$data_dim_1, y = p_dt$data_dim_2, size = .1, color = "gray") +
        geom_point() +
        annotate("segment", 
                 x= p_seg$source_prin_graph_dim_1, y = p_seg$source_prin_graph_dim_2, 
                 xend = p_seg$target_prin_graph_dim_1, yend = p_seg$target_prin_graph_dim_2, 
                 color = "black") +
        scale_color_viridis_c(option= "magma") +
        theme_classic() +
        labs(title = "Pseudotime of selection")
    ggsave(res_file("plot_pseudotime_of_selection.png"), p_pseudo, width = 6.85, height = 5.8)
    ggsave(res_file("plot_pseudotime_of_selection.pdf"), p_pseudo, width = 6.85, height = 5.8)
    #find gene that correlate with pseudotime
    
    
    sig_file = res_file(paste0("morans_test_significant.q_", max_q, ".morans_", min_morans, ".csv"))
    full_file = res_file("morans_test_full.csv")
    if(!file.exists(full_file)){
        pr_test_res.sub <- graph_test(sel_branch, neighbor_graph="principal_graph", cores=10)
        dt_test_res = as.data.table(pr_test_res.sub)
        dt_test_res = dt_test_res[order(q_value)]
        fwrite(dt_test_res, full_file)
    }else{
        dt_test_res = fread(full_file)
    }
    
    if(nrow(dt_test_res[morans_I >= min_morans]) < 50){
        old_morans = min_morans
        min_morans = round(dt_test_res[60,]$morans_I * .95, 3)
        # olq_q = max_q
        # max_q = round(max(dt_test_res[morans_I >= min_morans]$q_value)*1.05, 3)
        warning("reducing min_morans from ", old_morans, " to ", min_morans," to support module scoring.\nq_value will be ignored.")
        max_q = 1
    }
    
    if(!file.exists(sig_file)){
        dt_sig_res = dt_test_res[q_value <= max_q & morans_I >= min_morans]
        fwrite(dt_sig_res, sig_file)
    }else{
        dt_sig_res = fread(sig_file)
    }
    
    num_vars = colnames(dt_sig_res)[sapply(dt_sig_res, class) == "numeric"]
    DT::datatable(dt_sig_res) %>%
        DT::formatRound(columns=num_vars, digits=3)
    
    # calculate modules based on correlated genes
    module_file = res_file("module_assignment.csv")
    if(!file.exists(module_file)){
        tryCatch({
            gene_module_dt <- find_gene_modules(sel_branch[dt_sig_res$id,], 
                                                resolution=NULL, 
                                                cores = 20) %>% as.data.table
            message(gene_module_dt$module %>% unique %>% length, " modules found")
            fwrite(gene_module_dt, module_file)    
        }, 
        error = function(e){
            message("finding modules failed, considering all genes as single module.")
            #find_gene_modules UMAP method
            preprocess_mat <- sel_branch@preprocess_aux$gene_loadings
            if (is.null(sel_branch@preprocess_aux$beta) == FALSE) {
                preprocess_mat = preprocess_mat %*% sel_branch@preprocess_aux$beta
            }
            preprocess_mat = preprocess_mat[intersect(dt_sig_res$id, 
                                                      row.names(preprocess_mat)), ]
            umap_res = uwot::umap(as.matrix(preprocess_mat), n_components = 2, 
                                  metric = "cosine", min_dist = 0.1, n_neighbors = 15L, 
                                  fast_sgd = FALSE, n_threads = 20, verbose = FALSE, 
                                  nn_method = "annoy")
            row.names(umap_res) <- row.names(preprocess_mat)
            colnames(umap_res) <- paste0("dim_", 1:ncol(umap_res))
            reduced_dim_res <- umap_res
            gene_module_dt = data.table(id = rownames(reduced_dim_res), 
                                        module = 1, 
                                        supermodule = 1, 
                                        dim_1 = reduced_dim_res[, "dim_1"], 
                                        dim_2 = reduced_dim_res[, "dim_2"])
            fwrite(gene_module_dt, module_file)    
        })
        
    }else{
        gene_module_dt = fread(module_file)
    }
    
    # look at Seurat clusters in modules
    cell_group_df <- tibble::tibble(cell=row.names(colData(sel_branch)),
                                    cell_group=colData(sel_branch)$seurat_names)
    cell_group_df$cell_group = droplevels(cell_group_df$cell_group)
    cell_group_df = cell_group_df[order(cell_group_df$cell_group),]
    cell_group_df$cell_group %>% unique
    
    
    agg_mat <- my_aggregate_gene_expression(sel_branch, gene_module_dt, cell_group_df)
    row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
    lim = max(abs(range(agg_mat)))
    if(nrow(agg_mat) > 2){
        #this is ugly but it's the way I know how to do it
        pheat_parts = ssvRecipes::plot_hclust_heatmap(as.matrix(agg_mat), Colv = FALSE)
        
    }else{
        pheat_parts = ssvRecipes::plot_hclust_heatmap(as.matrix(agg_mat), Colv = FALSE, Rowv = FALSE)
    }
    pheat_parts$heatmap = pheat_parts$heatmap +
        scale_fill_gradientn(colours = c("blue", "white", "red"), limits = c(-lim, lim)) +
        labs(title = "Module expression in Seurat clusters") +
        theme(plot.title = element_text(size = 12))
    pheat_parts$heatmap
    p1 = ssvRecipes::plot_hclust_heatmap.assemble(pheat_parts)
    
    p_2_data = my_plot_cells(mon,
                             genes=gene_module_dt %>% dplyr::filter(module %in% c(1,2, 3,11, 6)),
                             label_cell_groups=FALSE,
                             show_trajectory_graph=FALSE, scale_to_range = TRUE, cell_size = .5, 
                             return_data = TRUE)
    p_2_dt = as.data.table(p_2_data$data_df)
    setnames(p_2_dt, c("data_dim_1", "data_dim_2"), c("UMAP_1", "UMAP_2"))
    p2 = ggplot(p_2_dt, aes(x = UMAP_1, y = UMAP_2, color = value)) +
        annotate("point", x = p_2_dt$UMAP_1, y = p_2_dt$UMAP_2, size = .1, color = "gray60") +
        geom_point(size = .4) +
        facet_grid(feature_label~seurat_names) +
        scale_color_viridis_c() +
        theme(panel.grid = element_blank(), panel.background = element_blank()) +
        labs(title = "Seurat clusters by pseudotime modules", color = "Module\nexpression")
    
    pg_seurat_clusters = cowplot::plot_grid(p1, p2, rel_widths = c(1, 1.5))
    
    ggsave(plot = pg_seurat_clusters, 
           res_file("plot_seurat_clusters.png"), 
           width = 14.4, 
           height = 5.4)
    
    p_sel_modules = my_plot_cells(sel_branch,
                               genes=gene_module_dt %>% dplyr::filter(module %in% c(1:4)),
                               label_cell_groups=FALSE,
                               show_trajectory_graph=FALSE, scale_to_range = TRUE, cell_size = .5) +
        labs(title = "Selected module aggregated expression in Monocle selected region") +
        theme(plot.title = element_text(size = 12))
    
    
    ggsave(plot = p_sel_modules, 
           res_file("plot_selected_module_expression.png"), 
           width = 5.8, 
           height = 5.4)
    
    #plot selected genes in pseudotime
    goi = dt_sig_res[order(morans_I, decreasing = TRUE)]$gene_short_name[1:12]
    p_pseudo_top12 = my_plot_genes_in_pseudotime(sel_branch, goi = goi, ncol = 4)
    ggsave(plot = p_pseudo_top12, 
           res_file("plot_morans_top12.png"), 
           width = 7.7, 
           height = 5.4)
    
    # plot module in pseudotime
    ps = pseudotime(sel_branch)
    if(max(ps) < 5){
        ps = ps * 4.999/max(ps)
        warning("pseudtotime has been inflated")
    }
    max(ps)
    cell_group_df <- tibble::tibble(cell=names(ps),
                                    cell_group=floor(ps))
    
    agg_mat <- my_aggregate_gene_expression(sel_branch, gene_module_dt, cell_group_df)
    row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
    set.seed(0)
    if(nrow(agg_mat) > 2){
        ph_res = pheatmap::pheatmap(agg_mat, cluster_cols = FALSE,
                                    scale="column", clustering_method="ward.D2")    
    }else{
        ph_res = list(tree_row = list(order = seq_len(nrow(agg_mat))))
    }
    
    
    my_melt = function(mat, rownames.var = "row", colnames.var = "column", value.var = "value", rownames.order = NULL, colnames.order = NULL){
        if(is.null(rownames.order)){
            rownames.order = rownames(mat)
        }
        if(is.null(colnames.order)){
            colnames.order = colnames(mat)
        }
        dt = melt(as.data.table(mat, keep.rownames = TRUE), id.vars = "rn")
        setnames(dt, c("rn", "variable", "value"), c(rownames.var, colnames.var, value.var)) 
        dt[[rownames.var]] = factor(dt[[rownames.var]], levels = rownames.order)
        dt[[colnames.var]] = factor(dt[[colnames.var]], levels = colnames.order)
        dt
    }
    
    agg_dt = my_melt(agg_mat, rownames.order = rev(rownames(agg_mat))[ph_res$tree_row$order])
    lim = max(abs(agg_dt$value))
    
    p_heatmap_modules = ggplot(agg_dt, aes(x = column, y = row, fill = value)) + geom_raster() +
        scale_fill_gradientn(colours = c("blue", "white", "red"), limits = c(-lim, lim)) +
        labs(x = "pseudotime", y = "gene module", fill = "z-score", 
             title = "modules along trajectory") +
        theme(panel.background = element_blank())
    p_heatmap_modules
    ggsave(res_file("plot_modules_heatmap.png"),p_heatmap_modules, width = 5.8, height = 3.8)
    
    tmp =  as.data.table(gene_module_dt)
    gm_dt = as.data.table(gene_module_dt)[, .(id, module)]
    gm_dt = merge(dt_sig_res, gm_dt, by = "id")
    gm_dt[, module_name := paste("Module", module)]
    gm_dt$module_name = factor(gm_dt$module_name, levels = levels(agg_dt$row))
    gl = split(gm_dt$gene_short_name, gm_dt$module_name)
    
    sel_module_genes = sapply(gl, function(g){
        dt_sig_res[gene_short_name %in% g, gene_short_name[which.max(morans_I)]]
    })
    lengths(gl)
    
    p_sel_module_genes = my_plot_genes_in_pseudotime(sel_branch, sel_module_genes, ncol = 3)
    p_sel_module_genes = p_sel_module_genes + labs(title = "highest correlated gene per module")
    ggsave(res_file("plot_pseudotime_best_per_module.png"), p_sel_module_genes, width = 5.8, height = 5.5)
    
    write_geneList_matrix = function(gl, file){
        if(is.null(names(gl))) stop("gene lists must be named")
        gmat = matrix("", nrow = max(lengths(gl)), ncol = length(gl))
        for(i in seq_along(gl)){
            x = gl[[i]]
            gmat[seq_along(x), i] = as.character(x)
        }
        colnames(gmat) = names(gl)
        write.csv(gmat, file, row.names = FALSE, quote = FALSE)
    }
    
    
    write_geneList_matrix(gl, heatmap_gl_file)
    
    bfc = BiocFileCache::BiocFileCache()
    gene_lists = lapply(gl, toupper)
    gene_lists = gl
    bg_genes = (rowData(mon)$gene_short_name)
    c_gl = gl
    names(c_gl) = sub("Module ", "", names(c_gl))
    cprof = monocle_clusterProfiler_fromGenes(c_gl, bg_genes = rowData(mon)$gene_short_name)
    saveRDS(cprof[[1]], res_file("module_clusterProfiler_result.Rds"))
    c_p = cprof[[2]]
    c_p = c_p + theme(axis.text.y = element_text(size = 8), 
                axis.text.x = element_text(size = 8)) +
        labs(x = "Module", title = "Gene Ontology per Gene Module")
    ggsave(res_file("module_clusterProfiler.pdf"), c_p, width = 5+2*length(c_gl), height = 2+2.5*length(c_gl))
    
}

monocle_clusterProfiler_fromGenes = function(gene_lists, bg_genes = NULL, force_overwrite = FALSE){
    if(is.null(bg_genes)){
        bg_genes = unique(unlist(gene_lists))
    }
    if(GEN == "mm10"){
        OrgDb = org.Mm.eg.db::org.Mm.eg.db
    }else if(GEN == "hg38"){
        OrgDb = org.Hs.eg.db::org.Hs.eg.db
    }else{
        stop("GEN value of ", GEN, " not implemented. Valid GEN are mm10 or hg38")
    }
    
    rname_go_dat = digest::digest(list(gene_lists, bg_genes, "BP", GEN))
    rname_go_plot = digest::digest(list(gene_lists, bg_genes, "BP", GEN, "plot"))
    bfc = BiocFileCache::BiocFileCache()
    res = bfcif(bfc, rname_go_plot, force_overwrite = force_overwrite, function(){
        message("calc go res")
        tryCatch(
            expr = {
                ck = bfcif(bfc, rname_go_dat, force_overwrite = force_overwrite, function(){
                    message("calc compareCluster")
                    clusterProfiler::compareCluster(geneCluster = gene_lists,
                                   universe      = bg_genes,
                                   fun = "enrichGO",
                                   OrgDb = OrgDb,
                                   keyType       = 'SYMBOL',
                                   ont           = "BP",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.01,
                                   qvalueCutoff  = 0.05)
                })
                # ck = my_clusterProfiler(qgr, clust_assign)
                p = clusterProfiler::dotplot(clusterProfiler::simplify(ck))
                list(ck, p)
            }, error = {
                function(e){
                    ck = NULL
                    p = ggplot() + annotate("text", x = 0, y = 0, label ="no GO results")
                    list(ck, p)
                }
            })
    })
    res
}
