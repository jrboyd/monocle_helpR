library(magrittr)
SEL_POSS = list(full = c(4, 13, 7, 0, 5, 3, 1, 16, 2, 10, 15, 6, 11, 14, 8, 9, 12),
                contiguous      = c(4, 13, 7, 0, 5, 3, 1, 16, 2, 10),
                "13reg" = c(4, 13, 7, 0),
                "10reg" = c(10, 3, 16, 2, 5, 1),
                "15_6reg" = c(15, 6))

SEL_POSS_NEW_CLUST = list(
    nc_contiguous78 = c(1:8),
    nc_2378 = c(2, 3, 7, 8),
    nc_not65 = c(1:4, 7:8),
    nc_5678 = c(5:8),
    nc_7812 = c(1, 2, 7, 8),
    nc_7813 = c(1, 3, 7, 8),
    nc_78 = c(7, 8),
    nc_56 = c(5, 6),
    nc_156 = c(1, 5, 6),
    nc_41 = c(4, 1),
    nc_21 = c(2, 1)
)

my_monocle_process = function(cds_input, 
                              seu_input,
                              sel_desc, 
                              align_apply,
                              force = FALSE,
                              new_clust = FALSE,
                              num_dim_pca = 30,
                              num_dim_umap = 30){
    if(new_clust){
        stopifnot(sel_desc %in% names(SEL_POSS_NEW_CLUST))
        sel_clust = SEL_POSS_NEW_CLUST[[sel_desc]]    
    }else{
        stopifnot(sel_desc %in% names(SEL_POSS))
        sel_clust = SEL_POSS[[sel_desc]]    
    }
    
    analysis_name = my_monocle_getwd(sel_desc, align_apply, new_clust = new_clust)
    dir.create(analysis_name, showWarnings = FALSE)
    
    
    res_file = function(file){
        file.path(analysis_name, file)
    }
    proc_data_file = res_file("data_processed.Rds")
    png_desc_file = res_file("UMAP_description.png")
    seu_meta_file = res_file("seurat_meta.Rds")
    digest_file = res_file("digest_data_processed.txt")
    
    if(file.exists(proc_data_file) & !force){
        cds = readRDS(proc_data_file)
        meta_dt = readRDS(seu_meta_file)
        # rname = read.table(digest_file)[1,1] %>% as.character
        # stopifnot(digest::digest(cds) == rname)
        
        
        if(new_clust){
            stopifnot(setequal(meta_dt[sel == TRUE]$meta_cluster %>% as.character, sel_clust))
        }else{
            stopifnot(setequal(meta_dt[sel == TRUE]$seurat_clusters %>% as.character, sel_clust))
        }
        
        plot_data = my_plot_cells(cds, return_data = TRUE)
    }else{
        # sel_clust = c(4, 13, 7, 0)
        meta_dt = get_meta_dt(seu_input)
        setkey(meta_dt, id)
        if(new_clust){
            meta_dt[, sel := meta_cluster %in% sel_clust]
        }else{
            meta_dt[, sel := seurat_clusters %in% sel_clust]    
        }
        
        meta_dt$sel = factor(meta_dt$sel, levels = c("TRUE", "FALSE"))
        
        saveRDS(meta_dt, seu_meta_file)
        
        lab_dt = meta_dt[, .(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)), .(seurat_clusters, sel)]
        # p = ggplot(meta_dt, aes(x = UMAP_1, y = UMAP_2, label = seurat_clusters, color = seurat_clusters)) +
        #     geom_point(size = .3) +
        #     geom_text(data = lab_dt, color= "black") +
        #     facet_wrap(~sel, drop = FALSE)
        # ggsave("clusters.png", p, width = 8.6, height = 4.2)
        
        cds = cds_input[, meta_dt[sel == TRUE]$id]
        
        pData(cds)$seurat_clusters = meta_dt[.(rownames(pData(cds))),]$seurat_clusters %>% factor
        pData(cds)$meta_cluster = meta_dt[.(rownames(pData(cds))),]$meta_cluster
        pData(cds)$genotype = meta_dt[.(rownames(pData(cds))),]$orig.ident
        
        table(pData(cds)$seurat_clusters)
        # cds <- preprocess_cds(cds, use_genes =  seu_input@assays$integrated %>% rownames,num_dim = 30)
        # ref = as.data.table(rowData(cds))
        # setkey(ref, gene_short_name)
        # seu_int_gs = pbmc@assays$integrated %>% rownames
        # sue_int_id = ref[.(seu_int_gs)]$id
        # sue_int_id = sue_int_id[!is.na(sue_int_id)]
        # seqsetvis::ssvFeatureVenn(list(cds = rownames(cds), seu = sue_int_id))
        # cds <- preprocess_cds(cds, use_genes =  sue_int_id,num_dim = 30)
        # browser()
        
        cds <- preprocess_cds(cds, num_dim = num_dim_pca)
        plot_pc_variance_explained(cds)
        if(align_apply){
            cds <- align_cds(cds, alignment_group = "genotype")    
        }
        
        
        cds = reduce_dimension(cds, reduction_method = "UMAP")
        cds = cluster_cells(cds)
        
        ## Step 5: Learn a graph
        cds <- learn_graph(cds, learn_graph_control = list(minimal_branch_len = 8))
        
        ## Step 6: Order cells
        # cds <- order_cells(cds)
        
        plot_data = my_plot_cells(cds, return_data = TRUE)
        p_dt = plot_data$data_df %>% as.data.table
        setnames(p_dt, c("data_dim_1", "data_dim_2", "sample_name"), c("UMAP_1", "UMAP_2", "id"))
        
        if(new_clust){
            lab_dt = p_dt[, .(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)), .(meta_cluster, genotype)]
            p_seu = ggplot(p_dt, aes(x = UMAP_1, y = UMAP_2, color = meta_cluster, label = meta_cluster)) + 
                geom_point(size = .2) +
                geom_label(data = lab_dt, show.legend = FALSE) +
                facet_wrap(~genotype) +
                cowplot::theme_cowplot()
        }else{
            lab_dt = p_dt[, .(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)), .(seurat_clusters, genotype)]
            p_seu = ggplot(p_dt, aes(x = UMAP_1, y = UMAP_2, color = seurat_clusters, label = seurat_clusters)) + 
                geom_point(size = .2) +
                geom_label(data = lab_dt, show.legend = FALSE) +
                facet_wrap(~genotype) +
                cowplot::theme_cowplot()    
        }
        
        
        pg = cowplot::plot_grid(plot_data$plot, p_seu, rel_widths = c(1, 2))
        ggsave(png_desc_file, pg, width = 13.6, height = 4.7)
        
        saveRDS(cds, file = proc_data_file)
        # write.table(digest::digest(cds), digest_file, 
        #             quote = FALSE, col.names = FALSE, row.names = FALSE)
        
    }
    invisible(list(cds = cds, meta_dt = meta_dt, plot_data = plot_data))
}

get_meta_dt.monocle = function(cds){
    plot_data = my_plot_cells(cds, return_data = TRUE)
    p_dt = plot_data$data_df %>% as.data.table
    setnames(p_dt, c("data_dim_1", "data_dim_2", "sample_name"), c("UMAP_1", "UMAP_2", "id"))
    p_dt
}

my_monocle_relearn = function(sel_desc, 
                              align_apply,
                              new_clust = FALSE,
                              ...){
    if(new_clust){
        stopifnot(sel_desc %in% names(SEL_POSS_NEW_CLUST))
        sel_clust = SEL_POSS_NEW_CLUST[[sel_desc]]    
    }else{
        stopifnot(sel_desc %in% names(SEL_POSS))
        sel_clust = SEL_POSS[[sel_desc]]    
    }
    analysis_name = my_monocle_getwd(sel_desc, align_apply, new_clust = new_clust)
    
    
    res_file = function(file){
        file.path(analysis_name, file)
    }
    proc_data_file = res_file("data_processed.Rds")
    png_desc_file = res_file("UMAP_description.png")
    seu_meta_file = res_file("seurat_meta.Rds")
    digest_file = res_file("digest_data_processed.txt")
    
    relearn_data_file =  res_file("relearn_data_processed.Rds")
    png_relearn_file = res_file("relearn_UMAP_description.png")
    digest_relearn_file = res_file("relearn_digest_data_processed.txt")
    if(file.exists(proc_data_file)){
        cds = readRDS(proc_data_file)
        meta_dt = readRDS(seu_meta_file)
        # rname = read.table(digest_file)[1,1] %>% as.character
        # stopifnot(digest::digest(cds) == rname)
        if(new_clust){
            meta_dt$seurat_clusters = NULL
            setnames(meta_dt, "meta_cluster", "seurat_clusters") 
        }
        stopifnot(setequal(meta_dt[sel == TRUE]$seurat_clusters %>% as.character, sel_clust))    
        
        
        
        cds@clusters$UMAP$partitions[cds@clusters$UMAP$partitions != "1"] <- "1"
        ## Step 5: Learn a graph
        cds <- learn_graph(cds, ...)
        
        ## Step 6: Order cells
        # cds <- order_cells(cds)
        
        plot_data = my_plot_cells(cds, return_data = TRUE)
        p_dt = plot_data$data_df %>% as.data.table
        setnames(p_dt, c("data_dim_1", "data_dim_2", "sample_name"), c("UMAP_1", "UMAP_2", "id"))
        lab_dt = p_dt[, .(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)), .(seurat_clusters, genotype)]
        p_seu = ggplot(p_dt, aes(x = UMAP_1, y = UMAP_2, color = seurat_clusters, label = seurat_clusters)) + 
            geom_point(size = .2) +
            geom_label(data = lab_dt, show.legend = FALSE) +
            facet_wrap(~genotype) +
            cowplot::theme_cowplot()
        
        pg = cowplot::plot_grid(plot_data$plot, p_seu, rel_widths = c(1, 2))
        ggsave(png_relearn_file, pg, width = 13.6, height = 4.7)
        saveRDS(cds, file = relearn_data_file)
        # write.table(digest::digest(cds), digest_relearn_file, 
        #             quote = FALSE, col.names = FALSE, row.names = FALSE)
        invisible(list(cds = cds, meta_dt = meta_dt, plot_data = plot_data))
    }else{
        stop("original data not found.")
    }
}

#fix for bug/unexpected cased when thereis a single gene group
my_aggregate_gene_expression = function (cds, gene_group_df = NULL, cell_group_df = NULL, norm_method = c("log", 
                                                                                                          "binary", "size_only"), pseudocount = 1, scale_agg_values = TRUE, 
                                         max_agg_value = 3, min_agg_value = -3, exclude.na = TRUE) 
{
    if (is.null(gene_group_df) && is.null(cell_group_df)) 
        stop("Error: one of either gene_group_df or cell_group_df must not be NULL")
    agg_mat <- normalized_counts(cds, norm_method = norm_method, 
                                 pseudocount = pseudocount)
    if (is.null(gene_group_df) == FALSE) {
        gene_group_df <- as.data.frame(gene_group_df)
        gene_group_df <- gene_group_df[gene_group_df[, 1] %in% 
                                           fData(cds)$gene_short_name | gene_group_df[, 1] %in% 
                                           row.names(fData(cds)), , drop = FALSE]
        agg_mat = as.matrix(monocle3:::my.aggregate.Matrix(agg_mat[gene_group_df[, 
                                                                                 1], ], as.factor(gene_group_df[, 2]), fun = "sum"))
        if (scale_agg_values) {
            agg_mat <- t(scale(t(agg_mat)))
            agg_mat[agg_mat < min_agg_value] <- min_agg_value
            agg_mat[agg_mat > max_agg_value] <- max_agg_value
        }
    }
    if (is.null(cell_group_df) == FALSE) {
        cell_group_df = as.data.frame(cell_group_df)
        cell_group_df = cell_group_df[cell_group_df[, 1] %in% 
                                          row.names(pData(cds)), , drop = FALSE]
        agg_mat = agg_mat[, cell_group_df[, 1], drop = FALSE]
        agg_mat = monocle3:::my.aggregate.Matrix(Matrix::t(agg_mat), as.factor(cell_group_df[, 
                                                                                             2]), fun = "mean")
        agg_mat = Matrix::t(agg_mat)
    }
    if (exclude.na) {
        agg_mat = agg_mat[row.names(agg_mat) != "NA", colnames(agg_mat) != 
                              "NA", drop = FALSE]
    }
    return(agg_mat)
}

#return data
my_plot_cells = function (cds, x = 1, y = 2, 
                          reduction_method = c("UMAP", "tSNE", "PCA", "LSI", "Aligned"), 
                          color_cells_by = "cluster", 
                          group_cells_by = c("cluster", "partition", colnames(colData(cds))), genes = NULL, show_trajectory_graph = TRUE, 
                          trajectory_graph_color = "grey28", trajectory_graph_segment_size = 0.75, 
                          norm_method = c("log", "size_only"), label_cell_groups = TRUE, 
                          label_groups_by_cluster = TRUE, group_label_size = 2, labels_per_group = 1, 
                          label_branch_points = TRUE, label_roots = TRUE, label_leaves = TRUE, 
                          graph_label_size = 2, cell_size = 0.35, cell_stroke = I(cell_size/2), 
                          alpha = 1, min_expr = 0.1, rasterize = FALSE, scale_to_range = FALSE, return_data = FALSE) 
{
    reduction_method <- match.arg(reduction_method)
    assertthat::assert_that(methods::is(cds, "cell_data_set"))
    assertthat::assert_that(!is.null(reducedDims(cds)[[reduction_method]]), 
                            msg = paste("No dimensionality reduction for", reduction_method, 
                                        "calculated.", "Please run reduce_dimensions with", 
                                        "reduction_method =", reduction_method, "before attempting to plot."))
    low_dim_coords <- reducedDims(cds)[[reduction_method]]
    assertthat::assert_that(ncol(low_dim_coords) >= max(x, y), 
                            msg = paste("x and/or y is too large. x and y must", 
                                        "be dimensions in reduced dimension", "space."))
    if (!is.null(color_cells_by)) {
        assertthat::assert_that(color_cells_by %in% c("cluster", 
                                                      "partition", "pseudotime") | color_cells_by %in% 
                                    names(colData(cds)), msg = paste("color_cells_by must one of", 
                                                                     "'cluster', 'partition', 'pseudotime,", "or a column in the colData table."))
        if (color_cells_by == "pseudotime") {
            tryCatch({
                pseudotime(cds, reduction_method = reduction_method)
            }, error = function(x) {
                stop(paste("No pseudotime for", reduction_method, 
                           "calculated. Please run order_cells with", 
                           "reduction_method =", reduction_method, "before attempting to color by pseudotime."))
            })
        }
    }
    assertthat::assert_that(!is.null(color_cells_by) || !is.null(markers), 
                            msg = paste("Either color_cells_by or markers must", 
                                        "be NULL, cannot color by both!"))
    norm_method = match.arg(norm_method)
    group_cells_by = match.arg(group_cells_by)
    assertthat::assert_that(!is.null(color_cells_by) || !is.null(genes), 
                            msg = paste("Either color_cells_by or genes must be", 
                                        "NULL, cannot color by both!"))
    if (show_trajectory_graph && is.null(principal_graph(cds)[[reduction_method]])) {
        message("No trajectory to plot. Has learn_graph() been called yet?")
        show_trajectory_graph = FALSE
    }
    gene_short_name <- NA
    sample_name <- NA
    data_dim_1 <- NA
    data_dim_2 <- NA
    if (rasterize) {
        plotting_func <- ggrastr::geom_point_rast
    }
    else {
        plotting_func <- ggplot2::geom_point
    }
    S_matrix <- reducedDims(cds)[[reduction_method]]
    data_df <- data.frame(S_matrix[, c(x, y)])
    colnames(data_df) <- c("data_dim_1", "data_dim_2")
    data_df$sample_name <- row.names(data_df)
    data_df <- as.data.frame(cbind(data_df, colData(cds)))
    if (group_cells_by == "cluster") {
        data_df$cell_group <- tryCatch({
            clusters(cds, reduction_method = reduction_method)[data_df$sample_name]
        }, error = function(e) {
            NULL
        })
    }
    else if (group_cells_by == "partition") {
        data_df$cell_group <- tryCatch({
            partitions(cds, reduction_method = reduction_method)[data_df$sample_name]
        }, error = function(e) {
            NULL
        })
    }
    else {
        stop("Error: unrecognized way of grouping cells.")
    }
    if (color_cells_by == "cluster") {
        data_df$cell_color <- tryCatch({
            clusters(cds, reduction_method = reduction_method)[data_df$sample_name]
        }, error = function(e) {
            NULL
        })
    }
    else if (color_cells_by == "partition") {
        data_df$cell_color <- tryCatch({
            partitions(cds, reduction_method = reduction_method)[data_df$sample_name]
        }, error = function(e) {
            NULL
        })
    }
    else if (color_cells_by == "pseudotime") {
        data_df$cell_color <- tryCatch({
            pseudotime(cds, reduction_method = reduction_method)[data_df$sample_name]
        }, error = function(e) {
            NULL
        })
    }
    else {
        data_df$cell_color <- colData(cds)[data_df$sample_name, 
                                           color_cells_by]
    }
    if (show_trajectory_graph) {
        ica_space_df <- t(cds@principal_graph_aux[[reduction_method]]$dp_mst) %>% 
            as.data.frame() %>% dplyr::select_(prin_graph_dim_1 = x, 
                                               prin_graph_dim_2 = y) %>% dplyr::mutate(sample_name = rownames(.), 
                                                                                       sample_state = rownames(.))
        dp_mst <- cds@principal_graph[[reduction_method]]
        edge_df <- dp_mst %>% igraph::as_data_frame() %>% dplyr::select_(source = "from", 
                                                                         target = "to") %>% dplyr::left_join(ica_space_df %>% 
                                                                                                                 dplyr::select_(source = "sample_name", source_prin_graph_dim_1 = "prin_graph_dim_1", 
                                                                                                                                source_prin_graph_dim_2 = "prin_graph_dim_2"), 
                                                                                                             by = "source") %>% dplyr::left_join(ica_space_df %>% 
                                                                                                                                                     dplyr::select_(target = "sample_name", target_prin_graph_dim_1 = "prin_graph_dim_1", 
                                                                                                                                                                    target_prin_graph_dim_2 = "prin_graph_dim_2"), 
                                                                                                                                                 by = "target")
    }
    markers_exprs <- NULL
    expression_legend_label <- NULL
    if (!is.null(genes)) {
        if (!is.null(dim(genes)) && dim(genes) >= 2) {
            markers = unlist(genes[, 1], use.names = FALSE)
        }
        else {
            markers = genes
        }
        markers_rowData <- as.data.frame(subset(rowData(cds), 
                                                gene_short_name %in% markers | row.names(rowData(cds)) %in% 
                                                    markers))
        if (nrow(markers_rowData) == 0) {
            stop("None of the provided genes were found in the cds")
        }
        if (nrow(markers_rowData) >= 1) {
            cds_exprs <- SingleCellExperiment::counts(cds)[row.names(markers_rowData), 
                                                           , drop = FALSE]
            cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/size_factors(cds))
            if (!is.null(dim(genes)) && dim(genes) >= 2) {
                genes = as.data.frame(genes)
                row.names(genes) = genes[, 1]
                genes = genes[row.names(cds_exprs), ]
                agg_mat = as.matrix(my_aggregate_gene_expression(cds, 
                                                              genes, norm_method = norm_method, scale_agg_values = FALSE))
                markers_exprs = agg_mat
                markers_exprs <- reshape2::melt(markers_exprs)
                colnames(markers_exprs)[1:2] <- c("feature_id", 
                                                  "cell_id")
                if (is.factor(genes[, 2])) 
                    markers_exprs$feature_id = factor(markers_exprs$feature_id, 
                                                      levels = levels(genes[, 2]))
                markers_exprs$feature_label <- markers_exprs$feature_id
                norm_method = "size_only"
                expression_legend_label = "Expression score"
            }
            else {
                cds_exprs@x = round(10000 * cds_exprs@x)/10000
                markers_exprs = matrix(cds_exprs, nrow = nrow(markers_rowData))
                colnames(markers_exprs) = colnames(SingleCellExperiment::counts(cds))
                row.names(markers_exprs) = row.names(markers_rowData)
                markers_exprs <- reshape2::melt(markers_exprs)
                colnames(markers_exprs)[1:2] <- c("feature_id", 
                                                  "cell_id")
                markers_exprs <- merge(markers_exprs, markers_rowData, 
                                       by.x = "feature_id", by.y = "row.names")
                if (is.null(markers_exprs$gene_short_name)) {
                    markers_exprs$feature_label <- as.character(markers_exprs$feature_id)
                }
                else {
                    markers_exprs$feature_label <- as.character(markers_exprs$gene_short_name)
                }
                markers_exprs$feature_label <- ifelse(is.na(markers_exprs$feature_label) | 
                                                          !as.character(markers_exprs$feature_label) %in% 
                                                          markers, as.character(markers_exprs$feature_id), 
                                                      as.character(markers_exprs$feature_label))
                markers_exprs$feature_label <- factor(markers_exprs$feature_label, 
                                                      levels = markers)
                if (norm_method == "size_only") 
                    expression_legend_label = "Expression"
                else expression_legend_label = "log10(Expression)"
            }
            if (scale_to_range) {
                markers_exprs = dplyr::group_by(markers_exprs, 
                                                feature_label) %>% dplyr::mutate(max_val_for_feature = max(value), 
                                                                                 min_val_for_feature = min(value)) %>% dplyr::mutate(value = 100 * 
                                                                                                                                         (value - min_val_for_feature)/(max_val_for_feature - 
                                                                                                                                                                            min_val_for_feature))
                expression_legend_label = "% Max"
            }
        }
    }
    if (label_cell_groups && is.null(color_cells_by) == FALSE) {
        if (is.null(data_df$cell_color)) {
            if (is.null(genes)) {
                message(paste(color_cells_by, "not found in colData(cds), cells will", 
                              "not be colored"))
            }
            text_df = NULL
            label_cell_groups = FALSE
        }
        else {
            if (is.character(data_df$cell_color) || is.factor(data_df$cell_color)) {
                if (label_groups_by_cluster && is.null(data_df$cell_group) == 
                    FALSE) {
                    text_df = data_df %>% dplyr::group_by(cell_group) %>% 
                        dplyr::mutate(cells_in_cluster = dplyr::n()) %>% 
                        dplyr::group_by(cell_color, add = TRUE) %>% 
                        dplyr::mutate(per = dplyr::n()/cells_in_cluster)
                    median_coord_df = text_df %>% dplyr::summarize(fraction_of_group = dplyr::n(), 
                                                                   text_x = stats::median(x = data_dim_1), text_y = stats::median(x = data_dim_2))
                    text_df = suppressMessages(text_df %>% dplyr::select(per) %>% 
                                                   dplyr::distinct())
                    text_df = suppressMessages(dplyr::inner_join(text_df, 
                                                                 median_coord_df))
                    text_df = text_df %>% dplyr::group_by(cell_group) %>% 
                        dplyr::top_n(labels_per_group, per)
                }
                else {
                    text_df = data_df %>% dplyr::group_by(cell_color) %>% 
                        dplyr::mutate(per = 1)
                    median_coord_df = text_df %>% dplyr::summarize(fraction_of_group = dplyr::n(), 
                                                                   text_x = stats::median(x = data_dim_1), text_y = stats::median(x = data_dim_2))
                    text_df = suppressMessages(text_df %>% dplyr::select(per) %>% 
                                                   dplyr::distinct())
                    text_df = suppressMessages(dplyr::inner_join(text_df, 
                                                                 median_coord_df))
                    text_df = text_df %>% dplyr::group_by(cell_color) %>% 
                        dplyr::top_n(labels_per_group, per)
                }
                text_df$label = as.character(text_df %>% dplyr::pull(cell_color))
            }
            else {
                message(paste("Cells aren't colored in a way that allows them to", 
                              "be grouped."))
                text_df = NULL
                label_cell_groups = FALSE
            }
        }
    }
    if (!is.null(markers_exprs) && nrow(markers_exprs) > 0) {
        data_df <- merge(data_df, markers_exprs, by.x = "sample_name", 
                         by.y = "cell_id")
        data_df$value <- with(data_df, ifelse(value >= min_expr, 
                                              value, NA))
        na_sub <- data_df[is.na(data_df$value), ]
        if (norm_method == "size_only") {
            g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2)) + 
                plotting_func(aes(data_dim_1, data_dim_2), size = I(cell_size), 
                              stroke = I(cell_stroke), color = "grey80", 
                              alpha = alpha, data = na_sub) + plotting_func(aes(color = value), 
                                                                            size = I(cell_size), stroke = I(cell_stroke), 
                                                                            na.rm = TRUE) + viridis::scale_color_viridis(option = "viridis", 
                                                                                                                         name = expression_legend_label, na.value = "grey80", 
                                                                                                                         end = 0.8, alpha = alpha) + guides(alpha = FALSE) + 
                facet_wrap(~feature_label)
        }
        else {
            g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2)) + 
                plotting_func(aes(data_dim_1, data_dim_2), size = I(cell_size), 
                              stroke = I(cell_stroke), color = "grey80", 
                              data = na_sub, alpha = alpha) + plotting_func(aes(color = log10(value + 
                                                                                                  min_expr)), size = I(cell_size), stroke = I(cell_stroke), 
                                                                            na.rm = TRUE, alpha = alpha) + viridis::scale_color_viridis(option = "viridis", 
                                                                                                                                        name = expression_legend_label, na.value = "grey80", 
                                                                                                                                        end = 0.8, alpha = alpha) + guides(alpha = FALSE) + 
                facet_wrap(~feature_label)
        }
    }
    else {
        g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2))
        if (color_cells_by %in% c("cluster", "partition")) {
            if (is.null(data_df$cell_color)) {
                g <- g + geom_point(color = I("gray"), size = I(cell_size), 
                                    stroke = I(cell_stroke), na.rm = TRUE, alpha = I(alpha))
                message(paste("cluster_cells() has not been called yet, can't", 
                              "color cells by cluster"))
            }
            else {
                g <- g + geom_point(aes(color = cell_color), 
                                    size = I(cell_size), stroke = I(cell_stroke), 
                                    na.rm = TRUE, alpha = alpha)
            }
            g <- g + guides(color = guide_legend(title = color_cells_by, 
                                                 override.aes = list(size = 4)))
        }
        else if (class(data_df$cell_color) == "numeric") {
            g <- g + geom_point(aes(color = cell_color), size = I(cell_size), 
                                stroke = I(cell_stroke), na.rm = TRUE, alpha = alpha)
            g <- g + viridis::scale_color_viridis(name = color_cells_by, 
                                                  option = "C")
        }
        else {
            g <- g + geom_point(aes(color = cell_color), size = I(cell_size), 
                                stroke = I(cell_stroke), na.rm = TRUE, alpha = alpha)
            g <- g + guides(color = guide_legend(title = color_cells_by, 
                                                 override.aes = list(size = 4)))
        }
    }
    if (show_trajectory_graph) {
        g <- g + geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                                         y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                                         yend = "target_prin_graph_dim_2"), size = trajectory_graph_segment_size, 
                              color = I(trajectory_graph_color), linetype = "solid", 
                              na.rm = TRUE, data = edge_df)
        if (label_branch_points) {
            mst_branch_nodes <- monocle3:::branch_nodes(cds)
            branch_point_df <- ica_space_df %>% dplyr::slice(match(names(mst_branch_nodes), 
                                                                   sample_name)) %>% dplyr::mutate(branch_point_idx = seq_len(dplyr::n()))
            g <- g + geom_point(aes_string(x = "prin_graph_dim_1", 
                                           y = "prin_graph_dim_2"), shape = 21, stroke = I(trajectory_graph_segment_size), 
                                color = "white", fill = "black", size = I(graph_label_size * 
                                                                              1.5), na.rm = TRUE, branch_point_df) + geom_text(aes_string(x = "prin_graph_dim_1", 
                                                                                                                                          y = "prin_graph_dim_2", label = "branch_point_idx"), 
                                                                                                                               size = I(graph_label_size), color = "white", 
                                                                                                                               na.rm = TRUE, branch_point_df)
        }
        if (label_leaves) {
            mst_leaf_nodes <- monocle3:::leaf_nodes(cds)
            leaf_df <- ica_space_df %>% dplyr::slice(match(names(mst_leaf_nodes), 
                                                           sample_name)) %>% dplyr::mutate(leaf_idx = seq_len(dplyr::n()))
            g <- g + geom_point(aes_string(x = "prin_graph_dim_1", 
                                           y = "prin_graph_dim_2"), shape = 21, stroke = I(trajectory_graph_segment_size), 
                                color = "black", fill = "lightgray", size = I(graph_label_size * 
                                                                                  1.5), na.rm = TRUE, leaf_df) + geom_text(aes_string(x = "prin_graph_dim_1", 
                                                                                                                                      y = "prin_graph_dim_2", label = "leaf_idx"), 
                                                                                                                           size = I(graph_label_size), color = "black", 
                                                                                                                           na.rm = TRUE, leaf_df)
        }
        if (label_roots) {
            mst_root_nodes <- monocle3:::root_nodes(cds)
            root_df <- ica_space_df %>% dplyr::slice(match(names(mst_root_nodes), 
                                                           sample_name)) %>% dplyr::mutate(root_idx = seq_len(dplyr::n()))
            g <- g + geom_point(aes_string(x = "prin_graph_dim_1", 
                                           y = "prin_graph_dim_2"), shape = 21, stroke = I(trajectory_graph_segment_size), 
                                color = "black", fill = "white", size = I(graph_label_size * 
                                                                              1.5), na.rm = TRUE, root_df) + geom_text(aes_string(x = "prin_graph_dim_1", 
                                                                                                                                  y = "prin_graph_dim_2", label = "root_idx"), 
                                                                                                                       size = I(graph_label_size), color = "black", 
                                                                                                                       na.rm = TRUE, root_df)
        }
    }
    if (label_cell_groups) {
        g <- g + ggrepel::geom_text_repel(data = text_df, mapping = aes_string(x = "text_x", 
                                                                               y = "text_y", label = "label"), size = I(group_label_size))
        if (is.null(markers_exprs)) 
            g <- g + theme(legend.position = "none")
    }
    g <- g + monocle3:::monocle_theme_opts() + xlab(paste(reduction_method, 
                                                          x)) + ylab(paste(reduction_method, y)) + theme(legend.key = element_blank()) + 
        theme(panel.background = element_rect(fill = "white"))
    if(!exists("text_df")) text_df = NULL
    if(!exists("edge_df")) edge_df = NULL
    if(!exists("root_df")) root_df = NULL
    if(!exists("leaf_df")) leaf_df = NULL
    if(!exists("branch_point_df")) branch_point_df = NULL
    if(return_data)return(list(plot = g, data_df = data_df, text_df = text_df, edge_df = edge_df, root_df = root_df, leaf_df = leaf_df, branch_point_df = branch_point_df))
    g
}

my_monocle_getwd = function(sel_desc, 
                            align_apply,
                            new_clust = FALSE){
    if(new_clust){
        stopifnot(sel_desc %in% names(SEL_POSS_NEW_CLUST))
        sel_clust = SEL_POSS_NEW_CLUST[[sel_desc]]
    }else{
        stopifnot(sel_desc %in% names(SEL_POSS))
        sel_clust = SEL_POSS[[sel_desc]]    
    }
    
    
    align_desc = ifelse(align_apply, "aligned", "unaligned")
    analysis_name = paste0("monocle_Bcell.", align_desc, ".", sel_desc)
    dir.create(analysis_name, showWarnings = FALSE)
    analysis_name
}

my_gene_xy = function(p_cds, a, b, point_size= .3){
    # a = "Pf4"
    # b = "Cd74"
    k_gene = rowData(p_cds)$gene_short_name %in% c(a, b)
    cds_exprs <- SingleCellExperiment::counts(p_cds[k_gene,])
    cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/size_factors(p_cds))
    cds_exprs <- reshape2::melt(round(as.matrix(cds_exprs)))
    colnames(cds_exprs) <- c("id", "Cell", "expression")
    cds_exprs %>% head
    
    # browser()
    cds_pseudo = data.table(Cell = colnames(p_cds), pseudotime = pseudotime(p_cds))
    
    cds_exprs = merge(cds_exprs, cds_pseudo, by = "Cell")
    
    p_dt = merge(cds_exprs, rowData(p_cds), by = "id") %>% as.data.table
    p_dt = dcast(p_dt, pseudotime+Cell~gene_short_name, value.var = "expression")
    setnames(p_dt, "Cell", "id")
    p_dt = merge(p_dt, meta_dt[, .(id, orig.ident, seurat_clusters)], by = "id")
    ggplot(p_dt, aes_string(x = a, y = b, color = "pseudotime")) + 
        geom_jitter(size = point_size) +
        facet_wrap(~seurat_clusters) +
        scale_color_viridis_c(option = "C") +
        theme(panel.grid = element_blank(),
              panel.background = element_rect(fill = "gray50"))
}

.plot_genes_in_pseudotime = function (cds_subset, min_expr = NULL, cell_size = 0.75, nrow = NULL, 
                                      ncol = 1, panel_order = NULL, color_cells_by = "pseudotime", 
                                      trend_formula = "~ splines::ns(pseudotime, df=3)", label_by_short_name = TRUE, 
                                      vertical_jitter = NULL, horizontal_jitter = NULL, return_data = FALSE) 
{
    assertthat::assert_that(methods::is(cds_subset, "cell_data_set"))
    tryCatch({
        pseudotime(cds_subset)
    }, error = function(x) {
        stop(paste("No pseudotime calculated. Must call order_cells first."))
    })
    colData(cds_subset)$pseudotime <- pseudotime(cds_subset)
    if (!is.null(min_expr)) {
        assertthat::assert_that(assertthat::is.number(min_expr))
    }
    assertthat::assert_that(assertthat::is.number(cell_size))
    if (!is.null(nrow)) {
        assertthat::assert_that(assertthat::is.count(nrow))
    }
    assertthat::assert_that(assertthat::is.count(ncol))
    assertthat::assert_that(is.logical(label_by_short_name))
    if (label_by_short_name) {
        assertthat::assert_that("gene_short_name" %in% names(rowData(cds_subset)), 
                                msg = paste("When label_by_short_name = TRUE,", "rowData must have a column of gene", 
                                            "names called gene_short_name."))
    }
    assertthat::assert_that(color_cells_by %in% c("cluster", 
                                                  "partition") | color_cells_by %in% names(colData(cds_subset)), 
                            msg = paste("color_cells_by must be a column in the", 
                                        "colData table."))
    if (!is.null(panel_order)) {
        if (label_by_short_name) {
            assertthat::assert_that(all(panel_order %in% rowData(cds_subset)$gene_short_name))
        }
        else {
            assertthat::assert_that(all(panel_order %in% row.names(rowData(cds_subset))))
        }
    }
    assertthat::assert_that(nrow(rowData(cds_subset)) <= 100, 
                            msg = paste("cds_subset has more than 100 genes -", "pass only the subset of the CDS to be", 
                                        "plotted."))
    assertthat::assert_that(methods::is(cds_subset, "cell_data_set"))
    assertthat::assert_that("pseudotime" %in% names(colData(cds_subset)), 
                            msg = paste("pseudotime must be a column in", "colData. Please run order_cells", 
                                        "before running", "plot_genes_in_pseudotime."))
    if (!is.null(min_expr)) {
        assertthat::assert_that(assertthat::is.number(min_expr))
    }
    assertthat::assert_that(assertthat::is.number(cell_size))
    assertthat::assert_that(!is.null(size_factors(cds_subset)))
    if (!is.null(nrow)) {
        assertthat::assert_that(assertthat::is.count(nrow))
    }
    assertthat::assert_that(assertthat::is.count(ncol))
    assertthat::assert_that(is.logical(label_by_short_name))
    if (label_by_short_name) {
        assertthat::assert_that("gene_short_name" %in% names(rowData(cds_subset)), 
                                msg = paste("When label_by_short_name = TRUE,", "rowData must have a column of gene", 
                                            "names called gene_short_name."))
    }
    assertthat::assert_that(color_cells_by %in% c("cluster", 
                                                  "partition") | color_cells_by %in% names(colData(cds_subset)), 
                            msg = paste("color_cells_by must be a column in the", 
                                        "colData table."))
    if (!is.null(panel_order)) {
        if (label_by_short_name) {
            assertthat::assert_that(all(panel_order %in% rowData(cds_subset)$gene_short_name))
        }
        else {
            assertthat::assert_that(all(panel_order %in% row.names(rowData(cds_subset))))
        }
    }
    assertthat::assert_that(nrow(rowData(cds_subset)) <= 100, 
                            msg = paste("cds_subset has more than 100 genes -", "pass only the subset of the CDS to be", 
                                        "plotted."))
    f_id <- NA
    Cell <- NA
    cds_subset = cds_subset[, is.finite(colData(cds_subset)$pseudotime)]
    cds_exprs <- SingleCellExperiment::counts(cds_subset)
    cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/size_factors(cds_subset))
    cds_exprs <- reshape2::melt(round(as.matrix(cds_exprs)))
    if (is.null(min_expr)) {
        min_expr <- 0
    }
    colnames(cds_exprs) <- c("f_id", "Cell", "expression")
    cds_colData <- colData(cds_subset)
    cds_rowData <- rowData(cds_subset)
    cds_exprs <- base::merge(cds_exprs, cds_rowData, by.x = "f_id", 
                             by.y = "row.names")
    cds_exprs <- base::merge(cds_exprs, cds_colData, by.x = "Cell", 
                             by.y = "row.names")
    cds_exprs$adjusted_expression <- cds_exprs$expression
    if (label_by_short_name == TRUE) {
        if (is.null(cds_exprs$gene_short_name) == FALSE) {
            cds_exprs$feature_label <- as.character(cds_exprs$gene_short_name)
            cds_exprs$feature_label[is.na(cds_exprs$feature_label)] <- cds_exprs$f_id
        }
        else {
            cds_exprs$feature_label <- cds_exprs$f_id
        }
    }
    else {
        cds_exprs$feature_label <- cds_exprs$f_id
    }
    cds_exprs$f_id <- as.character(cds_exprs$f_id)
    cds_exprs$feature_label <- factor(cds_exprs$feature_label)
    new_data <- data.frame(pseudotime = colData(cds_subset)$pseudotime)
    model_tbl = fit_models(cds_subset, model_formula_str = trend_formula)
    model_expectation <- model_predictions(model_tbl, new_data = colData(cds_subset))
    colnames(model_expectation) <- colnames(cds_subset)
    expectation <- plyr::ddply(cds_exprs, plyr::.(f_id, Cell), 
                               function(x) {
                                   data.frame(expectation = model_expectation[x$f_id, 
                                                                              x$Cell])
                               })
    cds_exprs <- merge(cds_exprs, expectation)
    cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
    cds_exprs$expectation[cds_exprs$expectation < min_expr] <- min_expr
    if (!is.null(panel_order)) {
        cds_exprs$feature_label <- factor(cds_exprs$feature_label, 
                                          levels = panel_order)
    }
    if(return_data) return(cds_exprs)
    q <- ggplot(aes(pseudotime, expression), data = cds_exprs)
    if (!is.null(color_cells_by)) {
        q <- q + geom_point(aes_string(color = color_cells_by), 
                            size = I(cell_size), position = position_jitter(horizontal_jitter, 
                                                                            vertical_jitter))
        if (class(colData(cds_subset)[, color_cells_by]) == "numeric") {
            q <- q + viridis::scale_color_viridis(option = "C")
        }
    }
    else {
        q <- q + geom_point(size = I(cell_size), position = position_jitter(horizontal_jitter, 
                                                                            vertical_jitter))
    }
    q <- q + geom_line(aes(x = pseudotime, y = expectation), 
                       data = cds_exprs)
    q <- q + scale_y_log10() + facet_wrap(~feature_label, nrow = nrow, 
                                          ncol = ncol, scales = "free_y")
    if (min_expr < 1) {
        q <- q + expand_limits(y = c(min_expr, 1))
    }
    q <- q + ylab("Expression")
    q <- q + xlab("pseudotime")
    q <- q + monocle_theme_opts()
    q
}


my_plot_genes_in_pseudotime = function(cds_full, goi, k_cell = TRUE, max_pseudotime = Inf, min_expr = .05, cell_size = 0.75, nrow = NULL, 
                                       ncol = 1, panel_order = NULL, color_cells_by = "pseudotime", 
                                       trend_formula = "~ splines::ns(pseudotime, df=3)", label_by_short_name = TRUE, 
                                       vertical_jitter = NULL, horizontal_jitter = NULL,
                                       use_rastr = FALSE, rastr.width = NULL, rastr.height = NULL){
    
    stopifnot(all(is.character(goi) | is.factor(goi)))
    dt = as.data.table(rowData(cds_full))
    setkey(dt, gene_short_name)
    k_gene = dt[.(goi)]$id
    cds_subset = cds_full[k_gene, k_cell]
    cds_subset = cds_subset[, pseudotime(cds_subset) <= max_pseudotime]
    
    cds_exprs = .plot_genes_in_pseudotime(cds_subset, 
                                          min_expr = min_expr, cell_size = cell_size, nrow = nrow, 
                                          ncol = ncol, panel_order = panel_order, color_cells_by = color_cells_by, 
                                          trend_formula = trend_formula, label_by_short_name = label_by_short_name, 
                                          vertical_jitter = vertical_jitter, horizontal_jitter = horizontal_jitter,
                                          return_data = TRUE)
    # browser()
    colData(cds_subset)$pseudotime <- pseudotime(cds_subset)
    if(all(is.character(goi))){
        cds_exprs$feature_label = factor(cds_exprs$feature_label, levels = goi)
    }else if(is.factor(goi)){
        cds_exprs$feature_label = factor(cds_exprs$feature_label, levels = goi %>% sort %>% as.character)
    }
    q <- ggplot(aes(pseudotime, expression), data = cds_exprs)
    if (!is.null(color_cells_by)) {
        if(use_rastr){
            q <- q + ggrastr::geom_point_rast(aes_string(color = color_cells_by), 
                                              size = I(cell_size), position = position_jitter(horizontal_jitter, 
                                                                                              vertical_jitter), 
                                              raster.width = rastr.width, raster.height = rastr.height)    
        }else{
            q <- q + geom_point(aes_string(color = color_cells_by), 
                                size = I(cell_size), position = position_jitter(horizontal_jitter, 
                                                                                vertical_jitter))    
        }
        
        if (class(colData(cds_subset)[, color_cells_by]) == "numeric") {
            q <- q + viridis::scale_color_viridis(option = "C")
        }
    }
    else {
        if(use_rastr){
            q <- q + ggrastr::geom_point_rast(size = I(cell_size), 
                                              position = position_jitter(horizontal_jitter, vertical_jitter), 
                                              raster.width = rastr.width, raster.height = rastr.height)    
        }else{
            q <- q + geom_point(size = I(cell_size), 
                                position = position_jitter(horizontal_jitter, vertical_jitter))
        }
        
    }
    q <- q + geom_line(aes(x = pseudotime, y = expectation), 
                       data = cds_exprs)
    q <- q + scale_y_log10() + facet_wrap(~feature_label, nrow = nrow, 
                                          ncol = ncol, scales = "free_y")
    if (min_expr < 1) {
        q <- q + expand_limits(y = c(min_expr, 1))
    }
    q <- q + ylab("Expression")
    q <- q + xlab("pseudotime")
    q <- q + monocle3:::monocle_theme_opts()
    q
    
}

my_plot_arbitrary_in_pseudotime = function(cds_obj, arbitrary_data, flat_size_factor = TRUE,
                                           max_pseudotime = Inf, min_expr = .05, cell_size = 0.75, nrow = NULL, 
                                           ncol = 1, panel_order = NULL, color_cells_by = "pseudotime", 
                                           trend_formula = "~ splines::ns(pseudotime, df=3)", label_by_short_name = TRUE, 
                                           vertical_jitter = NULL, horizontal_jitter = NULL,
                                           use_rastr = FALSE, rastr.width = NULL, rastr.height = NULL){
    stopifnot(colnames(arbitrary_data) %in% colnames(cds_obj))
    sel_cds = cds_obj[, colnames(arbitrary_data)]
    my_plot_genes_in_pseudotime(sel_cds, "Thbs1", max_pseudotime = 10, vertical_jitter = .1)
    M = arbitrary_data
    rn = rownames(arbitrary_data)
    DF = DataFrame(id = rn, gene_short_name = rn, row.names = rn)
    
    
    
    cds_toAdd = monocle3::new_cell_data_set(expression_data = as.matrix(M), 
                                            gene_metadata = DF, 
                                            cell_metadata = colData(sel_cds))
    cds_toAdd@int_colData = sel_cds@int_colData
    colData(cds_toAdd) = colData(sel_cds)
    
    cds_modded = rbind(sel_cds, cds_toAdd)#[rn,]
    if(flat_size_factor){
        colData(cds_modded)$Size_Factor = 1
    }
    my_plot_genes_in_pseudotime(cds_modded, 
                                rn, 
                                max_pseudotime = max_pseudotime, 
                                min_expr = min(arbitrary_data), cell_size = cell_size, nrow = nrow, 
                                ncol = ncol, panel_order = panel_order, color_cells_by = color_cells_by, 
                                trend_formula = trend_formula, label_by_short_name = label_by_short_name, 
                                vertical_jitter = vertical_jitter, horizontal_jitter = horizontal_jitter,
                                use_rastr = use_rastr, rastr.width = rastr.width, rastr.height = rastr.height)
}

default_choose_plot = function(cds, reduction_method, show.legend = FALSE){
    suppressMessages(plot_cells(cds, reduction_method = reduction_method, 
                                cell_size = 1, label_cell_groups = FALSE, rasterize = FALSE) + 
                         geom_point(alpha = colData(cds)$keep, show.legend = show.legend)) + 
        theme(legend.position = "none")
}

partition_choose_plot = function(mon, not_used, show.legend = FALSE, color_by = "seurat_names"){
    plot_data = my_plot_cells(mon, return_data = TRUE)
    
    p_dt = as.data.table(plot_data$data_df)
    setnames(p_dt, c("data_dim_1", "data_dim_2"), c("UMAP_1", "UMAP_2"))
    lab_dt = p_dt[, .(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)), c(color_by)]
    #combining p_dt with lab_dt allows labels to avoid points
    append_dt = p_dt[, .(COLOR_BY = "", UMAP_1, UMAP_2) ]
    setnames(append_dt, "COLOR_BY", color_by)
    lab_dt = rbind(lab_dt, append_dt)
    
    p_dt$partitions = partitions(mon)[p_dt$sample_name]
    ggplot(p_dt, aes(x = UMAP_1, y = UMAP_2, color = partitions)) +
        geom_point(size = .3, show.legend = show.legend) +
        geom_segment(data = plot_data$edge_df, aes_string(x = "source_prin_graph_dim_1", 
                                                          y = "source_prin_graph_dim_2", 
                                                          xend = "target_prin_graph_dim_1", 
                                                          yend = "target_prin_graph_dim_2"), size = 1, 
                     color = "black", linetype = "solid", 
                     na.rm = TRUE) +
        ggrepel::geom_text_repel(data = lab_dt, aes_string(label = color_by), color = "gray60", size = 5) +
        # scale_color_manual(values = get_clusters_colors()) +
        theme_classic() +
        guides(color = guide_legend(override.aes = list(size = 3)))
    
}

seurat_cluster_choose_plot.alpha = function(mon, not_used, show.legend = FALSE, col_scale = NULL, color_by = "seurat_names"){
    seurat_cluster_choose_plot(mon, not_used, show.legend, col_scale, color_by, point_alpha = .3) +
        geom_point(alpha = colData(mon)$keep, show.legend = FALSE)
}

seurat_cluster_choose_plot = function(mon, not_used, show.legend = FALSE, 
                                      col_scale = NULL, color_by = "seurat_names", 
                                      facet_by = character(),
                                      point_alpha = 1, point_size = .3){
    plot_data = my_plot_cells(mon, return_data = TRUE)
    
    p_dt = as.data.table(plot_data$data_df)
    setnames(p_dt, c("data_dim_1", "data_dim_2"), c("UMAP_1", "UMAP_2"))
    lab_dt = p_dt[, .(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)), c(color_by, facet_by)]
    #combining p_dt with lab_dt allows labels to avoid points
    
    if(length(facet_by) > 0){
        append_dt = p_dt[, .(COLOR_BY = "", UMAP_1, UMAP_2, FACET_BY = get(facet_by)) ]  
        setnames(append_dt, "FACET_BY", facet_by)
    }else{
        append_dt = p_dt[, .(COLOR_BY = "", UMAP_1, UMAP_2) ]    
    } 
    setnames(append_dt, "COLOR_BY", color_by)
    lab_dt = rbind(lab_dt, append_dt)
    
    p = ggplot(p_dt, aes_string(x = "UMAP_1", y = "UMAP_2", color = color_by)) +
        geom_point(size = point_size, show.legend = show.legend, alpha = point_alpha) +
        geom_segment(data = plot_data$edge_df, aes_string(x = "source_prin_graph_dim_1", 
                                                          y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                                                          yend = "target_prin_graph_dim_2"), size = 1, 
                     color = "black", linetype = "solid", 
                     na.rm = TRUE) +
        ggrepel::geom_text_repel(data = lab_dt, aes_string(label = color_by), color = "gray60", size = 5) +
        
        theme_classic()+
        guides(color = guide_legend(override.aes = list(size = 3)))
    if(length(facet_by) > 0) p = p + facet_wrap(paste0("~", facet_by))
    if(is.null(col_scale)){
        if(p_dt[[color_by]] %>% unique %>% length <= 8){
            p = p + scale_color_brewer(palette = "Dark2")    
        }else{
            p = p + scale_color_discrete()
        }
    }else{
        p = p + scale_color_manual(values = col_scale)
    }
    p
}

seurat_genotype_choose_plot = function(mon, not_used, show.legend = FALSE, 
                                       col_scale = NULL, color_by = "seurat_names"){
    plot_data = my_plot_cells(mon, return_data = TRUE)
    
    p_dt = as.data.table(plot_data$data_df)
    setnames(p_dt, c("data_dim_1", "data_dim_2"), c("UMAP_1", "UMAP_2"))
    lab_dt = p_dt[, .(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)), c(color_by)]
    #combining p_dt with lab_dt allows labels to avoid points
    append_dt = p_dt[, .(COLOR_BY = "", UMAP_1, UMAP_2) ]
    setnames(append_dt, "COLOR_BY", color_by)
    lab_dt = rbind(lab_dt, append_dt)
    p = ggplot(p_dt, aes(x = UMAP_1, y = UMAP_2, color = genotype)) +
        geom_point(size = .3, show.legend = show.legend) +
        geom_segment(data = plot_data$edge_df, aes_string(x = "source_prin_graph_dim_1", 
                                                          y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                                                          yend = "target_prin_graph_dim_2"), size = 1, 
                     color = "black", linetype = "solid", 
                     na.rm = TRUE) +
        ggrepel::geom_text_repel(data = lab_dt, aes_string(label = color_by), color = "gray60", size = 5) +
        theme_classic()+
        guides(color = guide_legend(override.aes = list(size = 3)))
    if(!is.null(col_scale)){
        p = p + scale_color_manual(values = col_scale)
    }else{
        p = p + scale_color_brewer(palette = "Set1")
    }
    p
}

seurat_cell_cycle_choose_plot = function(mon, not_used, show.legend = FALSE, plot_var = "Phase", show.trajectory = TRUE, color_by = "seurat_names"){
    plot_data = my_plot_cells(mon, return_data = TRUE)
    
    p_dt = as.data.table(plot_data$data_df)
    setnames(p_dt, c("data_dim_1", "data_dim_2"), c("UMAP_1", "UMAP_2"))
    lab_dt = p_dt[, .(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)), c(color_by)]
    #combining p_dt with lab_dt allows labels to avoid points
    append_dt = p_dt[, .(COLOR_BY = "", UMAP_1, UMAP_2) ]
    setnames(append_dt, "COLOR_BY", color_by)
    lab_dt = rbind(lab_dt, append_dt)
    if(is.numeric(p_dt[[plot_var]])){
        p = ggplot(p_dt, aes_string(x = "UMAP_1", y = "UMAP_2", color = plot_var)) +
            geom_point(size = .3, show.legend = show.legend) 
        if(show.trajectory){
            p = p + geom_segment(data = plot_data$edge_df, aes_string(x = "source_prin_graph_dim_1", 
                                                                      y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                                                                      yend = "target_prin_graph_dim_2"), size = 1, 
                                 color = "black", linetype = "solid", 
                                 na.rm = TRUE)
        }
        p = p + 
            ggrepel::geom_text_repel(data = lab_dt, aes_string(label = color_by), color = "gray60", size = 5) +
            scale_color_viridis_c(option = "magma") +
            theme_classic()
    }else{
        p = ggplot(p_dt, aes_string(x = "UMAP_1", y = "UMAP_2", color = plot_var)) +
            geom_point(size = .3, show.legend = show.legend) 
        if(show.trajectory){
            p = p + geom_segment(data = plot_data$edge_df, aes_string(x = "source_prin_graph_dim_1", 
                                                                      y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                                                                      yend = "target_prin_graph_dim_2"), size = 1, 
                                 color = "black", linetype = "solid", 
                                 na.rm = TRUE)
        }
        p = p + 
            ggrepel::geom_text_repel(data = lab_dt, aes_string(label = color_by), color = "gray60", size = 5) +
            scale_color_brewer(palette = "Dark2") +
            theme_classic()
    }
    p + guides(color = guide_legend(override.aes = list(size = 3)))
    
}

my_choose_cells = function (cds, reduction_method = c("UMAP", "tSNE", "PCA", "Aligned"), 
                            return_list = FALSE, plot_FUN = default_choose_plot) 
{
    reduction_method <- match.arg(reduction_method)
    assertthat::assert_that(methods::is(cds, "cell_data_set"))
    assertthat::assert_that(!is.null(reducedDims(cds)[[reduction_method]]), 
                            msg = paste0("No dimensionality reduction for ", reduction_method, 
                                         " calculated. ", "Please run reduce_dimensions with ", 
                                         "reduction_method = ", reduction_method, ", cluster_cells, and learn_graph ", 
                                         "before running choose_cells"))
    assertthat::assert_that(is.logical(return_list))
    assertthat::assert_that(interactive(), msg = paste("choose_cells only works in", 
                                                       "interactive mode."))
    reduced_dims <- as.data.frame(reducedDims(cds)[[reduction_method]])
    names(reduced_dims)[1:2] <- c("V1", "V2")
    ui <- shiny::fluidPage(shiny::titlePanel("Choose cells for a subset"), 
                           shiny::sidebarLayout(shiny::sidebarPanel(shiny::actionButton("choose_toggle", 
                                                                                        "Choose/unchoose"), shiny::actionButton("reset", 
                                                                                                                                "Clear"), shiny::actionButton("done", "Done"), shiny::h3("Instructions:"), 
                                                                    shiny::tags$ol(shiny::tags$li("Highlight points by clicking and dragging."), 
                                                                                   shiny::tags$li("Click the 'Choose/unchoose' button."), 
                                                                                   shiny::tags$li("Repeat until all of the desired cells are black."), 
                                                                                   shiny::tags$li("Click 'Done'.")), shiny::h4("Details:"), 
                                                                    shiny::tags$ul(shiny::tags$li("To start over, click 'Clear'"), 
                                                                                   shiny::tags$li(paste("You can also choose/unchoose specific cells", 
                                                                                                        "by clicking on them directly")))), shiny::mainPanel(shiny::plotOutput("plot1", 
                                                                                                                                                                               height = "auto", click = "plot1_click", brush = shiny::brushOpts(id = "plot1_brush")))))
    server <- function(input, output, session) {
        vals <- shiny::reactiveValues(keeprows = rep(FALSE, nrow(colData(cds))))
        output$plot1 <- shiny::renderPlot({
            colData(cds)$keep = vals$keeprows
            plot_FUN(cds, reduction_method)
        }, height = function() {
            session$clientData$output_plot1_width
        })
        shiny::observeEvent(input$plot1_click, {
            res <- shiny::nearPoints(reduced_dims, xvar = "V1", 
                                     yvar = "V2", input$plot1_click, allRows = TRUE)
            vals$keeprows <- vals$keeprows | res$selected_
        })
        shiny::observeEvent(input$choose_toggle, {
            res <- shiny::brushedPoints(reduced_dims, xvar = "V1", 
                                        yvar = "V2", input$plot1_brush, allRows = TRUE)
            vals$keeprows <- vals$keeprows | res$selected_
        })
        shiny::observeEvent(input$reset, {
            vals$keeprows <- rep(FALSE, nrow(colData(cds)))
        })
        shiny::observeEvent(input$done, {
            shiny::stopApp(vals$keeprows)
        })
    }
    sel <- suppressMessages(shiny::runApp(shiny::shinyApp(ui, 
                                                          server)))
    if (return_list) {
        return(row.names(colData(cds)[sel, ]))
    }
    else {
        return(cds[, sel])
    }
}    

barplot_cnt_per_seurat = function(full_mon, sel_mon, res_file = NULL, out_csv = "selection_count_per_cluster.csv"){
    meta_dt = as.data.table(pData(full_mon))
    meta_dt$id = colnames(full_mon)
    meta_dt[, in_selection := id %in% colnames(sel_mon)]
    cnt_dt = meta_dt[, .(.N), .(seurat_names, seurat_clusters, in_selection)][order(in_selection)][order(seurat_clusters)]
    cnt_dt = dcast(cnt_dt, seurat_names+seurat_clusters~in_selection, fill = 0, value.var = "N")
    if(is.null(cnt_dt[["FALSE"]])) cnt_dt[["FALSE"]] = 0
    setnames(cnt_dt, c("FALSE", "TRUE"), c("not_selected", "in_selection"))
    if(!is.null(res_file)){
        fwrite(cnt_dt, res_file(out_csv))    
    }
    
    p_sel_per_cluster = ggplot(meta_dt, aes(x = seurat_names, fill = in_selection)) +
        geom_bar() +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
        scale_fill_manual(values = c("FALSE" = "gray", 'TRUE' = "blue")) +
        labs(title = "Selection breakdown by Seurat cluster")
    p_sel_per_cluster
}

