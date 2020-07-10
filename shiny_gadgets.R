
library(shiny)
library(miniUI)
library(ggplot2)

shiny_overwrite <- function(msg = "Processed data found, overwrite?", yes_prompt = "Yes, overwrite.", no_prompt = "No, use existing.") {
    
    ui <- miniPage(
        gadgetTitleBar(msg, left = NULL),
        radioButtons("overwrite", label = c(""), 
                     selected = "no",
                     choiceValues = c("yes", "no"), 
                     choiceNames = c(yes_prompt, no_prompt))
    )
    
    server <- function(input, output, session) {
        # Handle the Done button being pressed.
        observeEvent(input$done, {
            # Return the brushed points. See ?shiny::brushedPoints.
            stopApp(input$overwrite == "yes")
        })
    }
    runGadget(ui, server, viewer = dialogViewer("overwrite"))    
}

shiny_choose_branch <- function(output_path = "output", 
                                processed_pattern = paste0("^", prefix_processed), 
                                anlaysis_pattern = paste0("^", prefix_subset),
                                allow_new = TRUE) {
    proc_dirs = dir(output_path, pattern = processed_pattern)
    
    add_proc_prefix = function(x){paste0(prefix_processed, ".", x)}
    rm_proc_prefix = function(x){sub(paste0(prefix_processed, "\\."), "", x)}
    
    
    add_subset_prefix = function(x){paste0(prefix_subset, ".", x)}
    rm_subset_prefix = function(x){sub(paste0(prefix_subset, "\\."), "", x)}
    
    ui <- miniPage(
        gadgetTitleBar("Select processed and subset", left = NULL, right = NULL),
        selectInput("select_processed", label = "Select Processed", choices = rm_proc_prefix(proc_dirs)),
        tags$br(),
        uiOutput("dynamic_subset"),
        actionButton("btnExisting", label = "Use Existing"),
        if(allow_new){
            tags$span(
                tags$h3("OR"),
                textInput("txt_new_subset", label = "New subset name", value = "full"),
                tags$h3("\"full\" is a special value and will skip the selection step."),
                actionButton("btnNew", label = "Use New")
            )
        }else{
            tags$span()
        }
        
    )
    
    server <- function(input, output, session) {
        
        output$dynamic_subset = renderUI({
            proc_dir = file.path(output_path, add_proc_prefix(input$select_processed))
            available_subsets = dir(proc_dir, pattern = anlaysis_pattern)
            selectInput("select_subset", label = "Select Subset", choices = rm_subset_prefix(available_subsets))
        })
        
        observeEvent(input$btnExisting, {
            # Return the brushed points. See ?shiny::brushedPoints.
            stopApp(file.path(output_path, add_proc_prefix(input$select_processed), add_subset_prefix(input$select_subset)))
        })
        
        # Handle the Done button being pressed.
        observeEvent(input$btnNew, {
            # Return the brushed points. See ?shiny::brushedPoints.
            stopApp(file.path(output_path, add_proc_prefix(input$select_processed), add_subset_prefix(tolower(input$txt_new_subset))))
        })
        
        
    }
    runGadget(ui, server, viewer = dialogViewer("subset"))    
}

shiny_new_name = function(starting_value = "", output_path = "output", processed_pattern = paste0("^", prefix_processed)){
    stopifnot(dir.exists(output_path))
    proc_dirs = dir(output_path, pattern = processed_pattern)
    proc_names = sub(paste0(prefix_processed, "."), "", proc_dirs)
    
    add_proc_prefix = function(x){paste0(prefix_processed, ".", x)}
    rm_proc_prefix = function(x){sub(paste0(prefix_processed, "\\."), "", x)}
    
    ui <- miniPage(
        shinyjs::useShinyjs(),
        gadgetTitleBar("Create new processed dataset", left = NULL, right = NULL),
        tags$br(),
        textInput("txt_new_processed", label = "New processed name", value = starting_value),
        shinyjs::hidden(tags$h3("dataset must not exist, change the name.", style="color:red;", id = "warn_exist")),
        actionButton("btnNew", label = "Use New")
        
    )
    
    server <- function(input, output, session) {
        
        observeEvent({
            input$txt_new_processed
        }, {
            if(input$txt_new_processed %in% proc_names){
                shinyjs::disable("btnNew")
                shinyjs::show("warn_exist")
                
            }else{
                shinyjs::enable("btnNew")
                shinyjs::hide("warn_exist")
            }
        })
        # Handle the Done button being pressed.
        observeEvent(input$btnNew, {
            if(!input$txt_new_processed %in% proc_names){
                stopApp(input$txt_new_processed)
            }
            # Return the brushed points. See ?shiny::brushedPoints.
            
        })
        
        
    }
    runGadget(ui, server, viewer = dialogViewer("new name"))    
}

processed_plots = function(mon, out_dir){
    #there are some annoying monocle3 plotting quirks, this is my way around them
    plot_data = my_plot_cells(mon, return_data = TRUE)
    p_dt = as.data.table(plot_data$data_df)
    seg_dt = as.data.table(plot_data$edge_df)
    setnames(p_dt, c("data_dim_1", "data_dim_2"), c("UMAP_1", "UMAP_2"))
    lab_dt = p_dt[, .(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)), .(seurat_names, genotype)]
    #combining p_dt with lab_dt allows labels to avoid points
    lab_dt = rbind(lab_dt, p_dt[, .(seurat_names = "", UMAP_1, UMAP_2, genotype) ])
    p_dt = p_dt[order(seurat_names)]
    p_dt$genotype = factor(p_dt$genotype, levels = names(get_treatments_colors()))
    p_seurat = ggplot(p_dt, aes(x = UMAP_1, y = UMAP_2, color = seurat_names)) +
        annotate("point", x = p_dt$UMAP_1, y = p_dt$UMAP_2, color = "gray80", size = .2) +
        geom_point(size = .6) +
        geom_segment(data = seg_dt, 
                     aes(x = source_prin_graph_dim_1, 
                         y = source_prin_graph_dim_2, 
                         xend = target_prin_graph_dim_1, 
                         yend = target_prin_graph_dim_2), 
                     color= "black", size = 1) +
        ggrepel::geom_text_repel(data = lab_dt, aes(label = seurat_names), color = "black", size = 5) +
        theme_classic() +
        guides(color = guide_legend(override.aes = list(size = 3))) +
        scale_color_manual(values = get_clusters_colors()) +
        labs(title = "Monocle processed data") +
        facet_wrap(~genotype)
    plot(p_seurat)
    
    ggsave(file.path(out_dir, "plot_processed_seurat_clusters.png"), p_seurat, 
           width = 10, height = 5)
    ggsave(file.path(out_dir, "plot_processed_seurat_clusters.pdf"), p_seurat, 
           width = 10, height = 5)
    
    
    
    plots = list(
        seurat_cluster_choose_plot(mon, show.legend = TRUE, col_scale = get_clusters_colors(), color_by = "seurat_names") + labs(title = "Seurat clusters"),
        seurat_genotype_choose_plot(mon, show.legend = TRUE, col_scale = get_treatments_colors()) + labs(title = "Genotype"),
        partition_choose_plot(mon, show.legend = TRUE) + labs(title = "Monocle partition"),
        seurat_cell_cycle_choose_plot(mon, show.legend = TRUE)  + labs(title = "Cell cycle")
    )
    pg = cowplot::plot_grid(plotlist = plots)
    
    ggsave(file.path(out_dir, "plot_procesed_overview.png"), pg, width = 8*1.5, height = 6*1.5)
    ggsave(file.path(out_dir, "plot_procesed_overview.pdf"), pg, width = 8*1.5, height = 6*1.5)
    
    pg_score = cowplot::plot_grid(ncol = 1,
                                  seurat_cell_cycle_choose_plot(mon, show.legend = TRUE, plot_var = "S.Score", show.trajectory = FALSE)  + 
                                      labs(title = "S") +
                                      facet_wrap(~genotype) +
                                      theme(panel.background = element_rect(fill = "gray30")),
                                  seurat_cell_cycle_choose_plot(mon, show.legend = TRUE, plot_var = "G2M.Score", show.trajectory = FALSE)  + 
                                      labs(title = "G2M") +
                                      facet_wrap(~genotype)+
                                      theme(panel.background = element_rect(fill = "gray30"))
    )
    ggsave(file.path(out_dir, "plot_processed_cell_cycle_scores.png"), pg_score, width = 6.5*1.5, height = 6*1.5)
    ggsave(file.path(out_dir, "plot_processed_cell_cycle_scores.pdf"), pg_score, width = 6.5*1.5, height = 6*1.5)
    
    red_line_y = .05
    p_var = plot_pc_variance_explained(mon) +
        annotate("line", x= c(0, num_pca), y = rep(red_line_y, 2), color = "red")
    
    ggsave(file.path(out_dir, "plot_processed_PCvariance.png"), p_var, width = 4, height = 4)
    ggsave(file.path(out_dir, "plot_processed_PCvariance.pdf"), p_var, width = 4, height = 4)
}
