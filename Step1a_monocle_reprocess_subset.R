#this takes a subsetted monocle object and initializes it as a new processed monocle object.
#this procedure increases local resolution of trajectories
source("Step0_setup.R")

select_dir = shiny_choose_branch(output_path = output_root_dir, allow_new = FALSE)
select_name = sub(paste0(prefix_subset, "."), "", basename(select_dir))
select_file = file.path(select_dir, file_subset)
stopifnot(file.exists(select_file))
Sys.sleep(1) #Sys.sleep is required to run multiple shiny gadgets in a row
new_proc_name = shiny_new_name(select_name, output_path = output_root_dir)
new_proc_dir = get_proc_dir(new_proc_name)
dir.create(new_proc_dir)
new_proc_file = file.path(new_proc_dir, file_processed)


run_process = TRUE
if(file.exists(default_monocle_processed_rds)){
    Sys.sleep(1)
    run_process = shiny_overwrite()
    Sys.sleep(1)
    message(run_process)
}
Sys.sleep(3)
if(run_process){
    message("Processing data")
    mon = readRDS(select_file)
    
    mon <- preprocess_cds(mon, num_dim = num_pca)
    if(align_apply){
        mon <- align_cds(mon, alignment_group = "genotype")    
    }
    mon = reduce_dimension(mon, reduction_method = "UMAP", cores = 20)
    mon = cluster_cells(mon)
    #minimal_branch_len is quite important and controls the detail of branches
    mon <- learn_graph(mon, learn_graph_control = list(minimal_branch_len = minimal_branch_len), use_partition = TRUE)
    plot_cells(mon)
    saveRDS(mon, new_proc_file)
}else{
    mon = readRDS(new_proc_file)
}
#same plots as normal step 1
processed_plots(mon, new_proc_dir)
