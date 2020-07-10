source("Step0_setup.R")
proc_dir = dir(output_root_dir, pattern = prefix_processed, full.names = TRUE)
sel_dir = sapply(proc_dir, function(x)dir(x, pattern = prefix_subset, full.names = TRUE)) %>% unlist
sel_file = file.path(sel_dir, file_subset)
sel_file = sel_file[file.exists(sel_file)]
# sel_file = sel_file[seq_along(sel_file) >= 9]
# sel_file = sel_file[11]

# undebug(run_step3)
# monocle_selected_rds = sel_file[1]
# force = TRUE
for(f in sel_file){
    message(f)
    run_step3(monocle_selected_rds = f, force = TRUE)
}

monocle_selected_rds = f
force = TRUE

message(which(sel_file == f), " / ", length(sel_file))
