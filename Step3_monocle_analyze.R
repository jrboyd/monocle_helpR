#after selecton, find genes correlated to pseudotime and generate modules
#should be run from monocle directory
source("Step0_setup.R")

use_loaded = FALSE
if(exists("monocle_selected_rds")){
   use_loaded = shiny_overwrite(msg = "Use loaded selection or select again?", yes_prompt = "Use loaded.", no_prompt = "Select again") 
   message(use_loaded)
}
Sys.sleep(1)
if(!use_loaded){
    Sys.sleep(1)
    message("select branch")
    monocle_selected_rds = file.path(shiny_choose_branch(output_path = output_root_dir, allow_new = FALSE), file_subset)
}else{
    message("skip select branch")
}

run_step3(monocle_selected_rds = monocle_selected_rds, 
          max_q = max_q, 
          min_morans = min_morans, 
          force = TRUE)


