
save_bw_rda <- function(bw_files_dir){
         
          library(usethis)
          
          bw_files <- list.files(bw_files_dir, pattern = paste(c("*.bw", "*.bedgraph"),collapse = '|'), recursive = T, full.names = T)
          
          names(bw_files) <- gsub(pattern = paste(c("_[[:upper:]]{6,}_.*_normalized.bw*", "_gencov_normalized.bedgraph*"),collapse = '|'),replacement = "", basename(bw_files))
          bw_files <- broom::tidy(bw_files)
          print(bw_files)
          
          xx <- bw_files %>%
                    dplyr::mutate(bw = purrr::map(x, function(ii) {
                              rtracklayer::import(ii)
                    }))
          
          
          mylist <- xx$bw
          
          names(mylist) <- xx$names
          
          purrr::walk2(mylist, names(mylist), function(obj, name) {
                    assign(name, obj)
                    do.call("use_data", list(as.name(name), overwrite=TRUE, compress = "bzip2"))
          })
          
}



