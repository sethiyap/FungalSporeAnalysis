
#' save bw files as GRanges in .RObject
#'
#' @param bw_files_dir Absolute path of bw or bdg files. Files should end with
#'   extension \code{ .bw or .bedgraph}
#' @param pattern A vector of pattern to replace from the bw or bedgraph file
#'   names. Suitable for long names, can be either a defined pattern or
#'   "default" pattern or NULL. \code{default:NULL}
#'
#' @return .rda for each bw/ bdg file input
#' @export
#' @import usethis
#' @examples
#' \dontrun{
#' 
#' save_bw_rda(bw_files_dir=".", pattern="default")
#' 
#' }
save_bw_rda <- function(bw_files_dir, pattern=NULL){
          
          bw_files <- list.files(bw_files_dir, pattern = paste(c("*.bw", "*.bedgraph"),collapse = '|'), recursive = T, full.names = T)
          
          if(is.null(pattern)==FALSE){
                    names(bw_files) <- gsub(pattern = pattern,replacement = "", basename(bw_files))
          }
          if(pattern=="default"){
                    names(bw_files) <- gsub(pattern = paste(c("_[[:upper:]]{6,}_.*_normalized.bw*", "_gencov_normalized.bedgraph*"),collapse = '|'),replacement = "", basename(bw_files))
          }
          else{
                    names(bw_files) <- basename(bw_files) 
          }
          
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



