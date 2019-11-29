

# genes   spore hypha input class
# cetA    125.   3.07  9.77 CET
# cetL    112.   2.82  5.38 CET
# tpsA     47.1  9.79  5.35 Trehalose_biosynthesis
# orlA     21.7  8.27  4.14 Trehalose_biosynthesis



#' heatmap of gene expressions by class
#' @author pooja sethiya
#' @description Heatmap of genes expression value categorized by classes
#' @param input_data A tibble, with four columns. First column name should be
#'   \code{genes} and second column name should be \code{class}. Other columns
#'   can have unique names to plot.
#' @param threshold Logical, to set maximum limit of the color
#'   key.\code{default:FALSE}. If \code{TRUE} on-run threshold value should be
#'   provided.
#' @param output_name A character vector, containing name of the output
#'   plot.\code{default:Sample}
#'
#' @return heatmap of expression of given genes categorized by class and output
#'   image file of same.
#' @export
#' @import ggplot2
#' @import magrittr
#' @importFrom tidyr gather
#' @importFrom dplyr mutate
#' @importFrom dplyr arrange
#' @importFrom dplyr if_else
#' @importFrom dplyr group_by
#'
#' @examples
#' \dontrun{
#'
#' dat <- readr::read_delim(system.file("extdata/genesets/an_spore_hypha_specificgenes.txt" , package = "FungalSporeAnalysis"), delim="\t", col_names = TRUE)
#' input_data <- dat %>% dplyr::filter(class=="spore_maturation")
#' ggplot_heatmap(input_data ,threshold = TRUE, output_name = "plots/An_spores_maturation_genes_exprsn")
#'
#' }
ggplot_heatmap <- function(input_data , threshold=FALSE, output_name="Sample"){


          dat_long <- input_data  %>%
                              tidyr::gather(columns, value, -genes, -class) %>%
                              dplyr::mutate(columns=factor(columns, levels = unique(columns))) %>%
                              dplyr::group_by(columns) %>%
                              dplyr::arrange(value) %>%
                              dplyr::mutate(genes=factor(genes, levels = unique(genes)))

          if(threshold==TRUE){

                    readnumber <- function()
                    {
                              n <- readline(prompt="Enter threshold to limit the color scale: ")
                              n <- as.double(n)
                              if (is.na(n)){
                                        n <- readnumber()
                              }
                              return(n)
                    }

                    read_num <- readnumber()
                    dat_long <- dat_long %>%
                                        dplyr::mutate(value=dplyr::if_else(value >= read_num,read_num , value))
          }
          else{
                    dat_long <- dat_long
          }

          gp <- dat_long %>% ggplot2::ggplot(ggplot2::aes(columns, genes, fill=value))+
                                        ggplot2::geom_tile(color="white", size=1)+
                                        ggplot2::scale_fill_gradient(low = "#ffeda0", high="#f03b20")+
                                        ggplot2::facet_wrap(~class, strip.position = "left", ncol = 1, scales = "free")


          # --- ggplot theme
          gg_theme <- function(gg){
                    gg_themed <- gg+ ggplot2::theme_bw()+
                                                  ggplot2::scale_y_discrete(position="right")+
                                                  ggplot2::scale_x_discrete(position="top")+
                                                  ggplot2::theme(axis.text.x= ggplot2::element_text(color="black",size=12,angle=0,vjust=0.8),
                                                        axis.text.y = ggplot2::element_text(color="black",size=12),
                                                        axis.title=ggplot2::element_blank(),
                                                        legend.position = "bottom",
                                                        legend.title=ggplot2::element_text(color="black",size=10),
                                                        legend.key.size = ggplot2::unit(1,"line"),
                                                        legend.box.spacing=ggplot2::unit(0.5,"mm"),
                                                        legend.text=ggplot2::element_text(color="black",size=8))+
                                                  ggplot2::guides(fill = ggplot2::guide_colourbar(title.position="top", title.hjust = 0.5,title ="RNAP Signal"))
                    return(gg_themed)
          }

          ggplot2::ggsave(gg_theme(gp),filename = paste(output_name, "ggplot_hm.png"), device="png", dpi = 300, width = 7, unit="cm",height=length(unique(input_data$class))*6,  pointsize=8)

          return(gg_theme(gp))


}
