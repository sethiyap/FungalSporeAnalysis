#' GO terms in diamond plot
#' @author pooja sethiya
#' @description Visualize GO enrichment result in a diamond plot
#' @param data A tibble, with following columns \itemize{\item term: character,
#'   GO term \item annotated: numeric, number of genes enriched in the GO-term
#'   \item total: numeric, total number of genes assigned to the GO-term \item
#'   pvalue: numeric, pvalue of enrichment \item category: character, a
#'   group-name to annotate the GO and color the diamonds \item class:
#'   character, optional if there are two categories, second category should be
#'   under column class.}
#' @param palette Character, to choose from available palette or define a
#'   palette.\code{default:TRUE} \itemize{ \item two_color: a palette of pink and green
#' \item three_color:  a palette of  orange, pink and green
#' \item NULL: user input of character vector of colors
#' \item TRUE: a palette of colors in rainbow
#'   }
#' @param output_name A character vector containing name of the output
#'   plot.\code{default:Sample}
#'
#' @return a diamond GO plot and image file of same.
#' @export
#'
#' @examples
#' \dontrun{
#'
#' data <- readr::read_delim("data/an_af_pm_spores_GO.txt", col_names = TRUE, delim="\t")
#' GO_diamondplot(data, output_name = "plots/An_Af_Pm_GO", palette="three_color")
#'
#' }
GO_diamondplot <- function(data, palette=TRUE, output_name){

          print(head(data))

          data_melt <- data %>%
                        dplyr::mutate(percent_enrichment=100*(data$total/data$annotated),
                                      term=factor(term, levels=rev(unique(term))),
                                      category=factor(category, levels=unique(category)))

          if(palette %in% "two_color"){
                    values=c("#c51b7d","#01665e")
          }
          if(palette %in% "three_color"){
                    values=c("#e41a1c","#1b9e77","#f0027f")
          }
          if(isTRUE(palette)){

                    values=rainbow(length(unique(data$class)))

          }
          if(is.null(palette)){
                    color_vector <- function()
                    {
                              n <- readline(prompt="Enter vector of colors (equivalent to categorical values): ")
                              n <- as.character(n)
                              if (is.na(n)){
                                        n <- color_vector()
                              }
                              return(n)
                    }
                    values=color_vector()
          }

          if("class" %in% colnames(data)){

                    data_melt <- data_melt %>%
                              dplyr::mutate(class=factor(class, levels=unique(class)))

                    gg= ggplot(data_melt , aes(y=term, x=class,ordered=TRUE))+xlab("")+ylab("")+
                              geom_point(aes(size=percent_enrichment),color="black",shape=23, stroke=0.8)+
                              geom_point(aes(size=percent_enrichment,alpha=-log10(pvalue), fill=category),shape=23)+
                              scale_fill_manual(values=values)
                              #scale_fill_manual(values=c("#e41a1c","#1b9e77","#f0027f"))


                    }
                    else{
                              gg= ggplot(data_melt , aes(y=term, x=category,ordered=TRUE))+
                                        geom_point(aes(size=percent_enrichment),color="black",shape=23, stroke=0.8)+
                                        geom_point(aes(size=percent_enrichment,alpha=-log10(pvalue), fill=category),shape=23)+
                                        scale_fill_manual(values=values)

                    }

                    # --- ggplot theme
                    gg_theme <- function(gg){
                              gg_themed <- gg+ ggplot2::theme_bw()+
                                                  ggplot2::labs(x="", y="")+
                                                  ggplot2::scale_size(range = c(1,15))+
                                                  ggplot2::scale_y_discrete(position="right")+
                                                  ggplot2::scale_x_discrete(position="top")+
                                                  ggplot2::theme(axis.text.x= ggplot2::element_text(color="black",size=10),
                                                                  axis.text.y = ggplot2::element_text(color="black",size=10),
                                                                  axis.title.x=ggplot2::element_text(color="black",size=10),
                                                                  legend.title=ggplot2::element_text(color="black",size=10),
                                                                  legend.key.size = ggplot2::unit(0.5,"line"),
                                                                  panel.spacing = ggplot2::unit(0.1, "lines"),
                                                                  legend.text=ggplot2::element_text(color="black",size=10))+
                                                  ggplot2::scale_color_manual(values ="black")+
                                                  ggplot2::guides(fill = ggplot2::guide_legend(title="",override.aes = list(size=8)),
                                                            size=ggplot2::guide_legend(title="percentage of \ngene over bgd",
                                                                                       override.aes = list(size=c(2,4,6,8))),
                                                            alpha=ggplot2::guide_legend(override.aes=list(size=5)))

                              return(gg_themed)
                    }



                    ggsave(plot = gg_theme(gg),filename = paste(output_name,"GoByColorPvalue.png",sep=""), height = 7, width=7, device = "png", units = "in", dpi = 300, pointsize=10)
                    return(gg_theme(gg))
}


