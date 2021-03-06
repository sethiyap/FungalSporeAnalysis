% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/save_bw_to_rdata.R
\name{save_bw_rda}
\alias{save_bw_rda}
\title{save bw files as GRanges in .RObject}
\usage{
save_bw_rda(bw_files_dir, pattern = NULL)
}
\arguments{
\item{bw_files_dir}{Absolute path of bw or bdg files. Files should end with
extension \code{ .bw or .bedgraph}}

\item{pattern}{A vector of pattern to replace from the bw or bedgraph file
names. Suitable for long names, can be either a defined pattern or
"default" pattern or NULL. \code{default:NULL}}
}
\value{
.rda for each bw/ bdg file input
}
\description{
save bw files as GRanges in .RObject
}
\examples{
\dontrun{

save_bw_rda(bw_files_dir=".", pattern="default")

}
}
