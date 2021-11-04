#' @title manhattan_xwas
#' @description  this is a function plot the p values of the XWAS results, analagous to a manhattan plot.
#'   since the y axis is in the -log scale, there may be issues with plotting if the p value is zero or very close to zero (taking the neg log of it will be infinite)

#' @author Yixuan He, \email{yixuan_he@@hms.harvard.edu}

#' @param xdff  matrix returned from XWAS function, rownames of matrix should be the X variables
#' @param pval column name of p value
#' @param thresh p value threshold for signficance

#' @export
#'
#' @importFrom magrittr %>%
#' @import naniar
#' @import glmnet
#' @import broom
#' @import ggplot2

manhattan_xwas = function(xdff, pval, thresh = 0.05) {

  xdff$manhp = -log10(unlist(xdff[which(colnames(xdff) == pval)]))
  gg <- ggplot2::ggplot(xdff, ggplot2::aes(x = row.names(xdff), y = manhp)) +
    ggplot2::geom_point(ggplot2::aes(color = 1:nrow(xdff)), size = 2) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::ylab('p-value (-log10)') +
    ggplot2::xlab('Exposure') +
    ggplot2::geom_hline(yintercept = -log10(thresh),
                        linetype = "dashed",
                        color = "blue") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ))
  return(gg)
}
