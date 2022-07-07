#' @title plot_coeff_xwas
#' @description  this is a function that plots the coefficients of the XWAS results

#' @author Yixuan He, \email{yixuan_he@@hms.harvard.edu}

#' @param xdff  matrix returned from XWAS function, rownames of matrix should be the X variables
#' @param pval column name of p value
#' @param coeff column name of coefficients
#' @param thresh p value threshold for signficance
#' @param all default is to plot only signficant associaitons, all=TRUE plots all associatons

#' @export
#'
#' @importFrom magrittr %>%
#' @import naniar
#' @import glmnet
#' @import broom
#' @import ggplot2

plot_coeff_xwas = function(xdff, pval, coeff, thresh = 0.05,all=FALSE,ang=90,va=0.5,ha=1) {

  if(all==TRUE){
    xdf=xdff
  }
  if(all==FALSE){
    xdf = xdff[which(xdff[, pval] < thresh), ]
  }
  colnames(xdf)[which(colnames(xdf) == coeff)] = 'COEFF'
  gg <- ggplot2::ggplot(xdf, ggplot2::aes(x = row.names(xdf), y = COEFF)) +
    ggplot2::geom_point(ggplot2::aes(color = 1:nrow(xdf)), size = 2) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::ylab('Effect Size') +
    ggplot2::xlab('Exposure') +
    ggplot2::theme(axis.text.x = ggplot2::element_text(
      angle = ang,
      vjust = va,
      hjust = ha
    ))
  return(gg)
}
