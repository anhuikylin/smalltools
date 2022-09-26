# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' Title
#' @title Row name of a col
#' @description Row name of a col
#' @param data    Numerical matrix
#' @param col    col
#' @import ggplot2
#' @import PerformanceAnalytics
#' @import gclus
#' @return data
#' @export col_to_rowname
#' @export qc_cor_plot
#'
#' @examples
#' mtcars <- data.frame(rownames(mtcars),mtcars)
#' col_to_rowname(mtcars,1)
col_to_rowname <- function(data,col) {
  data <- data.frame(data)
  rownames(data) <- data[,col]
  data <- data[,-col]
  return(data)
}
