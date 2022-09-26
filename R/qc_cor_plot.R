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

#' @title Correlation analysis
#' @description Correlation analysis
#' @param data Numerical matrix
#' @param col col
#' @import ggplot2
#' @import PerformanceAnalytics
#' @import gclus
#' @return data
#' @export one_to_rowname
#' @export qc_cor_plot
#'
#' @examples
#' qc_cor_plot(mtcars)
one_to_rowname <- function(data,col) {
  data <- data.frame(data)
  rownames(data) <- data[,col]
  data <- data[,-col]
  return(data)
}

panel.cor <- function(x, y, digits = 3, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor ,col = "#020260")
}

qc_cor_plot  <- function(data) {
  cpairs(scale(data),gap=1,
         pch = 21,#点的类型
         upper.panel = panel.cor,
         font.labels = 1,#字体类型
         cex.labels = 2,
         show.points = TRUE,
         border.color = 'black',
         cex = 2,
         row1attop=TRUE,
         col='#020260')#点的颜色)
}
