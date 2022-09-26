#' Title
#' @title Draw point line graph
#' @description Draw point line graph
#' @param x    a numerical vector.
#' @param y    a numerical vector.
#' @param xlab    x axis title
#' @param ylab    y axis title
#' @import ggplot2
#' @return data
#' @export line_plot
#'
#' @examples
#' line_plot(1:10,1:10)
line_plot <- function(x,y,xlab='x',ylab='y') {
  data <- data.frame(x,y)
  colnames(data) <- c("x","y")
  ggplot2::ggplot(data,aes(x = x,y = y))+
    geom_point(col = 'black')+
    geom_line(col = 'black')+
    labs(x=xlab,y=ylab)
}
