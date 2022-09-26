#' Title
#' @title Metabolic classification statistics
#' @description Metabolic classification statistics
#' @param classlist    a vector to class
#' @param fontsize    font size
#' @param labelsize     label size
#' @return b
#' @export class_bar_plot
#'
#' @examples
#' a <- rep(LETTERS[1:10],sample(1:10, 10, replace = TRUE))
#' class_bar_plot(a)
#' class_bar_plot(classlist = mtcars$cyl,fontsize = 20,labelsize = 7)
class_bar_plot <- function(classlist,fontsize = 20,labelsize = 7) {
  table(classlist)
  b <- data.frame(table(classlist))
  Class <- 'Class'
  colnames(b)[1] <- Class
  ggplot2::ggplot(b,aes(x = Freq ,y = Class))+
    geom_bar(stat = 'identity',
             col='black',
             fill='#159dd9',
             width = 0.4,
             show.legend = F)+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    theme(text = element_text(size = fontsize))+
    geom_text(aes(label=Freq,col='black'),size=labelsize,
              vjust=0.4,hjust=-0.5,col='black')+
    xlim(c(0,max(b$Freq)*1.25))
}
