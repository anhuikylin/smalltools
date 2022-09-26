
#' @title Metabonomics analysis
#' @description Metabonomics analysis
#' @param a    a matrix
#' @param group    group
#' @return b
#' @import ggplot2
#' @import ggrepel
#' @import RColorBrewer
#' @import ggpubr
#' @import ropls
#' @import FactoMineR
#' @import factoextra
#' @import PerformanceAnalytics
#' @import gclus
#' @import corrplot
#' @import eoffice
#' @export metabonomics_analysis
#'
#' @examples
#' metabonomics_analysis(a,group)
metabonomics_analysis <- function(a,group) {
  for (i in round(1:dim(group)[1])) {
    group_a <- rbind(a[which(a$Class_ID==as.character(group[i,2])),],
                     a[which(a$Class_ID==as.character(group[i,3])),])
    col_num <- dim(group_a)[1]
    row_num <- dim(group_a)[2]

    group_pca_data <- rbind(a[which(a$Class_ID==as.character(group[i,2])),],
                            a[which(a$Class_ID==as.character(group[i,3])),])
    group_pca_lable <- group_pca_data$Class_ID


    rownames(group_pca_data)=group_pca_data[,1]
    group_pca_data <- group_pca_data[,-c(1:2)]
    group_pca_data <- data.frame(group_pca_data)
    group_pca_variables <- PCA(group_pca_data,  graph = F)
    p3 <- factoextra::fviz_pca_ind(group_pca_variables,
                       geom.ind = c("point", "text"),
                       habillage = factor(group_pca_lable),
                       palette = c("#00AFBB", "#E7B800"),
                       addEllipses = TRUE,label = "all",
                       labelsize = 8, pointsize = 3,

    )+
      theme_bw()+
      labs(title = '',
           x=paste0('PC1','(',round(group_pca_variables[["eig"]][1,2],2),'%)'),
           y=paste0('PC2','(',round(group_pca_variables[["eig"]][2,2],2),'%)'))+
      theme(plot.title = element_text(hjust = 0.5),
            text = element_text(size = 25))

    ggsave(file = paste0(as.character(group[i,2]),'_vs_',
                      as.character(group[i,3]),
                      '_Contribution_rate_of_variance.pdf'),
        width = 12,height = 10)
    topptx(p3,paste0(as.character(group[i,2]),'_vs_',
                     as.character(group[i,3]),
                     '_Contribution_rate_of_variance.pptx'),
           width = 12,
           height = 10)

  }

}
