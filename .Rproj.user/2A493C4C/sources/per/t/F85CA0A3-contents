#' Title
#'
#' @param a    a matrix
#' @param group    group
#' @param all_compounds  all compounds
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
#' @import reshape2
#' @import eoffice
#' @import plyr
#' @export group_pca
#'
group_pca <- function(all_compounds,a,group) {
  #######################################all_pca
  for (i in round(1:dim(group)[1])) {
    group_a <- rbind(a[which(a$Class_ID==as.character(group[i,2])),],
                     a[which(a$Class_ID==as.character(group[i,3])),])
    col_num <- dim(group_a)[1]
    row_num <- dim(group_a)[2]

    group_pca_data <- rbind(a[which(a$Class_ID==as.character(group[i,2])),],
                            a[which(a$Class_ID==as.character(group[i,3])),])
    group_pca_lable <- group_pca_data$Class_ID

    group_pca_data <- data.frame(group_pca_data)
    rownames(group_pca_data)=group_pca_data[,1]
    group_pca_data <- group_pca_data[,-c(1:2)]
    group_pca_data <- data.frame(group_pca_data)
    group_pca_variables <- PCA(group_pca_data,  graph = F)
    pca_data2 <- prcomp(group_pca_data, scale = TRUE)
    #查看合适主成分个数
    #screeplot(pca_data2, type = "lines")
    summary(pca_data2)
    #查看行名，确认是否为样品的名称
    rownames(pca_data2$x)
    #提取PC1的百分比
    x_per <- round(summary(pca_data2)$importance[2, 1]*100, 1)
    #提取PC2的百分比
    y_per <- round(summary(pca_data2)$importance[2, 2]*100, 1)
    #按照样品名称添加组并且合并
    df_sample <- data.frame(samplenames=rownames(pca_data2$x), pca_data2$x)
    df_sample$group <- group_pca_lable
    pp2 <- ggplot(df_sample, aes(x = PC1, y = PC2,label=samplenames)) +
      geom_point(aes(colour=group),size=3,shape=16) +
      xlab(paste("PC1","(", x_per,"%)",sep=" ")) +
      ylab(paste("PC2","(", y_per,"%)",sep=" ")) +
      theme(legend.background = element_rect(colour="black", size=0.5))+
      scale_colour_manual(values = c("#E41A1C" ,"#377EB8"))+
      theme_bw()+
      geom_label_repel(aes(col=NULL),
                       min.segment.length = 0,
                       arrow = arrow(length=unit(0.01, "npc"),
                                     ends = 'last',
                                     type = "closed"),
                       box.padding = 1,
                       colour=rep(c('#F781BF',
                                    '#A1CAF1'),each=3),
                       segment.colour=rep(c('#F781BF',
                                            '#A1CAF1'),each=3),
                       point.padding=unit(1, "lines"))+
      theme(plot.title = element_text(hjust = 0.5),
            text = element_text(size = 15))
    #pca
    dir.create(paste0('pca/',as.character(group[i,2]),'_vs_',
                      as.character(group[i,3]),'/'), recursive = TRUE)
    ggsave(file = paste0('pca/',as.character(group[i,2]),'_vs_',
                         as.character(group[i,3]),'/',as.character(group[i,2]),'_vs_',
                         as.character(group[i,3]),
                         'pca.pdf'),
           width = 7,height = 5.5)
    topptx(pp2,
           paste0('pca/',as.character(group[i,2]),'_vs_',
                  as.character(group[i,3]),'/',as.character(group[i,2]),'_vs_',
                  as.character(group[i,3]),
                  'pca.pptx'),
           width = 7,height = 5.5)
    #################################################
    #################################################
    group_vc_rate <- data.frame(group_pca_variables[["eig"]])
    group_vc_rate$pc <- paste0('PC',1:dim(group_vc_rate)[1])
    group_vc_rate$pc <- factor(group_vc_rate$pc,levels=group_vc_rate$pc)
    group_vc_rate$group <- 'group'
    group_vc_rate <- group_vc_rate[,-1]
    group_vc_rate <- head(group_vc_rate,5)
    group_vc_rate <- reshape2::melt(group_vc_rate)
    # 这个是我们希望展示出来的标签名
    group_dose_labs <- c('Proportion of Variance',
                         'Cumulative Proportion')
    # 这个是我们希望隐藏的标签名
    names(group_dose_labs) <- unique(group_vc_rate$variable)
    ggplot(group_vc_rate,aes(pc,value,group=variable))+
      geom_point(size = 2.5)+
      geom_line(lty = 2,lwd = 0.8)+
      facet_wrap(~variable,
                 labeller = labeller(variable = group_dose_labs))+
      theme_bw()+
      labs(title = 'Variance explained by each principal component',
           x='Principal Component',
           y='variance')+
      theme(plot.title = element_text(hjust = 0.5))
    pp3 <- ggplot(group_vc_rate,aes(pc,value,group=variable))+
      geom_point(size = 2.5)+
      geom_line(lty = 2,lwd = 0.8)+
      facet_wrap(~variable,
                 labeller = labeller(variable = group_dose_labs))+
      theme_bw()+
      labs(title = 'Variance explained by each principal component',
           x='Principal Component',
           y='variance')+
      theme(plot.title = element_text(hjust = 0.5))
    ggsave(file = paste0('pca/',as.character(group[i,2]),'_vs_',
                         as.character(group[i,3]),'/',as.character(group[i,2]),'_vs_',
                         as.character(group[i,3]),
                         '_Contribution_rate_of_variance.pdf'),
           width = 8,height = 6)
    topptx(pp3,
           paste0('pca/',as.character(group[i,2]),'_vs_',
                  as.character(group[i,3]),'/',as.character(group[i,2]),'_vs_',
                  as.character(group[i,3]),
                  '_Contribution_rate_of_variance.pptx'),
           width = 8,height = 6)





    #输出方差贡献率和pca计算结果
    group_vc_rate <- data.frame(group_pca_variables[["eig"]])
    colnames(group_vc_rate)[2:3] <- c('Proportion of Variance',
                                      'Cumulative Proportion')
    rownames(group_vc_rate) <- paste0('PC',1:dim(group_vc_rate)[1])
    write.csv(group_vc_rate,paste0('pca/',as.character(group[i,2]),'_vs_',
                                   as.character(group[i,3]),'/',as.character(group[i,2]),'_vs_',
                                   as.character(group[i,3]),
                                   'Contribution_rate_of_variance.csv'))
    group_pca_result <- prcomp(group_pca_data, scale = TRUE)
    group_pca_result <- data.frame(group_pca_result$rotation)
    # ###########新增
    group_pca_result$index <- rownames(group_pca_result)
    group_pca_result <- join(group_pca_result,all_compounds,by='index')
    write.csv(group_pca_result,paste0('pca/',as.character(group[i,2]),'_vs_',
                                      as.character(group[i,3]),'/',as.character(group[i,2]),'_vs_',
                                      as.character(group[i,3]),
                                      'pca_result_rotation.csv'))
    ################################################
    ################################################
  }
}
