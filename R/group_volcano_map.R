#' @title Metabonomics analysis3
#' @description Metabonomics analysis3
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
#' @import reshape2
#' @import eoffice
#' @import plyr
#' @export group_volcano_map
#'
#' @examples
#' metabonomics_analysis(a,group)
group_volcano_map <- function(a,group) {
  for (i in round(1:dim(group)[1])) {
    group_a <- rbind(a[which(a$Class_ID==as.character(group[i,2])),],
                     a[which(a$Class_ID==as.character(group[i,3])),])
    col_num <- dim(group_a)[1]
    row_num <- dim(group_a)[2]
    #分组信息
    group_lable <- group_a$Class_ID
    #样品编号
    sample_id <- group_a$Secondary_ID
    rownames(group_a)=group_a[,1]
    group_a <- group_a[,-c(1:2)]
    group_a <- data.frame(group_a)
    group_opls <- group_a
    group_a <- data.frame(t(group_a))
    colnames(group_a) <- sample_id
    sample1 <- group_a[,1:(dim(group_a)[2]/2)]
    sample2 <- group_a[,(dim(group_a)[2]/2+1):dim(group_a)[2]]
    group_a$p_value <- sapply(1:nrow(sample1), function(x) t.test(sample1[x,], sample2[x,],
                                                                  alternative = 'two.sided',
                                                                  var.equal = FALSE)$p.value)

    group_a$fold_change <- apply(sample2,1,mean)/apply(sample1,1,mean)
    group_a$log2FC <- log(group_a$fold_change,2)
    group_lable <- group_lable
    ATR_oplsda <- opls(group_opls,
                       group_lable,predI = 1,
                       orthoI = 1,
                       crossvalI=6,
                       printL=T,
                       plotL=F,
                       log10L = TRUE,
                       permI=1)
    group_a$VIP <- ATR_oplsda@vipVn
    up <- which(group_a$VIP > 1 &  group_a$fold_change > 2)
    group_a$regulate <- 'insignificant'
    group_a$regulate[up] <- 'up'
    down <- which(group_a$VIP > 1 &  group_a$fold_change < 0.5)
    group_a$regulate[down] <- 'down'
    # # ###########新增
    # group_a$index <- rownames(group_a)
    # group_a <- plyr::join(group_a,all_compounds,by='index')
    # dir.create(paste0('group_table_of_vip_p_fc/',as.character(group[i,2]),'_vs_',
    #                   as.character(group[i,3])), recursive = TRUE)
    # write.csv(group_a,paste0('group_table_of_vip_p_fc/',as.character(group[i,2]),'_vs_',
    #                          as.character(group[i,3]),'/',as.character(group[i,2]),'_vs_',as.character(group[i,3]),'_all.csv'))
    # group_DEM_up <- group_a[which(group_a$VIP > 1 &  group_a$fold_change > 2),]
    # # # ###########新增
    # # group_DEM_up$index <- rownames(group_DEM_up)
    # # group_DEM_up <- plyr::join(group_DEM_up,all_compounds,by='index')
    # write.csv(group_DEM_up,paste0('group_table_of_vip_p_fc/',as.character(group[i,2]),'_vs_',
    #                               as.character(group[i,3]),'/',as.character(group[i,2]),'_vs_',as.character(group[i,3]),'_up.csv'))
    # group_DEM_down <- group_a[which(group_a$VIP > 1 &  group_a$fold_change < 0.5),]
    # # # ###########新增
    # # group_DEM_down$index <- rownames(group_DEM_down)
    # # group_DEM_down <- plyr::join(group_DEM_down,all_compounds,by='index')
    # write.csv(group_DEM_down,paste0('group_table_of_vip_p_fc/',as.character(group[i,2]),'_vs_',
    #                                 as.character(group[i,3]),'/',as.character(group[i,2]),'_vs_',as.character(group[i,3]),'_down.csv'))
    #
    # group_DEM_up_down <- rbind(group_DEM_up,group_DEM_down)
    # # # ###########新增
    # # group_DEM_up_down$index <- rownames(group_DEM_up_down)
    # # group_DEM_up_down <- plyr::join(group_DEM_up_down,all_compounds,by='index')
    # write.csv(group_DEM_up_down,paste0('group_table_of_vip_p_fc/',as.character(group[i,2]),'_vs_',
    #                                    as.character(group[i,3]),'/',as.character(group[i,2]),'_vs_',as.character(group[i,3]),'_up_and_down.csv'))
    #


    #####火山图
    max_x <- ceiling(max(group_a$log2FC))
    # max_x
    min_x <- floor(min(group_a$log2FC))
    # min_x
    breaks <- c(seq(min_x,max_x,by=4),1,-1)[order(c(seq(min_x,max_x, by = 4),1,-1))]
    group_volcano_map <- ggplot(group_a,aes(log2FC,VIP,col=regulate,size=p_value))+
      geom_point()+
      scale_color_manual(values = c('green',
                                    'grey',
                                    'orange'))+
      geom_vline(aes(xintercept=1),colour="gray55",
                 linetype="dashed",size=1)+
      geom_vline(aes(xintercept=-1),colour="gray55",
                 linetype="dashed",size=1)+
      geom_hline(aes(yintercept=1),colour="gray55",
                 linetype="dashed",size=1)+
      theme_bw()+
      theme(panel.grid=element_line(linetype = 'dashed'))+
      labs(title = '',x='',y='')+
      theme(plot.title = element_text(size = 20))+
      theme(text = element_text(size = 15))+
      theme(plot.title = element_text(hjust = 0.5))+
      labs(col='Statistics',size='p-value')+
      theme(legend.background = element_rect(colour="black", size=0.5))+
      theme_bw()+#删除此行图例加框
      labs(x=expression(paste("Log"["2"],"FC")),
           y='Variable Importance in Projection(VIP)')+
      scale_x_continuous(breaks = breaks)+
      theme(text = element_text(size = 13))
    dir.create(paste0('volcano_map/',as.character(group[i,2]),'_vs_',
                      as.character(group[i,3]),'/'), recursive = TRUE)
    ggsave(file = paste0('volcano_map/',as.character(group[i,2]),'_vs_',
                         as.character(group[i,3]),'/',as.character(group[i,2]),'_vs_',
                         as.character(group[i,3]),
                         '_Volcano_map.pdf'),
           width = 8,height = 7)
    eoffice::topptx(group_volcano_map,
                    paste0('volcano_map/',as.character(group[i,2]),'_vs_',
                           as.character(group[i,3]),'/',as.character(group[i,2]),'_vs_',
                           as.character(group[i,3]),
                           '_Volcano_map.pptx'),
                    width = 8,height = 7,top = 0,left = 0)

    print(paste0("第",i+1,"组数据分析中......"))
  }

}
