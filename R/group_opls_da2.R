
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
#' @export group_opls_da2
#'
#' @examples
#' metabonomics_analysis(a,group)
group_opls_da2 <- function(a,group) {
  #######################################all_pca


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
    group_lable <- group_lable
    ATR_oplsda <- opls(group_opls,
                       group_lable,predI = 1,
                       orthoI = 1,
                       crossvalI=6,
                       printL=T,
                       plotL=F,
                       log10L = TRUE,
                       permI=199)


    prem <- data.frame(ATR_oplsda@suppLs[["permMN"]],'location'=factor(1:200))
    prem <- prem[,c(2,3,7,8)]
    colnames(prem)
    colnames(prem)=c("R2Y","Q2Y","sim",'location')
    colnames(prem)
    prem2 <- reshape2::melt(prem[,c(1:2,4)])
    #prem <- melt(prem,
    #id.vars = "sim", #你不想改变的数据列
    #measure.vars = c("R2Y","Q2Y"), #你要melt的数据)
    prem <- data.frame(prem2,'sim'=rep(prem[,3],2))
    ATR_oplsda@summaryDF
    p_y <- round(ATR_oplsda@summaryDF[7],3)-round(ATR_oplsda@summaryDF[7],3)%%0.005
    q_y <- round(ATR_oplsda@summaryDF[8],3)-round(ATR_oplsda@summaryDF[8],3)%%0.005
    segment1 <- as.numeric(ATR_oplsda@summaryDF[1])
    segment2 <- as.numeric(ATR_oplsda@summaryDF[2])
    segment3 <- as.numeric(ATR_oplsda@summaryDF[3])

    Permutation_plot <- ggplot(data = prem,aes(x=value,fill=variable))+
      scale_fill_manual(values = c("#D95F02","#7570B3"),
                        name="",
                        breaks = c("R2Y","Q2Y"),
                        labels = c(substitute(paste("Perm ","R"^"2","Y")),
                                   substitute(paste("Perm ","Q"^"2"))
                        ))+
      geom_histogram(bins = 30,col='white')+
      labs(x='Permutations',y='Frequency',fill='')+
      theme_bw()+#主题主题主题主题主题主题主题主题主题主题主题主题主题主题
      # theme(legend.position = c(0, 1.03),
      #       legend.justification = c(0, 1))+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())+
      theme(legend.background = element_rect(colour=NA,
                                             size=0.5,fill=NA))+
      geom_segment(aes(x=segment1,
                       y=50,
                       xend=segment1,
                       yend=0),
                   arrow = arrow(length=unit(0.2, "cm"),type = 'closed'))+
      geom_segment(aes(x=segment2,
                       y=100,
                       xend=segment2,
                       yend=0),
                   arrow = arrow(length=unit(0.2, "cm"),type = 'closed'))+
      geom_segment(aes(x=segment3,
                       y=75,
                       xend=segment3,
                       yend=0),
                   arrow = arrow(length=unit(0.2, "cm"),type = 'closed'))+
      annotate("text",
               label = substitute(paste("R"^"2","Y = ",R2Y),list(R2Y=as.numeric(ATR_oplsda@summaryDF[2]))),
               x = as.numeric(ATR_oplsda@summaryDF[2]),#R2Y#R2Y#R2Y#R2Y#R2Y#R2Y#R2Y
               y = 109)+
      xlim(c(min(prem$value),1.2))+
      annotate("text",
               label = paste("p = ",
                             p_y,
                             '(',p_y*200,'/200)'),
               x = as.numeric(ATR_oplsda@summaryDF[2]),
               y = 104)+
      annotate("text",
               label = substitute(paste("Q"^"2"," = ",Q2),list(Q2=as.numeric(ATR_oplsda@summaryDF[3]))),
               x = as.numeric(ATR_oplsda@summaryDF[3]),#Q2#Q2#Q2#Q2#Q2#Q2#Q2#Q2
               y = 84)+
      annotate("text",
               label = paste("p = ",
                             q_y,
                             '(',q_y*200,'/200)'),
               x = as.numeric(ATR_oplsda@summaryDF[3]),
               y = 79)+
      annotate("text",
               label = substitute(paste("R"^"2","X = ",R2),list(R2=as.numeric(ATR_oplsda@summaryDF[1]))),
               x = as.numeric(ATR_oplsda@summaryDF[1]),#R2X#R2X#R2X#R2X#R2X#R2X#R2X
               y = 54)
    #opls-da Permutation
    dir.create(paste0('opls-da/',as.character(group[i,2]),'_vs_',
                      as.character(group[i,3])), recursive = TRUE)
    ggsave(file = paste0('opls-da/',as.character(group[i,2]),'_vs_',
                         as.character(group[i,3]),'/',as.character(group[i,2]),'_vs_',
                         as.character(group[i,3]),
                         '_Permutation_plot.pdf'),
           width = 10,height = 8)
    eoffice::topptx(Permutation_plot,
                    paste0('opls-da/',as.character(group[i,2]),'_vs_',
                           as.character(group[i,3]),'/',as.character(group[i,2]),'_vs_',
                           as.character(group[i,3]),
                           '_Permutation_plot.pptx'),
                    width = 10,height = 8)



    ######OPLS-DA得分图
    OPLS_defentu <- data.frame(ATR_oplsda@scoreMN,
                               ATR_oplsda@orthoScoreMN)
    OPLS_defentu$Group <- group_lable
    OPLS_defentu$label=rownames(OPLS_defentu)
    #library(ggplot2)
    #设置适合科研的背景色
    theme_set(theme_bw())
    Scores_OPLS_DA_Plot <- ggplot(OPLS_defentu,aes(p1,o1,label = label))+
      geom_point(aes(colour=Group),size=3,shape=16)+
      theme(legend.background = element_rect(colour="black", size=0.5))+
      scale_colour_manual(values =c('#2ea582','#d95f02'))+
      labs(x=paste0('T score[1] (',ATR_oplsda@modelDF$R2X[1]*100,'%)'),
           y=paste0('Orthogonal T score[1] (',ATR_oplsda@modelDF$R2X[2]*100,'%)'),
           title = 'Scores OPLS-DA Plot')+
      theme(plot.title = element_text(hjust = 0.5),
            text = element_text(size = 12))+
      stat_ellipse(aes(fill=Group),
                   type = "norm",geom = "polygon",
                   alpha=0.2,color=NA,
                   show.legend = FALSE)+
      scale_fill_manual(values =c('#2ea582','#d95f02'))+
      stat_ellipse(aes(fill=Group,color=Group),
                   type = "norm",
                   lwd=1.2,
                   linetype = 1,show.legend = FALSE)+
      geom_label_repel(min.segment.length = 0,
                       arrow = arrow(length=unit(0.01, "npc"),
                                     ends = 'last',
                                     type = "closed"),
                       box.padding = 1,
                       point.padding=unit(1, "lines"),size = 2#标签字体大小
      )
    ggsave(file = paste0('opls-da/',as.character(group[i,2]),'_vs_',
                         as.character(group[i,3]),'/',as.character(group[i,2]),'_vs_',
                         as.character(group[i,3]),
                         '_Scores_OPLS_DA_Plot.pdf'),
           width = 6.8,height = 5.5)
    eoffice::topptx(Scores_OPLS_DA_Plot,
                    paste0('opls-da/',as.character(group[i,2]),'_vs_',
                           as.character(group[i,3]),'/',as.character(group[i,2]),'_vs_',
                           as.character(group[i,3]),
                           '_Scores_OPLS_DA_Plot.pptx'),
                    width = 8,height = 6.5)
    print(paste0("已完成",i,"个Scores_OPLS_DA_Plot",
                 "和OPLS_DA_Permutation_plot"))

  }

}
