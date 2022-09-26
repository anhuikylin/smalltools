
#' @title Metabonomics analysis3
#' @description Metabonomics analysis3
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
#' @export metabonomics_analysis3_3
#'
#' @examples
#' metabonomics_analysis(a,group)
metabonomics_analysis3_3 <- function(all_compounds,a,group) {
  #######################################all_pca
  group_all_pca_data <- a
  group_all_pca_lable <- group_all_pca_data$Class_ID
  group_all_pca_data <- data.frame(group_all_pca_data)
  rownames(group_all_pca_data)=group_all_pca_data[,1]
  group_all_pca_data <- group_all_pca_data[,-c(1:2)]
  group_all_pca_data <- data.frame(group_all_pca_data)
  #iris.pca[["eig"]][1,2]
  #iris.pca[["eig"]][1,2]
  group_all_pca_variables <- PCA(group_all_pca_data,  graph = F)
  all_pca_no_lable <- fviz_pca_ind(group_all_pca_variables,
                                   geom.ind = c("point"),
                                   habillage = factor(group_all_pca_lable), # color by groups
                                   palette = rainbow(length(unique(a$`Class_ID`))),
                                   addEllipses = TRUE)+
    theme_bw()+
    labs(title = '',
         x=paste0('PC1','(',round(group_all_pca_variables[["eig"]][1,2],2),'%)'),
         y=paste0('PC2','(',round(group_all_pca_variables[["eig"]][2,2],2),'%)'))+
    theme(plot.title = element_text(hjust = 0.5),
          text = element_text(size = 20))
  dir.create('1.Data_assess/pca', recursive = TRUE)
  ggsave(file = paste0('1.Data_assess/pca/all_pca.pdf'),
         width = 10,height = 8)
  topptx(all_pca_no_lable,paste0('1.Data_assess/pca/all_pca.pptx'),
         width = 10,
         height = 8)
  #碎石图
  all_vc_rate <- data.frame(group_all_pca_variables[["eig"]])
  all_vc_rate$pc <- paste0('PC',1:dim(all_vc_rate)[1])
  all_vc_rate$pc <- factor(all_vc_rate$pc,levels=all_vc_rate$pc)
  all_vc_rate$group <- 'group'
  all_vc_rate <- all_vc_rate[,-1]
  all_vc_rate <- head(all_vc_rate,5)
  all_vc_rate <- reshape2::melt(all_vc_rate)
  # 这个是我们希望展示出来的标签名
  all_dose_labs <- c('Proportion of Variance',
                     'Cumulative Proportion')
  # 这个是我们希望隐藏的标签名
  names(all_dose_labs) <- unique(all_vc_rate$variable)
  ggplot(all_vc_rate,aes(pc,value,group=variable))+
    geom_point(size = 2.5)+
    geom_line(lty = 2,lwd = 0.8)+
    facet_wrap(~variable,
               labeller = labeller(variable = all_dose_labs))+
    theme_bw()+
    labs(title = 'Variance explained by each principal component',
         x='Principal Component',
         y='variance')+
    theme(plot.title = element_text(hjust = 0.5))
  pp1 <- ggplot(all_vc_rate,aes(pc,value,group=variable))+
    geom_point(size = 2.5)+
    geom_line(lty = 2,lwd = 0.8)+
    facet_wrap(~variable,
               labeller = labeller(variable = all_dose_labs))+
    theme_bw()+
    labs(title = 'Variance explained by each principal component',
         x='Principal Component',
         y='variance')+
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(file = paste0('1.Data_assess/pca/Contribution_rate_of_variance.pdf'),
         width = 8,height = 6)
  topptx(pp1,paste0('1.Data_assess/pca/Contribution_rate_of_variance.pptx'),
         width = 10,
         height = 8)



  #输出方差贡献率和pca计算结果
  all_vc_rate <- data.frame(group_all_pca_variables[["eig"]])
  colnames(all_vc_rate)[2:3] <- c('Proportion of Variance',
                                  'Cumulative Proportion')
  rownames(all_vc_rate) <- paste0('PC',1:dim(all_vc_rate)[1])
  write.csv(all_vc_rate,'1.Data_assess/pca/Contribution_rate_of_variance.csv')
  pca_result <- prcomp(group_all_pca_data, scale = TRUE)
  pca_result <- data.frame(pca_result$rotation)
  ###########新增
  pca_result$index <- rownames(pca_result)
  pca_result <- join(pca_result,all_compounds,by='index')
  write.csv(pca_result,'1.Data_assess/pca/pca_result_rotation.csv')

  #######################################
  #######################################all_cor,QC_cor
  qc_cor <- a[which(a$Class_ID=="QC"),]
  qc_cor <- data.frame(qc_cor)
  rownames(qc_cor) <- qc_cor[,1]
  col_qc_cor <- qc_cor$Secondary_ID
  qc_cor <- qc_cor[,-(1:2)]
  qc_cor <- data.frame(t(qc_cor))
  colnames(qc_cor) <- col_qc_cor


  dir.create('1.Data_assess/correlation_analysis', recursive = TRUE)
  pdf(file = '1.Data_assess/correlation_analysis/QC_cor.pdf',
      width = 8,height = 8)
  corrplot::corrplot(cor(qc_cor), type = 'lower',
           method = 'color',
           #order = 'AOE',
           diag = FALSE,
           tl.pos = 'ld',#字体位置
           cl.pos = 'r',#图例位置
           tl.col = 'black',#轴字体颜色
           cl.length = 5, #数字大小
           addCoef.col = 'white',#数字颜色
           col.lim = c(0,1),#范围必须覆盖矩阵
           col = colorRampPalette(c("#d53e4f",
                                    '#fdb366',
                                    "#6cc4a4",
                                    '#d6de96',
                                    "#d53e4f"))(100),
           outline = "white"
  )
  dev.off()

  qc_cor_result <- cor(qc_cor)
  qc_cor_result <- data.frame(qc_cor_result)
  colnames(qc_cor_result) <- rownames(qc_cor_result)
  write.csv(qc_cor_result,'1.Data_assess/correlation_analysis/QC_cor.csv')



  #所有样本相关性分析
  all_sample_cor <- a[which(a$Class_ID != "QC"),]
  #############热图分组数据
  all_sample_cor_group <- data.frame('group' = all_sample_cor$Class_ID)
  rownames(all_sample_cor_group) <- all_sample_cor$Secondary_ID
  ########################
  all_sample_cor <- data.frame(all_sample_cor)
  rownames(all_sample_cor) <- all_sample_cor[,1]
  col_all_sample_cor <- all_sample_cor$Secondary_ID
  all_sample_cor <- all_sample_cor[,-(1:2)]
  all_sample_cor <- data.frame(t(all_sample_cor))
  colnames(all_sample_cor) <- col_all_sample_cor
  library(pheatmap)
  library(RColorBrewer)
  display.brewer.all()
  brewer.pal(9,'Set1')
  brewer.pal(8,'Dark2')
  #暂时未修改颜色
  ann_colors_list <- c("#1B9E77" , "#D95F02" ,"#7570B3" ,"#E7298A" ,
                       "#66A61E" , "#E6AB02" ,"#A6761D" ,"#666666",
                       "#E41A1C" , "#377EB8" ,"#984EA3" ,"#FF7F00",
                       "#FFFF33", "#A65628" ,"#F781BF", "#999999" )
  #分组

  annotation_col = all_sample_cor_group

  ###所有样本相关性热图
  all_sample_cor_group_pheatmap <- pheatmap::pheatmap(cor(all_sample_cor),
                                            annotation_col=annotation_col,
                                            cluster_rows = FALSE,
                                            cluster_cols = FALSE,
                                            color = colorRampPalette(c('green','white','red'))(100),
                                            display_numbers = TRUE,
                                            cellwidth = 20,
                                            cellheight = 20)

  pdf(file = '1.Data_assess/correlation_analysis/all_sample_cor.pdf',
      width = length(all_sample_cor)/3,height = length(all_sample_cor)/3)
  print(all_sample_cor_group_pheatmap)
  dev.off()
  topptx(all_sample_cor_group_pheatmap,
         '1.Data_assess/correlation_analysis/all_sample_cor.pptx',
         width = length(all_sample_cor)/3,height = length(all_sample_cor)/3)

  #spearman
  cor_all_sample_cor <- cor(all_sample_cor,method = 'spearman')
  cor_all_sample_cor <- data.frame(cor_all_sample_cor)
  write.csv(cor_all_sample_cor,'1.Data_assess/correlation_analysis/all_sample_cor_spearman.csv')
  #pearson
  cor_all_sample_cor <- cor(all_sample_cor,method = 'pearson')
  cor_all_sample_cor <- data.frame(cor_all_sample_cor)
  write.csv(cor_all_sample_cor,'1.Data_assess/correlation_analysis/all_sample_cor_pearson.csv')

  all_sample_pheatmap <- pheatmap::pheatmap(all_sample_cor,
                                  scale = 'row',
                                  annotation_col=annotation_col,
                                  cluster_rows = TRUE,
                                  cluster_cols = FALSE,
                                  color = colorRampPalette(c('#1f9a51',
                                                             '#b0dd70',
                                                             '#fefcbb',
                                                             '#fdb96a',
                                                             '#d73127'))(100),
                                  show_rownames = FALSE,cellwidth = 20)
  dir.create('1.Data_assess/heatmap', recursive = TRUE)
  pdf(file = '1.Data_assess/heatmap/all_heatmap.pdf',
      width = length(all_sample_cor)/3,height = length(all_sample_cor)/4)
  print(all_sample_pheatmap)
  dev.off()
  topptx(all_sample_pheatmap,
         '1.Data_assess/heatmap/all_heatmap.pptx',
         width = length(all_sample_cor)/3,height = length(all_sample_cor)/4)

  for (i in round(1:dim(group)[1])) {
    group_a <- rbind(a[which(a$Class_ID==as.character(group[i,2])),],
                     a[which(a$Class_ID==as.character(group[i,3])),])
    col_num <- dim(group_a)[1]
    row_num <- dim(group_a)[2]

    group_pca_data <- rbind(a[which(a$Class_ID==as.character(group[i,2])),],
                            a[which(a$Class_ID==as.character(group[i,3])),])
    group_pca_lable <- group_pca_data$Class_ID

    group_pca_lable <- data.frame(group_pca_lable)
    group_pca_data <- data.frame(group_pca_data)
    rownames(group_pca_data)=group_pca_data[,1]
    group_pca_data <- group_pca_data[,-c(1:2)]
    group_pca_data <- data.frame(group_pca_data)
    group_pca_variables <- PCA(group_pca_data,  graph = F)
    factoextra::fviz_pca_ind(group_pca_variables,
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
    pp2 <- factoextra::fviz_pca_ind(group_pca_variables,
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
    #pca
    dir.create(paste0('2.Differential_analyse/',as.character(group[i,2]),'_vs_',
                      as.character(group[i,3]),'/pca/'), recursive = TRUE)
    ggsave(file = paste0('2.Differential_analyse/',as.character(group[i,2]),'_vs_',
                         as.character(group[i,3]),'/pca/',as.character(group[i,2]),'_vs_',
                      as.character(group[i,3]),
                      'pca.pdf'),
        width = 12,height = 12)
    topptx(pp2,
           paste0('2.Differential_analyse/',as.character(group[i,2]),'_vs_',
                  as.character(group[i,3]),'/pca/',as.character(group[i,2]),'_vs_',
                  as.character(group[i,3]),
                  'pca.pptx'),
           width = 12,height = 12)
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
    ggsave(file = paste0('2.Differential_analyse/',as.character(group[i,2]),'_vs_',
                         as.character(group[i,3]),'/pca/',as.character(group[i,2]),'_vs_',
                         as.character(group[i,3]),
                         '_Contribution_rate_of_variance.pdf'),
           width = 8,height = 6)
    topptx(pp3,
           paste0('2.Differential_analyse/',as.character(group[i,2]),'_vs_',
                  as.character(group[i,3]),'/pca/',as.character(group[i,2]),'_vs_',
                  as.character(group[i,3]),
                  '_Contribution_rate_of_variance.pptx'),
           width = 8,height = 6)





    #输出方差贡献率和pca计算结果
    group_vc_rate <- data.frame(group_pca_variables[["eig"]])
    colnames(group_vc_rate)[2:3] <- c('Proportion of Variance',
                                    'Cumulative Proportion')
    rownames(group_vc_rate) <- paste0('PC',1:dim(group_vc_rate)[1])
    write.csv(group_vc_rate,paste0('2.Differential_analyse/',as.character(group[i,2]),'_vs_',
                                   as.character(group[i,3]),'/pca/',as.character(group[i,2]),'_vs_',
                                   as.character(group[i,3]),
                                   'Contribution_rate_of_variance.csv'))
    group_pca_result <- prcomp(group_pca_data, scale = TRUE)
    group_pca_result <- data.frame(group_pca_result$rotation)
    # ###########新增
    group_pca_result$index <- rownames(group_pca_result)
    group_pca_result <- join(group_pca_result,all_compounds,by='index')
    write.csv(group_pca_result,paste0('2.Differential_analyse/',as.character(group[i,2]),'_vs_',
                                      as.character(group[i,3]),'/pca/',as.character(group[i,2]),'_vs_',
                                      as.character(group[i,3]),
                                      'pca_result_rotation.csv'))
    ################################################
    ################################################




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
                       permI=199)
    group_a$VIP <- ATR_oplsda@vipVn
    up <- which(group_a$VIP > 1 &  group_a$fold_change > 2)
    group_a$regulate <- 'insignificant'
    group_a$regulate[up] <- 'up'
    down <- which(group_a$VIP > 1 &  group_a$fold_change < 0.5)
    group_a$regulate[down] <- 'down'
    # ###########新增
    group_a$index <- rownames(group_a)
    group_a <- join(group_a,all_compounds,by='index')
    write.csv(group_a,paste0('2.Differential_analyse/',as.character(group[i,2]),'_vs_',
                             as.character(group[i,3]),'/',as.character(group[i,2]),'_vs_',as.character(group[i,3]),'_all.csv'))
    group_DEM_up <- group_a[which(group_a$VIP > 1 &  group_a$fold_change > 2),]
    # ###########新增
    group_DEM_up$index <- rownames(group_DEM_up)
    group_DEM_up <- join(group_DEM_up,all_compounds,by='index')
    write.csv(group_DEM_up,paste0('2.Differential_analyse/',as.character(group[i,2]),'_vs_',
                                  as.character(group[i,3]),'/',as.character(group[i,2]),'_vs_',as.character(group[i,3]),'_up.csv'))
    group_DEM_down <- group_a[which(group_a$VIP > 1 &  group_a$fold_change < 0.5),]
    # ###########新增
    group_DEM_down$index <- rownames(group_DEM_down)
    group_DEM_down <- join(group_DEM_down,all_compounds,by='index')
    write.csv(group_DEM_down,paste0('2.Differential_analyse/',as.character(group[i,2]),'_vs_',
                                    as.character(group[i,3]),'/',as.character(group[i,2]),'_vs_',as.character(group[i,3]),'_down.csv'))

    group_DEM_up_down <- rbind(group_DEM_up,group_DEM_down)
    # ###########新增
    group_DEM_up_down$index <- rownames(group_DEM_up_down)
    group_DEM_up_down <- join(group_DEM_up_down,all_compounds,by='index')
    write.csv(group_DEM_up_down,paste0('2.Differential_analyse/',as.character(group[i,2]),'_vs_',
                                       as.character(group[i,3]),'/',as.character(group[i,2]),'_vs_',as.character(group[i,3]),'_up_and_down.csv'))


    prem <- data.frame(ATR_oplsda@suppLs[["permMN"]],'location'=factor(1:200))
    prem <- prem[,c(2,3,7,8)]
    colnames(prem)
    colnames(prem)=c("R2Y","Q2Y","sim",'location')
    colnames(prem)
    library(reshape2)
    prem2 <- melt(prem[,c(1:2,4)])
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
      geom_histogram(bins = 30,col='black')+
      labs(x='Permutations',y='Frequency',fill='')+
      theme_bw()+#主题主题主题主题主题主题主题主题主题主题主题主题主题主题
      theme(legend.position = c(0, 1.03),
            legend.justification = c(0, 1))+
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
    dir.create(paste0('2.Differential_analyse/',as.character(group[i,2]),'_vs_',
                      as.character(group[i,3]),'/opls/'), recursive = TRUE)
    ggsave(file = paste0('2.Differential_analyse/',as.character(group[i,2]),'_vs_',
                         as.character(group[i,3]),'/opls/',as.character(group[i,2]),'_vs_',
                         as.character(group[i,3]),
                         '_Permutation_plot.pdf'),
           width = 8,height = 8)
    topptx(Permutation_plot,
           paste0('2.Differential_analyse/',as.character(group[i,2]),'_vs_',
                  as.character(group[i,3]),'/opls/',as.character(group[i,2]),'_vs_',
                  as.character(group[i,3]),
                  '_Permutation_plot.pptx'),
           width = 8,height = 6)


    #VIP值图
    group_vip_plot <- data.frame('id' = rownames(group_a),group_a)
    colnames(group_vip_plot)
    group_vip_plot <- group_vip_plot[which(group_vip_plot$regulate != 'insignificant'),]
    group_vip_plot <- group_vip_plot[order(group_vip_plot$VIP,decreasing = TRUE),]
    group_vip_plot <- head(group_vip_plot,20)
    group_vip_plot <- group_vip_plot[order(group_vip_plot$VIP),]
    group_vip_plot$id <- factor(group_vip_plot$id,levels = group_vip_plot$id)
    ggplot(group_vip_plot,aes(VIP,id,col = regulate))+
      geom_point(size = 6)+
      theme_bw()+
      labs(x='VIP scores',
           y='')+
      scale_color_manual(values = c('red','blue'))+
      theme(plot.title = element_text(hjust = 0.5),
            text = element_text(size = 18))
    pp4 <- ggplot(group_vip_plot,aes(VIP,id,col = regulate))+
      geom_point(size = 6)+
      theme_bw()+
      labs(x='VIP scores',
           y='')+
      scale_color_manual(values = c('red','blue'))+
      theme(plot.title = element_text(hjust = 0.5),
            text = element_text(size = 18))
    dir.create(paste0('2.Differential_analyse/',as.character(group[i,2]),'_vs_',
                      as.character(group[i,3]),'/vipscore/'), recursive = TRUE)
    ggsave(file = paste0('2.Differential_analyse/',as.character(group[i,2]),'_vs_',
                         as.character(group[i,3]),'/vipscore/',as.character(group[i,2]),'_vs_',
                         as.character(group[i,3]),
                         '_VIP_value_plot.pdf'),
           width = 8,height = 12)
    topptx(pp4,
           paste0('2.Differential_analyse/',as.character(group[i,2]),'_vs_',
                  as.character(group[i,3]),'/vipscore/',as.character(group[i,2]),'_vs_',
                  as.character(group[i,3]),
                  '_VIP_value_plot.pptx'),
           width = 8,height = 12)


    #############小提琴图
    vip_box <- data.frame(group_a)
    colnames(vip_box)
    vip_box <- vip_box[which(vip_box$regulate != 'insignificant'),]
    vip_box <- vip_box[order(vip_box$VIP,decreasing = TRUE),]
    vip_box <- head(vip_box,25)
    vip_box <- vip_box[,1:col_num]
    vip_box <- t(vip_box)
    vip_box <- data.frame(vip_box)
    vip_box$group <- rep(c(as.character(group[i,2]),as.character(group[i,3])),each = col_num/2)
    vip_box <- melt(vip_box)
    ggplot(vip_box,aes(group,value,fill=group))+
      geom_boxplot(width=0.2)+
      geom_violin(alpha = 0.1)+
      facet_wrap(~variable,scales = 'free')+
      scale_y_continuous(labels = function(x) format(x, scientific = TRUE))+
      labs(y = 'Raw Intensity',title = 'Violin Plot of Raw Values')+
      scale_fill_manual(values = c("#1B9E77","#D95F02"))+
      theme_bw()+
      theme(plot.title = element_text(hjust = 0.5))
    pp5 <- ggplot(vip_box,aes(group,value,fill=group))+
      geom_boxplot(width=0.2)+
      geom_violin(alpha = 0.1)+
      facet_wrap(~variable,scales = 'free')+
      scale_y_continuous(labels = function(x) format(x, scientific = TRUE))+
      labs(y = 'Raw Intensity',title = 'Violin Plot of Raw Values')+
      scale_fill_manual(values = c("#1B9E77","#D95F02"))+
      theme_bw()+
      theme(plot.title = element_text(hjust = 0.5))
    dir.create(paste0('2.Differential_analyse/',as.character(group[i,2]),'_vs_',
                      as.character(group[i,3]),'/violin/'), recursive = TRUE)
    ggsave(file = paste0('2.Differential_analyse/',as.character(group[i,2]),'_vs_',
                         as.character(group[i,3]),'/violin/',as.character(group[i,2]),'_vs_',
                         as.character(group[i,3]),
                         '_group_violin_plot.pdf'),
           width = 12,height = 12)
    topptx(pp5,
           paste0('2.Differential_analyse/',as.character(group[i,2]),'_vs_',
                  as.character(group[i,3]),'/violin/',as.character(group[i,2]),'_vs_',
                  as.character(group[i,3]),
                  '_group_violin_plot.pptx'),
           width = 8,height = 6)
    #####火山图
    max_x <- ceiling(max(group_a$log2FC))
    max_x
    min_x <- floor(min(group_a$log2FC))
    min_x
    breaks <- c(seq(min_x,max_x,by=4),1,-1)[order(c(seq(min_x,max_x, by = 4),1,-1))]
    ggplot(group_a,aes(log2FC,VIP,col=regulate,size=p_value))+
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
      theme(text = element_text(size = 20))+
      theme(plot.title = element_text(hjust = 0.5))+
      labs(col='Statistics',size='p-value')+
      theme(legend.background = element_rect(colour="black", size=0.5))+
      theme_bw()+#删除此行图例加框
      labs(x=expression(paste("Log"["2"],"FC")),
           y='Variable Importance in Projection(VIP)')+
      scale_x_continuous(breaks = breaks)+
      theme(text = element_text(size = 20))
    pp6 <- ggplot(group_a,aes(log2FC,VIP,col=regulate,size=p_value))+
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
      theme(text = element_text(size = 20))+
      theme(plot.title = element_text(hjust = 0.5))+
      labs(col='Statistics',size='p-value')+
      theme(legend.background = element_rect(colour="black", size=0.5))+
      theme_bw()+#删除此行图例加框
      labs(x=expression(paste("Log"["2"],"FC")),
           y='Variable Importance in Projection(VIP)')+
      scale_x_continuous(breaks = breaks)+
      theme(text = element_text(size = 20))
    dir.create(paste0('2.Differential_analyse/',as.character(group[i,2]),'_vs_',
                      as.character(group[i,3]),'/Vol/'), recursive = TRUE)
    ggsave(file = paste0('2.Differential_analyse/',as.character(group[i,2]),'_vs_',
                         as.character(group[i,3]),'/Vol/',as.character(group[i,2]),'_vs_',
                         as.character(group[i,3]),
                         '_Volcano_map.pdf'),
           width = 12,height = 12)
    topptx(pp6,
           paste0('2.Differential_analyse/',as.character(group[i,2]),'_vs_',
                  as.character(group[i,3]),'/Vol/',as.character(group[i,2]),'_vs_',
                  as.character(group[i,3]),
                  '_Volcano_map.pptx'),
           width = 12,height = 12)

    ######OPLS-DA得分图
    OPLS_defentu <- data.frame(ATR_oplsda@scoreMN,
                               ATR_oplsda@orthoScoreMN)
    OPLS_defentu$Group <- group_lable
    OPLS_defentu$label=rownames(OPLS_defentu)
    library(ggplot2)
    #设置适合科研的背景色
    theme_set(theme_bw())
    ggplot(OPLS_defentu,aes(p1,o1,label = label))+
      geom_point(aes(colour=Group),size=3,shape=16)+
      theme(legend.background = element_rect(colour="black", size=0.5))+
      scale_colour_manual(values =c('#2ea582','#d95f02'))+
      labs(x=paste0('T score[1] (',ATR_oplsda@modelDF$R2X[1]*100,'%)'),
           y=paste0('Orthogonal T score[1] (',ATR_oplsda@modelDF$R2X[2]*100,'%)'),
           title = 'Scores OPLS-DA Plot')+
      theme(plot.title = element_text(hjust = 0.5),
            text = element_text(size = 25))+
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
                       colour=rep(c('#2ea582','#d95f02'),each=col_num/2),
                       segment.colour=rep(c('#2ea582','#d95f02'),each=col_num[1]/2),
                       point.padding=unit(1, "lines"),size = 8#标签字体大小
      )
    pp7 <- ggplot(OPLS_defentu,aes(p1,o1,label = label))+
      geom_point(aes(colour=Group),size=3,shape=16)+
      theme(legend.background = element_rect(colour="black", size=0.5))+
      scale_colour_manual(values =c('#2ea582','#d95f02'))+
      labs(x=paste0('T score[1] (',ATR_oplsda@modelDF$R2X[1]*100,'%)'),
           y=paste0('Orthogonal T score[1] (',ATR_oplsda@modelDF$R2X[2]*100,'%)'),
           title = 'Scores OPLS-DA Plot')+
      theme(plot.title = element_text(hjust = 0.5),
            text = element_text(size = 25))+
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
                       colour=rep(c('#2ea582','#d95f02'),each=col_num/2),
                       segment.colour=rep(c('#2ea582','#d95f02'),each=col_num[1]/2),
                       point.padding=unit(1, "lines"),size = 8#标签字体大小
      )
    ggsave(file = paste0('2.Differential_analyse/',as.character(group[i,2]),'_vs_',
                         as.character(group[i,3]),'/opls/',as.character(group[i,2]),'_vs_',
                         as.character(group[i,3]),
                         '_Volcano_map.pdf'),
           width = 12,height = 12)
    topptx(pp7,
           paste0('2.Differential_analyse/',as.character(group[i,2]),'_vs_',
                  as.character(group[i,3]),'/opls/',as.character(group[i,2]),'_vs_',
                  as.character(group[i,3]),
                  '_Volcano_map.pptx'),
           width = 12,height = 12)








  }

}
