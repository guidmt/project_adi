plot_annotationGMT = function(annotated_regions, annotated_random, annotation_order = NULL,
                           plot_title, x_label, y_label, quiet = FALSE,perc=TRUE,select=T,list_selection) {

  # Tidy the GRanges into a tbl_df for use with dplyr functions
  annotated_regions = as.data.frame(annotated_regions, row.names = NULL)
  
  ########################################################################
  # Order and subset the annotations
  annotated_regions = annotatr:::subset_order_tbl(tbl = annotated_regions, col='annot.type', col_order=annotation_order)
  
  ########################################################################
  # If a region has multiple annotation types that are the same, count only one
  # from each type of annotation
  annotated_regions = dplyr::distinct_(
    dplyr::ungroup(annotated_regions),
    .dots=c('seqnames', 'start', 'end', 'annot.type'), .keep_all=TRUE)
  
  # Do particular things if annotated_random isn't NULL
  if(!missing(annotated_random)) {
    # Tidy the GRanges into a tbl_df for use with dplyr functions
    annotated_random = as.data.frame(annotated_random, row.names = NULL)
    
    # Order and subset the randomized annotations
    annotated_random = annotatr:::subset_order_tbl(tbl = annotated_random, col='annot.type', col_order=annotation_order)
    
    # If a region has multiple annotation types that are the same, count only one
    # from each type of annotation
    annotated_random = dplyr::distinct_(
      dplyr::ungroup(annotated_random),
      .dots=c('seqnames', 'start', 'end', 'annot.type'), .keep_all=TRUE)
    
    # Combine the tbl_dfs in preparation for visualization
    annotated_regions = dplyr::bind_rows("Data" = annotated_regions, "Random Regions" = annotated_random, .id = 'data_type')
  }
  
  ########################################################################
  # Construct the plot
  
  # Make the base ggplot
  # NOTE: binwidth may need to be a parameter
  annotated_regions$status<-as.factor(annotated_regions$status)
  
  tab_feat<- as.data.frame.matrix(table(annotated_regions[,c("annot.type","status")]))
  
  tab_feat<- tab_feat[order(tab_feat[,2],decreasing=F),]
  
  if(perc==TRUE){
  # 
  # tab_feat[,1]<-tab_feat[,1]/sum(tab_feat[,1])  
  # tab_feat[,2]<-tab_feat[,2]/sum(tab_feat[,2])  
  
    tab_feat<-(tab_feat/sum(rowSums(tab_feat)))*100
  }
  
  annotated_regions2<-reshape2::melt(as.matrix(tab_feat))
  colnames(annotated_regions2)<-c("annot.type","status","value")
  
  require(dplyr)
    
  annotated_regions2<-annotated_regions2 %>% mutate(value = ifelse(status == "DMRs-up",value,-1*value))
  
  breaks_values <-  pretty(annotated_regions2$value,n=10)
 
  if(select==TRUE){

	    annotated_regions2<-annotated_regions2[which(annotated_regions2$annot.type %in% list_selection),]

  } else {

	    annotated_regions2<-annotated_regions2

   }

  # https://www.onceupondata.com/2019/01/25/ggplot2-divergent-bars/
    
   p<- annotated_regions2 %>%
      ggplot(aes(x = annot.type, y = value, fill = status))+
      geom_hline(yintercept = 0) +
      geom_bar(stat = "identity",colour="black")+
      coord_flip()+
      scale_y_continuous(breaks = breaks_values,
                         labels = abs(breaks_values))+theme_classic()+scale_fill_manual(values = c("blue2", "red2"))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
   print(p)
  
}

#
# Define a function to create a piechart
#

plot_annotationPieChart_GMT = function(annotated_regions, annotated_random, annotation_order = NULL,
                              plot_title, x_label, y_label, quiet = FALSE) {
  
  # Tidy the GRanges into a tbl_df for use with dplyr functions
  annotated_regions = as.data.frame(annotated_regions, row.names = NULL)
  
  ########################################################################
  # Order and subset the annotations
  annotated_regions = annotatr:::subset_order_tbl(tbl = annotated_regions, col='annot.type', col_order=annotation_order)
  
  ########################################################################
  # If a region has multiple annotation types that are the same, count only one
  # from each type of annotation
  annotated_regions = dplyr::distinct_(
    dplyr::ungroup(annotated_regions),
    .dots=c('seqnames', 'start', 'end', 'annot.type'), .keep_all=TRUE)
  
  # Do particular things if annotated_random isn't NULL
  if(!missing(annotated_random)) {
    # Tidy the GRanges into a tbl_df for use with dplyr functions
    annotated_random = as.data.frame(annotated_random, row.names = NULL)
    
    # Order and subset the randomized annotations
    annotated_random = annotatr:::subset_order_tbl(tbl = annotated_random, col='annot.type', col_order=annotation_order)
    
    # If a region has multiple annotation types that are the same, count only one
    # from each type of annotation
    annotated_random = dplyr::distinct_(
      dplyr::ungroup(annotated_random),
      .dots=c('seqnames', 'start', 'end', 'annot.type'), .keep_all=TRUE)
    
    # Combine the tbl_dfs in preparation for visualization
    annotated_regions = dplyr::bind_rows("Data" = annotated_regions, "Random Regions" = annotated_random, .id = 'data_type')
  }
  
  ########################################################################
  # Construct the plot
  
  # Make the base ggplot
  # NOTE: binwidth may need to be a parameter
  annotated_regions$status<-as.factor(annotated_regions$status)
  
  tab_feat<- as.data.frame.matrix(table(annotated_regions[,c("annot.type","status")]))
  
  tab_feat<-tab_feat[order(tab_feat[,2],decreasing=F),]
  
   tab_feat[,1]<-tab_feat[,1]/sum(tab_feat[,1])  
   tab_feat[,2]<-tab_feat[,2]/sum(tab_feat[,2])  
  
  #tab_feat[,1]<-tab_feat[,1]/sum(rowSums(tab_feat))
  # tab_feat[,2]<-tab_feat[,2]/sum(rowSums(tab_feat))
  
  annotated_regions2<-reshape2::melt(as.matrix(tab_feat))
  colnames(annotated_regions2)<-c("annot.type","status","value")
  
  df2 <- annotated_regions2 %>% 
    arrange(status, desc(annot.type)) %>% 
    group_by(status) %>% 
    mutate(text_y = cumsum(value) - value/2)
  
  df2$value<-round(df2$value,2)
  # 
  # p1<-ggplot(data = df2, aes(x = "", y = value, fill = annot.type)) +
  #   geom_bar(stat = "identity") +scale_fill_manual(values = c("#6f6456", "#cddc49", "#cb7e94", "#e94b30", "#fee659", "#a1cfdd", "#749fd0", "#4eba06"))+
  #   geom_label_repel(aes(label = value, y = text_y)) + guides(fill = guide_legend(title = "Status"), nudge_x = 1) +
  #   facet_grid(. ~ status) +theme_void()+
  #   coord_polar(theta = "y")
  # 
  p1<-ggplot(df2, aes("", value, fill = annot.type)) +
    geom_bar(width = 1, size = 1, color = "white", stat = "identity") +
    coord_polar(theta = "y") +
    geom_label_repel(aes(label = value, y = text_y)) + guides(fill = guide_legend(title = "Status"), nudge_x = 1) +
    guides(fill = guide_legend(reverse = TRUE)) +
    scale_fill_manual(values = c("#6f6456", "#cddc49", "#cb7e94", "#e94b30", "#fee659", "#a1cfdd", "#749fd0", "#4eba06"))+
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())+facet_wrap(~ status)
  
  p2<-ggplot(annotated_regions2, aes("", value, fill = annot.type)) +
    geom_bar(width = 1, size = 1, color = "white", stat = "identity") +
    coord_polar("y") +
    guides(fill = guide_legend(reverse = TRUE)) +
    scale_fill_manual(values = c("#6f6456", "#cddc49", "#cb7e94", "#e94b30", "#fee659", "#a1cfdd", "#749fd0", "#4eba06"))+
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())+facet_wrap(~ status)
 
  #
  #single pie-chart
  #
  tab_feat<- as.data.frame.matrix(table(annotated_regions[,c("annot.type","status")]))
  
  tab_feat<-data.frame(rownames(tab_feat),rowSums(tab_feat))
  colnames(tab_feat)[1:2]<-c("annot.type","value")
  
  tab_feat<-tab_feat[order(tab_feat[,2],decreasing=F),]
  
  tab_feat$value<-round(tab_feat$value/sum(tab_feat$value),2)
  
  annotated_regions2<-tab_feat
  
  labs=as.character(round(annotated_regions2$value,3))
  
  p3<- ggpie(annotated_regions2, "value", label = labs,
             fill = "annot.type", color = "white",
             palette = c("#6f6456", "#cddc49", "#cb7e94", "#e94b30", "#fee659", "#a1cfdd", "#749fd0", "#4eba06"))
             
    print(p1)
    print(p2)
    print(p3)
}

#6f6456, #cddc49, #cb7e94, #e94b30, #fee659, #a1cfdd, #749fd0, #4eba06
