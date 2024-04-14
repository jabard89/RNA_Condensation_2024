# Main Style Sheet for Glauninger,Bard,Wong et al. 2023
theme_sedseq <- function(base_size=10) {
  theme_grey(base_size=base_size, base_family = "") %+replace%
    theme(panel.grid=element_blank(),
          axis.text = element_text(size = base_size),
          strip.text.x = element_text(size = base_size, margin = margin(t = 10)),
          strip.text.y = element_text(size = base_size, margin = margin(r = 10)),
          panel.background=element_blank(),
          axis.ticks=element_line(colour="grey20"),
          panel.border=element_rect(fill=NA),
          legend.background = element_blank(),
          legend.key.height = unit(3, "mm"),
          legend.key = element_blank(),
          strip.background = element_blank(),
          legend.box.spacing = unit(1, "mm"))
}


graycol <- "#333333cc"
orangecol <- "#cc5500cc"
bluecol <- "#0000aacc"
greencol <- "#22cc00cc"
purplecol <- "#cc22cccc"
cyancol <- "#2aa198cc"
redcol <- "#dc322fcc"
violetcol <- "#6c71c4cc"

label.levels <- c("other","glycolysis","ribosomal proteins","ribosome biogenesis","MSN2/4 targets","HSF1 targets")
label.cols <- c("other"="grey50",
                "ribosomal proteins"=palette.colors(palette="Paired")[2],
                "ribosome biogenesis"=palette.colors(palette="Paired")[10],
                "MSN2/4 targets"=palette.colors(palette="Paired")[8],
                "HSF1 targets"=palette.colors(palette="Paired")[6])

stress.cols <-  c("30C"="grey50",
                  "42C"=RColorBrewer::brewer.pal(n=7,"Oranges")[4],
                  "46C"=RColorBrewer::brewer.pal(n=7,"Oranges")[6],
                  "mock"="grey50",
                  "Azide\n0.5%"=RColorBrewer::brewer.pal(n=5,"Greens")[3],
                  "Azide\n0.8%"=RColorBrewer::brewer.pal(n=5,"Greens")[5],
                  "Ethanol\n5%"=RColorBrewer::brewer.pal(n=7,"Purples")[3],
                  "Ethanol\n7.5%"=RColorBrewer::brewer.pal(n=7,"Purples")[5],
                  "Ethanol\n10%"=RColorBrewer::brewer.pal(n=7,"Purples")[6],
                  "Ethanol\n15%"=RColorBrewer::brewer.pal(n=7,"Purples")[7])


plot_label_cross <- function(df,qlow=0.2,qhigh=0.8,xvar,yvar,cross_size=2){
  # example usage: plot_label_cross(df,ql=0.2,ql=0.8,"FC.vs.ctrl_RNA","FC.vs.ctrl_TE",cross_size=2)
  df_sum <- df %>% ungroup %>%
    select(ORF,label,Stress_label,sym(xvar),sym(yvar)) %>%
    pivot_longer(-c(ORF,label,Stress_label),names_to="Key",values_to="Value") %>%
    group_by(Stress_label,label,Key) %>%
    summarise(ql = quantile(Value,qlow,na.rm=T),
              qm = quantile(Value,0.5,na.rm=T),
              qh = quantile(Value,qhigh,na.rm=T)) %>%
    pivot_wider(values_from=c(ql,qm,qh),names_from="Key",names_sep=".")
  df_sum$label %>% unique %>% map(function(lab){
    list(
      geom_linerange(inherit.aes=F,data=df_sum %>% filter(label==!!lab),
                     aes(y=.data[[paste0("qm.",yvar)]],xmin=.data[[paste0("ql.",xvar)]],xmax=.data[[paste0("qh.",xvar)]],
                         color=label),
                     size=cross_size),
      geom_linerange(inherit.aes=F,data=df_sum %>% filter(label==!!lab),
                     aes(x=.data[[paste0("qm.",xvar)]],ymin=.data[[paste0("ql.",yvar)]],ymax=.data[[paste0("qh.",yvar)]],
                         color=label),
                     size=cross_size)
    )
  }) %>% list_flatten
}
