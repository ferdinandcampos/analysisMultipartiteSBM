#' Plot FC
#'
#' @param x à faire
#' @param y à faire
#' @return à faire
#' @examples
#' à faire
#' @import ggplot2
#' @export
plotFunctionalCartography=function(analyse,
                                   plotdim="2d",highlight.3d=FALSE,Indices=as.numeric(names(table(analyse$network_id)))){
  ############Functional Cartography 2d advanced : c_inter1,2 et PR en couleur
  scatter_type="p"
  metrique_x<-"c_inter1"
  metrique_y<-"c_inter2"
  metrique_color<-"PR_inter1_sym"
  x_lim<-range(analyse[[metrique_x]],na.rm=TRUE)
  y_lim<-range(analyse[[metrique_y]],na.rm=TRUE)
  # x_lim<-c(-0.5,2)
  # y_lim<-c(-0.5,2)

  FG_plotted<-"plants"
  data_plants<-analyse[analyse$FG==FG_plotted,]

  nb_divisions=6
  if(metrique_color=="cKnn_inter1" | metrique_color=="cKnn_inter2" | metrique_color=="cKnn_tot"){palette<-rev(RColorBrewer::brewer.pal(nb_divisions, "RdYlBu"))}
  if(metrique_color=="PR_inter1" | metrique_color=="PR_inter1_sym"){palette<-rev(RColorBrewer::brewer.pal(nb_divisions, "Spectral"))}
  values_colors<-rep("#000000",length(data_plants[[metrique_color]])) #NB si on modifie directement values colors en mettant des codes hexa au milieu des numeric ?a pose souci. faire deux objets
  m_values_colors<-min(data_plants[[metrique_color]])
  M_values_colors<-max(data_plants[[metrique_color]])
  #Pour forcer la sym?trie des couleurs et voir la s?paration >0 et <0 #sinon faire deux palettes diff?rentes..
  if(metrique_color=="cKnn_inter1" | metrique_color=="cKnn_inter2" | metrique_color=="cKnn_tot"){center_values_colors<-0}
  if(metrique_color=="PR_inter1" | metrique_color=="PR_inter1_sym"){center_values_colors<-0.5}
  abs_range_values_colors<-max(abs(m_values_colors-center_values_colors),abs(M_values_colors-center_values_colors))
  m_values_colors<-center_values_colors-abs_range_values_colors
  M_values_colors<-center_values_colors+abs_range_values_colors
  for (i in 1:nb_divisions){
    values_colors[data_plants[[metrique_color]] >= (m_values_colors+(i-1)*(M_values_colors-m_values_colors)/nb_divisions) & data_plants[[metrique_color]] < (m_values_colors+i*(M_values_colors-m_values_colors)/nb_divisions)]<-palette[i]
  }
  values_colors[data_plants[[metrique_color]]==M_values_colors]<-palette[nb_divisions]
  #pch_admissibles<-c(0,1,2,3,4,5,6,8,9,10,11,12,13,15,16,17,18)
  #pch_reseaux<-pch_admissibles[data_plants[data_plants$network_id %in% Indices, "network_id"]]
  pch_reseaux<-17

  x11(width = 11,height = 9)
  if(plotdim=="2d"){
    taille_blocks_markers<-1.6*exp(data_plants[["cluster_pi"]])
    lines_treshold1<-0
    lines_treshold2<-0

    plot(data_plants[[metrique_x]],data_plants[[metrique_y]],col=values_colors,pch=pch_reseaux,type=scatter_type,lwd=1,lty=2,cex=taille_blocks_markers,
         xlim=x_lim,ylim=y_lim,
         xlab=variables_to_Latex(metrique_x),ylab=variables_to_Latex(metrique_y),
         main="Functional Cartography")
    lines(c(lines_treshold1,lines_treshold1),y_lim+c(-1,1),lwd=1,lty=2,col="gray")
    lines(x_lim+c(-1,1),c(lines_treshold2,lines_treshold2),lwd=1,lty=2,col="gray")
    text(data_plants[[metrique_x]], data_plants[[metrique_y]],
         labels=as.character(data_plants[data_plants$network_id %in% Indices, "network_id"]),
         cex=0.7,col="black")
    legend("topright", legend = paste0("[",round(m_values_colors+(1:nb_divisions-1)*(M_values_colors-m_values_colors)/nb_divisions,2),",",round(m_values_colors+(1:nb_divisions)*(M_values_colors-m_values_colors)/nb_divisions,2),"]"),
           col = palette, pch = 17, bty = "n",title=variables_to_Latex(metrique_color),inset=0.02)

  } else if (plotdim=="3d"){
    metrique_z<-metrique_color
    if (highlight.3d){
      s3d<-scatterplot3d::scatterplot3d(data_plants[[metrique_x]],
                         data_plants[[metrique_y]],
                         data_plants[[metrique_z]],
                         main="Functionnal 3d cartography with highlighting3d colors",
                         xlab = metrique_x,
                         ylab = metrique_y,
                         zlab = metrique_z,
                         pch = pch_reseaux,
                         type="h",
                         cex.symbols=2.5,
                         highlight.3d=TRUE,
                         box=TRUE,
                         xlim=x_lim,ylim=y_lim,
                         lty.hplot=1,
                         lwd=1)
      text(s3d$xyz.convert(data_plants[[metrique_x]],
                           data_plants[[metrique_y]],
                           data_plants[[metrique_z]]), labels=as.character(data_plants[data_plants$network_id %in% Indices, "network_id"]),
           cex=0.6,col="white")
    } else {
      s3d<-scatterplot3d::scatterplot3d(data_plants[[metrique_x]],
                         data_plants[[metrique_y]],
                         data_plants[[metrique_z]],
                         main="3d Functional Cartography",
                         xlab = metrique_x,
                         ylab = metrique_y,
                         zlab = metrique_z,
                         pch = pch_reseaux,
                         type="h",
                         cex.symbols=2.5,
                         color=values_colors,
                         box=TRUE,
                         xlim=x_lim,ylim=y_lim,
                         lty.hplot=1,
                         lwd=1
      )
      legend("topright", legend = paste0("[",round((1:nb_divisions-1)/nb_divisions,2),",",round((1:nb_divisions)/nb_divisions,2),"]"),
             col = palette, pch = 17, bty = "n",title=variables_to_Latex(metrique_color),inset=0.02)
      text(s3d$xyz.convert(data_plants[[metrique_x]],
                           data_plants[[metrique_y]],
                           data_plants[[metrique_z]]), labels=as.character(data_plants[data_plants$network_id %in% Indices, "network_id"]),
           cex=0.6,col="white")

    }
  }
}
