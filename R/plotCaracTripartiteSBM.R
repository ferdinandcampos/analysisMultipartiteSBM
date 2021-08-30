#' Plot carac
#'
#' @param x à faire
#' @param y à faire
#' @return à faire
#' @examples
#' à faire
#'
#' @export
plotCaracTripartiteSBM=function(caracterisation,
                                plotIndicesblocs=1:nrow(caracterisation$data$analysePlants)){

  km_centers<-caracterisation$centers
  km_clusters<-caracterisation$memberships
  analysePlants<-caracterisation$data$analysePlants
  km_centering<-caracterisation$data$km_centering
  km_scaling<-caracterisation$data$km_scaling
  res.pca<-caracterisation$result.PCA
  fit.MDS <- caracterisation$result.MDS

  metrique_color<-"PR_inter1_sym"
  scatter_type="p"
  lines_treshold1<-0
  lines_treshold2<-0

  pch_reseaux<-17

  #Palette PR
  nb_divisions=6
  if(metrique_color=="cKnn_inter1" | metrique_color=="cKnn_inter2" | metrique_color=="cKnn_tot"){palette<-rev(RColorBrewer::brewer.pal(nb_divisions, "RdYlBu"))}
  if(metrique_color=="PR_inter1" | metrique_color=="PR_inter1_sym"){palette<-rev(RColorBrewer::brewer.pal(nb_divisions, "Spectral"))}
  values_colors<-rep("#000000",length(analysePlants[[metrique_color]])) #NB si on modifie directement values colors en mettant des codes hexa au milieu des numeric ?a pose souci. faire deux objets
  m_values_colors<-min(analysePlants[[metrique_color]])
  M_values_colors<-max(analysePlants[[metrique_color]])
  #Pour forcer la sym?trie des couleurs et voir la s?paration >0 et <0 #sinon faire deux palettes diff?rentes..
  if(metrique_color=="cKnn_inter1" | metrique_color=="cKnn_inter2" | metrique_color=="cKnn_tot"){center_values_colors<-0}
  if(metrique_color=="PR_inter1" | metrique_color=="PR_inter1_sym"){center_values_colors<-0.5}
  abs_range_values_colors<-max(abs(m_values_colors-center_values_colors),abs(M_values_colors-center_values_colors))
  m_values_colors<-center_values_colors-abs_range_values_colors
  M_values_colors<-center_values_colors+abs_range_values_colors
  for (i in 1:nb_divisions){
    values_colors[analysePlants[[metrique_color]] >= (m_values_colors+(i-1)*(M_values_colors-m_values_colors)/nb_divisions) & analysePlants[[metrique_color]] < (m_values_colors+i*(M_values_colors-m_values_colors)/nb_divisions)]<-palette[i]
  }
  values_colors[analysePlants[[metrique_color]]==M_values_colors]<-palette[nb_divisions]

  #Palette k-means
  if (!study_insectes_only){palette2<-colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(nb_clusters+2)
  }else{palette2<-colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(6+2)
  palette2<-palette2[-c(5,6)]}
  values_colors_reseaux<-palette2[km_clusters]


  ##################################MAIN PLOT

  x11()
  par(mfrow=c(2,2),mar = c(3.8, 4.3, 0.3, 0.3))
  taille_lab = 1.3
  taille_axes = 1
  taille_axesmetriques = 1.5
  taille_points = 3
  taille_pointstext = 1
  taille_legende = 1.2
  taille_soustitres = 1.2

  ###PCA
  for (metriques_xy in list(c("PC1","PC2"),c("PC2","PC3"))){
    metrique_x<-metriques_xy[1]
    metrique_y<-metriques_xy[2]

    x_lim<-1*range(res.pca$x[,metrique_x],na.rm=TRUE)
    y_lim<-1*range(res.pca$x[,metrique_y],na.rm=TRUE)
    scalefactor<-0.95*min(abs(c(x_lim,y_lim)))/max(abs(res.pca$rotation[,c(metrique_x,metrique_y)]))
    x_lim<-c(-max(abs(x_lim)),max(abs(x_lim)))
    y_lim<-c(-max(abs(y_lim)),max(abs(y_lim)))
    # x_lim<-range(res.pca$x[,metrique_x],na.rm=TRUE)
    # y_lim<-range(res.pca$x[,metrique_y],na.rm=TRUE)
    #Palette k-means
    if (!study_insectes_only){palette2<-colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(nb_clusters+2)
    }else{palette2<-colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(6+2)
    palette2<-palette2[-c(5,6)]}
    values_colors_reseaux<-palette2[km_clusters]

    ###Map with k-means clusters
    plot(res.pca$x[,metrique_x],res.pca$x[,metrique_y],col=values_colors_reseaux,pch=pch_reseaux,type=scatter_type,lwd=1,lty=2,cex=taille_points,
         xlim=x_lim,ylim=y_lim,xlab=variables_to_Latex(metrique_x),ylab=variables_to_Latex(metrique_y),
         #main=paste("PCA with",nb_clusters,"kmeans clusters,",weighted_by_pi_text),
         cex.lab = taille_lab, cex.axis = taille_axes)
    lines(c(lines_treshold1,lines_treshold1),y_lim+c(-1,1),lwd=1,lty=2,col="gray")
    lines(x_lim+c(-1,1),c(lines_treshold2,lines_treshold2),lwd=1,lty=2,col="gray")

    km_centers_normalised<-t(apply(km_centers,1,function(x){x*km_scaling})+km_centering) #weighted k-means unscaling
    km_centers_normalised<-t(apply(km_centers_normalised,1,function(x){x-res.pca$center})/res.pca$scale) #PCA scaling
    km_centers_projected<-t(solve(res.pca$rotation)%*%t(km_centers_normalised)) #rotating
    #(km_centers_projected)<-colnames(res.pca$x)#Ne fonctionne
    #km_centers_projected<-t((solve(res.pca$rotation))%*%(apply(km_centers,1,function(x){x-res.pca$center})/res.pca$scale)) #Ne fonctionne pas
    points(km_centers_projected[,metrique_x],km_centers_projected[,metrique_y],col=palette2[1:nb_clusters],pch=19,cex=1.5)
    points(km_centers_projected[,metrique_x],km_centers_projected[,metrique_y],col="black",pch=8,cex=1)

    text(res.pca$x[,metrique_x], res.pca$x[,metrique_y]+0,
         labels=plotIndicesblocs,
         cex=taille_pointstext,col="white",font=2)
    legend("topright", legend = paste0("C",1:nb_clusters), col = palette2,text.col=palette2,title.col="black", pch = 17,
           bty = "n",title="Clusters",cex=taille_legende,y.intersp=0.8,x.intersp=0.4)

    for (i in 1:nrow(res.pca$rotation)){
      arrows(0, 0, scalefactor*res.pca$rotation[i,metrique_x], scalefactor*res.pca$rotation[i,metrique_y],
             code=2, lwd=1, length = 0.08)
      text(1.02*scalefactor*res.pca$rotation[i,metrique_x], 1.02*scalefactor*res.pca$rotation[i,metrique_y]-0.2,
           labels=variables_to_Latex(rownames(res.pca$rotation)[i]),cex=taille_axesmetriques,col="black")
    }
  }

  ###FC Map with k-means clusters
  metrique_x<-"c_inter1"
  metrique_y<-"c_inter2"
  x_lim<-range(analysePlants[[metrique_x]],na.rm=TRUE)
  y_lim<-range(analysePlants[[metrique_y]],na.rm=TRUE)
  # x_lim<-c(-0.5,2)
  # y_lim<-c(-0.5,2)


  plot(analysePlants[[metrique_x]],analysePlants[[metrique_y]],col=values_colors_reseaux,pch=pch_reseaux,type=scatter_type,lwd=1,lty=2,cex=taille_points,
       xlim=x_lim,ylim=y_lim,xlab=variables_to_Latex(metrique_x),ylab=variables_to_Latex(metrique_y),
       #main=paste("F. Cart with",nb_clusters,"kmeans clusters,",weighted_by_pi_text),
       cex.lab = taille_lab, cex.axis = taille_axes)
  lines(c(lines_treshold1,lines_treshold1),y_lim+c(-1,1),lwd=1,lty=2,col="gray")
  lines(x_lim+c(-1,1),c(lines_treshold2,lines_treshold2),lwd=1,lty=2,col="gray")
  # points(km_centers[,metrique_x],km_centers[,metrique_y],col=palette2[1:nb_clusters],pch=19,cex=1.5)
  # points(km_centers[,metrique_x],km_centers[,metrique_y],col="black",pch=8,cex=1)
  text(analysePlants[[metrique_x]], analysePlants[[metrique_y]],
       labels=plotIndicesblocs,
       cex=taille_pointstext,col="white",font=2)
  legend("topright", legend = paste0("C",1:nb_clusters), col = palette2,text.col=palette2,title.col="black", pch = 17,
         bty = "n",title="Clusters",cex=taille_legende,y.intersp=0.8,x.intersp=0.4)

  ###FC Map with PR colors
  plot(analysePlants[[metrique_x]],analysePlants[[metrique_y]],col=values_colors,pch=pch_reseaux,type=scatter_type,lwd=1,lty=2,cex=taille_points,
       xlim=x_lim,ylim=y_lim,xlab=variables_to_Latex(metrique_x),ylab=variables_to_Latex(metrique_y),
       #main=paste("Functional Cartography"),
       cex.lab = taille_lab, cex.axis = taille_axes)
  lines(c(lines_treshold1,lines_treshold1),y_lim+c(-1,1),lwd=1,lty=2,col="gray")
  lines(x_lim+c(-1,1),c(lines_treshold2,lines_treshold2),lwd=1,lty=2,col="gray")
  text(analysePlants[[metrique_x]]+0, analysePlants[[metrique_y]]+0,
       labels=as.character(analysePlants[analysePlants$network_id %in% Indices, "network_id"]),
       cex=taille_pointstext,col="black",font=2)
  legend("topright", legend = paste0("[",round(m_values_colors+(1:nb_divisions-1)*(M_values_colors-m_values_colors)/nb_divisions,2),",",round(m_values_colors+(1:nb_divisions)*(M_values_colors-m_values_colors)/nb_divisions,2),"]"),
         col = palette, pch = 17, bty = "n",title=variables_to_Latex(metrique_color),cex=taille_legende,y.intersp=0.8,x.intersp=0.4,title.adj=0.25)

  ###################################FIG ANNEXE

  ###Affichage multiple, par r?seaux
  plot_multiple=TRUE
  if(plot_multiple){
    x11(width = 18,height = 9)
    par(mfrow=c(3,6),mar = c(3.8, 4.3, 0.3, 0.3))
    #par(mfrow=c(1,2))
    # plot(analysePlants[[metrique_x]],analysePlants[[metrique_y]],col=values_colors_reseaux,pch=pch_reseaux,type=scatter_type,lwd=1,lty=2,cex=taille_blocks_markers,
    #      xlim=x_lim,ylim=y_lim,xlab=variables_to_Latex(metrique_x),ylab=variables_to_Latex(metrique_y),
    #      main=paste("F. Cart with",nb_clusters,"kmeans clusters,",weighted_by_pi_text))
    for (index in Indices){
      data_plants_index<-analysePlants[analysePlants$network_id==index,]
      taille_blocks_markers_index<-3*exp(analysePlants[analysePlants$network_id==index,"cluster_pi"])
      values_colors_index<-palette2[km_clusters[names(km_clusters) %in% rownames(data_plants_index)]]
      # x_lim<-range(data_plants_index[[metrique_x]],na.rm=TRUE)
      # y_lim<-range(data_plants_index[[metrique_y]],na.rm=TRUE)
      x_lim<-c(-0.5,2)
      y_lim<-c(-0.5,2)
      plot(data_plants_index[[metrique_x]],data_plants_index[[metrique_y]],col=values_colors_index,pch=pch_reseaux,type="p",lwd=1,lty=2,cex=taille_blocks_markers_index,
           xlim=x_lim,ylim=y_lim,xlab=variables_to_Latex(metrique_x),ylab=variables_to_Latex(metrique_y),
           #main=paste("Id. =",index),
           cex.lab = taille_lab, cex.axis = taille_axes,cex.main=taille_soustitres)
      lines(c(lines_treshold1,lines_treshold1),y_lim+c(-1,1),lwd=1,lty=2,col="gray")
      lines(x_lim+c(-1,1),c(lines_treshold2,lines_treshold2),lwd=1,lty=2,col="gray")
      text(data_plants_index[[metrique_x]], data_plants_index[[metrique_y]],
           labels=plotIndicesblocs[sapply(rownames(data_plants_index),function(x){which(rownames(analysePlants)==x)})],
           cex=taille_pointstext,col="white",font=2)
      legend("topright", legend = paste("Id. =",index), bty = "n",cex=taille_legende,text.font=2)
    }
  }

  ###Autre affichage : MDS
  x11()
  par(mar = c(3.8, 4.3, 1.5, 1))
  x <- fit.MDS $points[,1]
  y <- fit.MDS $points[,2]
  plot(x, y, xlab=variables_to_Latex("MDS1"), ylab=variables_to_Latex("MDS2"),
       main="Multidimensional Scaling (MDS)", type="p",col=values_colors_reseaux,pch=17,
       cex=taille_points,cex.lab = taille_lab, cex.axis = taille_axes)
  text(x, y, labels = plotIndicesblocs, cex=taille_pointstext,col="white",font=2)
  legend("topright", legend = paste0("C",1:nb_clusters), col = palette2,text.col=palette2,title.col="black", pch = 17,
         bty = "n",title="Clusters",cex=taille_legende,y.intersp=0.8,x.intersp=0.4)
}

# source("toolsFunctions.R")
