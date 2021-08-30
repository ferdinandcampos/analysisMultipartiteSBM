#' Plot deg-deg inversion
#'
#' @param x à faire
#' @param y à faire
#' @return à faire
#' @examples
#' à faire
#'
#' @export
plotDegInversion=function(analyse,
                          Indices=as.numeric(names(table(analyse$network_id))),seuilPRb_hidecluster=0){
  x11(width = 18,height = 9)
  par(mfrow=c(1,2),mar = c(4, 4.3, 1, 1))
  taille_lab = 1.3
  taille_axes = 1
  taille_pointstext = 1
  taille_legende = 1.2
  taille_soustitres = 1
  for (FG in c("plants","pollination")){
    if(FG=="plants"){
      FG_name="plantes"
    }else if(FG=="pollination"){
      FG_name="pollinisateurs"
    }
    scatter_type="b"
    metrique_x<-"cKnn_inter1"
    metrique_y<-"c_inter1"
    x_lim<-range(analyse[[metrique_x]],na.rm=TRUE)
    y_lim<-range(analyse[[metrique_y]],na.rm=TRUE)
    # x_lim<-c(-0.4,0.4)
    # y_lim<-c(-0.5,2)
    palette<-colorRampPalette(brewer.pal(9, "Set1"))(17)

    index=Indices[1]
    data_plants_index<-analyse[analyse$network_id==index & analyse$FG==FG,]
    data_plants_index<-data_plants_index[order(data_plants_index[[metrique_x]]), ]
    plot(data_plants_index[[metrique_x]],data_plants_index[[metrique_y]],col=alpha(palette[1],data_plants_index[["PR_inter1_sym"]]*(data_plants_index[["PR_inter1_sym"]]>seuilPRb_hidecluster)),
         pch=17,type="p",lwd=1,lty=1,cex=exp(data_plants_index[["cluster_pi"]]),
         xlim=x_lim,ylim=y_lim,xlab=variables_to_Latex(metrique_x),ylab=variables_to_Latex(metrique_y),
         main=paste("Inversion degr?-degr? des voisins :",FG_name),
         cex.lab = taille_lab, cex.axis = taille_axes,cex.main=taille_soustitres)
    lines(data_plants_index[[metrique_x]],data_plants_index[[metrique_y]],col=alpha(palette[index],0.3),type="l",lwd=1,lty=1)

    for (index in Indices[-c(Indices[1])]){
      data_plants_index<-analyse[analyse$network_id==index & analyse$FG==FG,]
      data_plants_index<-data_plants_index[order(data_plants_index[[metrique_x]]), ] #tri par ordre croissant
      #lines(data_plants_index[[metrique_x]],data_plants_index[[metrique_y]],col=alpha(palette[index],data_plants_index[["PR_inter1_sym"]]),pch=17,type=scatter_type,lwd=1,lty=2,cex=exp(data_plants_index[["cluster_pi"]]))
      lines(data_plants_index[[metrique_x]],data_plants_index[[metrique_y]],col=alpha(palette[index],data_plants_index[["PR_inter1_sym"]]*(data_plants_index[["PR_inter1_sym"]]>seuilPRb_hidecluster)),
              cex=exp(data_plants_index[["cluster_pi"]]),pch=17,type="p",lwd=1,lty=1)
      lines(data_plants_index[[metrique_x]],data_plants_index[[metrique_y]],col=alpha(palette[index],0.3),type="l",lwd=1,lty=1)
    }
    if(FG=="plants"){
      legend("topright", legend = c(0.1,0.25,0.5,0.75,1), col = alpha(palette[1],c(0.1,0.25,0.5,0.75,1)),
             cex=taille_legende, pch = 17, bty = "n", title=variables_to_Latex("PR_inter1_sym"),inset=0.02,y.intersp=0.8,x.intersp=0.4)}
  }
}

# library(RColorBrewer)
# library(ggplot2)
