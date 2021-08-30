#' Plot network with blocks
#'
#' @param x à faire
#' @param y à faire
#' @return à faire
#' @examples
#' à faire
#' @import igraph sbm
#' @export
plotNetworkTripartiteSBM=function(myMSBM,seuil_hesitation=0.9,plot_types=1,
                                  plot_graph_no_partites=TRUE,show_titles=TRUE,
                                  nodes_size=3,nodes_pie=TRUE){
  Incidence_matrices<-formatSBMtoIncidenceMatrices(myMSBM$networkData)
  partites<-myMSBM$dimLabels
  nb_parties=length(partites)
  linking_set<-partites[1]
  interactions<-list(partites[2],partites[3])
  #Les noms pourraient ?tre r?cup?r?s depuis mySMBM directement mais visiblement disparaissent quand un seul cluster pour l'instant
  linking_names<-rownames(Incidence_matrices[[1]])
  interaction_1_names<-colnames(Incidence_matrices[[eval(interactions[[1]])]])
  interaction_2_names<-colnames(Incidence_matrices[[eval(interactions[[2]])]])
  #### Il semblerait qu'on puisse en fait les r?cup?rer non pas dans membership mais l? dedans sans passer par l'incidence
  # linking_names<-rownames(myMSBM2$networkData[[1]]$networkData)
  # interaction_1_names<-colnames(myMSBM2$networkData[[1]]$networkData)
  # interaction_2_names<-colnames(myMSBM2$networkData[[2]]$networkData)

  #Cr?ation du graph
  # g<-define_graph_from_indicence_matrice(Incidence_matrices) #? remplacer plus tard, car need des matrices m_pi et m_hp en dessous
  m_pi = as.matrix(Incidence_matrices[[1]])
  dimnames(m_pi)<-list(
    consuming=rownames(Incidence_matrices[[1]]),
    consumed=colnames(Incidence_matrices[[1]]))
  m_hp =  as.matrix(Incidence_matrices[[2]])
  dimnames(m_hp)<-list(
    consumed=rownames(Incidence_matrices[[2]]),
    consuming=colnames(Incidence_matrices[[2]]))
  m_hp<-t(m_hp)

  el_pi <- as.data.frame(as.table(m_pi), stringsAsFactors = FALSE)
  el_hp <- as.data.frame(as.table(m_hp), stringsAsFactors = FALSE)
  el <- rbind(el_pi, el_hp)
  g <- graph.data.frame( el[el$Freq != 0 , ]  )
  V(g)$type <- substr(V(g)$name, 1, 1)

  ##D?finition d'un layout multipartite personnalis?
  nb_c_int2<- myMSBM$nbBlocks[[3]]
  nb_c_link<-myMSBM$nbBlocks[[1]]
  nb_c_int1<-myMSBM$nbBlocks[[2]]
  all_nodes_clusters=do.call(c,list(myMSBM$memberships[[eval(interactions[[2]])]],myMSBM$memberships[[eval(linking_set)]]+nb_c_int2,myMSBM$memberships[[eval(interactions[[1]])]]+nb_c_link+nb_c_int2))
  all_nodes_probMemberships<-unlist(c(apply(myMSBM$probMemberships[[eval(interactions[[2]])]],1, list),
                                      apply(myMSBM$probMemberships[[eval(linking_set)]],1, list),
                                      apply(myMSBM$probMemberships[[eval(interactions[[1]])]],1, list)),
                                    recursive = FALSE)
  names(all_nodes_clusters)<-do.call(c,list(interaction_2_names,linking_names,interaction_1_names))
  names(all_nodes_probMemberships)<-do.call(c,list(interaction_2_names,linking_names,interaction_1_names))
  #all_nodes_clusters=do.call(c,list(myMSBM$memberships[[3]],myMSBM$memberships[[1]]+myMSBM$nbBlocks[[3]],myMSBM$memberships[[2]]+myMSBM$nbBlocks[[1]]+myMSBM$nbBlocks[[3]]))
  permutation_selon_clusters<-Matrix::invPerm(order(all_nodes_clusters[names(V(g))]))
  g=permute(g,permutation_selon_clusters)

  #Couleurs
  couleurs_clusters<-all_nodes_clusters[names(V(g))]
  names(couleurs_clusters)<-NULL
  couleurs_clusters_pie<-as.list(all_nodes_clusters[names(V(g))])
  if (nb_c_int2<=8){
    palette_int2<-RColorBrewer::brewer.pal(max(nb_c_int2+1,3), "BuPu")
  } else {palette_int2<-grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "BuPu"))(nb_c_int2+1)}

  if (nb_c_link<=8){
    palette_link<-RColorBrewer::brewer.pal(max(nb_c_link+1,3), "PiYG")
  } else{palette_link<-grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "PiYG"))(nb_c_link+1)}

  if (nb_c_int1<=8){
    palette_int1<-RColorBrewer::brewer.pal(max(nb_c_int1+1,3), "OrRd")
  } else{palette_int1<-grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "OrRd"))(nb_c_int1+1)}

  palette_perso=c(palette_int2[2:(nb_c_int2+1)],palette_link[2:(nb_c_link+1)],palette_int1[2:(nb_c_int1+1)])
  for (i in 1:(nb_c_int2+nb_c_link+nb_c_int1)){
    couleurs_clusters[couleurs_clusters==i]<-palette_perso[i]
  }

  couleurs_clusters_pie[couleurs_clusters_pie %in% 1:nb_c_int2]<-list(palette_perso[1:nb_c_int2])
  couleurs_clusters_pie[couleurs_clusters_pie %in% (nb_c_int2+1):(nb_c_int2+nb_c_link)]<-list(palette_perso[(nb_c_int2+1):(nb_c_int2+nb_c_link)])
  couleurs_clusters_pie[couleurs_clusters_pie %in% (nb_c_int2+nb_c_link+1):(nb_c_int2+nb_c_link+nb_c_int1)]<-list(palette_perso[(nb_c_int2+nb_c_link+1):(nb_c_int2+nb_c_link+nb_c_int1)])
  #Pie
  probMemberships_pie<-all_nodes_probMemberships[names(V(g))]
  #Nodes hesitating
  nodes_hesitation<-(rapply(probMemberships_pie,max)<seuil_hesitation)
  # couleurs_borders_hesitation<-as.list(all_nodes_clusters[names(V(g))])
  couleurs_borders_hesitation<-all_nodes_clusters[names(V(g))]
  couleurs_borders_hesitation[]<-"black"
  couleurs_borders_hesitation[nodes_hesitation]<-"red"

  font_hesitation<-all_nodes_clusters[names(V(g))]
  font_hesitation[]<-1
  font_hesitation[nodes_hesitation]<-2

  # label_dist_hesitation<-all_nodes_clusters[names(V(g))]
  # label_dist_hesitation[]<-0
  # label_dist_hesitation[nodes_hesitation]<-1

  #Positions
  positions_FG <- match(V(g)$type, c(species_type_letter(partites[3]), species_type_letter(partites[1]), species_type_letter(partites[2])) )
  #positions_clusters<-do.call(c,list(1:nrow(m_hp),1:nrow(m_pi),1:ncol(m_pi))) #all vertices
  positions_clusters<-do.call(c,list(which(!apply(m_hp == 0, 1, all)),which(!apply(cbind(m_pi,t(m_hp)) == 0, 1, all)),which(!apply(m_pi == 0, 2, all)))) #only vertices with at least 1 interaction
  names(positions_clusters)<-NULL
  layout_multipartite <- layout_with_fr(g, minx=positions_FG, maxx=positions_FG,miny=positions_clusters,maxy=positions_clusters)
  layout_multipartite2 <- layout_with_fr(g, minx=positions_FG, maxx=positions_FG) #SANS CONFIG MANUEL SUR AXE y des CLUSTERS


  if (1 %in% plot_types){
    #Affichage (layout personnalis? en x et y)
    x11()
    if (plot_graph_no_partites){
      par(mfrow=c(1,2))
      if (nodes_pie){plot(g, edge.width=E(g)$Freq, edge.arrow.size=0,vertex.size=6,vertex.label.cex=0.8,vertex.label.degree=0,
                          vertex.shape="pie",vertex.pie=probMemberships_pie, vertex.pie.color=couleurs_clusters_pie,
                          vertex.frame.color=couleurs_borders_hesitation,vertex.label.font=font_hesitation,vertex.label.color=couleurs_borders_hesitation)
      }
      else{plot(g, vertex.color=couleurs_clusters, edge.width=E(g)$Freq, edge.arrow.size=0,vertex.size=5,vertex.label.cex=0.8)}
      if(show_titles){title("Network (no partites)")}
    }

    if (nodes_pie){
      plot(g, vertex.color=couleurs_clusters, layout=layout_multipartite, edge.width=E(g)$Freq, edge.arrow.size=0,vertex.size=nodes_size,vertex.label.dist=1,vertex.label.degree=pi,vertex.label.cex=0.8,
           vertex.frame.color=couleurs_borders_hesitation,vertex.label.color=couleurs_borders_hesitation)
    }
    else{plot(g, vertex.color=couleurs_clusters, layout=layout_multipartite, edge.width=E(g)$Freq, edge.arrow.size=0,vertex.size=nodes_size,vertex.label.dist=1,vertex.label.degree=pi,vertex.label.cex=0.8)}
    if(show_titles){title("Tripartite network (1)")}
  }

  if (2 %in% plot_types){
    #Affichage (layout personnalis? en x, auto en y)
    x11()
    if (plot_graph_no_partites){
      par(mfrow=c(1,2))
      if (nodes_pie){plot(g, edge.width=E(g)$Freq, edge.arrow.size=0,vertex.size=6,vertex.label.cex=0.8,vertex.label.degree=0,
                          vertex.shape="pie",vertex.pie=probMemberships_pie, vertex.pie.color=couleurs_clusters_pie,
                          vertex.frame.color=couleurs_borders_hesitation,vertex.label.font=font_hesitation,vertex.label.color=couleurs_borders_hesitation)
      }
      else{plot(g, vertex.color=couleurs_clusters, edge.width=E(g)$Freq, edge.arrow.size=0,vertex.size=5,vertex.label.cex=0.8)}
      if(show_titles){title("Network (no partites)")}
    }


    if (nodes_pie){
      plot(g, vertex.color=couleurs_clusters, layout=layout_multipartite2, edge.width=E(g)$Freq, edge.arrow.size=0,vertex.size=nodes_size,vertex.label.dist=1,vertex.label.degree=pi,vertex.label.cex=0.8,
           vertex.frame.color=couleurs_borders_hesitation,vertex.label.color=couleurs_borders_hesitation)
    }
    else{plot(g, vertex.color=couleurs_clusters, layout=layout_multipartite2, edge.width=E(g)$Freq, edge.arrow.size=0,vertex.size=nodes_size,vertex.label.dist=1,vertex.label.degree=pi,vertex.label.cex=0.8)}
    if(show_titles){title("Tripartite network (2)")}

  }


  # if (nodes_pie){
  #   x11()
  #   plot(g, edge.width=E(g)$Freq, edge.arrow.size=0,vertex.size=6,vertex.label.cex=0.8,vertex.label.degree=0,
  #        vertex.shape="pie",vertex.pie=probMemberships_pie, vertex.pie.color=couleurs_clusters_pie,
  #        vertex.frame.color=couleurs_borders_hesitation,vertex.label.font=font_hesitation,vertex.label.color=couleurs_borders_hesitation,vertex.label.dist=label_dist_hesitation)
  #
  # }
  #
}
