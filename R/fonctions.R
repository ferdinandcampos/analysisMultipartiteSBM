library(igraph)
library(Matrix)
library(RColorBrewer)



####################################PREPARE DATA

find_partites_type=function(dataset_names,index){
  
  if (startsWith( dataset_name, 'Ibanez')){
    interactions = list("pollination", "herbivory")
    linking_set='plants'
  } else if (startsWith( dataset_name, 'Sinohara')){
    interactions = list("pollination", "herbivory")
    linking_set='plants'
  } else if (startsWith( dataset_name, 'Melian')){
    if (str_detect( dataset_name, 'HSD')) {
      interactions = list("dispersion", "herbivory")
      linking_set='plants'
    } else if (str_detect( dataset_name, 'PH')) {
      interactions = list("pollination", "herbivory")
      linking_set='plants'
    } else if (str_detect(dataset_name, 'PSD')) {
      interactions = list("pollination","dispersion")
      linking_set='plants'
    }
  } else if (startsWith( dataset_name, 'Hackett')){
    if (str_detect(dataset_name,'_PH')){
      interactions = list("pollination", "herbivory")
      linking_set='plants'
    } else if (str_detect(dataset_name,'_HPa')) {
      interactions = list('herbivory', 'parasitism')
      linking_set='hosts'
    } else if (str_detect(dataset_name,'_SHPa')) {
      interactions = list('herbivory', 'parasitism')
      linking_set='hosts'
    } else if (str_detect(dataset_name,'_LHSH')) {
      interactions = list('leaf_herbivory', 'seed_herbivory')
      linking_set='plants'
    }
  } else if (startsWith( dataset_name, 'Pocock')){
    if (str_detect(dataset_name,'_PH')){
      interactions = list("pollination", "herbivory")
      linking_set='plants'
    } else if (str_detect(dataset_name,'_SHH')){
      interactions = list("seed_herbivory", "herbivory")
      linking_set='plants'
    }
  } else if (startsWith( dataset_name, 'Genrich')){
    interactions = list("dispersion","seed_herbivory")
    linking_set='plants'
  } else if (startsWith( dataset_name, 'Melian_DH')){
    interactions = list("dispersion", "herbivory")
    linking_set='plants'
  } else if (startsWith( dataset_name,  'Mello')){
    interactions = list("frugivory", "nectarivory")
    linking_set='plants'
  } else if (str_detect( dataset_name,  '_PSD')) {
    interactions = list("dispersion", "pollination")
    linking_set='plants'
  } else if (str_detect( dataset_name,  '_PA')){
    interactions = list("ant", "pollination")
    linking_set='plants'
  } else if (str_detect( dataset_name,  '_SDA')){
    interactions = list("ant", "dispersion")
    linking_set='plants'
  } else if (startsWith( dataset_name,  'McFayden')){
    interactions = list('herbivory', 'parasitism')
    linking_set='hosts'
  } else if (startsWith( dataset_name,  'Dattilo')){
    if (endsWith( dataset_name, 'PSD')) {
      interactions = list("pollination","dispersion")
      linking_set='plants'
    }
    else if (str_detect( dataset_name, 'PA')){
      interactions = list("pollination","ant")
      linking_set='plants'
    }
  } else {
    interactions = list("pollination", "herbivory")
    linking_set='plants'
  }
  return(list("interactions"=interactions,"linking_set"=linking_set))
}


species_type_letter=function(species_type){
  if (species_type=="plants"){
    species_letter="p"
  } else if(species_type=="hosts"){ #à revoir
    species_letter="p"
  } else if(species_type=="pollination"){
    species_letter="i"
  } else if(species_type=="herbivory"){
    species_letter="h"
  } else if(species_type=="dispersion"){
    species_letter="d"
  } else if(species_type=="ant"){
    species_letter="a"
  } else if(species_type=="parasitism"){
    species_letter="K"
  } else if(species_type=="leaf_herbivory"){
    species_letter="L"
  } else if(species_type=="seed_herbivory"){
    species_letter="S"
  } else if(species_type=="frugivory"){
    species_letter="F"
  } else if(species_type=="nectarivory"){
    species_letter="N"
  } else{
    species_letter="X"
  }
  #Améliorer la gestion des lettres (en autoriser deux)
  return(species_letter)
}

rename_species=function(species_names,species_type,rename_short=FALSE){
  species_letter=species_type_letter(species_type)

  for (i in 1:length(species_names)){
    if (rename_short==TRUE){
      species_names[i]<-paste0(species_letter,i)
    } else{
      species_names[i]<-sprintf(paste0(species_letter,"%s.",species_names[i]),i)}
  }
  #Refaire avec syntaxe de vecteur sans boucles : paste0("p",1:40)
  return(species_names)
}

format_sbm_networks_to_incidence_matrix=function(Networks){
  #Création d'un objet Incidence_matrices avec même format que d'habitude
  #à partir d'un objet networks SBM (chemin inverse)
  Incidence_matrices=list()
  interactions=list()
  for(network in Networks)
  {
    linking_set<-network$dimLabels[1]
    linking_species_letter=species_type_letter(linking_set)
    nb_linking_species<-network$nbNodes[1]
    
    interaction<-network$dimLabels[2]
    interaction_species_letter<-species_type_letter(interaction)
    nb_interaction_species<-network$nbNodes[2]
    
    Adj<-network$networkData
    rownames(Adj)=paste0(linking_species_letter,1:nb_linking_species)
    colnames(Adj)=paste0(interaction_species_letter,1:nb_interaction_species)
    Incidence_matrices=list.append(Incidence_matrices, Adj)
    interactions=list.append(interactions, interaction)
  }
  names(Incidence_matrices) <- interactions
  return(Incidence_matrices)
}

####################################SIMULATION OF GRAPHS

sample_tripartite_graph=function(nbNodes,blockProp,connectParam){
  ### SAMPLE TRIPARTITE SBM  : 2 networks between 3 Functional Groups
  ## Graph parameters
  #Functional Groups (FG)
  names(blockProp)<-c("plants","pollination","herbivory")
  #Interactions between the FG
  archiMultipartite  <-  rbind(c(1,2),c(1,3)) #
  model <- c('bernoulli','bernoulli') # type of distribution in each network
  # for each network : directed or not (not required for an interaction between two different FG)
  directed <- c(NA,NA)
  dimLabels <- c("plants","pollination","herbivory")
  ## Graph Sampling
  mySampleMBM <- sampleMultipartiteSBM(nbNodes, blockProp,
                                       archiMultipartite,
                                       connectParam, model, directed,
                                       dimLabels)
  Networks <- mySampleMBM$listSBM
  memberships <- mySampleMBM$memberships
  Incidence_matrices<-format_sbm_networks_to_incidence_matrix(Networks)
  return(list("Networks"=Networks,"Incidence_matrices"=Incidence_matrices,"memberships"=memberships))
}


####################################ANALYSE DE GRAPH

define_graph_from_indicence_matrice=function(Incidence_matrices){
  #Création du graph (objet igraph)
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
  return(g)
}

find_nodes_hesitation=function(myMSBM,Incidence_matrices,seuil_hesitation=0.9){
  partites<-myMSBM$dimLabels
  linking_set<-partites[1]
  interactions<-list(partites[2],partites[3])
  linking_names<-rownames(Incidence_matrices[[1]])
  interaction_1_names<-colnames(Incidence_matrices[[eval(interactions[[1]])]])
  interaction_2_names<-colnames(Incidence_matrices[[eval(interactions[[2]])]])
  
  all_nodes_probMemberships<-unlist(c(apply(myMSBM$probMemberships[[eval(interactions[[2]])]],1, list),
                                      apply(myMSBM$probMemberships[[eval(linking_set)]],1, list),
                                      apply(myMSBM$probMemberships[[eval(interactions[[1]])]],1, list)),
                                    recursive = FALSE)
  names(all_nodes_probMemberships)<-do.call(c,list(interaction_2_names,linking_names,interaction_1_names))
  nodes_hesitation<-(rapply(all_nodes_probMemberships,max)<seuil_hesitation)
  return(all_nodes_probMemberships[nodes_hesitation])
}


####################################FONCTIONS DE PLOT

plot_graph_tripartite_with_clusters=function(myMSBM,Incidence_matrices,
                                             newfigures=TRUE,plot_types=c(1,2),
                                             plot_graph_no_partites=TRUE,show_titles=TRUE,
                                             nodes_size=3,nodes_pie=TRUE,seuil_hesitation=0.9){
  #Possibilité de rajouter le calcul de la matrice d'incidence à partir de myMSMB dans cette fonction.
  partites<-myMSBM$dimLabels
  nb_parties=length(partites)
  linking_set<-partites[1]
  interactions<-list(partites[2],partites[3])
  #Les noms pourraient être récupérés depuis mySMBM directement mais visiblement disparaissent quand un seul cluster pour l'instant
  linking_names<-rownames(Incidence_matrices[[1]])
  interaction_1_names<-colnames(Incidence_matrices[[eval(interactions[[1]])]]) 
  interaction_2_names<-colnames(Incidence_matrices[[eval(interactions[[2]])]])

  #Création du graph
  # g<-define_graph_from_indicence_matrice(Incidence_matrices) #à remplacer plus tard, car need des matrices m_pi et m_hp en dessous
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
  
  ##Définition d'un layout multipartite personnalisé
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
    palette_int2<-brewer.pal(max(nb_c_int2+1,3), "BuPu")
  } else {palette_int2<-colorRampPalette(brewer.pal(9, "BuPu"))(nb_c_int2+1)}
  
  if (nb_c_link<=8){
    palette_link<-brewer.pal(max(nb_c_link+1,3), "PiYG")
  } else{palette_link<-colorRampPalette(brewer.pal(9, "PiYG"))(nb_c_link+1)}
  
  if (nb_c_int1<=8){
    palette_int1<-brewer.pal(max(nb_c_int1+1,3), "OrRd")
  } else{palette_int1<-colorRampPalette(brewer.pal(9, "OrRd"))(nb_c_int1+1)}
  
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
    #Affichage (layout personnalisé en x et y)
    if (newfigures){
      x11()
      if (plot_graph_no_partites){par(mfrow=c(1,2))}
      }
    if (plot_graph_no_partites){
      if (nodes_pie){plot(g, edge.width=E(g)$Freq, edge.arrow.size=0,vertex.size=6,vertex.label.cex=0.8,vertex.label.degree=0,
                          vertex.shape="pie",vertex.pie=probMemberships_pie, vertex.pie.color=couleurs_clusters_pie,
                          vertex.frame.color=couleurs_borders_hesitation,vertex.label.font=font_hesitation,vertex.label.color=couleurs_borders_hesitation)
      }
      else{plot(g, vertex.color=couleurs_clusters, edge.width=E(g)$Freq, edge.arrow.size=0,vertex.size=5,vertex.label.cex=0.8)}
      if(show_titles){title("Graph (no partites)")}
    }
    
    if (nodes_pie){
      plot(g, vertex.color=couleurs_clusters, layout=layout_multipartite, edge.width=E(g)$Freq, edge.arrow.size=0,vertex.size=nodes_size,vertex.label.dist=1,vertex.label.degree=pi,vertex.label.cex=0.8,
         vertex.frame.color=couleurs_borders_hesitation,vertex.label.color=couleurs_borders_hesitation)
    }
    else{plot(g, vertex.color=couleurs_clusters, layout=layout_multipartite, edge.width=E(g)$Freq, edge.arrow.size=0,vertex.size=nodes_size,vertex.label.dist=1,vertex.label.degree=pi,vertex.label.cex=0.8)}
    if(show_titles){title("Tripartite graph (1)")}
  }
  
  if (2 %in% plot_types){
    #Affichage (layout personnalisé en x, auto en y)
    if (newfigures){
      x11()
      if (plot_graph_no_partites){par(mfrow=c(1,2))}
      }
    if (plot_graph_no_partites){
      if (nodes_pie){plot(g, edge.width=E(g)$Freq, edge.arrow.size=0,vertex.size=6,vertex.label.cex=0.8,vertex.label.degree=0,
                          vertex.shape="pie",vertex.pie=probMemberships_pie, vertex.pie.color=couleurs_clusters_pie,
                          vertex.frame.color=couleurs_borders_hesitation,vertex.label.font=font_hesitation,vertex.label.color=couleurs_borders_hesitation)
      }
      else{plot(g, vertex.color=couleurs_clusters, edge.width=E(g)$Freq, edge.arrow.size=0,vertex.size=5,vertex.label.cex=0.8)}
      if(show_titles){title("Graph (no partites)")}
    }
    
    
    if (nodes_pie){
      plot(g, vertex.color=couleurs_clusters, layout=layout_multipartite2, edge.width=E(g)$Freq, edge.arrow.size=0,vertex.size=nodes_size,vertex.label.dist=1,vertex.label.degree=pi,vertex.label.cex=0.8,
         vertex.frame.color=couleurs_borders_hesitation,vertex.label.color=couleurs_borders_hesitation)
    }
    else{plot(g, vertex.color=couleurs_clusters, layout=layout_multipartite2, edge.width=E(g)$Freq, edge.arrow.size=0,vertex.size=nodes_size,vertex.label.dist=1,vertex.label.degree=pi,vertex.label.cex=0.8)}
    if(show_titles){title("Tripartite graph (2)")}

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


plot_graph_tripartite_without_clusters=function(myMSBM,Incidence_matrices,nodes_size=3){
  #Possibilité de rajouter le calcul de la matrice d'incidence à partir de myMSMB dans cette fonction.
  partites<-myMSBM$dimLabels
  nb_parties=length(partites)
  linking_set<-partites[1]
  interactions<-list(partites[2],partites[3])
  #Les noms pourraient être récupérés depuis mySMBM directement mais visiblement disparaissent quand un seul cluster pour l'instant
  linking_names<-rownames(Incidence_matrices[[1]])
  interaction_1_names<-colnames(Incidence_matrices[[eval(interactions[[1]])]]) 
  interaction_2_names<-colnames(Incidence_matrices[[eval(interactions[[2]])]])
  
  #Création du graph
  # g<-define_graph_from_indicence_matrice(Incidence_matrices) #à remplacer plus tard, car need des matrices m_pi et m_hp en dessous
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
  
  permutation<-Matrix::invPerm(order(V(g)$type))
  g=permute(g,permutation)
  
  
  positions_FG <- match(V(g)$type, c(species_type_letter(partites[3]), species_type_letter(partites[1]), species_type_letter(partites[2])) )
  positions_clusters<-do.call(c,list(which(!apply(m_hp == 0, 1, all)),which(!apply(cbind(m_pi,t(m_hp)) == 0, 1, all)),which(!apply(m_pi == 0, 2, all)))) #only vertices with at least 1 interaction
  names(positions_clusters)<-NULL
  # color_herbi<-brewer.pal(3, "OrRd")[3]
  # color_plants<-brewer.pal(3, "PiYG")[3]
  # color_polli<-brewer.pal(3, "BuPu")[2]
  color_plants<-"forestgreen"
  color_polli<-"firebrick1"
  color_herbi<-"deepskyblue3"
  couleurs_clusters=c(rep(color_herbi,sum(!apply(m_hp == 0, 1, all))),rep(color_polli,sum(!apply(m_pi == 0, 2, all))),rep(color_plants,sum(!apply(cbind(m_pi,t(m_hp)) == 0, 1, all))))
  layout_multipartite <- layout_with_fr(g, minx=positions_FG, maxx=positions_FG,miny=positions_clusters,maxy=positions_clusters)
  x11()
  plot(g, vertex.color=couleurs_clusters, layout=layout_multipartite, edge.width=E(g)$Freq, edge.arrow.size=0,vertex.size=nodes_size,vertex.label.dist=1,vertex.label.degree=pi,vertex.label.cex=0.8)
}