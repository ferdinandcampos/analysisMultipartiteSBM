#' Compute structural metrics on a tripartite network model MBM
#'
#' @param myMSBM the model to analyze
#' @param normalise enable metrics' normalization if TRUE
#' @return a data.frame including for each cluster of each FG of the network the metrics values and the relative sizes of clusters \code{cluster_pi}.
#' @examples
#' getMetricsTripartiteSBM(myMSBM)
#'
#' @export
getMetricsTripartiteSBM=function(myMSBM,normalise=TRUE){
  metric_degrees<-data.frame(FG=character(),cluster_id=integer(),cluster_pi=numeric(),
                             c_inter1=numeric(),c_inter2=numeric(),c_tot=numeric(),c_tot_sym=numeric(),PR_inter1=numeric(),PR_inter1_sym=numeric(),
                             cKnn_inter1=numeric(),cKnn_inter2=numeric(),cKnn_tot=numeric())
  linking_set<-myMSBM$dimLabels[1]
  inter1<-myMSBM$dimLabels[2]
  inter2<-myMSBM$dimLabels[3]
  interactions<-list(inter1,inter2)
  inter_percents<-as.numeric(myMSBM$nbNodes[2:3]/sum(myMSBM$nbNodes[2:3]))
  connect_alphas<-myMSBM$connectParam
  names(connect_alphas)<-unlist(interactions)

  #Calcul de connectances et PR
  for (FG in c(linking_set,inter1,inter2)){
    nb_FG_clust<-myMSBM$nbBlocks[[eval(FG)]]
    if(FG==linking_set){
      for (cluster_id in (1:nb_FG_clust)){
        cluster_pi<-myMSBM$blockProp[[eval(FG)]][cluster_id]
        c_inter1<-weighted.mean(connect_alphas[[inter1]]$mean[cluster_id,],myMSBM$blockProp[[eval(inter1)]])
        c_inter2<-weighted.mean(connect_alphas[[inter2]]$mean[cluster_id,],myMSBM$blockProp[[eval(inter2)]])
        c_tot<-weighted.mean(c(c_inter1,c_inter2),inter_percents) #DISCUTABLE
        c_tot_sym<-mean(c(c_inter1,c_inter2))
        PR_inter1<-c_inter1*inter_percents[1]/sum(c(c_inter1,c_inter2)*inter_percents) #DISCUTABLE
        PR_inter1_sym<-c_inter1/(c_inter1+c_inter2)

        metric_degrees[nrow(metric_degrees)+1,]<-c(FG,cluster_id,cluster_pi,c_inter1,c_inter2,c_tot,c_tot_sym,PR_inter1,PR_inter1_sym,rep(NA,3))
      }
    }else{
      for (cluster_id in (1:nb_FG_clust)){
        cluster_pi<-myMSBM$blockProp[[eval(FG)]][cluster_id]
        c_inter<-weighted.mean(connect_alphas[[eval(FG)]]$mean[,cluster_id],myMSBM$blockProp[[eval(linking_set)]])
        c_tot<-c_inter
        c_tot_sym<-c_inter

        metric_degrees[nrow(metric_degrees)+1,]<-c(FG,cluster_id,cluster_pi,rep(NA,2),c_tot,c_tot_sym,rep(NA,2),rep(NA,3))
        metric_degrees[nrow(metric_degrees),2+which(myMSBM$dimLabels==FG)]<-c_inter
      }
    }

  }
  metric_degrees[,3:ncol(metric_degrees)]<-lapply(metric_degrees[,3:ncol(metric_degrees)], function(x) as.numeric(x)) #force numerics pour les calculs suivants... #pourquoi ils sont chang?s en character?

  #Calcul des Knn  - Attention d?pend ordre colonnes
  for (FG in c(linking_set,inter1,inter2)){
    nb_FG_clust<-myMSBM$nbBlocks[[eval(FG)]]
    if(FG==linking_set){
      for (cluster_id in (1:nb_FG_clust)){
        cKnn_inter1<-sum(metric_degrees[metric_degrees$FG==inter1,"c_inter1"]*connect_alphas[[inter1]]$mean[cluster_id,]*myMSBM$blockProp[[eval(inter1)]])
        cKnn_inter2<-sum(metric_degrees[metric_degrees$FG==inter2,"c_inter2"]*connect_alphas[[inter2]]$mean[cluster_id,]*myMSBM$blockProp[[eval(inter2)]])
        cKnn_inter1<-cKnn_inter1/metric_degrees[metric_degrees$FG==FG & metric_degrees$cluster_id==cluster_id,"c_inter1"]
        cKnn_inter2<-cKnn_inter2/metric_degrees[metric_degrees$FG==FG & metric_degrees$cluster_id==cluster_id,"c_inter2"]
        PR_inter1_cluster<-metric_degrees[metric_degrees$FG==FG & metric_degrees$cluster_id==cluster_id,"PR_inter1"]
        cKnn_tot<-weighted.mean(c(cKnn_inter1,cKnn_inter2),c(PR_inter1_cluster,1-PR_inter1_cluster)) #DISCUTABLE

        metric_degrees[metric_degrees$FG==FG & metric_degrees$cluster_id==cluster_id,c("cKnn_inter1","cKnn_inter2","cKnn_tot")]<-c(cKnn_inter1,cKnn_inter2,cKnn_tot)
      }
    }else{
      for (cluster_id in (1:nb_FG_clust)){

        cKnn_linking<-sum(metric_degrees[metric_degrees$FG==linking_set,2+which(myMSBM$dimLabels==FG)]*connect_alphas[[eval(FG)]]$mean[,cluster_id]*myMSBM$blockProp[[eval(linking_set)]])
        cKnn_linking<-cKnn_linking/metric_degrees[metric_degrees$FG==FG & metric_degrees$cluster_id==cluster_id,2+which(myMSBM$dimLabels==FG)]
        cKnn_tot<-cKnn_linking

        if(FG==inter1){ #peut ?tre + propre de faire 3 elsif que des which car ?a d?pend de l'ordre des colonnes...
          metric_degrees[metric_degrees$FG==FG & metric_degrees$cluster_id==cluster_id,c("cKnn_inter1","cKnn_tot")]<-c(cKnn_linking,cKnn_tot)
        }else{
            metric_degrees[metric_degrees$FG==FG & metric_degrees$cluster_id==cluster_id,c("cKnn_inter2","cKnn_tot")]<-c(cKnn_linking,cKnn_tot)}
      }
    }

  }

  if (normalise){ #NORMALISATION

    #Vecteurs des moyennes de c pour les 3 FG
    c_interactions_mean<-list(weighted.mean(connect_alphas[[inter1]]$mean,myMSBM$blockProp[[eval(linking_set)]]%*%t((myMSBM$blockProp[[eval(inter1)]])))
                              ,weighted.mean(connect_alphas[[inter2]]$mean,myMSBM$blockProp[[eval(linking_set)]]%*%t((myMSBM$blockProp[[eval(inter2)]]))))
    names(c_interactions_mean)<-c(inter1,inter2)
    #Vecteurs des moyennes et std de Knn pour les 3 FG
    #Quantit?s ql
    produits_pis_c_inter1<-myMSBM$blockProp[[eval(linking_set)]]*metric_degrees[metric_degrees$FG==linking_set,4] #refaire propre le "4" et "5"
    produits_pis_c_inter1<-produits_pis_c_inter1%*%t((myMSBM$blockProp[[eval(inter1)]]))
    produits_pis_c_inter2<-myMSBM$blockProp[[eval(linking_set)]]*metric_degrees[metric_degrees$FG==linking_set,5]
    produits_pis_c_inter2<-produits_pis_c_inter2%*%t((myMSBM$blockProp[[eval(inter2)]]))
    cKnn_interactions_mean<-list(sum(connect_alphas[[inter1]]$mean*produits_pis_c_inter1)/c_interactions_mean[[1]],sum(connect_alphas[[inter2]]$mean*produits_pis_c_inter2)/c_interactions_mean[[2]])
    names(cKnn_interactions_mean)<-c(inter1,inter2)
    cKnn_interactions_std<-list(sqrt(cKnn_interactions_mean[[1]]*(1-cKnn_interactions_mean[[1]])),
                                sqrt(cKnn_interactions_mean[[2]]*(1-cKnn_interactions_mean[[2]])))
    names(cKnn_interactions_std)<-c(inter1,inter2)
    #Quantit?s lq
    produits_pis_c_inter1<-myMSBM$blockProp[[eval(inter1)]]*metric_degrees[metric_degrees$FG==inter1,4]
    produits_pis_c_inter1<-produits_pis_c_inter1%*%t((myMSBM$blockProp[[eval(linking_set)]]))
    produits_pis_c_inter2<-myMSBM$blockProp[[eval(inter2)]]*metric_degrees[metric_degrees$FG==inter2,5]
    produits_pis_c_inter2<-produits_pis_c_inter2%*%t((myMSBM$blockProp[[eval(linking_set)]]))
    cKnn_linking_set_mean<-list(sum(connect_alphas[[inter1]]$mean*t(produits_pis_c_inter1))/c_interactions_mean[[1]],sum(connect_alphas[[inter2]]$mean*t(produits_pis_c_inter2))/c_interactions_mean[[2]])
    names(cKnn_linking_set_mean)<-c(inter1,inter2)
    cKnn_linking_set_std<-list(sqrt(cKnn_linking_set_mean[[1]]*(1-cKnn_linking_set_mean[[1]])),
                               sqrt(cKnn_linking_set_mean[[2]]*(1-cKnn_linking_set_mean[[2]])))
    names(cKnn_linking_set_std)<-c(inter1,inter2)
    #Cas limite o? un seul block dans un FG
    if (myMSBM$nbBlocks[[eval(inter1)]]==1){ #on force std=1 pour ?viter NaN quand (cas K_a = 1)
      cKnn_interactions_std[[1]]<-1
    }
    if (myMSBM$nbBlocks[[eval(inter2)]]==1){ #on force std=1 pour ?viter NaN quand (cas K_b = 1)
      cKnn_interactions_std[[2]]<-1
    }
    if (myMSBM$nbBlocks[[eval(linking_set)]]==1){ #on force std=1 pour ?viter NaN quand (cas K_l = 1)
      cKnn_linking_set_std[[1]]<-1
      cKnn_linking_set_std[[2]]<-1
    }

    #Calculs de normalisation
    for (FG in c(linking_set,inter1,inter2)){
      nb_FG_clust<-myMSBM$nbBlocks[[eval(FG)]]
      if(FG==linking_set){
        for (cluster_id in (1:nb_FG_clust)){
          c_inter1<-metric_degrees[metric_degrees$FG==FG & metric_degrees$cluster_id==cluster_id,"c_inter1"]
          c_inter2<-metric_degrees[metric_degrees$FG==FG & metric_degrees$cluster_id==cluster_id,"c_inter2"]
          cKnn_inter1<-metric_degrees[metric_degrees$FG==FG & metric_degrees$cluster_id==cluster_id,"cKnn_inter1"]
          cKnn_inter2<-metric_degrees[metric_degrees$FG==FG & metric_degrees$cluster_id==cluster_id,"cKnn_inter2"]

          c_inter1<-(c_inter1-c_interactions_mean[[1]])/sqrt(c_interactions_mean[[1]]*(1-c_interactions_mean[[1]]))
          c_inter2<-(c_inter2-c_interactions_mean[[2]])/sqrt(c_interactions_mean[[2]]*(1-c_interactions_mean[[2]]))
          c_tot<-weighted.mean(c(c_inter1,c_inter2),inter_percents)
          c_tot_sym<-mean(c(c_inter1,c_inter2))
          cKnn_inter1<- (cKnn_inter1-cKnn_linking_set_mean[[1]])/cKnn_linking_set_std[[1]]
          cKnn_inter2<- (cKnn_inter2-cKnn_linking_set_mean[[2]])/cKnn_linking_set_std[[2]]
          cKnn_tot<-weighted.mean(c(cKnn_inter1,cKnn_inter2),c(PR_inter1_cluster,1-PR_inter1_cluster))

          metric_degrees[metric_degrees$FG==FG & metric_degrees$cluster_id==cluster_id,c("c_inter1","c_inter2","c_tot","c_tot_sym","cKnn_inter1","cKnn_inter2","cKnn_tot")]<-c(c_inter1,c_inter2,c_tot,c_tot_sym,cKnn_inter1,cKnn_inter2,cKnn_tot)
        }
      }else{
        for (cluster_id in (1:nb_FG_clust)){
          if(FG==inter1){ #peut ?tre + propre de faire 3 elsif que des which car ?a d?pend de l'ordre des colonnes...
            c_inter<-metric_degrees[metric_degrees$FG==FG & metric_degrees$cluster_id==cluster_id,"c_inter1"]
            cKnn_linking<-metric_degrees[metric_degrees$FG==FG & metric_degrees$cluster_id==cluster_id,"cKnn_inter1"]
          }else{
            c_inter<-metric_degrees[metric_degrees$FG==FG & metric_degrees$cluster_id==cluster_id,"c_inter2"]
            cKnn_linking<-metric_degrees[metric_degrees$FG==FG & metric_degrees$cluster_id==cluster_id,"cKnn_inter2"]
          }
          c_inter<-(c_inter-c_interactions_mean[[eval(FG)]])/sqrt(c_interactions_mean[[eval(FG)]]*(1-c_interactions_mean[[eval(FG)]]))
          c_tot<-c_inter
          c_tot_sym<-c_inter
          cKnn_linking<-(cKnn_linking-cKnn_interactions_mean[[eval(FG)]])/cKnn_interactions_std[[eval(FG)]]
          cKnn_tot<-cKnn_linking

          if(FG==inter1){
            metric_degrees[metric_degrees$FG==FG & metric_degrees$cluster_id==cluster_id,c("c_inter1","c_tot","c_tot_sym","cKnn_inter1","cKnn_tot")]<-c(c_inter,c_tot,c_tot_sym,cKnn_linking,cKnn_tot)
          }else{
            metric_degrees[metric_degrees$FG==FG & metric_degrees$cluster_id==cluster_id,c("c_inter2","c_tot","c_tot_sym","cKnn_inter2","cKnn_tot")]<-c(c_inter,c_tot,c_tot_sym,cKnn_linking,cKnn_tot)
          }
        }
      }
  }
  }
  #Types et arrondis
  metric_degrees[,3:ncol(metric_degrees)]<-lapply(metric_degrees[,3:ncol(metric_degrees)], function(x) as.numeric(x)) #pourquoi ils sont en character?
  #is.num <- sapply(metric_degrees, is.numeric)
  metric_degrees[,3:ncol(metric_degrees)] <- lapply(metric_degrees[,3:ncol(metric_degrees)], round, 4)
  return(metric_degrees)
}


