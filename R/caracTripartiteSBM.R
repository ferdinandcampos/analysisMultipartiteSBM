#' Structural caracterisation of a collection of tripartite networks.
#'
#' @param x à faire
#' @param y à faire
#' @return à faire
#' @examples
#' à faire
#'
#' @export

caracTripartiteSBM=function(analysePlants,
                            metrics=c("c_inter1","c_inter2","PR_inter1_sym","cKnn_inter1","cKnn_inter2"),
                            nb_clusters=NULL,order_clusters=NULL,weighted_by_pi=TRUE,init_centers=NULL){


  ###weighted k-means
  df <- analysePlants[,metrics]
  df <- na.omit(df)
  df_nodoubles<-df
  ###WEIGHT A FAIRE VIS A VIS DE TAILLE CLUSTER
  if (weighted_by_pi){
    #Discr?tisation/pr?cision au 1%, prendre int(100*pi) fois l'observation a une taille pi
    df<-df[rep(seq_len(nrow(df)),round(100*analysePlants[,"cluster_pi"],0)),]
  }

  if (is.null(nb_clusters) & is.null(order_clusters)){ ##Nombre de clusters ? d?terminer
    x11()
    elbowrule<-fviz_nbclust(df, kmeans, method = "wss",nstart = 50)
    print(elbowrule)
    print ("Enter the number of clusters")
    nb_clusters <-scan(n=1,quiet=TRUE)
  } else {if(is.null(nb_clusters)){nb_clusters<-length(order_clusters)}}

  #df_nodoubles<-df[rownames(analysePlants),]

  #SCALING
  df <- scale(df) #scale each variable to have a mean of 0 and sd of 1
  df_centering<-attributes(df)$`scaled:center`
  df_scaling<-attributes(df)$`scaled:scale`

  if (is.null(init_centers)){init_centers<-nb_clusters}

  km <- kmeans(df, centers = init_centers, nstart = 100)
  km_clusters<-km$cluster[rownames(analysePlants)]
  km_centers<-km$centers

  #Sorting clusters to set non-random colors
  clusters_order<-order(table(km_clusters),decreasing = TRUE)#Unique sorting : by size
  if(!is.null(order_clusters)){
    clusters_order<-clusters_order[order_clusters]
  }
  km_centers<-km_centers[clusters_order,]
  rownames(km_centers)<-1:nb_clusters
  km_clusters_copy<-km_clusters
  for (k in 1:nb_clusters){
    km_clusters[km_clusters_copy==clusters_order[k]]<-k
  }


  ##MEMBERSHIPS NETWORKS DISTRIBUTION
  km_memberships=data.frame(rep(NA,length(km_clusters)))
  for (km_cluster_id in 1:nb_clusters){
    km_memberships[,km_cluster_id]=1*(km_clusters==km_cluster_id)
  }
  names(km_memberships)<-paste0("cluster",1:nb_clusters)
  rownames(km_memberships)<-names(km_clusters)
  nb_blocks_per_kmclusters<-aggregate(km_memberships, by=list(network_id=analysePlants$network_id), FUN=sum)
  rownames(nb_blocks_per_kmclusters)<-as.numeric(names(table(analysePlants$network_id)))

  ##PCA
  res.pca<-prcomp(df_nodoubles,scale=TRUE)

  ##MDS
  d <- dist(scale(df_nodoubles))
  fit.MDS <- cmdscale(d,eig=TRUE, k=2)

  return(list(centers = km_centers, memberships = km_clusters, distrib.networks = nb_blocks_per_kmclusters,
              data=list(analysePlants=analysePlants,km_centering=df_centering,km_scaling=df_scaling),
              result.kmeans = km ,result.PCA = res.pca, result.MDS=fit.MDS))
}

# library(ggplot2)
# library(stats)
# library(factoextra)
