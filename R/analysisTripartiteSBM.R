#' Compute structural analysis on a collection of tripartite networks.
#' First it calibrates MBM on each network then computes structural metrics.
#'
#' @param x à faire
#' @param y à faire
#' @return à faire
#' @examples
#' à faire
#' @import sbm rlist stringr
#'
#' @export

analysisTripartiteSBM=function(dataset_path,dataset_names,Indices=as.numeric(names(table(analyse$network_id))),
                               weight=FALSE,must_rename_species=TRUE,rename_short=TRUE,
                               force_run_MSBM=FALSE,must_save=FALSE, quiet=FALSE,
                               guess_K_clusters_type="free",range_clusters_constrained=list(c(2,10),c(2,10),c(2,10)),
                               must_normalise_metrics=TRUE){


  analyse_reseaux_reels<-data.frame(network_id=integer(),network_name=character(),nb_nodes=integer(),nb_nodes_plants=integer(),nb_nodes_polli=integer(),nb_nodes_herbi=integer(),nb_plants_clust=integer(),nb_polli_clust=integer(),nb_herbi_clust=integer(),
                                    FG=character(),cluster_id=integer(),cluster_pi=numeric(),
                                    c_inter1=numeric(),c_inter2=numeric(),c_tot=numeric(),c_tot_sym=numeric(),PR_inter1=numeric(),PR_inter1_sym=numeric(),
                                    cKnn_inter1=numeric(),cKnn_inter2=numeric(),cKnn_tot=numeric())
  for (index in Indices){
    dataset_name=dataset_names[index]

    ####PREPARATION DES DONNEES
    if (weight){ #AVEC WEIGHT
      dataset_filename=file.path(dataset_path,sprintf("Mnet_%s_w.csv",dataset_name))
    } else{ #SANS WEIGHT
      dataset_filename=file.path(dataset_path,sprintf("Mnet_%s.csv",dataset_name))
      }

    #name of interactions
    partites=find_partites_type(dataset_name)
    interactions<-partites$interactions
    linking_set<-partites$linking_set
    n_interactions=length(interactions)
    if(guess_K_clusters_type=="constrained"){
      names(range_clusters_constrained)<-c(linking_set,interactions[[1]],interactions[[2]])
    }

    #store linking species names (usually plants)
    net1 = read.csv(file = dataset_filename,header = TRUE,skip = 1,check.names = FALSE)
    if (must_rename_species==TRUE){
      linking_names=rename_species(net1[,2],linking_set,rename_short)
    } else{linking_names = net1[,2]}

    #store other species names (usually animals)
    net2 = read.csv(file = dataset_filename,header = FALSE)
    type_animals = net2[1,]
    type_animals = as.factor(as.character(type_animals[-c(1,2)]))

    Incidence_matrices=list()
    for(interaction in interactions)
    {
      ind_col=which(net2[1,] == interaction)
      Adj=net1[,ind_col]
      rownames(Adj)=linking_names
      Adj=as.matrix(Adj)
      if (must_rename_species==TRUE){
        colnames(Adj)=rename_species(colnames(Adj),interaction,rename_short)
      }
      Incidence_matrices=list.append(Incidence_matrices, Adj)
    }
    names(Incidence_matrices) <- interactions

    Networks<-list()
    for (interaction in interactions)
    {
      Net=defineSBM(Incidence_matrices[[eval(interaction)]],'bernoulli','bipartite',FALSE,dimLabels=c(linking_set,interaction))
      Networks=list.append(Networks,Net)
    }


    ####Nom du fichier mod?le
    if (weight)
    {
      outputfile=sprintf("./MSBM_models_output/resMSBM_%s_%s_w.Rdata",dataset_name,guess_K_clusters_type)
      outputfile_free=sprintf("./MSBM_models_output/resMSBM_%s_free_w.Rdata",dataset_name)
    }  else {
      outputfile=sprintf("./MSBM_models_output/resMSBM_%s_%s.Rdata",dataset_name,guess_K_clusters_type)
      outputfile_free=sprintf("./MSBM_models_output/resMSBM_%s_free.Rdata",dataset_name)
    }

    ####CHARGEMENT OU RUN DE MSBM
    if (!file.exists(outputfile) | force_run_MSBM) #ou a) RUN DE MSBM
    {
      if (guess_K_clusters_type=="free"){
        estimOptions = list(initBM = TRUE)
        myMSBM <- estimateMultipartiteSBM(Networks,list.append(estimOptions,"verbosity"=0))
      } else if (guess_K_clusters_type=="constrained"){
        estimOptions = list(initBM = TRUE,nbBlocksRange=range_clusters_constrained)
        myMSBM <- estimateMultipartiteSBM(Networks,list.append(estimOptions,"verbosity"=0))
      } else if (guess_K_clusters_type=="forcedsplit" | guess_K_clusters_type=="forcedsplitsoft"){

        if (file.exists(outputfile_free)){
          load(outputfile_free) #v?rifier que cela charge bien l'objet sous le nom "myMSBM" ?
          nb_blocks<-myMSBM$nbBlocks
          if(sum(nb_blocks==1)){
            nb_blocks_forcedsplit<-as.numeric(nb_blocks+(nb_blocks==1))
            if (guess_K_clusters_type=="forcedsplit"){
              range_clusters_forcedsplit<-list(rep(nb_blocks_forcedsplit[1],2),rep(nb_blocks_forcedsplit[2],2),rep(nb_blocks_forcedsplit[3],2))
            } else if (guess_K_clusters_type=="forcedsplitsoft"){
              range_clusters_forcedsplit<-list(c(nb_blocks_forcedsplit[1],10),c(nb_blocks_forcedsplit[2],10),c(nb_blocks_forcedsplit[3],10))}

            names(range_clusters_forcedsplit)<-c(linking_set,interactions[[1]],interactions[[2]])
            estimOptions = list(initBM = TRUE,nbBlocksRange=range_clusters_forcedsplit,maxiterVEM=1000) #Si initBM=TRUE et contrainte large --> ERREUR.. ICL=-inf et pas de cv
            myMSBM <- estimateMultipartiteSBM(Networks,list.append(estimOptions,"verbosity"=1))
            newnb_blocks<-myMSBM$nbBlocks
            if(!quiet){print(paste(guess_K_clusters_type,"| Split needed for index=",index,
                        "with nb of blocks from =", nb_blocks[1],nb_blocks[2],nb_blocks[3],
                        "to =",newnb_blocks[1],newnb_blocks[2],newnb_blocks[3]))}
          } else{if(!quiet){print(paste("Split not needed for index=",index,"with nb of blocks =", nb_blocks[1],nb_blocks[2],nb_blocks[3]))}}
        } else{if(!quiet){print(paste("! Must run free MSBM first on index=",index))}}
      }

      if(must_save){save(myMSBM,file=outputfile)
        if(!quiet){print(paste("MSBM computed and saved index =",index))}
      } else{if(!quiet){print(paste("MSBM computed index =",index))}}

    } else{
      load(outputfile)
      if(!quiet){print(paste("MSBM loaded index =",index))}
      } #ou b) CHARGEMENT DU MODELE

    ####ANALYSE
    metric_degrees<-getMetricsTripartiteSBM(myMSBM,normalise = must_normalise_metrics)
    for (FG in c(linking_set,interactions[[1]],interactions[[2]])){
      nb_FG_clust<-myMSBM$nbBlocks[[eval(FG)]]
      for (cluster_id in (1:nb_FG_clust)){
        analyse_reseaux_reels[nrow(analyse_reseaux_reels)+1,]<-c(index,dataset_name, sum(myMSBM$nbNodes), as.numeric(myMSBM$nbNodes),as.numeric(myMSBM$nbBlocks),
                                                                 metric_degrees[metric_degrees$FG==FG & metric_degrees$cluster_id==cluster_id,])
      }
    }


  if(!quiet){print(paste("Done analysis index =",index))}
  }
  if(must_save){
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "analyse_r?seaux")
    openxlsx::writeData(wb, "analyse_r?seaux", analyse_reseaux_reels, rowNames=FALSE)
    openxlsx::saveWorkbook(wb, file = "analyse_r?seaux_r?els.xlsx", overwrite = TRUE)
    if(!quiet){print("Analysis saved in file .xlsx")}
  }
  return(analyse_reseaux_reels)
}




## ---------------------------------------------------------
####
# weight
# guess_K_clusters_type #free (1:10 per default), constrained ou forcedsplit ou forcedsplitsoft
# #range_clusters_constrained #if constrained, list(c(K_min,K_max)) for each FG
# must_rename_species
# rename_short
# must_normalise_metrics
# force_run_MSBM #force run even if a file already exists
# must_save
####
