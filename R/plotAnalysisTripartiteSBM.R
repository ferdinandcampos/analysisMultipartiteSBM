source("toolsFunctions.R")
library(knitr)
library(rlist)
library(stringr)
library(rgl)
library(car)
library(carData) #apparemment il faut ça aussi

plotAnalysisTripartiteSBM=function(analyse,
                                   Indices=as.numeric(names(table(analyse$network_id)))){
  ############Histogrammes des nombres de blocks
  x11(width = 19,height = 12)
  par(mfrow=c(2,3),mar = c(4, 5, 1, 0.3))
  distribution_K_all<-analyse[analyse$network_id %in% Indices & analyse$FG=="plants" & analyse$cluster_id==1,c("nb_herbi_clust","nb_plants_clust","nb_polli_clust")]
  rownames(distribution_K_all)<-Indices
  FG_names<-c("espèces périphériques a","espèces connectrices l","espèces périphériques b")
  FG_colors<-c("#E0ECF4","#A1D76A","#E34A33")
  K_max<-max(distribution_K_all)
  for (i in 1:3){
    distribution_K<-table(factor(distribution_K_all[,i],levels=1:K_max))/length(Indices)
    bp<-barplot(distribution_K,col=FG_colors[i],border="white",
                xlab=variables_to_Latex(names(distribution_K_all)[i]),ylab="Fréquence (%)",main=FG_names[i],
                cex.lab = 1.7, cex.axis = 1.3,cex=1.3)
    text(bp, 0, round(distribution_K, 2), cex=1.3, pos=3,col="black")
  }
  ############Histogrammes des tailles des blocks
  for (FG in c("herbivory","plants","mutualism")){
    if(FG=="plants"){
      nb_clust_names<-"nb_plants_clust"
      pi_clust_names<-"cluster_pi_plants"
      test2<-c("#F1B6DA", "#B8E186", "#4DAC26","#DE77AE", "#FDE0EF", "#E6F5D0","#F7F7F7")
    }else if(FG=="herbivory"){
      nb_clust_names<-"nb_herbi_clust"
      pi_clust_names<-"cluster_pi_herbi"
      test2<-c("#9EBCDA", "#8856A7","#8C96C6")
    }else if(FG=="mutualism"){
      FG<-c("pollination","dispersion")
      nb_clust_names<-"nb_polli_clust"
      pi_clust_names<-"cluster_pi_polli"
      test2<-c("#FDBB84", "#E34A33","#FC8D59")
    }
    distribution_pi_all<-matrix(0,K_max,length(Indices))
    colnames(distribution_pi_all)<-Indices
    #test<-matrix("red",K_max,length(Indices))
    #test2<-brewer.pal(max(3+1,3), "PiYG")[2:(3+1)]
    #test2<-brewer.pal(4, "PiYG")
    # palette_link_all<-c() #fonctionne PAS
    for (index in Indices){
      distribution_pi<-analyse[analyse$network_id==index & analyse$FG%in%FG,c("cluster_pi")]
      nb_c_link<-analyse[analyse$network_id==index ,nb_clust_names][1]
      # palette_link<-brewer.pal(max(nb_c_link+1,3), "PiYG")
      # palette_link_all<-c(palette_link_all,palette_link[2:(nb_c_link+1)])
      #test[1:length(distribution_pi),index]<-palette_link[2:(nb_c_link+1)]
      distribution_pi_all[1:length(distribution_pi),c(as.character(index))]<-sort(analyse[analyse$network_id==index & analyse$FG%in%FG,c("cluster_pi")])
    }
    bp<-barplot(distribution_pi_all,beside=FALSE,col=test2,border="white",
                xlab="Indice du réseau",ylab=variables_to_Latex(pi_clust_names),
                axisnames=FALSE,
                #main=paste("Répartition des tailles de blocks",FG[1],"par réseau"),
                #main=variables_to_Latex("cluster_pi_polli"),
                #cex.main=2,
                cex.lab = 1.7, cex.axis = 1.3,cex=1.3)
    text(bp,0, Indices, cex=1.3, pos=3,col="black")
  }
}  
