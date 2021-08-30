


#################################### MSBM DATA PREPARATION, MANIPULATION
#' tools letter
#'
#' @param x à faire
#' @param y à faire
#' @return à faire
#' @examples
#' à faire
#'
#' @export
species_type_letter=function(species_type){
  if (species_type=="plants"){
    species_letter="p"
  } else if(species_type=="hosts"){
    species_letter="H"
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
  #Am?liorer la gestion des lettres (en autoriser deux)
  return(species_letter)
}
#' tools rename
#'
#' @param x à faire
#' @param y à faire
#' @return à faire
#' @examples
#' à faire
#'
#' @export
rename_species=function(species_names,species_type,rename_short=TRUE){
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
#' tools incidence matrix
#'
#' @param x à faire
#' @param y à faire
#' @return à faire
#' @examples
#' à faire
#'
#' @import rlist sbm
#' @export
formatSBMtoIncidenceMatrices=function(Networks){
  #Cr?ation d'un objet Incidence_matrices avec m?me format que d'habitude
  #? partir d'un objet networks SBM (chemin inverse)
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
#################################### PLOTTING
#' tools variables printing in plot
#'
#' @param x à faire
#' @param y à faire
#' @return à faire
#' @examples
#' à faire
#'
#' @export
variables_to_Latex=function(variable){
  if (variable=="c_inter1"){latex_string=latex2exp::TeX(r'($\bar{c}_b$)')
  }else if (variable=="c_inter2"){latex_string=latex2exp::TeX(r'($\bar{c}_a$)')
  }else if (variable=="c_tot"){latex_string=latex2exp::TeX(r'($\bar{c}_{tot}$)')
  }else if (variable=="c_tot_sym"){latex_string=latex2exp::TeX(r'($\bar{c}_{tot}^{sym}$)')
  }else if (variable=="PR_inter1"){latex_string=latex2exp::TeX(r'($PR_b$)')
  }else if (variable=="PR_inter1_sym"){latex_string=latex2exp::TeX(r'($PR_b^{sym}$)')
  }else if (variable=="cKnn_inter1"){latex_string=latex2exp::TeX(r'($\bar{c}_b^{N}$)')
  }else if (variable=="cKnn_inter2"){latex_string=latex2exp::TeX(r'($\bar{c}_a^{N}$)')
  }else if (variable=="cKnn_tot"){latex_string=latex2exp::TeX(r'($\bar{c}_{tot}^{N}$)')
  }else if (variable=="PC1"){latex_string=latex2exp::TeX(r'($PCA_1$)')
  }else if (variable=="PC2"){latex_string=latex2exp::TeX(r'($PCA_2$)')
  }else if (variable=="PC3"){latex_string=latex2exp::TeX(r'($PCA_3$)')
  }else if (variable=="MDS1"){latex_string=latex2exp::TeX(r'($MDS_1$)')
  }else if (variable=="MDS2"){latex_string=latex2exp::TeX(r'($MDS_2$)')
  }else if (variable=="nb_herbi_clust"){latex_string=latex2exp::TeX(r'($K_a$)')
  }else if (variable=="nb_plants_clust"){latex_string=latex2exp::TeX(r'($K_l$)')
  }else if (variable=="nb_polli_clust"){latex_string=latex2exp::TeX(r'($K_b$)')
  }else if (variable=="cluster_pi_herbi"){latex_string=latex2exp::TeX(r'($\pi_k^a$)')
  }else if (variable=="cluster_pi_plants"){latex_string=latex2exp::TeX(r'($\pi_{k'}^l$)')
  }else if (variable=="cluster_pi_polli"){latex_string=latex2exp::TeX(r'($\pi_k^b$)')
  }else{latex_string="X"}
  return(latex_string)
}
