---
title: "Structural analysis of ecological tripartite networks"
subtitle: "An application of methods to a new collection of ecological tripartite networks"
author: "Ferdinand Campos"
output: rmarkdown::html_vignette
bibliography: biblio.bib
link-citations: yes
vignette: >
  %\VignetteIndexEntry{Application_EcologicalNetworksDataset}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



## Preliminaries

### Requirements
The packages required for the analysis are **analysisMultipartiteSBM** and some others for data manipulation or representation:

```r
library(analysisMultipartiteSBM)
```

### A new collection of networks


Notre étude est motivée par une nouvelle collection de 42 réseaux tripartites récemment rassemblée et préparée par Virginia Domínguez-García dans @dominguez2021notes à partir des principales études écologiques multipartites existant dans la
littérature. Les réseaux collectés sont issus de zones géographiques - et donc de type d’espèces, climat, nature
des sols - différentes compte tenu de la grande variété des pays d’origine des études sources (@melian2009diversity, @shinohara2019contrasting, @pocock2012robustness, @hackett2019reshaping, @macfadyen2009differences et @dattilo2016unravelling).
Il s’agit néanmoins du premier jeu de réseaux écologiques multipartites réels unifiant différents travaux
écologiques et le seul de cette ampleur.
Les données comprennent 44 réseaux tripartites structurés comme ci-dessus, dont les deux interactions
bipartites sont de trois natures différentes : mutualiste–antagoniste (MA), antagoniste–antagoniste (AA)
et mutualiste–mutualiste (MM).


Les données utilisées sont les 17 réseaux de type mutualiste-antagoniste (MA) de ce jeu de données tripartites.
Les interactions de ces 17 réseaux sont de deux types, herbivorie–pollinisation (i.e. a sont des
herbivores, l des plantes et b des pollinisateurs) pour 16 d'entre eux et herbivorie–dispersion de graines (i.e. a sont des herbivores, l des plantes et b des espèces qui disséminent des graines) pour un d'entre eux. Par souci de clarté, on utilisera l’appellation commune d’espèces mutualistes pour l’ensemble des 17 réseaux.




```r
data(ecologicalNetworksDataset)
```

### Preparation of the data


```r
dataset_path=""
dataset_names_MA <- c('Sinohara_1_ALL_PH', 'Sinohara_2_ALL_PH', 'Sinohara_3_ALL_PH', 'Sinohara_4_ALL_PH', 'Sinohara_ALL_A_PH', 
                      'Sinohara_ALL_E_PH', 'Sinohara_ALL_I_PH', 'Sinohara_2_E_PH', 'Sinohara_3_E_PH', 'Sinohara_4_I_PH', 
                      'Melian_OO_OO_PH', 'Hackett_1_ALL_PH', 'Hackett_2_ALL_PH', 'Hackett_1_S_PH', 'Hackett_1_GL_PH',
                      'Pocock_OO_OO_PH','Melian_OO_OO_HSD')


dataset_names<-dataset_names_MA
Indices=1:17
study_insectes_only=FALSE
if (study_insectes_only){dataset_names[16]<-"PocockInsects_OO_OO_PH"
                         Indices<-Indices[-c(11,17)]}
```

### Multipartite Block Model (MBM) 
Le modèle MBM introduit dans @Barhen2020 étend pour des réseaux multipartites l'approche des modèles à blocs. C'est un modèle de mélange probabiliste qui généralise conjointement le célèbre \textit{Stochastic Block Model (SBM)} sur réseau simple et le \textit{Latent Block Model (LBM)}, son extension aux réseaux bipartites. Il identifie des blocs d'espèces à l'intérieur des groupes d'espèces partageant des propriétés de connectivité similaires au sein de la collection de réseaux dans laquelle elles sont impliquées.

Pour tout couple de groupe d'espèces (GE) en interaction $(q,q^\prime)\in\mathcal E=\left\{(a,l),(b,l)\right\}$ on définit la matrice $X^{qq^\prime}$ de sorte que pour tout couple d'espèces $(i,i^\prime)\in 1,\dots,n_q \times  1,\dots, n_{q^\prime}$ des GE $q$ et $q^\prime$ respectivement :


$$
X^{qq^\prime}_{ii^\prime} = \left\{
    \begin{array}{ll}
        1 & \mbox{si $i$ et $i^\prime$ sont en interaction,} \\ 
        0 & \mbox{sinon.}
    \end{array}
\right.
$$


On suppose que chaque GE $q\in\{a,l,b\}$ est divisé en $K_q\in\mathbb N^*$. Pour toute espèce $i\in 1,\dots, n_q$ de $q$, soit $Z_i^q$ la variable latente de sorte que $Z_i^q=k$ si l'espèce $i$ du GE $q$ appartient au bloc $k\in 1,\dots,K_q$. On suppose alors que :

Les $Z_i^q$ sont indépendantes et $\forall q\in\{a,l,b\}, \forall k \in 1,\dots,K_q,
    \forall i \in 1,\dots,n_q,$ 
$$
    \mathbb P(Z_{i}^q=k)=\pi_k^q \quad\quad\quad \textrm{avec} \quad\quad\quad
    \pi_k^q\in[0,1], \ \sum_{k=1}^{K_q}\pi_k^q=1,
$$
et puis $\forall (q,q^\prime)\in\mathcal E, 
    \forall (k,k^\prime)\in\mathcal A^{qq^\prime},
    \forall (i,i^\prime)\in\mathcal S^{qq^\prime},$
    
$$
    X_{ii^\prime}^{qq^\prime}| \{Z_{i}^q=k, Z_{i^\prime}^{q^\prime}=k^\prime\} \sim_{ind}  \mathcal Bern(\alpha_{kk^\prime}^{qq^\prime}) \quad\quad\quad \textrm{avec} \quad\quad\quad \alpha_{kk^\prime}^{qq^\prime}\in [0,1].
$$

La quantité $\theta_K=(\alpha,\pi)$, où $\alpha=\left(\alpha_{kk^\prime}^{qq^\prime}\right)_{(k,k^\prime)\in\mathcal A^{qq^\prime}, (q,q^\prime)\in\mathcal E}$ sont les probabilités de connexion entre blocs et $\pi=\left(\pi_k^q\right)_{k\in 1,\dots,K_q, q\in\{a,l,b\}}$ les tailles relatives des blocs, donne les paramètres inconnus du modèle à estimer à $K=(K_a,K_l,K_b)$ fixé. Le vecteur des nombres de blocs $K$ est aussi à déterminer.

### Calibrate MBM and compute metrics
Pour chaque réseau tripartite individuel, on calibre MBM et puis on s'intéresse aux blocs inférés dans les trois groupes d'espèces (GE).
On fait le choix des 5 métriques suivantes pour classifier les blocs du groupe des connecteurs $l$, définies pour chaque bloc $k^\prime\in1,\dots,K_l$ d'un groupe de connecteurs $l$ donné : $$\left(\bar{c}_a, \bar{c}_b, \textrm{PR}_b^{sym}, \bar{c}_a^{\textrm{N}}, \bar{c}_b^{\textrm{N}}\right).$$


```r
#analyse<-analysisTripartiteSBM(dataset_path,dataset_names,Indices,quiet=TRUE)

if (!study_insectes_only){
  data(analyseEcologicalNetworksDataset)
} else{
  data(analyseEcologicalNetworksDataset_insects)
}
#> Warning in data(analyseEcologicalNetworksDataset): data set 'analyseEcologicalNetworksDataset' not found
```

## Caracterisation of the ecological structure

### Blocks analysis
Les données issues des réseaux représentent $N=53$ blocs $k^\prime$ de plantes à classifier, en interactions avec $21$ blocs d'herbivores et $28$ blocs de mutualistes. 


```r
plotAnalysisTripartiteSBM(analyse)
```

### Preparation of caracterisation

```r
FG="plants"
if (!study_insectes_only){#For our analysis 1
                          analysePlants<-analyse[analyse$FG==FG,]
                          metrics<-c("c_inter1","c_inter2","PR_inter1_sym","cKnn_inter1","cKnn_inter2")
                          order_clusters<-c(1,3,2,5,4,6)
                          nb_clusters<-6
                          plotIndicesblocs<-1:nrow(analysePlants)
}else{#For our analysis 2 (insects only)
      analysePlants<-analyse[analyse$FG==FG,]
      metrics<-c("c_inter1","c_inter2","PR_inter1_sym","cKnn_inter1")
      order_clusters<-c(2,1,3,5,4)
      nb_clusters<-5
      plotIndicesblocs<-1:nrow(analysePlants)
      plotIndicesblocs<-plotIndicesblocs[-c(27,28,29,30)] #on enl?ve les 4 indices propres au r?seau 11 removed
      last_indice_blocs<-plotIndicesblocs[length(plotIndicesblocs)]
      plotIndicesblocs<-c(plotIndicesblocs,last_indice_blocs+1,last_indice_blocs+2,last_indice_blocs+3,last_indice_blocs+4) #on ajoute 4 derniers num?ros
}
```


### Clustering analysis to caracterise blocks of plants


```r
caracterisation<-caracTripartiteSBM(analysePlants,metrics,order_clusters=order_clusters)

caracterisation$centers
```


```r
plotCaracTripartiteSBM(caracterisation,plotIndicesblocs)
```

## Additionnal analyses

### Second application, insects only



```r
study_insectes_only=TRUE
```

### Appendix 1 : functional cartography
Une façon d'illustrer les connectances partielles et taux de participation réside dans l'outil de la cartographie fonctionnelle ($\textit{functional cartography}$) introduite en écologie dans le célèbre @Guimera2005A, qui trace et classe des espèces dans un plan du type $(\textit{taux de participation}, \textit{degré})$.

```r
plotFunctionalCartography(analyse)
```

### Appendix 2 :degree-degree of neighbors inversion
Dans les écosystèmes plantes--pollinisateurs, une inversion en moyenne du caractère spécialiste/généraliste entre les espèces et leurs voisins est reportée dans la littérature (@bascompte2003nested, @jordano1987patterns). Ici on propose de retrouver ce signal à l'échelle des blocs en traçant pour les plantes  et les pollinisateurs les courbes de connectances partielles $\bar c_b$ en fonction des connectances partielles moyennes des voisins $\bar c_b^{\textrm{N}}$ associées aux blocs des 16 réseaux de type herbivorie--pollinisation à l'étude.

```r
Indices_plantspollination<-Indices[-c(5,8,11,13,14,15,17)]
plotDegInversion(analyse,Indices_plantspollination)
```
---
