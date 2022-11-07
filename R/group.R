#' @import dplyr
#' @import ape
#' @import phangorn
#' @import seqinr
#' @import ggplot2
#' @import parallel
#' @import foreach
#' @title Group of cells more informative
#' @description Group of cells more informative in a dataset
#' @details Computes from a table cell counts of RNASeq the group of cells more informative
#' @param WantedData Matrix data. \strong{NECESSARY}
#' @param n Number of subsampling. Optional. By default = 100.
#' @param Filter_value Cut-off value. Optional. By default = 1.
#' @param num_cells Cells to select. Optional. By default = 20 (maximum value allowed. Minimal value allowed is 3).
#' @param num_set_seed Set seed value. Optional. By default = 44.
#' @param num_col Column number in which the data matrix begins. Optional. By default = 1.
#' @examples no examples yet
#' @author Jose Jordan, JuanVi
#' @export
group <- function(WantedData,n,Filter_value,num_cells,num_set_seed,num_col){

  message('PLEASE, BE PATIENT WITH n > 1000 o.O
It is highly recommended to use an n < 50,000
/////////////////////////////////////////////')

  if (missing(n)) {
    n=100
    message('
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!++++++++++++++++++++++++++++++++++++++++++!!
!!+  ALERT!: n (subsampl.) = 100 (DEFAULT) +!!
!!++++++++++++++++++++++++++++++++++++++++++!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ')

  }

  if (missing(Filter_value)) {
    Filter_value=1
    message('
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!++++++++++++++++++++++++++++++++++++++++++!!
!!+   ALERT!: Filter_value = 1 (DEFAULT)   +!!
!!++++++++++++++++++++++++++++++++++++++++++!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ')
  }


  if (missing(num_cells)) {
    num_cells=20
    message('
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!++++++++++++++++++++++++++++++++++++++++++!!
!!+    ALERT!: num_cells = 20 (DEFAULT)    +!!
!!++++++++++++++++++++++++++++++++++++++++++!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ')

  }

  if (missing(num_set_seed)) {
    num_set_seed=44
    message('
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!++++++++++++++++++++++++++++++++++++++++++!!
!!+   ALERT!: num_set_seed = 44 (DEFAULT)  +!!
!!++++++++++++++++++++++++++++++++++++++++++!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ')

  }

  if (missing(num_col)) {
    num_col=1
    message('
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!++++++++++++++++++++++++++++++++++++++++++!!
!!+     ALERT!: num_col = 1 (DEFAULT)      +!!
!!++++++++++++++++++++++++++++++++++++++++++!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ')

  }
  x=num_col # columna en la que empieza el contaje de tus datos
  ncol(WantedData)->num_cols # numero total de columnas
  WantedData[,num_col:num_cols] -> WantedData2
  #Binarizamos la matriz de TPMs:

  WantedData2[WantedData2 < Filter_value] <- 0
  WantedData2[WantedData2 >= Filter_value] <- 1

  dist(t(WantedData2),method="euclidean",diag = TRUE, upper = TRUE) -> distancia_usando_logaritmo_base2_plus1
  as.data.frame(as.matrix(distancia_usando_logaritmo_base2_plus1)) -> jj

  which(jj==max(jj))->jj2
  #Buscamos posicion max.value:
  ceiling(jj2[1]/ncol(jj))->aaj
  ceiling(jj2[2]/ncol(jj))->bbj
  # View(jj[-c(aaj,bbj)])
  # View(jj[c(aaj,bbj),-c(aaj,bbj)])
  apply(jj[c(aaj,bbj),-c(aaj,bbj)], MARGIN=2, FUN=min)->jj3
  ##
  max(jj3)
  which(jj3==max(jj3))->hjkk
  hjkk
  grep(names(hjkk),colnames(jj))->hjkk
  hjkk
  ##
  apply(jj[c(aaj,bbj,hjkk),-c(aaj,bbj,hjkk)], MARGIN=2, FUN=min)->jj4
  max(jj4)
  which(jj4==max(jj4))->hjkk2
  hjkk2
  grep(names(hjkk2),colnames(jj))->hjkk2
  hjkk2
  ##
  apply(jj[c(aaj,bbj,hjkk,hjkk2),-c(aaj,bbj,hjkk,hjkk2)], MARGIN=2, FUN=min)->jj5
  max(jj5)
  which(jj5==max(jj5))->hjkk3
  hjkk3
  grep(names(hjkk3),colnames(jj))->hjkk3
  hjkk3
  ##
  apply(jj[c(aaj,bbj,hjkk,hjkk2,hjkk3),-c(aaj,bbj,hjkk,hjkk2,hjkk3)], MARGIN=2, FUN=min)->jj6
  max(jj6)
  which(jj6==max(jj6))->hjkk4
  hjkk4
  grep(names(hjkk4),colnames(jj))->hjkk4
  hjkk4
  ##
  apply(jj[c(aaj,bbj,hjkk,hjkk2,hjkk3,hjkk4),-c(aaj,bbj,hjkk,hjkk2,hjkk3,hjkk4)], MARGIN=2, FUN=min)->jj7
  max(jj7)
  which(jj7==max(jj7))->hjkk5
  hjkk5
  grep(names(hjkk5),colnames(jj))->hjkk5
  hjkk5
  ##
  apply(jj[c(aaj,bbj,hjkk,hjkk2,hjkk3,hjkk4,hjkk5),-c(aaj,bbj,hjkk,hjkk2,hjkk3,hjkk4,hjkk5)], MARGIN=2, FUN=min)->jj8
  max(jj8)
  which(jj8==max(jj8))->hjkk6
  hjkk6
  grep(names(hjkk6),colnames(jj))->hjkk6
  hjkk6
  ##
  apply(jj[c(aaj,bbj,hjkk,hjkk2,hjkk3,hjkk4,hjkk5,hjkk6),-c(aaj,bbj,hjkk,hjkk2,hjkk3,hjkk4,hjkk5,hjkk6)], MARGIN=2, FUN=min)->jj9
  max(jj9)
  which(jj9==max(jj9))->hjkk7
  hjkk7
  grep(names(hjkk7),colnames(jj))->hjkk7
  hjkk7
  ##
  apply(jj[c(aaj,bbj,hjkk,hjkk2,hjkk3,hjkk4,hjkk5,hjkk6,hjkk7),-c(aaj,bbj,hjkk,hjkk2,hjkk3,hjkk4,hjkk5,hjkk6,hjkk7)], MARGIN=2, FUN=min)->jj10
  max(jj10)
  which(jj10==max(jj10))->hjkk8
  hjkk8
  grep(names(hjkk8),colnames(jj))->hjkk8
  hjkk8
  ##
  apply(jj[c(aaj,bbj,hjkk,hjkk2,hjkk3,hjkk4,hjkk5,hjkk6,hjkk7,hjkk8),-c(aaj,bbj,hjkk,hjkk2,hjkk3,hjkk4,hjkk5,hjkk6,hjkk7,hjkk8)], MARGIN=2, FUN=min)->jj11
  max(jj11)
  which(jj11==max(jj11))->hjkk9
  hjkk9
  grep(names(hjkk9),colnames(jj))->hjkk9
  hjkk9
  ##
  apply(jj[c(aaj,bbj,hjkk,hjkk2,hjkk3,hjkk4,hjkk5,hjkk6,hjkk7,hjkk8,hjkk9),-c(aaj,bbj,hjkk,hjkk2,hjkk3,hjkk4,hjkk5,hjkk6,hjkk7,hjkk8,hjkk9)], MARGIN=2, FUN=min)->jj12
  max(jj12)
  which(jj12==max(jj12))->hjkk10
  hjkk10
  grep(names(hjkk10),colnames(jj))->hjkk10
  hjkk10
  ##
  apply(jj[c(aaj,bbj,hjkk,hjkk2,hjkk3,hjkk4,hjkk5,hjkk6,hjkk7,hjkk8,hjkk9,hjkk10),-c(aaj,bbj,hjkk,hjkk2,hjkk3,hjkk4,hjkk5,hjkk6,hjkk7,hjkk8,hjkk9,hjkk10)], MARGIN=2, FUN=min)->jj13
  max(jj13)
  which(jj13==max(jj13))->hjkk11
  hjkk11
  grep(names(hjkk11),colnames(jj))->hjkk11
  hjkk11
  ##
  apply(jj[c(aaj,bbj,hjkk,hjkk2,hjkk3,hjkk4,hjkk5,hjkk6,hjkk7,hjkk8,hjkk9,hjkk10,hjkk11),-c(aaj,bbj,hjkk,hjkk2,hjkk3,hjkk4,hjkk5,hjkk6,hjkk7,hjkk8,hjkk9,hjkk10,hjkk11)], MARGIN=2, FUN=min)->jj14
  max(jj14)
  which(jj14==max(jj14))->hjkk12
  hjkk12
  grep(names(hjkk12),colnames(jj))->hjkk12
  hjkk12
  ##
  apply(jj[c(aaj,bbj,hjkk,hjkk2,hjkk3,hjkk4,hjkk5,hjkk6,hjkk7,hjkk8,hjkk9,hjkk10,hjkk11,hjkk12),-c(aaj,bbj,hjkk,hjkk2,hjkk3,hjkk4,hjkk5,hjkk6,hjkk7,hjkk8,hjkk9,hjkk10,hjkk11,hjkk12)], MARGIN=2, FUN=min)->jj15
  max(jj15)
  which(jj15==max(jj15))->hjkk13
  hjkk13
  grep(names(hjkk13),colnames(jj))->hjkk13
  hjkk13
  ##
  apply(jj[c(aaj,bbj,hjkk,hjkk2,hjkk3,hjkk4,hjkk5,hjkk6,hjkk7,hjkk8,hjkk9,hjkk10,hjkk11,hjkk12,hjkk13),-c(aaj,bbj,hjkk,hjkk2,hjkk3,hjkk4,hjkk5,hjkk6,hjkk7,hjkk8,hjkk9,hjkk10,hjkk11,hjkk12,hjkk13)], MARGIN=2, FUN=min)->jj16
  max(jj16)
  which(jj16==max(jj16))->hjkk14
  hjkk14
  grep(names(hjkk14),colnames(jj))->hjkk14
  hjkk14
  ##
  apply(jj[c(aaj,bbj,hjkk,hjkk2,hjkk3,hjkk4,hjkk5,hjkk6,hjkk7,hjkk8,hjkk9,hjkk10,hjkk11,hjkk12,hjkk13,hjkk14),-c(aaj,bbj,hjkk,hjkk2,hjkk3,hjkk4,hjkk5,hjkk6,hjkk7,hjkk8,hjkk9,hjkk10,hjkk11,hjkk12,hjkk13,hjkk14)], MARGIN=2, FUN=min)->jj17
  max(jj17)
  which(jj17==max(jj17))->hjkk15
  hjkk15
  grep(names(hjkk15),colnames(jj))->hjkk15
  hjkk15
  ##
  apply(jj[c(aaj,bbj,hjkk,hjkk2,hjkk3,hjkk4,hjkk5,hjkk6,hjkk7,hjkk8,hjkk9,hjkk10,hjkk11,hjkk12,hjkk13,hjkk14,hjkk15),-c(aaj,bbj,hjkk,hjkk2,hjkk3,hjkk4,hjkk5,hjkk6,hjkk7,hjkk8,hjkk9,hjkk10,hjkk11,hjkk12,hjkk13,hjkk14,hjkk15)], MARGIN=2, FUN=min)->jj18
  max(jj18)
  which(jj18==max(jj18))->hjkk16
  hjkk16
  grep(names(hjkk16),colnames(jj))->hjkk16
  hjkk16
  ##
  apply(jj[c(aaj,bbj,hjkk,hjkk2,hjkk3,hjkk4,hjkk5,hjkk6,hjkk7,hjkk8,hjkk9,hjkk10,hjkk11,hjkk12,hjkk13,hjkk14,hjkk15,hjkk16),-c(aaj,bbj,hjkk,hjkk2,hjkk3,hjkk4,hjkk5,hjkk6,hjkk7,hjkk8,hjkk9,hjkk10,hjkk11,hjkk12,hjkk13,hjkk14,hjkk15,hjkk16)], MARGIN=2, FUN=min)->jj19
  max(jj19)
  which(jj19==max(jj19))->hjkk17
  hjkk17
  grep(names(hjkk17),colnames(jj))->hjkk17
  hjkk17
  ##
  apply(jj[c(aaj,bbj,hjkk,hjkk2,hjkk3,hjkk4,hjkk5,hjkk6,hjkk7,hjkk8,hjkk9,hjkk10,hjkk11,hjkk12,hjkk13,hjkk14,hjkk15,hjkk16,hjkk17),-c(aaj,bbj,hjkk,hjkk2,hjkk3,hjkk4,hjkk5,hjkk6,hjkk7,hjkk8,hjkk9,hjkk10,hjkk11,hjkk12,hjkk13,hjkk14,hjkk15,hjkk16,hjkk17)], MARGIN=2, FUN=min)->jj20
  max(jj20)
  which(jj20==max(jj20))->hjkk18
  hjkk18
  grep(names(hjkk18),colnames(jj))->hjkk18
  hjkk18

  select(WantedData, c(names(jj[aaj]),
                       names(jj[bbj]),
                       names(jj[hjkk]),
                       names(jj[hjkk2]),
                       names(jj[hjkk3]),
                       names(jj[hjkk4]),
                       names(jj[hjkk5]),
                       names(jj[hjkk6]),
                       names(jj[hjkk7]),
                       names(jj[hjkk8]),
                       names(jj[hjkk9]),
                       names(jj[hjkk10]),
                       names(jj[hjkk11]),
                       names(jj[hjkk12]),
                       names(jj[hjkk13]),
                       names(jj[hjkk14]),
                       names(jj[hjkk15]),
                       names(jj[hjkk16]),
                       names(jj[hjkk17]),
                       names(jj[hjkk18])
  ))->WantedData3

  WantedData3[,1:num_cells]->WantedData3

  m_WantedData <- as.data.frame(WantedData3)
  m_WantedData[m_WantedData < Filter_value] <- 0
  m_WantedData[m_WantedData >= Filter_value] <- 1

  message('||||||||||||||........................... 1/3
      MOST INFORMATIVE GROUP OF CELLS')

  WantedData[,x:num_cols]->WantedData4
  m_WantedData2 <- as.data.frame(WantedData4)
  m_WantedData2[m_WantedData2 < Filter_value] <- 0
  m_WantedData2[m_WantedData2 >= Filter_value] <- 1
  dist(t(m_WantedData),method =  "euclidean",diag = TRUE, upper = TRUE) -> distanciaa
  arbolito <- NJ(distanciaa)

  set.seed(num_set_seed)
  selected_cells=num_cells
  sample_cells=list()
  n=n # num. replicas

  for (i in 1:n) {
    sample_cells=c(sample_cells, list(sample(1:ncol(m_WantedData2),selected_cells,replace=F)))
  }

  message('|||||||||||||||||||||||||||.............. 2/3
             SUBSAMPLING CELLS')

  distancias=vector()
  evolucion_total=list()
  porc_surfaceoma=list()
  #mean_sd_porc_surfaceoma=list()
  #data_frame_mean_sd_porc_surfaceoma=list()
  table_collapsador=list()
  colapsador=function(x) {paste0(x, collapse="")}

  parallel::detectCores()
  n.cores <- parallel::detectCores() - 1
  #create the cluster:
  my.cluster <- parallel::makeCluster(
    n.cores,
    type = "PSOCK"
  )
  print(my.cluster) # checkeamos el cluster, opcional
  #register it to be used by %dopar%
  doParallel::registerDoParallel(cl = my.cluster)
  #check if it is registered (optional). Must be TRUE!
  foreach::getDoParRegistered()
  foreach::getDoParWorkers()

  evolucion_total<-foreach(j = 1:n, .packages=c("ape","phangorn","seqinr")) %dopar% {
    dist(t(m_WantedData2[,sample_cells[[j]]]),method =  "euclidean",diag = TRUE, upper = TRUE) -> distancias
    tree <- NJ(distancias)
    sum(tree$edge.length)
  }

  porc_surfaceoma<-foreach(j = 1:n, .packages=c("base")) %dopar% {
    1-(length(which(rowSums(m_WantedData2[,sample_cells[[j]]])==0))/nrow(m_WantedData2))->porc_surfaceoma[[j]]
  }

  table_collapsador<-foreach(j = 1:n, .packages=c("base")) %dopar% {
    table(table(apply(m_WantedData2[,sample_cells[[j]]],1, colapsador)))->table_collapsador[[j]]
  }


  #for (j in 1:n) {
  #  dist(t(m_WantedData2[,sample_cells[[j]]]),method =  "euclidean",diag = TRUE, upper = TRUE) -> distancias
  #  tree <- NJ(distancias)
  #  sum(tree$edge.length)->evolucion_total[j]
#
 #   1-(length(which(rowSums(m_WantedData2[,sample_cells[[j]]])==0))/nrow(m_WantedData2))->porc_surfaceoma[[j]]
#
 #   table(table(apply(m_WantedData2[,sample_cells[[j]]],1, colapsador)))->table_collapsador[[j]]
#
 # }

  parallel::stopCluster(cl = my.cluster)

  message('||||||||||||||||||||||||||||||||||||||||| 3/3
                 PLOTTING')

  qplot(unlist(evolucion_total),xlab = paste("Sum value lenght of each brach tree simulation (n= ", n, ")",  sep = ""),ylab = "Frecuency",binwidth = 1,fill=I("#1ac6ff"),col=I("black"),
        breaks=seq(min(unlist(evolucion_total)), max(unlist(evolucion_total)), by = 2)) +
    geom_vline(xintercept = sum(arbolito$edge.length), linetype="solid", color = "red", size=1)+theme_bw()+ggtitle(paste("Sum branch distance of each tree simulation (n= ", n, ")",  sep = ""))->>qgraph

  qgraph

  my_vector1 <- unlist(porc_surfaceoma)
  as.data.frame(my_vector1)->valores

  mean(my_vector1)->mean_porc_surfaceoma
  sd(my_vector1)->sd_porc_surfaceoma

  # Obtenemos el valor porcentaje surfaceoma en la matriz metodo secuencial:
  1-(length(which(rowSums(m_WantedData)==0))/nrow(m_WantedData))->porc_surfaceoma_metodo_sec

  x1=paste("Simulation data (n= ", n, ")",  sep = "")

  ggplot(data=valores, aes(x=x1,y=my_vector1))+
    stat_boxplot(geom ='errorbar')+geom_boxplot(fatten = 3)+
    geom_hline(yintercept = porc_surfaceoma_metodo_sec, linetype="solid", color = "red", size=1)+ labs(x = "", y="% Genes")+ggtitle(paste("Boxplot % Genes: \nSequential method vs. No sequential method (n= ", n, ")",  sep = ""))+theme_bw()+geom_point(aes(y=porc_surfaceoma_metodo_sec), colour="red",size=I(3))->>plot_surfaceoma
  plot_surfaceoma

  lapply(table_collapsador, `[[`, 1)->primer_table_collapsador
  lapply(table_collapsador, `[[`, 2)->segon_table_collapsador
  lapply(table_collapsador, `[[`, 3)->tercer_table_collapsador
  ### solo primer###
  mean(unlist(primer_table_collapsador))->mean_primer_table_collapsador
  sd(unlist(primer_table_collapsador))->sd_primer_table_collapsador


  #mean_primer_table_collapsador+sd_primer_table_collapsador->valor_superior
  #mean_primer_table_collapsador+sd_primer_table_collapsador->valor_inferior
  mean_primer_table_collapsador/nrow(m_WantedData2)->coc_genes
  #as.data.frame(coc_genes)->cociente_genes
  #valor_superior/nrow(m_WantedData2)->coc_genes_superior
  #valor_inferior/nrow(m_WantedData2)->coc_genes_inferior
  #as.data.frame(coc_genes_superior)->cociente_genes_superior
  #as.data.frame(coc_genes_inferior)->cociente_genes_inferior


  x2="Data simulation (mean and sd (upper and lower))"

  table(table(apply(m_WantedData2,1, colapsador)))->table_collapsador2

  qplot(x2,coc_genes,size=I(3))+geom_errorbar(aes(x=x2, ymin=coc_genes-(sd_primer_table_collapsador/nrow(m_WantedData2)), ymax=coc_genes+(sd_primer_table_collapsador/nrow(m_WantedData2))), width=0.12)+geom_point(aes(y=(table_collapsador2[1]/nrow(m_WantedData2))), colour="red",size=I(3))+theme_bw()+xlab("") + ylab("% unique genes, expression profile")+ ggtitle(paste("Value unique genes \nExpression profile (n= ", n, ")",  sep = ""))->>qplot_unique_expression_profile
  plot(qplot_unique_expression_profile)
  plot(plot_surfaceoma)
  plot(qgraph)
  names(WantedData3)->>names_group_cells
  message('Plots saved in these R objects: qplot_unique_expression_profile, plot_surfaceoma and qgraph. In names_group_cells are the names of cells')
}
