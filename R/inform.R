#' @import dplyr
#' @title More informatives
#' @description More informatives in a dataset
#' @details Computes a binary table from a table counts using a cut-off value (= 1 BY DEFAULT
#' if you don't specify a value).
#' If Filter_value is = or > will be computed as 1. If is < will be computed as 0.
#'
#' Then, the function computes the 20 most informative columns / samples / cells of the dataset,
#' showing you on the screen the percentage increase of information and the list
#' of the 20 most informative cells/groups (ordered from highest to lowest).
#'
#' \strong{IMPORTANT!: YOU MUST USE A MATRIX ONLY WITH DATA}
#' @param x Object in R that contains the table (raw counts, tpm counts, pkm counts...) to apply the filter.
#' @param Filter_value Cut-off value. Optional. By default = 1.
#' @examples
#' addinfo(data_RNASeq_raw_counts) # Internal dataset
#' hg38 # When you use the function extract(), you can load a gene data file from Rsubread package
#' tpm_method(raw_counts_addinfo,gene_lenght_file)
#' # To use inform() you MUST use only the matrix with the dataset
#' tpm_sample[,2:6] -> new_sample
#' inform(new_sample,10000)
#' # The most informative columns are S1,S2,S3 and S5 (S4 is not informative).
#' # S1.1, ... appears because we have less than 20 columns in the dataset.
#' # Don't worry, ignore.
#' @author Jose Jordan
#' @export
inform <- function(x,Filter_value) {

  m_WantedData <- as.data.frame(x)

  if (missing(Filter_value)) {
    Filter_value=1
    message('
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!++++++++++++++++++++++++++++++++++++++!!
!!+ ALERT!: Filter_value = 1 (DEFAULT) +!!
!!++++++++++++++++++++++++++++++++++++++!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ')

  }

  if (!is.numeric(Filter_value)) {
    Filter_value=1
    message('
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!++++++++++++++++++++++++++++++++++++++!!
!!+ ALERT!: Filter_value = 1 (DEFAULT) +!!
!!++++++++++++++++++++++++++++++++++++++!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ')

  }


  message('WAIT... BE PATIENT... O.o o.O
          ')


  m_WantedData[m_WantedData < Filter_value] <- 0
  m_WantedData[m_WantedData >= Filter_value] <- 1
  m_WantedData[as.logical(rowSums(m_WantedData != 0)), ]->m_WantedData_final
  m=m_WantedData_final
  s=nrow(m)
  c=combn(1:ncol(m),2)
  result=vector()
  for(i in 1:ncol(c)){
    result[i]=sum(apply(X = cbind(m[,c[1,i]],m[,c[2,i]]),MARGIN = 1, FUN = sum)>0)/s

  }
  c[,which.max(result)]->maxresult
  m[,c(maxresult[1],maxresult[2])]->informativos
  rowSums(informativos)->infor2
  infor2[infor2 < 1] <- 0
  infor2[infor2 >= 1] <- 1
  as.matrix(infor2)->infor3
  sum(infor3)/nrow(infor3)->val_porc
  val=vector()
  for (j in 1:ncol(m)) {
    infor3+m[,j]->mm
    mm[mm < 1] <- 0
    mm[mm >= 1] <- 1
    as.matrix(mm)->mmm
    sum(mmm)/nrow(mmm)->val[j]
    val[j]
  }
  val
  which.max(val)->maxresultval
  infor3+m[,which.max(val)]->mm
  mm[mm < 1] <- 0
  mm[mm >= 1] <- 1
  as.matrix(mm)->mmm
  sum(mmm)/nrow(mmm)->val
  val

  val2=vector()
  for (j in 1:ncol(m)) {
    mmm++m[,j]->jj
    jj[jj < 1] <- 0
    jj[jj >= 1] <- 1
    as.matrix(jj)->jjj
    sum(jjj)/nrow(jjj)->val2[j]
    val2[j]
  }
  val2
  which.max(val2)->maxresultval2
  mmm++m[,which.max(val2)]->jj
  jj[jj < 1] <- 0
  jj[jj >= 1] <- 1
  as.matrix(jj)->jjj
  sum(jjj)/nrow(jjj)->val2
  val2

  val3=vector()
  for (j in 1:ncol(m)) {
    jjj+m[,j]->kk
    kk[kk < 1] <- 0
    kk[kk >= 1] <- 1
    as.matrix(kk)->kkk
    sum(kkk)/nrow(kkk)->val3[j]
    val3[j]
  }
  val3
  which.max(val3)->maxresultval3
  jjj+m[,which.max(val3)]->kk
  kk[kk < 1] <- 0
  kk[kk >= 1] <- 1
  as.matrix(kk)->kkk
  sum(kkk)/nrow(kkk)->val3
  val3

  val4=vector()
  for (j in 1:ncol(m)) {
    kkk+m[,j]->ll
    ll[ll < 1] <- 0
    ll[ll >= 1] <- 1
    as.matrix(ll)->lll
    sum(lll)/nrow(lll)->val4[j]
    val4[j]
  }
  val4
  which.max(val4)->maxresultval4
  kkk+m[,which.max(val4)]->ll
  ll[ll < 1] <- 0
  ll[ll >= 1] <- 1
  as.matrix(ll)->lll
  sum(lll)/nrow(lll)->val4
  val4


  val5=vector()
  for (j in 1:ncol(m)) {
    lll+m[,j]->nn
    nn[nn < 1] <- 0
    nn[nn >= 1] <- 1
    as.matrix(nn)->nnn
    sum(nnn)/nrow(nnn)->val5[j]
    val5[j]
  }
  val5
  which.max(val5)->maxresultval5
  lll+m[,which.max(val5)]->nn
  nn[nn < 1] <- 0
  nn[nn >= 1] <- 1
  as.matrix(nn)->nnn
  sum(nnn)/nrow(nnn)->val5
  val5


  val6=vector()
  for (j in 1:ncol(m)) {
    nnn+m[,j]->bb
    bb[bb < 1] <- 0
    bb[bb >= 1] <- 1
    as.matrix(bb)->bbb
    sum(bbb)/nrow(bbb)->val6[j]
    val6[j]
  }
  val6
  which.max(val6)->maxresultval6
  nnn+m[,which.max(val6)]->bb
  bb[bb < 1] <- 0
  bb[bb >= 1] <- 1
  as.matrix(bb)->bbb
  sum(bbb)/nrow(bbb)->val6
  val6


  val7=vector()
  for (j in 1:ncol(m)) {
    bbb+m[,j]->pp
    pp[pp < 1] <- 0
    pp[pp >= 1] <- 1
    as.matrix(pp)->ppp
    sum(ppp)/nrow(ppp)->val7[j]
    val7[j]
  }
  val7
  which.max(val7)->maxresultval7
  bbb+m[,which.max(val7)]->pp
  pp[pp < 1] <- 0
  pp[pp >= 1] <- 1
  as.matrix(pp)->ppp
  sum(ppp)/nrow(ppp)->val7
  val7


  val8=vector()
  for (j in 1:ncol(m)) {
    ppp+m[,j]->qq
    qq[qq < 1] <- 0
    qq[qq >= 1] <- 1
    as.matrix(qq)->qqq
    sum(qqq)/nrow(qqq)->val8[j]
    val8[j]
  }
  val8
  which.max(val8)->maxresultval8
  ppp+m[,which.max(val8)]->qq
  qq[qq < 1] <- 0
  qq[qq >= 1] <- 1
  as.matrix(qq)->qqq
  sum(qqq)/nrow(qqq)->val8
  val8




  val9=vector()
  for (j in 1:ncol(m)) {
    qqq+m[,j]->rr
    rr[rr < 1] <- 0
    rr[rr >= 1] <- 1
    as.matrix(rr)->rrr
    sum(rrr)/nrow(rrr)->val9[j]
    val9[j]
  }
  val9
  which.max(val9)->maxresultval9
  qqq+m[,which.max(val9)]->rr
  rr[rr < 1] <- 0
  rr[rr >= 1] <- 1
  as.matrix(rr)->rrr
  sum(rrr)/nrow(rrr)->val9
  val9



  val10=vector()
  for (j in 1:ncol(m)) {
    rrr+m[,j]->ss
    ss[ss < 1] <- 0
    ss[ss >= 1] <- 1
    as.matrix(ss)->sss
    sum(sss)/nrow(sss)->val10[j]
    val10[j]
  }
  val10
  which.max(val10)->maxresultval10
  rrr+m[,which.max(val10)]->ss
  ss[ss < 1] <- 0
  ss[ss >= 1] <- 1
  as.matrix(ss)->sss
  sum(sss)/nrow(sss)->val10
  val10


  val11=vector()
  for (j in 1:ncol(m)) {
    sss+m[,j]->tt
    tt[tt < 1] <- 0
    tt[tt >= 1] <- 1
    as.matrix(tt)->ttt
    sum(ttt)/nrow(ttt)->val11[j]
    val11[j]
  }
  val11
  which.max(val11)->maxresultval11
  sss+m[,which.max(val11)]->tt
  tt[tt < 1] <- 0
  tt[tt >= 1] <- 1
  as.matrix(tt)->ttt
  sum(ttt)/nrow(ttt)->val11
  val11

  val12=vector()
  for (j in 1:ncol(m)) {
    ttt+m[,j]->uu
    uu[uu < 1] <- 0
    uu[uu >= 1] <- 1
    as.matrix(uu)->uuu
    sum(uuu)/nrow(uuu)->val12[j]
    val12[j]
  }
  val12
  which.max(val12)->maxresultval12
  ttt+m[,which.max(val12)]->uu
  uu[uu < 1] <- 0
  uu[uu >= 1] <- 1
  as.matrix(uu)->uuu
  sum(uuu)/nrow(uuu)->val12
  val12


  val13=vector()
  for (j in 1:ncol(m)) {
    uuu+m[,j]->vv
    vv[vv < 1] <- 0
    vv[vv >= 1] <- 1
    as.matrix(vv)->vvv
    sum(vvv)/nrow(vvv)->val13[j]
    val13[j]
  }
  val13
  which.max(val13)->maxresultval13
  uuu+m[,which.max(val13)]->vv
  vv[vv < 1] <- 0
  vv[vv >= 1] <- 1
  as.matrix(vv)->vvv
  sum(vvv)/nrow(vvv)->val13
  val13


  val14=vector()
  for (j in 1:ncol(m)) {
    vvv+m[,j]->ww
    ww[ww < 1] <- 0
    ww[ww >= 1] <- 1
    as.matrix(ww)->www
    sum(www)/nrow(www)->val14[j]
    val14[j]
  }
  val14
  which.max(val14)->maxresultval14
  vvv+m[,which.max(val14)]->ww
  ww[ww < 1] <- 0
  ww[ww >= 1] <- 1
  as.matrix(ww)->www
  sum(www)/nrow(www)->val14
  val14



  val15=vector()
  for (j in 1:ncol(m)) {
    www+m[,j]->xx
    xx[xx < 1] <- 0
    xx[xx >= 1] <- 1
    as.matrix(xx)->xxx
    sum(xxx)/nrow(xxx)->val15[j]
    val15[j]
  }
  val15
  which.max(val15)->maxresultval15
  www+m[,which.max(val15)]->xx
  xx[xx < 1] <- 0
  xx[xx >= 1] <- 1
  as.matrix(xx)->xxx
  sum(xxx)/nrow(xxx)->val15
  val15


  val16=vector()
  for (j in 1:ncol(m)) {
    xxx+m[,j]->yi
    yi[yi < 1] <- 0
    yi[yi >= 1] <- 1
    as.matrix(yi)->yiy
    sum(yiy)/nrow(yiy)->val16[j]
    val16[j]
  }
  val16
  which.max(val16)->maxresultval16
  xxx+m[,which.max(val16)]->yi
  yi[yi < 1] <- 0
  yi[yi >= 1] <- 1
  as.matrix(yi)->yiy
  sum(yiy)/nrow(yiy)->val16
  val16


  val17=vector()
  for (j in 1:ncol(m)) {
    yiy+m[,j]->zz
    zz[zz < 1] <- 0
    zz[zz >= 1] <- 1
    as.matrix(zz)->zzz
    sum(zzz)/nrow(zzz)->val17[j]
    val17[j]
  }
  val17
  which.max(val17)->maxresultval17
  yiy+m[,which.max(val17)]->zz
  zz[zz < 1] <- 0
  zz[zz >= 1] <- 1
  as.matrix(zz)->zzz
  sum(zzz)/nrow(zzz)->val17
  val17


  val18=vector()
  for (j in 1:ncol(m)) {
    zzz+m[,j]->zaa
    zaa[zaa < 1] <- 0
    zaa[zaa >= 1] <- 1
    as.matrix(zaa)->qwe
    sum(qwe)/nrow(qwe)->val18[j]
    val18[j]
  }
  val18
  which.max(val18)->maxresultval18
  zzz+m[,which.max(val18)]->zaa
  zaa[zaa < 1] <- 0
  zaa[zaa >= 1] <- 1
  as.matrix(zaa)->qwe
  sum(qwe)/nrow(qwe)->val18
  val18


  resultado_valor_porcentual <- c(val_porc,val,val2,val3,val4,val5,val6,val7,val8,val9,val10,val11,val12,val13,val14,val15,val16,val17,val18)
  celulas <- names(m[,c(maxresult,maxresultval,maxresultval2,maxresultval3,maxresultval4,maxresultval5,maxresultval6,maxresultval7,maxresultval8,
                        maxresultval9,maxresultval10,maxresultval11,maxresultval12,maxresultval13,maxresultval14,maxresultval15,maxresultval16,
                        maxresultval17,maxresultval18)])
  message('PERCENTAGE VALUES FOR EACH CELL ADDED:')
  print(resultado_valor_porcentual)
  message('SELECTED CELLS (from + to - informative in the dataset):')
  print(celulas)

  message('PLOTING :)')
  plot(y=resultado_valor_porcentual,x=c(2:20),xlab = "# cells used",ylab = "% value sum",type = "b", main = "% value variation (20 + informative cells)")


  rm(resultado_valor_porcentual,celulas)
}
