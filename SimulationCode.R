# Filename: SimulationCode.R
# Authors: Francesca Graziano, Valsecchi Maria Grazia, Paola Rebora
# Article: Sampling strategies to evaluate the prognostic value of a new biomarker on a time-to-event end-point

#  @param lambda scale parameter
#  @param k shape parameter
#  @param beta list of beta coefficents 
#  @param X matrix of beta values
#  @param rate rate of censoring
#  @param follow follow-up 
#  @param N sample size of full cohort
#  @param n sample size of two-phase
#  @param ssed = Random Number Generation
#  Returns list of B results in each different sampling design   

# All rights reserved.

simul<-function(lambda, k,  beta,X, rate,follow,N,n,B,ssed){
  i<-1;Full_Cohort<-NULL;Rsimple<-NULL;CC_ev<- NULL; CC_intEV_elfin<-NULL;CC_intEV_conf<-NULL;CC_intEV_surr<-NULL;
  CC_stra_elfin<-NULL;CC_stra_conf<-NULL;CC_stra_surr<-NULL;
  PPS_EVe<-NULL; PPS_EVe_elf<-NULL; PPS_EVsurr<-NULL;PPS_EVconf<-NULL;CC_allCases<-NULL;
  NCC_1<-NULL;NCCsur<-NULL; CM<-NULL;
  
  set.seed(ssed)
 genSAMPLEexpWeiC<-function(N, lambda, k, beta, X, rate, follow) {
  u <- runif(N, min=0, max=1)
  
  if (k=="1") 
   t<- -1/(lambda*exp(X %*% beta))*log(1-u)
   if (k!="1")
   t <- (-log(1-u)/(lambda^k*exp( X %*% beta)))^(1/k)
  
  cens.time <- rexp(length(t), rate = rate)                 
  cens.time <- ifelse(cens.time > follow, follow, cens.time)
  cens <- ifelse (t > cens.time, 0 , 1)
  time <- ifelse(t < cens.time, t, cens.time)
  id<-1:N
  return(data.frame(id, X, time, cens))
   }
  for (i in 1:B) {
   
    #GENERATE 1st PHASE - Exponential/Weibull
    dati<-genSAMPLEexpWeiC(N, lambda, k, beta, X, rate, follow)
    Fit<-(coxph(Surv(time,cens)~x1 + x4, dati))
    Summ<-summary(Fit)
    sample.est<- c(Summ$coefficients[1,],
                   Summ$conf.int[1,3:4],
                   Summ$coefficients[2,c(1,3)])
    
    surv.prob <- summary(survfit(Surv(time=time,event=cens)~x1, dati),time=2)$surv; 
    if (length(surv.prob)==0) {
      surv.prob<-c(NA,NA) }
    if (length(surv.prob)==1)  surv.prob<- c(surv.prob,surv.prob)
    if (length(surv.prob)==2)   surv.prob<- surv.prob
    
    pX<-length(dati$x1[dati$x1=="1"])/N
    
    T1<-table(dati$x3,dati$x1)
    SENS<-T1[2,2]/(T1[2,2]+T1[1,2]);SENS
    SPEC<-T1[1,1]/(T1[1,1]+T1[2,1]);SPEC
    
    Full_Cohort<-rbind(Full_Cohort,c(sample=i,N,pX, surv.prob, sample.est, NA, SENS, SPEC))
    colnames(Full_Cohort)<-c("sample", "n", "p_x1", "surv_x0", "surv_x1", "coef", "exp coef", "se.coef", "z", "pvalue", "lower95HR", "upper95HR","BConf","seConf","DEFF", "sens", "spec")
    
    ############TYPE OF SAMPLING: ################### 
    #1. SRS: Ogni campione ha la stessa probabilità di essere estratto
    
    ev=0
    idsampled<-sample(dati$id, n, replace = FALSE, prob = NULL) ;   
    campione<-dati[dati$id %in% idsampled,]  # CAMPIONE SELEZIONATO
    ev=sum(campione$cens==1)
    dati$incl<-ifelse(dati$id %in% idsampled,1,0); #table(dati$incl)
    des<-twophase(id=list(~id,~id),subset=~incl==1,strata=list(NULL,NULL),data=dati)
    
    Fit<-(svycoxph(Surv(time,cens)~x1 + x4, des))
    Summ<-summary(Fit)
    sample.est<- c(Summ$coefficients[1,],
                   Summ$conf.int[1,3:4],
                   Summ$coefficients[2,c(1,3)])
    
    DATse<-dati[dati$incl==1,]; dim(DATse) #table(DATse$x1) BIOMARKER
    
    T1<-table(DATse$x3,DATse$x1)
    SENS<-T1[2,2]/(T1[2,2]+T1[1,2]);SENS
    SPEC<-T1[1,1]/(T1[1,1]+T1[2,1]);SPEC
    
     surv.prob <- summary(survfit(Surv(time=time,event=cens)~x1, DATse),time=2)$surv; 
    if (length(surv.prob)==0) {
      surv.prob<-c(NA,NA) }
    if (length(surv.prob)==1)  surv.prob<- c(surv.prob,surv.prob)
    if (length(surv.prob)==2)  surv.prob<- surv.prob
    
    pX<-length(DATse$x1[DATse$x1=="1"])/length(DATse$x1)
    
    VAR_srs<-sample.est[3]^2
    DEFF<- 1
    
    Rsimple<-rbind(Rsimple,c(sample=i,round(length(DATse$x1),1),pX, surv.prob, sample.est,DEFF, SENS, SPEC) )
    colnames(Rsimple)<-c("sample", "n", "p_x1", "survx0", "survx1", "coef", "exp coef", "se.coef", "z", "pvalue", "lower95HR", "upper95HR","BConf","seConf","DEFF", "sens","spec")
    
    #######2. CASO-CONTROLLO - CC
    
    cases<-NULL;cntl<-NULL #length(dati$id[dati$cens==1])
    if(length(dati$id[dati$cens==1])>=n/2){
      cases<-sample(dati$id[dati$cens==1], n/2, replace = FALSE, prob = NULL)
      cntl<-sample(dati$id[!dati$cens==1],n-length(cases), replace = FALSE, prob = NULL)
    }
    if (length(dati$id[dati$cens==1])<n/2) cases<-dati$id[dati$cens==1]
    cntl<-sample(dati$id[!dati$cens==1],n/2, replace = FALSE, prob = NULL)
    
    idsampled<-c(cases,cntl); #length(idsampled)
   
    campione<-dati[dati$id %in% idsampled,]  #save selected ids 
    dati$incl<-ifelse(dati$id %in% idsampled,1,0); #table(dati$incl)
    
    DATse<-dati[dati$incl==1,]
    surv.prob <- summary(survfit(Surv(time=time,event=cens)~x1, DATse),time=2)$surv; 
    if (length(surv.prob)==0) {
      surv.prob<-c(NA,NA) }
    if (length(surv.prob)==1)  surv.prob<- c(surv.prob,surv.prob)
    if (length(surv.prob)==2)   surv.prob<- surv.prob
    
    T1<-table(DATse$x3,DATse$x1)
    SENS<-T1[2,2]/(T1[2,2]+T1[1,2]);SENS
    SPEC<-T1[1,1]/(T1[1,1]+T1[2,1]);SPEC
    
    pX<-length(DATse$x1[DATse$x1=="1"])/length(DATse$x1)

    ##2.1 CC EVENTO
    des<-twophase(id=list(~id,~id),subset=~incl==1,strata=list(NULL,~cens),data=dati)
    
    Fit<-svycoxph(Surv(time,cens)~x1 + x4, des)
    Summ<-summary(Fit)
    sample.est<-c(Summ$coefficients[1,],
                   Summ$conf.int[1,3:4],
                   Summ$coefficients[2,c(1,3)])
    
    VAR_d<-sample.est[3]^2
    DEFF<- VAR_srs/VAR_d
    
    CC_ev<-rbind(CC_ev,c(sample=i,length(DATse$x1),pX, surv.prob, sample.est,DEFF, SENS, SPEC))
    colnames(CC_ev)<-c("sample", "n", "p_x1", "surv_x0", "surv_x1", "coef", "exp coef", "se.coef", "z", "pvalue", "lower95HR", "upper95HR","BConf","seConf", "DEFF", "sens","spec")
    
    ####2.2 CC - Evento Post stratified Elfin 
    des<-twophase(id=list(~id,~id),subset=~incl==1,strata=list(NULL,~interaction(x2, cens)),data=dati)
    
    Fit<-svycoxph(Surv(time,cens)~x1 + x4, des)
    Summ<-summary(Fit)
    sample.est<-c(Summ$coefficients[1,],
                  Summ$conf.int[1,3:4],
                  Summ$coefficients[2,c(1,3)])
    
    VAR_d<-sample.est[3]^2
    DEFF<- VAR_srs/VAR_d
    
    CC_intEV_elfin<-rbind( CC_intEV_elfin,c(sample=i,length(DATse$x1),pX, surv.prob, sample.est, DEFF, SENS, SPEC) )
    colnames( CC_intEV_elfin)<-c("sample", "n", "p_x1", "surv_x0", "surv_x1", "coef", "exp coef", "se.coef", "z", "pvalue", "lower95HR", "upper95HR","BConf","seConf", "DEFF", "sens","spec")
    
    #######2.3 CC - Evento Post stratified Confounder
    
    des<-twophase(id=list(~id,~id),subset=~incl==1,strata=list(NULL,~interaction(x4, cens)),data=dati)
    
    Fit<-svycoxph(Surv(time,cens)~x1 + x4, des)
    Summ<-summary(Fit)
    sample.est<-c(Summ$coefficients[1,],
                  Summ$conf.int[1,3:4],
                  Summ$coefficients[2,c(1,3)])
    
    VAR_d<-sample.est[3]^2
    DEFF<- VAR_srs/VAR_d
    
    CC_intEV_conf<-rbind( CC_intEV_conf,c(sample=i,length(DATse$x1),pX, surv.prob, sample.est, DEFF, SENS, SPEC) )
    colnames( CC_intEV_conf)<-c("sample","n", "p_x1", "surv_x0", "surv_x1",  "coef", "exp coef", "se.coef", "z", "pvalue", "lower95HR", "upper95HR","BConf","seConf", "DEFF", "sens","spec")
    
    #######2.4 CC - Evento Post stratified Surrogate/Auxiliary
    
    des<-twophase(id=list(~id,~id),subset=~incl==1,strata=list(NULL,~interaction(x3, cens)),data=dati)
    
    Fit<-svycoxph(Surv(time,cens)~x1 + x4, des)
    Summ<-summary(Fit)
    sample.est<-c(Summ$coefficients[1,],
                  Summ$conf.int[1,3:4],
                  Summ$coefficients[2,c(1,3)])
    
    VAR_d<-sample.est[3]^2
    DEFF<- VAR_srs/VAR_d
    
    CC_intEV_surr<-rbind(CC_intEV_surr,c(sample=i,length(DATse$x1),pX, surv.prob, sample.est, DEFF, SENS, SPEC) )
    colnames(CC_intEV_surr)<-c("sample", "n", "p_x1", "surv_x0", "surv_x1", "coef", "exp coef", "se.coef", "z", "pvalue", "lower95HR", "upper95HR","BConf","seConf","DEFF", "sens","spec")
    
    
    ##2.5 CC - Evento stra_elfin
    cases0<-NULL;cntl0<-NULL; cases1<-NULL;cntl1<-NULL #length(dati$id[dati$cens==1])
  
    
    if(length(dati$id[dati$cens==1 & dati$x2==0])>=n/4){
      cases0<-sample(dati$id[dati$cens==1 & dati$x2==0], n/4 ,replace = FALSE, prob = NULL) } 
    if(length(dati$id[dati$cens==1 & dati$x2==0])<n/4) cases0<-dati$id[dati$cens==1 & dati$x2==0]
    
    cntl0<-sample(dati$id[!dati$cens==1 & dati$x2==0], n/4 ,replace = FALSE , prob = NULL)
    
    if(length(dati$id[dati$cens==1 & dati$x2==1])>=n/4){
      cases1<-sample(dati$id[dati$cens==1 & dati$x2==1], n/4,replace = FALSE , prob = NULL)
      }
    if(length(dati$id[dati$cens==1 & dati$x2==1])<n/4) cases1<-dati$id[dati$cens==1 & dati$x2==1]

    cntl1<-sample(dati$id[!dati$cens==1 & dati$x2==1], n/4 , replace = FALSE, prob = NULL)
    
    idsampled<-c(cases0,cases1,cntl0,cntl1); length(idsampled)
    
    campione<-dati[dati$id %in% idsampled,]  #save selected ids 
    dati$incl<-ifelse(dati$id %in% idsampled,1,0); #table(dati$incl)
    
    DATse<-dati[dati$incl==1,]
    surv.prob <- summary(survfit(Surv(time=time,event=cens)~x1, DATse),time=2)$surv; 
    if (length(surv.prob)==0) {
      surv.prob<-c(NA,NA) }
    if (length(surv.prob)==1)  surv.prob<- c(surv.prob,surv.prob)
    if (length(surv.prob)==2)   surv.prob<- surv.prob
    
    T1<-table(DATse$x3,DATse$x1)
    SENS<-T1[2,2]/(T1[2,2]+T1[1,2]);SENS
    SPEC<-T1[1,1]/(T1[1,1]+T1[2,1]);SPEC
    
    pX<-length(DATse$x1[DATse$x1=="1"])/length(DATse$x1)
    
    des<-twophase(id=list(~id,~id),subset=~incl==1,strata=list(NULL,~interaction(x2, cens)),data=dati)
    
    Fit<-svycoxph(Surv(time,cens)~x1 + x4, des)
    Summ<-summary(Fit)
    sample.est<-c(Summ$coefficients[1,],
                  Summ$conf.int[1,3:4],
                  Summ$coefficients[2,c(1,3)])
    
    VAR_d<-sample.est[3]^2
    DEFF<- VAR_srs/VAR_d
    
    CC_stra_elfin<-rbind(CC_stra_elfin,c(sample=i,length(DATse$x1),pX, surv.prob, sample.est,DEFF, SENS, SPEC))
    colnames(CC_stra_elfin)<-c("sample", "n", "p_x1", "surv_x0", "surv_x1", "coef", "exp coef", "se.coef", "z", "pvalue", "lower95HR", "upper95HR","BConf","seConf", "DEFF", "sens","spec")
    
    ##2.6 CC - Evento stra_conf
    cases0<-NULL;cntl0<-NULL; cases1<-NULL;cntl1<-NULL #length(dati$id[dati$cens==1])
    
    if(length(dati$id[dati$cens==1 & dati$x4==0])>=n/4){
      cases0<-sample(dati$id[dati$cens==1 & dati$x4==0], n/4 ,replace = FALSE, prob = NULL)} 
    if(length(dati$id[dati$cens==1 & dati$x4==0])<n/4) cases0<-dati$id[dati$cens==1 & dati$x4==0]
    
    cntl0<-sample(dati$id[!dati$cens==1 & dati$x4==0], n/4 ,replace = FALSE , prob = NULL)
    
    if(length(dati$id[dati$cens==1 & dati$x4==1])>=n/4){
      cases1<-sample(dati$id[dati$cens==1 & dati$x4==1], n/4,replace = FALSE , prob = NULL)}
    if(length(dati$id[dati$cens==1 & dati$x4==1])<n/4) cases1<-dati$id[dati$cens==1 & dati$x4==1]
    
    cntl1<-sample(dati$id[!dati$cens==1 & dati$x4==1], n/4 , replace = FALSE, prob = NULL)
    
    idsampled<-c(cases0,cases1,cntl0,cntl1); length(idsampled)
    
    campione<-dati[dati$id %in% idsampled,]  #save selected ids 
    dati$incl<-ifelse(dati$id %in% idsampled,1,0); #table(dati$incl)
    
    DATse<-dati[dati$incl==1,]
    surv.prob <- summary(survfit(Surv(time=time,event=cens)~x1, DATse),time=2)$surv; 
    if (length(surv.prob)==0) {
      surv.prob<-c(NA,NA) }
    if (length(surv.prob)==1)  surv.prob<- c(surv.prob,surv.prob)
    if (length(surv.prob)==2)   surv.prob<- surv.prob
    T1<-table(DATse$x3,DATse$x1)
    SENS<-T1[2,2]/(T1[2,2]+T1[1,2]);SENS
    pX<-length(DATse$x1[DATse$x1=="1"])/length(DATse$x1)
    des<-twophase(id=list(~id,~id),subset=~incl==1,strata=list(NULL,~interaction(x4, cens)),data=dati)
    Fit<-svycoxph(Surv(time,cens)~x1 + x4, des)
    Summ<-summary(Fit)
    sample.est<-c(Summ$coefficients[1,],
                  Summ$conf.int[1,3:4],
                  Summ$coefficients[2,c(1,3)])
    VAR_d<-sample.est[3]^2
    DEFF<- VAR_srs/VAR_d
    CC_stra_conf<-rbind(CC_stra_conf,c(sample=i,length(DATse$x1),pX, surv.prob, sample.est,DEFF, SENS, SPEC))
    colnames(CC_stra_conf)<-c("sample", "n", "p_x1", "surv_x0", "surv_x1", "coef", "exp coef", "se.coef", "z", "pvalue", "lower95HR", "upper95HR","BConf","seConf", "DEFF", "sens","spec")
    
    ###################################
    ##2.7 CC - Evento stra_surr
   cases0<-NULL;cntl0<-NULL; cases1<-NULL;cntl1<-NULL #length(dati$id[dati$cens==1])
    if(length(dati$id[dati$cens==1 & dati$x3==0])>=n/4){
    cases0<-sample(dati$id[dati$cens==1 & dati$x3==0], n/4 ,replace = FALSE, prob = NULL) }
    if(length(dati$id[dati$cens==1 & dati$x3==0])<n/4) cases0<-dati$id[dati$cens==1 & dati$x3==0]
      
    cntl0<-sample(dati$id[!dati$cens==1 & dati$x3==0], n/4 ,replace = FALSE , prob = NULL)
    
    if(length(dati$id[dati$cens==1 & dati$x3==1])>=n/4){
    cases1<-sample(dati$id[dati$cens==1 & dati$x3==1], n/4,replace = FALSE , prob = NULL)}
    if(length(dati$id[dati$cens==1 & dati$x3==1])<n/4) cases1<-dati$id[dati$cens==1 & dati$x3==1]
    
    if(length(dati$id[!dati$cens==1 & dati$x3==1])>=n/4){
    cntl1<-sample(dati$id[!dati$cens==1 & dati$x3==1], n/4 , replace = FALSE, prob = NULL)}
    if(length(dati$id[!dati$cens==1 & dati$x3==1])<n/4) cntl1<-dati$id[!dati$cens==1 & dati$x3==1]
    #modificato con le linee sopra
    #cntl1<-sample(dati$id[!dati$cens==1 & dati$x3==1], n/4 , replace = FALSE, prob = NULL)
    
    idsampled<-c(cases0,cases1,cntl0,cntl1); length(idsampled)
    
    campione<-dati[dati$id %in% idsampled,]  #save selected ids 
    dati$incl<-ifelse(dati$id %in% idsampled,1,0); #table(dati$incl)
    
    DATse<-dati[dati$incl==1,]
    surv.prob <- summary(survfit(Surv(time=time,event=cens)~x1, DATse),time=2)$surv; 
    if (length(surv.prob)==0) {
      surv.prob<-c(NA,NA) }
    if (length(surv.prob)==1)  surv.prob<- c(surv.prob,surv.prob)
    if (length(surv.prob)==2)   surv.prob<- surv.prob
    
    T1<-table(DATse$x3,DATse$x1)
    SENS<-T1[2,2]/(T1[2,2]+T1[1,2]);SENS
    SPEC<-T1[1,1]/(T1[1,1]+T1[2,1]);SPEC
    
    pX<-length(DATse$x1[DATse$x1=="1"])/length(DATse$x1)
    
    des<-twophase(id=list(~id,~id),subset=~incl==1,strata=list(NULL,~interaction(x3, cens)),data=dati)
    
    Fit<-svycoxph(Surv(time,cens)~x1 + x4, des)
    Summ<-summary(Fit)
    sample.est<-c(Summ$coefficients[1,],
                  Summ$conf.int[1,3:4],
                  Summ$coefficients[2,c(1,3)])
    
    VAR_d<-sample.est[3]^2
    DEFF<- VAR_srs/VAR_d
    
    CC_stra_surr<-rbind(CC_stra_surr,c(sample=i,length(DATse$x1),pX, surv.prob, sample.est,DEFF, SENS, SPEC))
    colnames(CC_stra_surr)<-c("sample", "n", "p_x1", "surv_x0", "surv_x1", "coef", "exp coef", "se.coef", "z", "pvalue", "lower95HR", "upper95HR","BConf","seConf", "DEFF", "sens","spec")
    
    ###3. probability proportional to size (PPS) for 2 strata - EVENT
    
    T1<-table(dati$cens)
    cases<-sample(dati$id[dati$cens==1], size= T1[2]*n/N, replace = FALSE)
    cntl<-sample(dati$id[!dati$cens==1], size= T1[1]*n/N, replace = FALSE)
    idsampled<-c(cases,cntl) ; #length(idsampled)
    
    campione<-dati[dati$id %in% idsampled,]  # CAMPIONE SELEZIONATO
    dati$incl<-ifelse(dati$id %in% idsampled,1,0); #table(dati$incl)
    
    des<-twophase(id=list(~id,~id),subset=~incl==1,strata=list(NULL,~cens), data=dati)
    
    Fit<-svycoxph(Surv(time,cens)~x1 + x4, des)
    Summ<-summary(Fit)
    sample.est<-c(Summ$coefficients[1,],
                  Summ$conf.int[1,3:4],
                  Summ$coefficients[2,c(1,3)])
    
    DATse<-dati[dati$incl==1,]
   
    surv.prob <- summary(survfit(Surv(time=time,event=cens)~x1, DATse),time=2)$surv; 
    if (length(surv.prob)==0) {
      surv.prob<-c(NA,NA) }
    if (length(surv.prob)==1)  surv.prob<- c(surv.prob,surv.prob)
    if (length(surv.prob)==2)   surv.prob<- surv.prob
    
    pX<-length(DATse$x1[DATse$x1=="1"])/length(DATse$x1)
    
    T1<-table(DATse$x3,DATse$x1)
    SENS<-T1[2,2]/(T1[2,2]+T1[1,2]);SENS
    SPEC<-T1[1,1]/(T1[1,1]+T1[2,1]);SPEC
    
    VAR_d<-sample.est[3]^2
    DEFF<- VAR_srs/VAR_d   
    
    PPS_EVe<-rbind( PPS_EVe,c(sample=i,length(DATse$x1),pX, surv.prob, sample.est, DEFF, SENS, SPEC) )
    colnames(PPS_EVe)<-c("sample", "n", "p_x1", "surv_x0", "surv_x1", "coef", "exp coef", "se.coef", "z", "pvalue", "lower95HR", "upper95HR","BConf","seConf", "DEFF", "sens","spec")
    
    #####4. probability proportional to size (PPS) for 4 strata - ELFIN
    cases0<-NULL; cases1<-NULL; cntl0<- NULL ; cntl1<-NULL
    
    T1<-table(dati$cens,dati$x2); #str(T1)
    cases0<-sample(dati$id[dati$cens==1 & dati$x2==0], T1[2,1]*n/N ,replace = FALSE, prob = NULL) 
    cntl0<-sample(dati$id[!dati$cens==1 & dati$x2==0], T1[1,1]*n/N ,replace = FALSE , prob = NULL)
    cases1<-sample(dati$id[dati$cens==1 & dati$x2==1], T1[2,2]*n/N,replace = FALSE , prob = NULL)
    cntl1<-sample(dati$id[!dati$cens==1 & dati$x2==1],T1[1,2]*n/N , replace = FALSE, prob = NULL)
    
    idsampled<-c(cases0,cases1,cntl0,cntl1)
    #campione<-dati[dati$id %in% idsampled,]  # CAMPIONE SELEZIONATO
    dati$incl<-ifelse(dati$id %in% idsampled,1,0); #table(dati$incl)
    
    des<-twophase(id=list(~id,~id),subset=~incl==1,strata=list(NULL,~interaction(x2, cens)), data=dati)
    
    Fit<-svycoxph(Surv(time,cens)~x1 + x4, des)
    Summ<-summary(Fit)
    sample.est<-c(Summ$coefficients[1,],
                  Summ$conf.int[1,3:4],
                  Summ$coefficients[2,c(1,3)])
    
    DATse<-dati[dati$incl==1,]
    #fit <- survfit(Surv(time=time,event=cens)~x1, DATse)
    surv.prob <- summary(survfit(Surv(time=time,event=cens)~x1, DATse),time=2)$surv; 
    if (length(surv.prob)==0) {
      surv.prob<-c(NA,NA) }
    if (length(surv.prob)==1)  surv.prob<- c(surv.prob,surv.prob)
    if (length(surv.prob)==2)   surv.prob<- surv.prob
    
    pX<-length(DATse$x1[DATse$x1=="1"])/length(DATse$x1)
    
    T1<-table(DATse$x3,DATse$x1)
    SENS<-T1[2,2]/(T1[2,2]+T1[1,2]);SENS
    SPEC<-T1[1,1]/(T1[1,1]+T1[2,1]);SPEC
    
    VAR_d<-sample.est[3]^2
    DEFF<- VAR_srs/VAR_d
    
    PPS_EVe_elf<-rbind(PPS_EVe_elf,c(sample=i,length(DATse$x1),pX, surv.prob, sample.est, DEFF, SENS, SPEC))
    colnames(PPS_EVe_elf)<-c("sample", "n", "p_x1", "surv_x0", "surv_x1", "coef", "exp coef", "se.coef", "z", "pvalue", "lower95HR", "upper95HR","BConf","seConf", "DEFF", "sens","spec")
    
    ##### 5. probability proportional to size (PPS) for 4 strata SURROGATE/AUXILIARY
    
    T1<-table(dati$cens,dati$x3); 
    cases0<-sample(dati$id[dati$cens==1 & dati$x3==0], T1[2,1]*n/N ,replace = FALSE, prob = NULL) 
    cntl0<-sample(dati$id[!dati$cens==1 & dati$x3==0], T1[1,1]*n/N,replace = FALSE , prob = NULL)
    cases1<-sample(dati$id[dati$cens==1 & dati$x3==1], T1[2,2]*n/N,replace = FALSE , prob = NULL)
    cntl1<-sample(dati$id[!dati$cens==1 & dati$x3==1],T1[1,2]*n/N, replace = FALSE, prob = NULL)
    
    idsampled<-c(cases0,cases1,cntl0,cntl1)
    
    dati$incl<-ifelse(dati$id %in% idsampled,1,0);#table(dati$incl)
    
    des<-twophase(id=list(~id,~id),subset=~incl==1,strata=list(NULL,~interaction(x3, cens)), data=dati)
  
    Fit<-svycoxph(Surv(time,cens)~x1 + x4, des)
    Summ<-summary(Fit)
    sample.est<-c(Summ$coefficients[1,],
                  Summ$conf.int[1,3:4],
                  Summ$coefficients[2,c(1,3)])
    
    DATse<-dati[dati$incl==1,]
    #fit <- survfit(Surv(time=time,event=cens)~x1, DATse)
    surv.prob <- summary(survfit(Surv(time=time,event=cens)~x1, DATse),time=2,extend=TRUE)$surv; 
    if (length(surv.prob)==0) {
      surv.prob<-c(NA,NA) }
    if (length(surv.prob)==1)  surv.prob<- c(surv.prob,surv.prob)
    if (length(surv.prob)==2)   surv.prob<- surv.prob
    
    pX<-length(DATse$x1[DATse$x1=="1"])/length(DATse$x1)
    
    T1<-table(DATse$x3,DATse$x1)
    SENS<-T1[2,2]/(T1[2,2]+T1[1,2]);SENS
    SPEC<-T1[1,1]/(T1[1,1]+T1[2,1]);SPEC
    
    VAR_d<-sample.est[3]^2
    DEFF<- VAR_srs/VAR_d
    
    PPS_EVsurr<-rbind(PPS_EVsurr,c(sample=i,length(DATse$x1),pX, surv.prob, sample.est, DEFF, SENS,SPEC) )
    colnames( PPS_EVsurr)<-c("sample", "n", "p_x1", "surv_x0", "surv_x1", "coef", "exp coef", "se.coef", "z", "pvalue", "lower95HR", "upper95HR","BConf","seConf", "DEFF", "sens","spec")
  
    ##########6. probability proportional to size (PPS) for 4 strata - CONFOUNDER
    
    T1<-table(dati$cens,dati$x4);
    cases0<-sample(dati$id[dati$cens==1 & dati$x4==0], T1[2,1]*n/N ,replace = FALSE, prob = NULL) 
    cntl0<-sample(dati$id[!dati$cens==1 & dati$x4==0], T1[1,1]*n/N ,replace = FALSE , prob = NULL)
    cases1<-sample(dati$id[dati$cens==1 & dati$x4==1], T1[2,2]*n/N,replace = FALSE , prob = NULL)
    cntl1<-sample(dati$id[!dati$cens==1 & dati$x4==1], T1[1,2]*n/N , replace = FALSE, prob = NULL)
    
    idsampled<-c(cases0,cases1,cntl0,cntl1)
    
    dati$incl<-ifelse(dati$id %in% idsampled,1,0);
    
    des<-twophase(id=list(~id,~id),subset=~incl==1,strata=list(NULL,~interaction(x4, cens)), data=dati)
    
    Fit<-svycoxph(Surv(time,cens)~x1 + x4, des)
    Summ<-summary(Fit)
    sample.est<-c(Summ$coefficients[1,],
                  Summ$conf.int[1,3:4],
                  Summ$coefficients[2,c(1,3)])
    
    DATse<-dati[dati$incl==1,]
    surv.prob <- summary(survfit(Surv(time=time,event=cens)~x1, DATse),time=2, extend=TRUE)$surv; 
    if (length(surv.prob)==0) {
      surv.prob<-c(NA,NA) }
    if (length(surv.prob)==1)  surv.prob<- c(surv.prob,surv.prob)
    if (length(surv.prob)==2)   surv.prob<- surv.prob
    
    pX<-length(DATse$x1[DATse$x1=="1"])/length(DATse$x1)
    
    T1<-table(DATse$x3,DATse$x1)
    SENS<-T1[2,2]/(T1[2,2]+T1[1,2]);SENS
    SPEC<-T1[1,1]/(T1[1,1]+T1[2,1]);SPEC
    
    VAR_d<-sample.est[3]^2
    DEFF<- VAR_srs/VAR_d
    
    PPS_EVconf<-rbind(PPS_EVconf,c(sample=i,length(DATse$x1),pX, surv.prob, sample.est, DEFF, SENS,SPEC) )
    colnames( PPS_EVconf)<-c("sample", "n", "p_x1", "surv_x0", "surv_x1", "coef", "exp coef", "se.coef", "z", "pvalue", "lower95HR", "upper95HR","BConf","seConf", "DEFF", "sens","spec")
    
    #########################################################################################################
    ######7. CASO - CONTROLLO CON INclUSIONE DI TUTTI I CASI 
  
    cases<-NULL; cntl<-NULL
    cases<-dati$id[dati$cens==1]
    if (length(cases)<n) {
      cntl<-sample(dati$id[!dati$cens==1],round(n-length(cases)), replace = FALSE, prob = NULL)
      }
    if (length(cases)>=n)  cntl<-0

    idsampled<-c(cases,cntl) ;
    dati$incl<-ifelse(dati$id %in% idsampled,1,0);table(dati$incl)
    
    DATse<-dati[dati$incl==1,]; 
    pX<-length(DATse$x1[DATse$x1=="1"])/length(DATse$x1)

    if (dim(table(DATse$cens))==1) {
      surv.prob<-c(NA,NA) 
      sample.est<-rep(NA,9)
    }   else 
 
      if (dim(table(DATse$cens))>1)  {
      surv.prob<- summary(survfit(Surv(time=time,event=cens)~x1, DATse),time=2,extend=TRUE)$surv
      des<-twophase(id=list(~id,~id),subset=~incl==1,strata=list(NULL,~cens),data=dati)
      
      Fit<-svycoxph(Surv(time,cens)~x1 + x4, des)
      Summ<-summary(Fit)
      sample.est<-c(Summ$coefficients[1,],
                    Summ$conf.int[1,3:4],
                    Summ$coefficients[2,c(1,3)])
      }
  
    T1<-table(DATse$x3,DATse$x1)
    SENS<-T1[2,2]/(T1[2,2]+T1[1,2]);SENS
    SPEC<-T1[1,1]/(T1[1,1]+T1[2,1]);SPEC
    
    VAR_d<-sample.est[3]^2
    DEFF<- VAR_srs/VAR_d

    CC_allCases<-rbind(CC_allCases,c(sample=i,length(DATse$x1),pX, surv.prob, sample.est,DEFF, SENS, SPEC) )
    colnames(CC_allCases)<-c("sample", "n", "p_x1", "surv_x0", "surv_x1", "coef", "exp coef", "se.coef", "z", 
                             "pvalue", "lower95HR", "upper95HR","BConf","seConf", "DEFF", "sens","spec")
 
    ###### Nested Case control 1:1 - two PHASE METHODS
    
    campione<-suppressWarnings(ccwc(entry=0,exit=time,fail=cens,controls=1,data=dati,include=id,silent=TRUE))
    if((length(campione$Set)/2)>=n/2){
      SetId<-sample(1:((length(campione$Set))/2), n/2, replace = FALSE, prob = NULL)
    }
    if ((length(campione$Set)/2)<n/2) SetId<-unique(campione$Set)
    
    IDCASE<-campione$id[campione$Set %in% SetId & campione$Fail=="1"]; #length(IDsampl)
    IDContr<-campione$id[campione$Set %in% SetId & campione$Fail=="0"]; #length(IDsampl)
    
    dati$sampleidstat<-ifelse(dati$id %in% IDCASE,2,ifelse(dati$id %in% IDContr,1,0)) 
    
    dati$pro<-KMprob(dati$time,dati$sampleidstat,m=1); 
    pCas<-length(IDCASE)/(length(dati$id[dati$cens=="1"])); #pCas
    dati$pro<-ifelse(dati$cens=="1", pCas,dati$pro)
    #table(dati$pro)
    dati$incl<-ifelse(dati$id %in% c(IDCASE,IDContr),1,0); #table(dati$incl)
    
   des<-twophase(id=list(~id,~id),subset=~incl==1,probs=list(NULL,~pro),data=dati)

    Fit<-svycoxph(Surv(time,cens)~x1 + x4, des)
    Summ<-summary(Fit)
    sample.est<-c(Summ$coefficients[1,],
                  Summ$conf.int[1,3:4],
                  Summ$coefficients[2,c(1,3)])
      
     DATse<-dati[dati$incl==1,]
    surv.prob <- summary(survfit(Surv(time=time,event=cens)~x1, DATse),time=2)$surv; 
    if (length(surv.prob)==0) {
    surv.prob<-c(NA,NA) }
    if (length(surv.prob)==1)  surv.prob<- c(surv.prob,surv.prob)
    if (length(surv.prob)==2)   surv.prob<- surv.prob
    
    pX<-length(DATse$x1[DATse$x1=="1"])/length(DATse$x1)

    T1<-table(DATse$x3,DATse$x1)
    SENS<-T1[2,2]/(T1[2,2]+T1[1,2]);SENS
    SPEC<-T1[1,1]/(T1[1,1]+T1[2,1]);SPEC
    
    VAR_d<-sample.est[3]^2    
    DEFF<- VAR_srs/VAR_d
    
    NCC_1<-rbind(NCC_1,c(sample=i,round(length(DATse$x1),1),pX, surv.prob, sample.est,DEFF, SENS,SPEC) )
    colnames(NCC_1)<-c("sample", "n", "p_x1", "survx0", "survx1", "coef", "exp coef", "se.coef", "z", "pvalue", "lower95HR", "upper95HR","BConf","seConf","DEFF","sens","spec")
    
    #####################################################
    ###### Nested Case control with matching SURR 1:1 - two PHASE METHODS

    campione<-suppressWarnings(ccwc(entry=0,exit=time,fail=cens,controls=1,data=dati,include=id,silent=TRUE,match=x3))
    if((length(campione$Set)/2)>=n/2){
      SetId<-sample(1:((length(campione$Set))/2), n/2, replace = FALSE, prob = NULL)
    }
    if ((length(campione$Set)/2)<n/2) SetId<-unique(campione$Set)
   
    IDCASE<-campione$id[campione$Set %in% SetId & campione$Fail=="1"]; #length(IDsampl)
    IDContr<-campione$id[campione$Set %in% SetId & campione$Fail=="0"]; #length(IDsampl)
    dati$incl<-ifelse(dati$id %in% c(IDCASE,IDContr),1,0); #table(dati$incl)
    dati$sampleidstat<-ifelse(dati$id %in% IDCASE,2,ifelse(dati$id %in% IDContr,1,0)) 
    dati$pro<-KMprob(dati$time,dati$sampleidstat,m=1,match.var=cbind(dati$x3),match.int=c(0,0)); 
     pCas<-length(IDCASE)/(length(dati$id[dati$cens=="1"])); #pCas
    dati$pro<-ifelse(dati$cens=="1", pCas,dati$pro)
    dati$incl<-ifelse(dati$id %in% c(IDCASE,IDContr),1,0); #table(dati$incl)
    des<-twophase(id=list(~id,~id),subset=~incl==1,probs = list (NULL, ~pro),data=dati)

    Fit<-svycoxph(Surv(time,cens)~x1 + x4, des)
    Summ<-summary(Fit)
    sample.est<- c(Summ$coefficients[1,],
                   Summ$conf.int[1,3:4],
                   Summ$coefficients[2,c(1,3)])
    
    DATse<-dati[dati$incl==1,]
    surv.prob <- summary(survfit(Surv(time=time,event=cens)~x1, DATse),time=2)$surv; 
    if (length(surv.prob)==0) {
      surv.prob<-c(NA,NA) }
    if (length(surv.prob)==1)  surv.prob<- c(surv.prob,surv.prob)
    if (length(surv.prob)==2)   surv.prob<- surv.prob
    
    #PROp of BM in selected two-phase dataset 
    pX<-length(DATse$x1[DATse$x1=="1"])/length(DATse$x1)
    
    T1<-table(DATse$x3,DATse$x1)
    SENS<-T1[2,2]/(T1[2,2]+T1[1,2]);SENS
    SPEC<-T1[1,1]/(T1[1,1]+T1[2,1]);SPEC
    
    VAR_d<-sample.est[3]^2    
    DEFF<- VAR_srs/VAR_d
    
    NCCsur<-rbind( NCCsur,c(sample=i,round(length(DATse$x1),1),pX, surv.prob, sample.est, DEFF, SENS,SPEC) )
    colnames( NCCsur)<-c("sample", "n", "p_x1", "survx0", "survx1", "coef", "exp coef", "se.coef", "z", "pvalue", "lower95HR", "upper95HR","BConf","seConf","DEFF","sens","spec")
  
    ###################
    ###COUNTER MATCHING

    dati<-dati[order(dati$time),] #ORDER Time
    dati$surr<-dati$x3+1 #Ricode Surr 1-2
    colstrat<-which( colnames(dati)=="surr")
    coltime<-which( colnames(dati)=="time")
    colcensor<-which( colnames(dati)=="cens")
    cmsample1<- function(data, m, colstrat) {
    data$strat<-data[, colstrat]
    iter<-nrow(data)
     mat<-numeric()# This matrix will save cases and control, it is done for using
    for(j  in (1:iter)[data$cens==1]){
    wj<-numeric()
    st_j<-as.vector(table(data$strat[j:iter])); L<-length(st_j); stat<-1:L  # Does not take into account strata names
    c_j<-data$strat[j];  dl<-numeric(); m_j<-m ;m_ja<-m_j; m_j[c_j]<- m_j[c_j]-1; w<- st_j/m_ja;  
    sj<-numeric(); staj<-numeric()
    for(l in 1:L){
      saux<-(((j+1):iter)[data$strat[(j+1):iter]==stat[l]  ]); si<- as.numeric(sample(as.character(saux), m_j[l])); sj<-c(sj,si); w_jl<-rep(st_j[l]/m_ja[l], nrow(data.frame(si)))#}; 
      staj<-c(staj,data$strat[si]); sl<-data[si,]; dl<-rbind(dl,sl); wj<-c(wj, w_jl)
    }
    mat<-rbind(mat, cbind(t(cbind(t(data[j,]), t(dl ) ) ) ,  c(st_j[c_j]/(m_ja[c_j] ), wj )))} 
  mat
     }
    
    ccount<-cmsample1(data=dati, m=c(1,1), colstrat) #c(1,1) un controllo per strato!
 
    rownames(ccount)<-1:nrow(ccount); 
    ccount<-data.frame(ccount);
    
    risize<-rep(0,sum(ccount$cens=="1")); 
    id<-1:nrow(ccount)
    
    risize<-rep(2, nrow(ccount)/2)
    ccount$riskset<-rep((1:sum(ccount$x3)), (risize));

    if((length(ccount$riskset)/2)>=n/2){
      SetId<-sample(1:((length(ccount$riskset))/2), n/2, replace = FALSE, prob = NULL)
    }
    if ((length(ccount$riskset)/2)<n/2) SetId<-unique(ccount$riskset)
    
    Prob<-function(m, data, coltime,colcensor, colstrat, t0=NULL,...){
    Pkt_f<-function(t, dat,iter, t0=NULL,m, colstrat){
  dat$strat<-dat[, colstrat]
  if(is.null(t0)){t0=rep(min(data$t), iter)} # If t0 is null, then it assumes all subjects enter the study at the same time
  p_t<-numeric(iter)
  if(dat$d[t]==1) {
    c_t<- dat$strat[t]; c_k<- dat$strat;  L<- length(unique(dat$strat));  n_t<-numeric(L);  len<- numeric(L)
    for(r in 1:L)   {n_t[r]<-sum((rep(1, iter-t+1))[dat$strat[t:iter]==r &  t0[t:iter]<=dat$t[t]  ]  )
    if( n_t[r]>0) len[r]<-0    else  len[r]<-1    } ;   m_t<-m
    
    if( n_t[c_t]>=m_t[c_t]) {p_t[c_t==c_k ]<- (m_t[c_t]-1)/(n_t[c_t]-1)}
    if(n_t[c_t]>=m_t[c_t]){ p_t[c_t==c_k]<- (m_t[c_t]-1)/(n_t[c_t]-1) }
    
    p_t[c_t==c_k & n_t[c_k]<m_t[c_k]& n_t[c_k]>1]<-1
    p_t[c_t==c_k & n_t[c_k]<=1]<-0
    
    a<-as.numeric(c_t!=c_k & n_t[c_k]>=m_t[c_k])
    b<-c_k[a==1]
    p_t[c_t!=c_k & n_t[c_k]>=m_t[c_k]]<- (m_t[b])/(n_t[b])
    p_t[c_t!=c_k & n_t[c_k]<m_t[c_k]& n_t[c_k]>=1]<-1
    p_t[c_t!=c_k & n_t[c_k]==0]<-0
    p_t[1:t]<- 0
  }
  p_t}

  data$strat<-data[, colstrat]
  iter<-dim(data)[1]; colnames(data)[coltime]<-'t'; colnames(data)[colcensor]<-'d';dat<-data[order(data$t),];
  strat<-data$strat;k<-1:iter;t<-1:iter; 
  if(is.null(t0)){t0=rep(min(data$t), iter)} 
  output = apply(as.matrix(1:iter,iter,1), 1,function(x) Pkt_f(x,dat,iter, t0, m,  colstrat))
  Pkt<-matrix(output,iter,iter); p<-numeric(iter)
  daux<-dat$d
  for(j  in 1:iter){
    if(dat$d[j]==1){  p[j]<-1}
    if(dat$d[j]==0){  taux<- min(data$t[t0[j]<= data$t]);bj<-(1:iter)[data$t==taux]
    pl<-((1-Pkt[j,])[bj:(j)]);  pl[j]<-1 ; p[j]<- 1-prod(pl)}}
  list(p=p )}
    
    dati$pro<-Prob(m=c(1,1), dati, coltime,colcensor, colstrat)$p
    dati$pro<-ifelse(dati$cens=="0", dati$pro, (n/2)/(length(dati$id[dati$cens==1])))

    IDsample<-ccount$id[ccount$riskset %in% SetId]; length(IDsample) 
    dati$incl<-ifelse(dati$id %in% IDsample,1,0); #table(dati$incl)
    
    des<-twophase(id=list(~id,~id),subset=~incl==1, probs = list (NULL, ~pro),data=dati)
    Fit<-svycoxph(Surv(time,cens)~x1 + x4, des)
    Summ<-summary(Fit)
    sample.est<- c(Summ$coefficients[1,],
                   Summ$conf.int[1,3:4],
                   Summ$coefficients[2,c(1,3)])
    
    DATse<-dati[dati$incl==1,]
    surv.prob <- summary(survfit(Surv(time=time,event=cens)~x1, DATse),time=2)$surv; 
    if (length(surv.prob)==0) {
      surv.prob<-c(NA,NA) }
    if (length(surv.prob)==1)  surv.prob<- c(surv.prob,surv.prob)
    if (length(surv.prob)==2)   surv.prob<- surv.prob
    
    #PROp of BM in selected two-phase dataset 
    pX<-length(DATse$x1[DATse$x1=="1"])/length(DATse$x1)
    
    T1<-table(DATse$x3,DATse$x1)
    SENS<-T1[2,2]/(T1[2,2]+T1[1,2]);SENS
    SPEC<-T1[1,1]/(T1[1,1]+T1[2,1]);SPEC
    
    VAR_d<-sample.est[3]^2    
    DEFF<- VAR_srs/VAR_d
    
    CM<-rbind(CM,c(sample=i,round(length(DATse$x1),1),pX, surv.prob, sample.est, DEFF, SENS,SPEC) )
    colnames(CM)<-c("sample", "n", "p_x1", "survx0", "survx1", "coef", "exp coef", "se.coef", "z", "pvalue", "lower95HR", "upper95HR","BConf","seConf","DEFF","sens","spec")
    
    ###################################
    print(i)
  }
  return(list(Full_Cohort=data.frame(Full_Cohort), Rsimple=data.frame(Rsimple),CC_ev=data.frame(CC_ev), CC_intEV_elfin= data.frame(CC_intEV_elfin),CC_intEV_conf=data.frame(CC_intEV_conf),
              CC_intEV_surr=data.frame(CC_intEV_surr),CC_stra_elfin=data.frame(CC_stra_elfin),CC_stra_conf=data.frame(CC_stra_conf),CC_stra_surr=data.frame(CC_stra_surr),
              PPS_EVe=data.frame(PPS_EVe), PPS_EVe_elf= data.frame(PPS_EVe_elf), PPS_EVsurr= data.frame(PPS_EVsurr),
              PPS_EVconf=data.frame(PPS_EVconf),CC_allCases=data.frame(CC_allCases),NCC_1=data.frame(NCC_1), NCCsur=data.frame(NCCsur), CM=data.frame(CM)))
}
