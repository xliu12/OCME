
library(mvtnorm)
library(tidyverse)

#load("ocmeZ_sensi.RData")

sens.tab <- function(
  R.obs, 
  covariates.MYboth=NULL, 
  covariates.Monly=NULL, 
  covariates.Yonly=NULL, 
  n, 
  alpha=0.05,
  input.re=NULL, 
  r_mc_value=.2, 
  r_yc_value=.2
){
  source("fun_ocmeZ.R")
  n=n; alpha=alpha; R.obs=R.obs
  cols.Z1=covariates.MYboth
  cols.Z2=covariates.Monly
  cols.Z3=covariates.Yonly
  
  dim.Z1=ifelse(length(cols.Z1)==0, 0, length(cols.Z1)); 
  dim.Z2=ifelse(length(cols.Z2)==0, 0, length(cols.Z2)); 
  dim.Z3=ifelse(length(cols.Z3)==0, 0, length(cols.Z3)); 
  dim.Z=dim.Z1+dim.Z2+dim.Z3

  # out_m for sensitivity plot ----
  
  if(length(input.re)!=0){

    if(  length(input.re) %% (ncol(R.obs)-1) != 0  ){
      stop("input.re does not match R.obs: length(input.re) %% (ncol(R.obs)-1) != 0 ")
    }
    if( length(input.re) %% (2+dim.Z) != 0 ){
      stop("input.re does not match dimension of observed covariates: length(input.re) %% (2+dim.Z) != 0 ")
    }

    input.re=matrix(input.re, ncol = (2+dim.Z) )
  }

  if(length(input.re)==0){ # default reliability
    r_yy_seq=c(.7,1); r_mm_seq=c(.7,1)
    re.ym=expand.grid(r_yy=r_yy_seq, r_mm=r_mm_seq)
    re.z=matrix(1, nrow = nrow(re.ym), ncol = dim.Z )
    input.re=cbind(re.ym, re.z)
  }

  if(is.null(colnames(input.re))){
    if(dim.Z==0){
      colnames(input.re)=c("r_yy","r_mm")
    }
    if(dim.Z>0){
      colnames(input.re)=c("r_yy","r_mm",paste("re.ObsCovariate_",1:dim.Z,sep = ""))
    }
  }
  colnames(input.re)[1:2]=c("r_yy","r_mm")
  
  tab_m = NULL
  
  for(i in 1:nrow(input.re)){
    tab_m1=AdjOCME(R.obs=R.obs, 
                   input.re=input.re[i,], 
                   r_mc=r_mc_value, 
                   r_yc=r_yc_value, 
                   cols.Z1=covariates.MYboth,
                   cols.Z2=covariates.Monly,
                   cols.Z3=covariates.Yonly,
                   n=n, 
                   alpha=alpha)
    tab_m=rbind(unlist(tab_m1),tab_m)
  }
  
  df_m=as.data.frame(tab_m)
  sens_tab=cbind(r_mc_value, r_yc_value, input.re, 
                 df_m[,c("hat.a.std", "hat.b.std","hat.ab.std","F.a","F.b","IfSig.a", "IfSig.b")])
  
  return(sens_tab)
  
  
  
}




sens.plot <- function(
  R.obs, 
  covariates.MYboth=NULL, 
  covariates.Monly=NULL, 
  covariates.Yonly=NULL, 
  n, 
  alpha=0.05,
  input.re=NULL, 
  subplot.layout="",
  r_mc_range=c(0,1), 
  r_mc_value=.2, 
  r_yc_value=.2
){
  source("fun_ocmeZ.R")
  n=n; alpha=alpha; R.obs=R.obs
  cols.Z1=covariates.MYboth
  cols.Z2=covariates.Monly
  cols.Z3=covariates.Yonly
  
  dim.Z1=ifelse(length(cols.Z1)==0, 0, length(cols.Z1)); 
  dim.Z2=ifelse(length(cols.Z2)==0, 0, length(cols.Z2)); 
  dim.Z3=ifelse(length(cols.Z3)==0, 0, length(cols.Z3)); 
  dim.Z=dim.Z1+dim.Z2+dim.Z3
  
  if(length(input.re)!=0){
    
    if(  length(input.re) %% (ncol(R.obs)-1) != 0  ){
      stop("input.re does not match R.obs: length(input.re) %% (ncol(R.obs)-1) != 0 ")
    }
    if( length(input.re) %% (2+dim.Z) != 0 ){
      stop("input.re does not match dimension of observed covariates: length(input.re) %% (2+dim.Z) != 0 ")
    }
    
    input.re=as.matrix(input.re)
  }
  
  
  if(length(input.re)==0){ # default reliability
    r_yy_seq=c(.7,1); r_mm_seq=c(.7,1)
    re.ym=expand.grid(r_yy=r_yy_seq, r_mm=r_mm_seq)
    re.z=matrix(1, nrow = nrow(re.ym), ncol = dim.Z )
    input.re=cbind(re.ym, re.z)
    
    facetV=colnames(input.re)[which(apply(input.re, 2, function(x){length(unique(x))>1}))]
    facetZ=facetV[ !facetV%in%c("r_yy","r_mm") ]
    subplot.layout = paste(
      paste( c("r_yy",facetZ),collapse = "+" ),"~" ,
      paste( c("r_mm"),collapse = "+" ) 
      )
  }
  
  if(length(colnames(input.re)) == 0){
    if(dim.Z==0){
      colnames(input.re)=c("r_yy","r_mm")
    }
    if(dim.Z>0){
      colnames(input.re)=c("r_yy","r_mm",paste("re.ObsCovariate_",1:dim.Z,sep = ""))
    }
  }
  
  # colnames(input.re)[1:2]=c("r_yy","r_mm")
 
  
  if(subplot.layout==""){
    facetV=colnames(input.re)[which(apply(input.re, 2, function(x){length(unique(x))>1}))]
    # facetZ=facetV[ !facetV%in%c("r_yy","r_mm") ]
    subplot.layout=paste(paste( facetV,collapse = "+" ),"~ .")
  }
  
  r_mc_seq=seq(r_mc_range[1], r_mc_range[2], length.out = 1000)
  input_m=cbind(r_mc=r_mc_seq, r_yc=r_yc_value,
                input.re[rep(seq_len(nrow(input.re)), each=length(r_mc_seq)),])
  
  # out_m for sensitivity plot ----
  out_m = NULL
  
  for(i in 1:nrow(input_m)){
    out_m1=AdjOCME(R.obs=R.obs, 
                   input.re=input_m[i,-c(1,2)], 
                   r_mc=input_m[i,1], r_yc=input_m[i,2], 
                   cols.Z1=covariates.MYboth,cols.Z2=covariates.Monly,cols.Z3=covariates.Yonly,
                   n=n, alpha=alpha)
    out_m=rbind(unlist(out_m1),out_m)
  }
  out_mz=cbind(input_m[,-c(1:4)], out_m)
  
  # results from model M1
  R2.M1=R.obs[ c(3,cols.Z1,cols.Z2),c(3,cols.Z1,cols.Z2) ]
  r2.M1=as.numeric(R.obs[ 2,c(3,cols.Z1,cols.Z2) ])
  Beta2.M1=solve(R2.M1) %*% r2.M1
  Rsq2.M1=t(Beta2.M1)%*%r2.M1
  df2.M1=n-1-length(Beta2.M1)
  hat.a.M1=Beta2.M1[1]; 
  F.a.M1=( hat.a.M1^2/solve(R2.M1)[1,1] ) / ( (1-Rsq2.M1)/df2.M1 )
  
  R3.M1=R.obs[ c(2, 3,cols.Z1,cols.Z3),c(2,3,cols.Z1,cols.Z3) ]
  r3.M1=as.numeric(R.obs[ 1,c(2, 3,cols.Z1,cols.Z3) ])
  Beta3.M1=solve(R3.M1) %*% r3.M1
  Rsq3.M1=t(Beta3.M1)%*%r3.M1
  df3.M1=n-1-length(Beta3.M1)
  hat.b.M1=Beta3.M1[1]; 
  F.b.M1=( hat.b.M1^2/solve(R3.M1)[1,1] ) / ( (1-Rsq3.M1)/df3.M1 )
  
  IfSig.a.M1= (F.a.M1 > qf(1-alpha,1, df2.M1 ))
  IfSig.b.M1= (F.b.M1 > qf(1-alpha,1, df3.M1 ))
  
  # nonsignificant mediation from M1 ----
  
  if( IfSig.a.M1*IfSig.b.M1==0 ){
    
    out.tbl.s=as_tibble(out_mz) %>%
      filter(r_mc > r_mc_range_l & r_mc < r_mc_range_u) %>%
      filter(IfSig.a==1) %>% 
      filter( (r_yc_zero >= r_yc_range_l)&(r_yc_zero <= r_yc_range_u) ) 
    
    sens_plot = ggplot(  ) 
    
    if( nrow(filter(out.tbl.s, r_mc>0))>0 ){  # r_mc>0
      sens_plot = sens_plot +
        geom_ribbon( data = filter(out.tbl.s
                                   , r_mc>0
                                   , IfSig.a==1
                                   ,(r_yc_lower >= r_yc_range_l)&(r_yc_upper <= r_yc_range_u)
                                   ),
                     aes( x = r_mc, ymin = r_yc_range_l, ymax = r_yc_lower
                          , fill="Mediation becomes significant")
                     , alpha=0.8) +
        geom_ribbon( data = filter(out.tbl.s
                                   , r_mc>0
                                   , IfSig.a==1
                                   ,(r_yc_lower >= r_yc_range_l)&(r_yc_upper <= r_yc_range_u)
        ),
                     aes( x = r_mc, ymin = r_yc_upper, ymax = r_yc_range_u
                          , fill="Mediation becomes significant" )
                     , alpha=0.8) +
        geom_path( data = filter(out.tbl.s
                                 , r_mc>0
                                 # , IfSig.a==1
                                 ,(r_yc_zero >= r_yc_range_l)&(r_yc_zero <= r_yc_range_u)
        ),
                   aes( x = r_mc, y = r_yc_zero
                        , linetype="Point estimate of path b equals 0"
                   ), size=1) +
        
        geom_path( data = filter(out.tbl.s
                                 , r_mc>0
                                 # , IfSig.a==1
                                 ,(r_yc_lower >= r_yc_range_l)&(r_yc_lower <= r_yc_range_u)
        ),
                   aes( x = r_mc, y = r_yc_lower
                        , linetype="F-statistic of testing path b equals the critical value"
                   ), size=0.5) +
        geom_path( data = filter(out.tbl.s
                                 , r_mc>0
                                 # , IfSig.a==1
                                 ,(r_yc_upper >= r_yc_range_l)&(r_yc_upper <= r_yc_range_u)
        ),
                   aes( x = r_mc, y = r_yc_upper
                        , linetype="F-statistic of testing path b equals the critical value"
                   ), size=0.5) +
        
        geom_path( data = filter(out.tbl.s, r_mc>0),
                   aes( x = r_mc, y = r_yc_range_l
                        , linetype="Range limit of r_yc given r_mc"
                   ), color="black", size=.5) +
        geom_path( data = filter(out.tbl.s, r_mc>0),
                   aes( x = r_mc, y = r_yc_range_u
                        , linetype="Range limit of r_yc given r_mc"
                   ), color="black", size=.5) 
      
    } 
    
    if( nrow(filter(out.tbl.s, r_mc<0))>0 ){
        # r_mc<=0
      sens_plot = sens_plot +
        geom_ribbon(data = filter(out.tbl.s
                                  , r_mc<=0
                                  , IfSig.a==1
                                  ,(r_yc_lower >= r_yc_range_l)&(r_yc_upper <= r_yc_range_u)
        ),
                    aes( x = r_mc, ymin = r_yc_range_l, ymax = r_yc_lower
                         , fill="Mediation becomes significant")
                    , alpha=0.8) +
        geom_ribbon( data = filter(out.tbl.s
                                   , r_mc<=0
                                   , IfSig.a==1
                                   ,(r_yc_lower >= r_yc_range_l)&(r_yc_upper <= r_yc_range_u)
        ),
                     aes( x = r_mc, ymin = r_yc_upper, ymax = r_yc_range_u
                          , fill="Mediation becomes significant" )
                     , alpha=0.8) +
        geom_path( data = filter(out.tbl.s
                                 , r_mc<=0
                                 # , IfSig.a==1
                                 ,(r_yc_zero >= r_yc_range_l)&(r_yc_zero <= r_yc_range_u)
        ),
                   aes( x = r_mc, y = r_yc_zero
                        , linetype="Point estimate of path b equals 0"
                   ), size=1) +
        
        geom_path( data = filter(out.tbl.s
                                 , r_mc<=0
                                 # , IfSig.a==1
                                 ,(r_yc_lower >= r_yc_range_l)&(r_yc_lower <= r_yc_range_u)
        ),
                   aes( x = r_mc, y = r_yc_lower
                        , linetype="F-statistic of testing path b equals the critical value"
                   ), size=0.5) +
        geom_path(data = filter(out.tbl.s
                                , r_mc<=0
                                # , IfSig.a==1
                                ,(r_yc_upper >= r_yc_range_l)&(r_yc_upper <= r_yc_range_u)
        ),
                  aes( x = r_mc, y = r_yc_upper
                       , linetype="F-statistic of testing path b equals the critical value"
                  ), size=0.5) +
        
        geom_path(data = filter(out.tbl.s, r_mc<=0),
                  aes( x = r_mc, y = r_yc_range_l
                       , linetype="Range limit of r_yc given r_mc"
                  ), color="black", size=.5) +
        geom_path(data = filter(out.tbl.s, r_mc<=0),
                  aes( x = r_mc, y = r_yc_range_u
                       , linetype="Range limit of r_yc given r_mc"
                  ), color="black", size=.5)
    }
    
    sens_plot = sens_plot +
      # theme 
      scale_fill_manual(guide = "legend",name="After accounting for confounding and measurement error",
                        values = c("Mediation becomes significant"="green4")) +
      scale_linetype_manual(name="After accounting for confounding and measurement error"
                            ,guide = "legend"
                            ,values=c("Point estimate of path b equals 0"=2
                                      ,"F-statistic of testing path b equals the critical value"=1
                                      ,"Range limit of r_yc given r_mc"=3)
      ) +
      facet_grid( as.formula( subplot.layout ), labeller = label_both ) +
      labs( y = "r_yc" ) +
      xlim(r_mc_range[1],r_mc_range[2]) +
      ylim(-1,1) + 
      theme_bw() + 
      theme( axis.title = element_text( size = 15 )
             , legend.text = element_text( size = 15 )
             , legend.title = element_text( size = 15 )
             , legend.position = "bottom"
             , legend.box = "vertical"
             , legend.direction = "vertical"
             , legend.key.width = unit(1.5,"cm")
             , legend.key = element_rect(fill = NULL)
             , strip.text = element_text(size = 15)) 
    
  }

  # significant mediation from M1 ----
  
  if( IfSig.a.M1*IfSig.b.M1!=0 ){
    out.tbl.s=as_tibble(out_mz) %>%
      # filter(r_mc > r_mc_range_l & r_mc < r_mc_range_u) %>%
      filter(IfSig.a==1) %>%
      filter( (r_yc_zero >= r_yc_range_l)&(r_yc_zero <= r_yc_range_u) ) 
    out.tbl.s=as_tibble(out_mz)
    
    sens_plot = ggplot(  ) 
    if( nrow(filter(out.tbl.s, r_mc>0))>0 ){
      # r_mc>0
      sens_plot = sens_plot +
        geom_ribbon(data = filter(out.tbl.s
                                  , r_mc>0
                                  , IfSig.a==1
                                  , (r_yc_lower >= r_yc_range_l)&(r_yc_upper <= r_yc_range_u) 
                                  ) ,
                    aes( x = r_mc, ymin = r_yc_lower, ymax = r_yc_upper
                         , fill="Mediation becomes non-significant (point estimate of path b in the same direction)" 
                    )
                    , alpha=1) +
        geom_ribbon(data = filter(out.tbl.s
                                  , r_mc>0
                                  , IfSig.a==1
                                  , (r_yc_zero >= r_yc_range_l)&(r_yc_zero <= r_yc_range_u) 
                    ) ,
                    aes( x = r_mc, 
                         ymin = (r_yc_zero)*( abs(r_yc_upper*r_mc)>abs(r_yc_zero*r_mc) )+
                           (r_yc_lower)*( abs(r_yc_upper*r_mc)<=abs(r_yc_zero*r_mc) ), 
                         ymax = (r_yc_upper)*( abs(r_yc_upper*r_mc)>abs(r_yc_zero*r_mc) )+
                           (r_yc_zero)*( abs(r_yc_upper*r_mc)<=abs(r_yc_zero*r_mc) )
                         , fill="Mediation becomes non-significant (point estimate of path b in the opposite direction)" )
                    , alpha=1) +
        geom_ribbon(data = filter(out.tbl.s
                                  , r_mc>0
                                  , IfSig.a==1
                                  , (r_yc_zero >= r_yc_range_l)&(r_yc_zero <= r_yc_range_u) 
                    ) ,
                    aes( x = r_mc, 
                         ymin = (r_yc_upper)*( abs(r_yc_upper*r_mc)>abs(r_yc_zero*r_mc) )+
                           (r_yc_range_l)*( abs(r_yc_upper*r_mc)<=abs(r_yc_zero*r_mc) ), 
                         ymax = (r_yc_range_u)*( abs(r_yc_upper*r_mc)>abs(r_yc_zero*r_mc) )+
                           (r_yc_lower)*( abs(r_yc_upper*r_mc)<=abs(r_yc_zero*r_mc) )
                         , fill="Mediation remains significant but point estimate of path b in the opposite direction" )
                    , alpha=1) +
        
        geom_path( data = filter(out.tbl.s
                                 , r_mc>0
                                 , (r_yc_zero >= r_yc_range_l)&(r_yc_zero <= r_yc_range_u) 
                   ),
                   aes( x = r_mc, y = r_yc_zero
                        , linetype="Point estimate of path b equals 0"
                   ), size=1) +
        
        geom_path( data = filter(out.tbl.s
                                 , r_mc>0
                                 , (r_yc_lower >= r_yc_range_l)&(r_yc_lower <= r_yc_range_u) 
        ),
                   aes( x = r_mc, y = r_yc_lower
                        , linetype="F-statistic of testing path b equals the critical value"
                   ), size=0.5) +
        geom_path( data = filter(out.tbl.s
                                 , r_mc>0
                                 , (r_yc_upper >= r_yc_range_l)&(r_yc_upper <= r_yc_range_u) 
        ),
                   aes( x = r_mc, y = r_yc_upper
                        , linetype="F-statistic of testing path b equals the critical value"
                   ), size=0.5) +
        
        geom_path( data = filter(out.tbl.s, r_mc>0),
                   aes( x = r_mc, y = r_yc_range_l
                        , linetype="Range limit of r_yc given r_mc"
                   ), color="black", size=.5) +
        geom_path( data = filter(out.tbl.s, r_mc>0),
                   aes( x = r_mc, y = r_yc_range_u
                        , linetype="Range limit of r_yc given r_mc"
                   ), color="black", size=.5) 
    }
    
    if( nrow(filter(out.tbl.s, r_mc < 0))>0 ){
      sens_plot = sens_plot +
        # r_mc < 0
        geom_ribbon( data = filter(out.tbl.s
                                   , r_mc <= 0
                                   , IfSig.a==1
                                   , (r_yc_lower >= r_yc_range_l)&(r_yc_upper <= r_yc_range_u) 
        ),
                     aes( x = r_mc, ymin = r_yc_lower, ymax = r_yc_upper
                          , fill="Mediation becomes non-significant (point estimate of path b in the same direction)" 
                     )
                     , alpha=1) +
        geom_ribbon( data = filter(out.tbl.s
                                   , r_mc <= 0
                                   , IfSig.a==1
                                   , (r_yc_zero >= r_yc_range_l)&(r_yc_zero <= r_yc_range_u) 
        ),
                     aes( x = r_mc, 
                          ymin = (r_yc_zero)*( abs(r_yc_upper*r_mc)>abs(r_yc_zero*r_mc) )+
                            (r_yc_lower)*( abs(r_yc_upper*r_mc)<=abs(r_yc_zero*r_mc) ), 
                          ymax = (r_yc_upper)*( abs(r_yc_upper*r_mc)>abs(r_yc_zero*r_mc) )+
                            (r_yc_zero)*( abs(r_yc_upper*r_mc)<=abs(r_yc_zero*r_mc) )
                          , fill="Mediation becomes non-significant (point estimate of path b in the opposite direction)" )
                     , alpha=1) +
        geom_ribbon( data = filter(out.tbl.s
                                   , r_mc <= 0
                                   , IfSig.a==1
                                   , (r_yc_zero >= r_yc_range_l)&(r_yc_zero <= r_yc_range_u) 
        ),
                     aes( x = r_mc, 
                          ymin = (r_yc_upper)*( abs(r_yc_upper*r_mc)>abs(r_yc_zero*r_mc) )+
                            (r_yc_range_l)*( abs(r_yc_upper*r_mc)<=abs(r_yc_zero*r_mc) ), 
                          ymax = (r_yc_range_u)*( abs(r_yc_upper*r_mc)>abs(r_yc_zero*r_mc) )+
                            (r_yc_lower)*( abs(r_yc_upper*r_mc)<=abs(r_yc_zero*r_mc) )
                          , fill="Mediation remains significant but point estimate of path b in the opposite direction" )
                     , alpha=1) +
        
        geom_path( data = filter(out.tbl.s
                                 , r_mc <= 0
                                 , (r_yc_zero >= r_yc_range_l)&(r_yc_zero <= r_yc_range_u) 
        ),
                   aes( x = r_mc, y = r_yc_zero
                        , linetype="Point estimate of path b equals 0"
                   ), size=1) +
        
        geom_path( data = filter(out.tbl.s
                                 , r_mc <= 0
                                 , (r_yc_lower >= r_yc_range_l)&(r_yc_lower <= r_yc_range_u) 
        ),
                   aes( x = r_mc, y = r_yc_lower
                        , linetype="F-statistic of testing path b equals the critical value"
                   ), size=0.5) +
        geom_path( data = filter(out.tbl.s
                                 , r_mc <= 0
                                 , (r_yc_upper >= r_yc_range_l)&(r_yc_upper <= r_yc_range_u) 
        ),
                   aes( x = r_mc, y = r_yc_upper
                        , linetype="F-statistic of testing path b equals the critical value"
                   ), size=0.5) +
        
        geom_path( data = filter(out.tbl.s, r_mc <= 0),
                   aes( x = r_mc, y = r_yc_range_l
                        , linetype="Range limit of r_yc given r_mc"
                   ), color="black", size=.5) +
        geom_path( data = filter(out.tbl.s, r_mc <= 0),
                   aes( x = r_mc, y = r_yc_range_u
                        , linetype="Range limit of r_yc given r_mc"
                   ), color="black", size=.5)
    }
    
    sens_plot = sens_plot +
      # theme
      scale_fill_manual(name="After accounting for confounding and measurement error",
                        guide = "legend"
                        ,values =c("Mediation becomes non-significant (point estimate of path b in the same direction)"="red"
                                   ,"Mediation remains significant but point estimate of path b in the opposite direction"="yellow2"
                                   ,"Mediation becomes non-significant (point estimate of path b in the opposite direction)"="orange") 
      ) +
      scale_linetype_manual(name="After accounting for confounding and measurement error"
                            ,guide = "legend"
                            ,values=c("Point estimate of path b equals 0"=2
                                      ,"F-statistic of testing path b equals the critical value"=1
                                      ,"Range limit of r_yc given r_mc"=3)
      ) +
      facet_grid( formula(subplot.layout), labeller = label_both ) +
      labs( y = "r_yc" ) +
      xlim(r_mc_range[1],r_mc_range[2]) +
      ylim(-1,1) +
      theme_bw() + 
      theme( axis.title = element_text( size = 15 )
             , legend.text = element_text( size = 15 )
             , legend.title = element_text( size = 15 )
             , legend.position = "bottom"
             , legend.box = "vertical"
             , legend.direction = "vertical"
             , legend.key.width = unit(1.5,"cm")
             , legend.key = element_rect(fill = NULL)
             , strip.text = element_text(size = 15)) 
    
  }
  
  print(sens_plot)
  
  return(sens_plot)
  
  
  
}




