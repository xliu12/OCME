

AdjRe <- function( input.re, R.obs ){ # input.re: length of (2+dim.Z)
  input.re = as.numeric(as.matrix(input.re))
  r_yy = input.re[1];r_mm = input.re[2];r_zz = input.re[-c(1:2)];
  # if(length(input.re)==2){
  #   r_zz=NULL
  # }
  # if(length(input.re)>2){
  #   r_zz = input.re[-c(1:2)];
  # }
  
  r_xx=1
  D=diag( 1/sqrt(c(r_yy, r_mm, r_xx, r_zz)) )
  Rstar.obs = D%*%R.obs%*%D
  diag(Rstar.obs)=1
  colnames(Rstar.obs)=rownames(Rstar.obs)=colnames(R.obs)
  return(Rstar.obs)
}


AdjOC <- function(r_mc, r_yc, input.R,cols.Z1,cols.Z2,cols.Z3,
                  n, alpha=0.05){
  
  r_mc=as.numeric(r_mc); r_yc=as.numeric(r_yc)
  dim.Z1=length(cols.Z1); dim.Z2=length(cols.Z2); dim.Z3=length(cols.Z3)
  # for Eq2 
  R2=diag(1,(2+dim.Z1+dim.Z2)) # X C Z1 Z2 # sample corr(X,Z) not zero
  R2[1,2]=R2[2,1]=0 # r_xc = 0 by assumption
  R2[-2,-2]=input.R[ c(3,cols.Z1,cols.Z2),c(3,cols.Z1,cols.Z2) ]
  r2=rep(0, nrow(R2))
  r2[2]=r_mc
  r2[-2]=input.R[ 2,c(3,cols.Z1,cols.Z2) ]
  invR2=solve(R2)
  # estimation Eq2
  Beta2=invR2%*%r2
  Rsq2=t(Beta2)%*%r2
  df2=n-1-length(Beta2)
  hat.a=Beta2[1]; 
  F.a=( hat.a^2/invR2[1,1] ) / ( (1-Rsq2)/df2 )
  
  # for Eq3
  R3=diag(1,(3+dim.Z1+dim.Z3)) # M X C Z1 Z2
  R3[1,3]=R3[3,1]=r_mc
  R3[-3,-3]=input.R[ c(2,3,cols.Z1,cols.Z3 ),c(2,3,cols.Z1,cols.Z3 )]
  r3=rep(0,nrow(R3))
  r3[3]=r_yc
  r3[-3]=input.R[ 1, c(2,3,cols.Z1,cols.Z3 ) ]
  invR3=solve(R3)
  # estimation Eq3
  Beta3=invR3%*%r3
  Rsq3=t(Beta3)%*%r3
  df3=n-1-length(Beta3)
  hat.b=Beta3[1]
  F.b=( hat.b^2/invR3[1,1] ) / ( (1-Rsq3)/df3 )
  
  
  # solve r_mc based on Eq1 ----
  
  # hat.a not changed by r_mc: invR2[1,2] = 0
  # r_mc_zero = - invR2[1,-2]%*%r2[-2] / invR2[1,2] 
  
  # solve for r_mc that makes F.a = qF(1-alpha, 1, df2)
  #alpha=.05 # input: significance level
  qF2=qf(1-alpha, 1, df2)
  Q2 = qF2*invR2[1,1]/df2
  A2 = invR2[1,2]^2 + Q2*invR2[2,2]
  B2 = 2 * (invR2[1,2]*invR2[1,-2]%*%r2[-2] + Q2*invR2[2,-2]%*%r2[-2])
  C2 = (invR2[1,-2]%*%r2[-2])^2 + Q2*t(r2[-2])%*%invR2[-2,-2]%*%r2[-2] - Q2
  Delta2 = B2^2-4*A2*C2
  if(abs(Delta2)<1e-5) {Delta2=0}
  
  r_mc_s1 = (-B2+sqrt(Delta2))/(2*A2) # NaNs produced if Delta<0: is.nan()
  r_mc_s2 = (-B2-sqrt(Delta2))/(2*A2)
  r_mc_lower = min(r_mc_s1, r_mc_s2)
  r_mc_upper = max(r_mc_s1, r_mc_s2)
  
  # solve r_yc based on Eq3 ----
  # solve for r_yc that makes hat.b = 0
  r_yc_zero = - invR3[1,-3]%*%r3[-3] / invR3[1,3]
  
  # solve for r_yc that makes F.b = qF(1-alpha, 1, df3)
  #alpha=.05 # input: significance level
  qF3=qf(1-alpha, 1, df3)
  Q3 = qF3*invR3[1,1]/df3
  A3 = invR3[1,3]^2 + Q3*invR3[3,3]
  B3 = 2 * (invR3[1,3]*invR3[1,-3]%*%r3[-3] + Q3*invR3[3,-3]%*%r3[-3])
  C3 = (invR3[1,-3]%*%r3[-3])^2 + Q3*t(r3[-3])%*%invR3[-3,-3]%*%r3[-3] - Q3
  Delta3 = B3^2-4*A3*C3
  if(abs(Delta3)<1e-5) {Delta3=0}
  
  r_yc_s1 = (-B3+sqrt(Delta3))/(2*A3)
  r_yc_s2 = (-B3-sqrt(Delta3))/(2*A3)
  r_yc_lower = min(r_yc_s1, r_yc_s2)
  r_yc_upper = max(r_yc_s1, r_yc_s2)
  
  
  # possible range for r_yc, r_mc: ----
  # Eq2: r_mc
  Ro2=input.R[ c(3,cols.Z1,cols.Z2),c(3,cols.Z1,cols.Z2) ]
  invRo2=solve(Ro2)
  rmo2=input.R[ c(3,cols.Z1,cols.Z2),c(2) ]
  rco2=rep(0, length(rmo2))
  r_mc2_range_u= t(rco2)%*%invRo2%*%rmo2 + sqrt( (1-t(rco2)%*%invRo2%*%rco2)*(1-t(rmo2)%*%invRo2%*%rmo2) )
  r_mc2_range_l= t(rco2)%*%invRo2%*%rmo2 - sqrt( (1-t(rco2)%*%invRo2%*%rco2)*(1-t(rmo2)%*%invRo2%*%rmo2) )

  # Eq3:
  # r_mc
  Ro3m=input.R[c(1,3,cols.Z1,cols.Z3),c(1,3,cols.Z1,cols.Z3)]
  invRo3m=solve(Ro3m)
  rmo3=input.R[ c(1,3,cols.Z1,cols.Z3),c(2) ]
  rco3m=c(r_yc=0,rep(0,(length(rmo3)-1)) )
  r_mc3_range_u= t(rco3m)%*%invRo3m%*%rmo3 + sqrt( (1-t(rco3m)%*%invRo3m%*%rco3m)*(1-t(rmo3)%*%invRo3m%*%rmo3) )
  r_mc3_range_l= t(rco3m)%*%invRo3m%*%rmo3 - sqrt( (1-t(rco3m)%*%invRo3m%*%rco3m)*(1-t(rmo3)%*%invRo3m%*%rmo3) )
  
  # r_mc_range_l = max(r_mc3_range_l, r_mc2_range_l)
  # r_mc_range_u = min(r_mc3_range_u, r_mc2_range_u)
  r_mc_range_l = r_mc2_range_l
  r_mc_range_u = r_mc2_range_u

  # r_yc
  Ro3y=input.R[c(2,3,cols.Z1,cols.Z3),c(2,3,cols.Z1,cols.Z3)] 
  invRo3y=solve(Ro3y)
  ryo3=input.R[ c(2,3,cols.Z1,cols.Z3),c(1) ]
  rco3y=c(r_mc,rep(0,(length(ryo3)-1)) )
  r_yc3_range_u= t(rco3y)%*%invRo3y%*%ryo3 + sqrt( (1-t(rco3y)%*%invRo3y%*%rco3y)*(1-t(ryo3)%*%invRo3y%*%ryo3) )
  r_yc3_range_l= t(rco3y)%*%invRo3y%*%ryo3 - sqrt( (1-t(rco3y)%*%invRo3y%*%rco3y)*(1-t(ryo3)%*%invRo3y%*%ryo3) )
  
  r_yc_range_u=r_yc3_range_u
  r_yc_range_l=r_yc3_range_l
  
  # return the following for M2 (if input.R=R.M1=AdjRe( input.re = rep(1,(2+dim.Z)) )
  # or M4 (if input.R=R.M3=AdjRe( input.re = input.re ) )
  
  IfSig.a= (F.a>qF2)
  IfSig.b= (F.b>qF3)
  
  names_out=c("hat.a", "hat.b", "F.a", "F.b", "IfSig.a", "IfSig.b",
              "r_mc_lower","r_mc_upper","r_mc_range_l","r_mc_range_u",
              "r_yc_lower","r_yc_upper","r_yc_range_l","r_yc_range_u","r_yc_zero")
  out=unlist(mget(names_out))
  
  return(out)
}



AdjOCME <- function(R.obs, input.re, r_mc, r_yc, 
                    cols.Z1,cols.Z2,cols.Z3, n, alpha=0.05){
  # adjust for reliability 
  R.star=AdjRe( input.re, R.obs = R.obs ) # for Y M X Z1 Z2 Z3
  r_yy = as.numeric(input.re[1]) ;r_mm = as.numeric(input.re[2])
  rstar_yc=r_yc/sqrt(r_yy)
  rstar_mc=r_mc/sqrt(r_mm)
  # adjust for OC
  out1 = AdjOC(r_mc = rstar_mc, r_yc = rstar_yc, input.R=R.star,cols.Z1,cols.Z2,cols.Z3,n,alpha)
  hat.a.std=out1["hat.a"]*sqrt(r_mm)
  hat.b.std=out1["hat.b"]*sqrt(r_yy)/sqrt(r_mm)
  hat.ab.std=hat.a.std*hat.b.std
  
  out2 = unlist(c(r_mc, r_yc, r_yy, r_mm, 
                hat.a.std, hat.b.std,hat.ab.std, out1))
  out2[14:17]=
    sqrt(r_mm)*out2[14:17]
  out2[18:22]=
    sqrt(r_yy)*out2[18:22]

  names(out2)[1:7]=c("r_mc", "r_yc","r_yy","r_mm","hat.a.std", "hat.b.std","hat.ab.std")
  return(out2)
}











