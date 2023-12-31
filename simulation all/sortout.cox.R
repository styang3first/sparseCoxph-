load("E:/Google ���ݵw��/program/yang/group penalized regression 2013/group cox/�T�wweightver2/simulation all/OK/allresult.RData")
library(survival)
# AGLASSO0.sel.400[,1]$sel1$info$penalty
# AGLASSO0.sel.400[,1]$sel1$info$gamma
# AGLASSO0.sel.400[,1]$sel2$info$penalty
# AGLASSO0.sel.400[,1]$sel2$info$gamma
# AGLASSO0.sel.800[,1]$sel1$info$penalty
# AGLASSO0.sel.800[,1]$sel1$info$gamma
# AGLASSO0.sel.800[,1]$sel2$info$penalty
# AGLASSO0.sel.800[,1]$sel2$info$gamma
# AGLASSO0.sel.1200[,1]$sel1$info$penalty
# AGLASSO0.sel.1200[,1]$sel1$info$gamma
# AGLASSO0.sel.1200[,1]$sel2$info$penalty
# AGLASSO0.sel.1200[,1]$sel2$info$gamma
# 
# AGLASSO1.sel.400[,1]$sel1$info$penalty
# AGLASSO1.sel.400[,1]$sel1$info$gamma
# AGLASSO1.sel.400[,1]$sel2$info$penalty
# AGLASSO1.sel.400[,1]$sel2$info$gamma
# AGLASSO1.sel.800[,1]$sel1$info$penalty
# AGLASSO1.sel.800[,1]$sel1$info$gamma
# AGLASSO1.sel.800[,1]$sel2$info$penalty
# AGLASSO1.sel.800[,1]$sel2$info$gamma
# AGLASSO1.sel.1200[,1]$sel1$info$penalty
# AGLASSO1.sel.1200[,1]$sel1$info$gamma
# AGLASSO1.sel.1200[,1]$sel2$info$penalty
# AGLASSO1.sel.1200[,1]$sel2$info$gamma
# 
# GSCAD.sel.400[,1]$sel1$info$penalty
# GSCAD.sel.400[,1]$sel1$info$gamma
# GSCAD.sel.400[,1]$sel2$info$penalty
# GSCAD.sel.400[,1]$sel2$info$gamma
# GSCAD.sel.800[,1]$sel1$info$penalty
# GSCAD.sel.800[,1]$sel1$info$gamma
# GSCAD.sel.800[,1]$sel2$info$penalty
# GSCAD.sel.800[,1]$sel2$info$gamma
# GSCAD.sel.1200[,1]$sel1$info$penalty
# GSCAD.sel.1200[,1]$sel1$info$gamma
# GSCAD.sel.1200[,1]$sel2$info$penalty
# GSCAD.sel.1200[,1]$sel2$info$gamma
# 
# GMCP.sel.400[,1]$sel1$info$penalty
# GMCP.sel.400[,1]$sel1$info$gamma
# GMCP.sel.400[,1]$sel2$info$penalty
# GMCP.sel.400[,1]$sel2$info$gamma
# GMCP.sel.800[,1]$sel1$info$penalty
# GMCP.sel.800[,1]$sel1$info$gamma
# GMCP.sel.800[,1]$sel2$info$penalty
# GMCP.sel.800[,1]$sel2$info$gamma
# GMCP.sel.1200[,1]$sel1$info$penalty
# GMCP.sel.1200[,1]$sel1$info$gamma
# GMCP.sel.1200[,1]$sel2$info$penalty
# GMCP.sel.1200[,1]$sel2$info$gamma


posi = c(1, 2, 5, 6, 11);
beta.true = rep(0, 15)
beta.true[posi] = c(3, 3, -3.5, -3.5, 2)/10

####### MLE, MLE_true
### 400
cox1.400 = apply(data400, 2, function(x) coxph(Surv(y[,2], y[,3])~x.m, data=x)$coef)
cox2.400 = apply(data400, 2, function(x) coxph(Surv(y[,1], y[,2], y[,3])~x.m, data=x)$coef)
cox1.t.400 = apply(data400, 2, function(x) coxph(Surv(y[,2], y[,3])~x.m[,posi], data=x)$coef)
cox2.t.400 = apply(data400, 2, function(x) coxph(Surv(y[,1], y[,2], y[,3])~x.m[,posi], data=x)$coef)
ASE.cox1.t.400 = apply(data400, 2, function(x) sqrt(diag(coxph(Surv(y[,2], y[,3])~x.m[,posi], data=x)$var)))
ASE.cox2.t.400 = apply(data400, 2, function(x) sqrt(diag(coxph(Surv(y[,1], y[,2], y[,3])~x.m[,posi], data=x)$var)))
### 800
cox1.800 = apply(data800, 2, function(x) coxph(Surv(y[,2], y[,3])~x.m, data=x)$coef)
cox2.800 = apply(data800, 2, function(x) coxph(Surv(y[,1], y[,2], y[,3])~x.m, data=x)$coef)
cox1.t.800 = apply(data800, 2, function(x) coxph(Surv(y[,2], y[,3])~x.m[,posi], data=x)$coef)
cox2.t.800 = apply(data800, 2, function(x) coxph(Surv(y[,1], y[,2], y[,3])~x.m[,posi], data=x)$coef)
ASE.cox1.t.800 = apply(data800, 2, function(x) sqrt(diag(coxph(Surv(y[,2], y[,3])~x.m[,posi], data=x)$var)))
ASE.cox2.t.800 = apply(data800, 2, function(x) sqrt(diag(coxph(Surv(y[,1], y[,2], y[,3])~x.m[,posi], data=x)$var)))
### 1200
cox1.1200 = apply(data1200, 2, function(x) coxph(Surv(y[,2], y[,3])~x.m, data=x)$coef)
cox2.1200 = apply(data1200, 2, function(x) coxph(Surv(y[,1], y[,2], y[,3])~x.m, data=x)$coef)
cox1.t.1200 = apply(data1200, 2, function(x) coxph(Surv(y[,2], y[,3])~x.m[,posi], data=x)$coef)
cox2.t.1200 = apply(data1200, 2, function(x) coxph(Surv(y[,1], y[,2], y[,3])~x.m[,posi], data=x)$coef)
ASE.cox1.t.1200 = apply(data1200, 2, function(x) sqrt(diag(coxph(Surv(y[,2], y[,3])~x.m[,posi], data=x)$var)))
ASE.cox2.t.1200 = apply(data1200, 2, function(x) sqrt(diag(coxph(Surv(y[,1], y[,2], y[,3])~x.m[,posi], data=x)$var)))
#### MSE cox
{
  MSE.mle1.400 = sapply(1:100, function(i){
    V = cov(data400[,i]$x.m)
    beta = cox1.400[,i]
    output = t(beta-beta.true)%*%V%*%(beta-beta.true)
    return(output)
  })
  
  MSE.mle2.400 = sapply(1:100, function(i){
    V = cov(data400[,i]$x.m)
    beta = cox2.400[,i]
    output = t(beta-beta.true)%*%V%*%(beta-beta.true)
    return(output)
  })
  
  MSE.mle1.800 = sapply(1:100, function(i){
    V = cov(data800[,i]$x.m)
    beta = cox1.800[,i]
    output = t(beta-beta.true)%*%V%*%(beta-beta.true)
    return(output)
  })
  
  MSE.mle2.800 = sapply(1:100, function(i){
    V = cov(data1200[,i]$x.m)
    beta = cox2.800[,i]
    output = t(beta-beta.true)%*%V%*%(beta-beta.true)
    return(output)
  })
  
  MSE.mle1.1200 = sapply(1:100, function(i){
    V = cov(data1200[,i]$x.m)
    beta = cox1.1200[,i]
    output = t(beta-beta.true)%*%V%*%(beta-beta.true)
    return(output)
  })
  
  MSE.mle2.1200 = sapply(1:100, function(i){
    V = cov(data1200[,i]$x.m)
    beta = cox2.1200[,i]
    output = t(beta-beta.true)%*%V%*%(beta-beta.true)
    return(output)
  })
}

######  beta1.record
### 5-fold CV
{
  beta1.CV.400 = list(
    GLASSO = apply(AGLASSO0.CV.400, 2, function(x) beta.cox(x$CV1, "CV")),
    GSCAD = apply(GSCAD.CV.400, 2, function(x) beta.cox(x$CV1, "CV")),
    GMCP = apply(GMCP.CV.400, 2, function(x) beta.cox(x$CV1, "CV")),
    AGLASSO = apply(AGLASSO1.CV.400, 2, function(x) beta.cox(x$CV1, "CV"))
  )  
  beta1.CV.800 = list(
    GLASSO = apply(AGLASSO0.CV.800, 2, function(x) beta.cox(x$CV1, "CV")),
    GSCAD = apply(GSCAD.CV.800, 2, function(x) beta.cox(x$CV1, "CV")),
    GMCP = apply(GMCP.CV.800, 2, function(x) beta.cox(x$CV1, "CV")),
    AGLASSO = apply(AGLASSO1.CV.800, 2, function(x) beta.cox(x$CV1, "CV"))
  )
  beta1.CV.1200 = list(
    GLASSO = apply(AGLASSO0.CV.1200, 2, function(x) beta.cox(x$CV1, "CV")),
    GSCAD = apply(GSCAD.CV.1200, 2, function(x) beta.cox(x$CV1, "CV")),
    GMCP = apply(GMCP.CV.1200, 2, function(x) beta.cox(x$CV1, "CV")),
    AGLASSO = apply(AGLASSO1.CV.1200, 2, function(x) beta.cox(x$CV1, "CV"))
  )
}
### GCV
{
  beta1.GCV.400 = list(
    GLASSO = apply(AGLASSO0.sel.400, 2, function(x) beta.cox(x$sel1, "GCV")),
    GSCAD = apply(GSCAD.sel.400, 2, function(x) beta.cox(x$sel1, "GCV")),
    GMCP = apply(GMCP.sel.400, 2, function(x) beta.cox(x$sel1, "GCV")),
    AGLASSO = apply(AGLASSO1.sel.400, 2, function(x) beta.cox(x$sel1, "GCV"))
  )  
  beta1.GCV.800 = list(
    GLASSO = apply(AGLASSO0.sel.800, 2, function(x) beta.cox(x$sel1, "GCV")),
    GSCAD = apply(GSCAD.sel.800, 2, function(x) beta.cox(x$sel1, "GCV")),
    GMCP = apply(GMCP.sel.800, 2, function(x) beta.cox(x$sel1, "GCV")),
    AGLASSO = apply(AGLASSO1.sel.800, 2, function(x) beta.cox(x$sel1, "GCV"))
  )
  beta1.GCV.1200 = list(
    GLASSO = apply(AGLASSO0.sel.1200, 2, function(x) beta.cox(x$sel1, "GCV")),
    GSCAD = apply(GSCAD.sel.1200, 2, function(x) beta.cox(x$sel1, "GCV")),
    GMCP = apply(GMCP.sel.1200, 2, function(x) beta.cox(x$sel1, "GCV")),
    AGLASSO = apply(AGLASSO1.sel.1200, 2, function(x) beta.cox(x$sel1, "GCV"))
  )
}
### BIC
{
  beta1.BIC.400 = list(
    GLASSO = apply(AGLASSO0.sel.400, 2, function(x) beta.cox(x$sel1, "BIC")),
    GSCAD = apply(GSCAD.sel.400, 2, function(x) beta.cox(x$sel1, "BIC")),
    GMCP = apply(GMCP.sel.400, 2, function(x) beta.cox(x$sel1, "BIC")),
    AGLASSO = apply(AGLASSO1.sel.400, 2, function(x) beta.cox(x$sel1, "BIC"))
  )  
  beta1.BIC.800 = list(
    GLASSO = apply(AGLASSO0.sel.800, 2, function(x) beta.cox(x$sel1, "BIC")),
    GSCAD = apply(GSCAD.sel.800, 2, function(x) beta.cox(x$sel1, "BIC")),
    GMCP = apply(GMCP.sel.800, 2, function(x) beta.cox(x$sel1, "BIC")),
    AGLASSO = apply(AGLASSO1.sel.800, 2, function(x) beta.cox(x$sel1, "BIC"))
  )
  beta1.BIC.1200 = list(
    GLASSO = apply(AGLASSO0.sel.1200, 2, function(x) beta.cox(x$sel1, "BIC")),
    GSCAD = apply(GSCAD.sel.1200, 2, function(x) beta.cox(x$sel1, "BIC")),
    GMCP = apply(GMCP.sel.1200, 2, function(x) beta.cox(x$sel1, "BIC")),
    AGLASSO = apply(AGLASSO1.sel.1200, 2, function(x) beta.cox(x$sel1, "BIC"))
  )
}
### AIC
{
  beta1.AIC.400 = list(
    GLASSO = apply(AGLASSO0.sel.400, 2, function(x) beta.cox(x$sel1, "AIC")),
    GSCAD = apply(GSCAD.sel.400, 2, function(x) beta.cox(x$sel1, "AIC")),
    GMCP = apply(GMCP.sel.400, 2, function(x) beta.cox(x$sel1, "AIC")),
    AGLASSO = apply(AGLASSO1.sel.400, 2, function(x) beta.cox(x$sel1, "AIC"))
  )  
  beta1.AIC.800 = list(
    GLASSO = apply(AGLASSO0.sel.800, 2, function(x) beta.cox(x$sel1, "AIC")),
    GSCAD = apply(GSCAD.sel.800, 2, function(x) beta.cox(x$sel1, "AIC")),
    GMCP = apply(GMCP.sel.800, 2, function(x) beta.cox(x$sel1, "AIC")),
    AGLASSO = apply(AGLASSO1.sel.800, 2, function(x) beta.cox(x$sel1, "AIC"))
  )
  beta1.AIC.1200 = list(
    GLASSO = apply(AGLASSO0.sel.1200, 2, function(x) beta.cox(x$sel1, "AIC")),
    GSCAD = apply(GSCAD.sel.1200, 2, function(x) beta.cox(x$sel1, "AIC")),
    GMCP = apply(GMCP.sel.1200, 2, function(x) beta.cox(x$sel1, "AIC")),
    AGLASSO = apply(AGLASSO1.sel.1200, 2, function(x) beta.cox(x$sel1, "AIC"))
  )
}

######  beta2.coxrecord
### 5-fold CV
{
  beta2.CV.400 = list(
    GLASSO = apply(AGLASSO0.CV.400, 2, function(x) beta.cox(x$CV2, "CV")),
    GSCAD = apply(GSCAD.CV.400, 2, function(x) beta.cox(x$CV2, "CV")),
    GMCP = apply(GMCP.CV.400, 2, function(x) beta.cox(x$CV2, "CV")),
    AGLASSO = apply(AGLASSO1.CV.400, 2, function(x) beta.cox(x$CV2, "CV"))
  )  
  beta2.CV.800 = list(
    GLASSO = apply(AGLASSO0.CV.800, 2, function(x) beta.cox(x$CV2, "CV")),
    GSCAD = apply(GSCAD.CV.800, 2, function(x) beta.cox(x$CV2, "CV")),
    GMCP = apply(GMCP.CV.800, 2, function(x) beta.cox(x$CV2, "CV")),
    AGLASSO = apply(AGLASSO1.CV.800, 2, function(x) beta.cox(x$CV2, "CV"))
  )
  beta2.CV.1200 = list(
    GLASSO = apply(AGLASSO0.CV.1200, 2, function(x) beta.cox(x$CV2, "CV")),
    GSCAD = apply(GSCAD.CV.1200, 2, function(x) beta.cox(x$CV2, "CV")),
    GMCP = apply(GMCP.CV.1200, 2, function(x) beta.cox(x$CV2, "CV")),
    AGLASSO = apply(AGLASSO1.CV.1200, 2, function(x) beta.cox(x$CV2, "CV"))
  )
}
### GCV
{
  beta2.GCV.400 = list(
    GLASSO = apply(AGLASSO0.sel.400, 2, function(x) beta.cox(x$sel2, "GCV")),
    GSCAD = apply(GSCAD.sel.400, 2, function(x) beta.cox(x$sel2, "GCV")),
    GMCP = apply(GMCP.sel.400, 2, function(x) beta.cox(x$sel2, "GCV")),
    AGLASSO = apply(AGLASSO1.sel.400, 2, function(x) beta.cox(x$sel2, "GCV"))
  )  
  beta2.GCV.800 = list(
    GLASSO = apply(AGLASSO0.sel.800, 2, function(x) beta.cox(x$sel2, "GCV")),
    GSCAD = apply(GSCAD.sel.800, 2, function(x) beta.cox(x$sel2, "GCV")),
    GMCP = apply(GMCP.sel.800, 2, function(x) beta.cox(x$sel2, "GCV")),
    AGLASSO = apply(AGLASSO1.sel.800, 2, function(x) beta.cox(x$sel2, "GCV"))
  )
  beta2.GCV.1200 = list(
    GLASSO = apply(AGLASSO0.sel.1200, 2, function(x) beta.cox(x$sel2, "GCV")),
    GSCAD = apply(GSCAD.sel.1200, 2, function(x) beta.cox(x$sel2, "GCV")),
    GMCP = apply(GMCP.sel.1200, 2, function(x) beta.cox(x$sel2, "GCV")),
    AGLASSO = apply(AGLASSO1.sel.1200, 2, function(x) beta.cox(x$sel2, "GCV"))
  )
}
### BIC
{
  beta2.BIC.400 = list(
    GLASSO = apply(AGLASSO0.sel.400, 2, function(x) beta.cox(x$sel2, "BIC")),
    GSCAD = apply(GSCAD.sel.400, 2, function(x) beta.cox(x$sel2, "BIC")),
    GMCP = apply(GMCP.sel.400, 2, function(x) beta.cox(x$sel2, "BIC")),
    AGLASSO = apply(AGLASSO1.sel.400, 2, function(x) beta.cox(x$sel2, "BIC"))
  )  
  beta2.BIC.800 = list(
    GLASSO = apply(AGLASSO0.sel.800, 2, function(x) beta.cox(x$sel2, "BIC")),
    GSCAD = apply(GSCAD.sel.800, 2, function(x) beta.cox(x$sel2, "BIC")),
    GMCP = apply(GMCP.sel.800, 2, function(x) beta.cox(x$sel2, "BIC")),
    AGLASSO = apply(AGLASSO1.sel.800, 2, function(x) beta.cox(x$sel2, "BIC"))
  )
  beta2.BIC.1200 = list(
    GLASSO = apply(AGLASSO0.sel.1200, 2, function(x) beta.cox(x$sel2, "BIC")),
    GSCAD = apply(GSCAD.sel.1200, 2, function(x) beta.cox(x$sel2, "BIC")),
    GMCP = apply(GMCP.sel.1200, 2, function(x) beta.cox(x$sel2, "BIC")),
    AGLASSO = apply(AGLASSO1.sel.1200, 2, function(x) beta.cox(x$sel2, "BIC"))
  )
}
### AIC
{
  beta2.AIC.400 = list(
    GLASSO = apply(AGLASSO0.sel.400, 2, function(x) beta.cox(x$sel2, "AIC")),
    GSCAD = apply(GSCAD.sel.400, 2, function(x) beta.cox(x$sel2, "AIC")),
    GMCP = apply(GMCP.sel.400, 2, function(x) beta.cox(x$sel2, "AIC")),
    AGLASSO = apply(AGLASSO1.sel.400, 2, function(x) beta.cox(x$sel2, "AIC"))
  )  
  beta2.AIC.800 = list(
    GLASSO = apply(AGLASSO0.sel.800, 2, function(x) beta.cox(x$sel2, "AIC")),
    GSCAD = apply(GSCAD.sel.800, 2, function(x) beta.cox(x$sel2, "AIC")),
    GMCP = apply(GMCP.sel.800, 2, function(x) beta.cox(x$sel2, "AIC")),
    AGLASSO = apply(AGLASSO1.sel.800, 2, function(x) beta.cox(x$sel2, "AIC"))
  )
  beta2.AIC.1200 = list(
    GLASSO = apply(AGLASSO0.sel.1200, 2, function(x) beta.cox(x$sel2, "AIC")),
    GSCAD = apply(GSCAD.sel.1200, 2, function(x) beta.cox(x$sel2, "AIC")),
    GMCP = apply(GMCP.sel.1200, 2, function(x) beta.cox(x$sel2, "AIC")),
    AGLASSO = apply(AGLASSO1.sel.1200, 2, function(x) beta.cox(x$sel2, "AIC"))
  )
}



######  ASE2.recored
#### 5-fold CV
{
  ASE2.CV.400 = list(
    GLASSO = apply(AGLASSO0.CV.400, 2, function(x) ASE.cox(x$CV2, type="CV")),
    GSCAD = apply(GSCAD.CV.400, 2, function(x)ASE.cox(x$CV2, type="CV")),
    GMCP = apply(GMCP.CV.400, 2, function(x) ASE.cox(x$CV2, type="CV")),
    AGLASSO = apply(AGLASSO1.CV.400, 2, function(x) ASE.cox(x$CV2, type="CV"))
  )
  
  
  ASE2.CV.800 = list(
    GLASSO = apply(AGLASSO0.CV.800, 2, function(x) ASE.cox(x$CV2, type="CV")),
    GSCAD = apply(GSCAD.CV.800, 2, function(x)ASE.cox(x$CV2, type="CV")),
    GMCP = apply(GMCP.CV.800, 2, function(x) ASE.cox(x$CV2, type="CV")),
    AGLASSO = apply(AGLASSO1.CV.800, 2, function(x) ASE.cox(x$CV2, type="CV"))
  )
  
  ASE2.CV.1200 = list(
    GLASSO = apply(AGLASSO0.CV.1200, 2, function(x) ASE.cox(x$CV2, type="CV")),
    GSCAD = apply(GSCAD.CV.1200, 2, function(x)ASE.cox(x$CV2, type="CV")),
    GMCP = apply(GMCP.CV.1200, 2, function(x) ASE.cox(x$CV2, type="CV")),
    AGLASSO = apply(AGLASSO1.CV.1200, 2, function(x) ASE.cox(x$CV2, type="CV"))
  )
}
#### AIC
{
  ASE2.AIC.400 = list(
    GLASSO = apply(AGLASSO0.sel.400, 2, function(x) ASE.cox(x$sel2, type="AIC")),
    GSCAD = apply(GSCAD.sel.400, 2, function(x)ASE.cox(x$sel2, type="AIC")),
    GMCP = apply(GMCP.sel.400, 2, function(x) ASE.cox(x$sel2, type="AIC")),
    AGLASSO = apply(AGLASSO1.sel.400, 2, function(x) ASE.cox(x$sel2, type="AIC"))
  )
  
  ASE2.AIC.800 = list(
    GLASSO = apply(AGLASSO0.sel.800, 2, function(x) ASE.cox(x$sel2, type="AIC")),
    GSCAD = apply(GSCAD.sel.800, 2, function(x)ASE.cox(x$sel2, type="AIC")),
    GMCP = apply(GMCP.sel.800, 2, function(x) ASE.cox(x$sel2, type="AIC")),
    AGLASSO = apply(AGLASSO1.sel.800, 2, function(x) ASE.cox(x$sel2, type="AIC"))
  )
  
  ASE2.AIC.1200 = list(
    GLASSO = apply(AGLASSO0.sel.1200, 2, function(x) ASE.cox(x$sel2, type="AIC")),
    GSCAD = apply(GSCAD.sel.1200, 2, function(x)ASE.cox(x$sel2, type="AIC")),
    GMCP = apply(GMCP.sel.1200, 2, function(x) ASE.cox(x$sel2, type="AIC")),
    AGLASSO = apply(AGLASSO1.sel.1200, 2, function(x) ASE.cox(x$sel2, type="AIC"))
  )  
}
#### BIC
{
  ASE2.BIC.400 = list(
    GLASSO = apply(AGLASSO0.sel.400, 2, function(x) ASE.cox(x$sel2, type="BIC")),
    GSCAD = apply(GSCAD.sel.400, 2, function(x)ASE.cox(x$sel2, type="BIC")),
    GMCP = apply(GMCP.sel.400, 2, function(x) ASE.cox(x$sel2, type="BIC")),
    AGLASSO = apply(AGLASSO1.sel.400, 2, function(x) ASE.cox(x$sel2, type="BIC"))
  )
  
  ASE2.BIC.800 = list(
    GLASSO = apply(AGLASSO0.sel.800, 2, function(x) ASE.cox(x$sel2, type="BIC")),
    GSCAD = apply(GSCAD.sel.800, 2, function(x)ASE.cox(x$sel2, type="BIC")),
    GMCP = apply(GMCP.sel.800, 2, function(x) ASE.cox(x$sel2, type="BIC")),
    AGLASSO = apply(AGLASSO1.sel.800, 2, function(x) ASE.cox(x$sel2, type="BIC"))
  )
  
  ASE2.BIC.1200 = list(
    GLASSO = apply(AGLASSO0.sel.1200, 2, function(x) ASE.cox(x$sel2, type="BIC")),
    GSCAD = apply(GSCAD.sel.1200, 2, function(x)ASE.cox(x$sel2, type="BIC")),
    GMCP = apply(GMCP.sel.1200, 2, function(x) ASE.cox(x$sel2, type="BIC")),
    AGLASSO = apply(AGLASSO1.sel.1200, 2, function(x) ASE.cox(x$sel2, type="BIC"))
  )  
}
#### GCV
{
  ASE2.GCV.400 = list(
    GLASSO = apply(AGLASSO0.sel.400, 2, function(x) ASE.cox(x$sel2, type="GCV")),
    GSCAD = apply(GSCAD.sel.400, 2, function(x)ASE.cox(x$sel2, type="GCV")),
    GMCP = apply(GMCP.sel.400, 2, function(x) ASE.cox(x$sel2, type="GCV")),
    AGLASSO = apply(AGLASSO1.sel.400, 2, function(x) ASE.cox(x$sel2, type="GCV"))
  )
  
  ASE2.GCV.800 = list(
    GLASSO = apply(AGLASSO0.sel.800, 2, function(x) ASE.cox(x$sel2, type="GCV")),
    GSCAD = apply(GSCAD.sel.800, 2, function(x)ASE.cox(x$sel2, type="GCV")),
    GMCP = apply(GMCP.sel.800, 2, function(x) ASE.cox(x$sel2, type="GCV")),
    AGLASSO = apply(AGLASSO1.sel.800, 2, function(x) ASE.cox(x$sel2, type="GCV"))
  )
  
  ASE2.GCV.1200 = list(
    GLASSO = apply(AGLASSO0.sel.1200, 2, function(x) ASE.cox(x$sel2, type="GCV")),
    GSCAD = apply(GSCAD.sel.1200, 2, function(x)ASE.cox(x$sel2, type="GCV")),
    GMCP = apply(GMCP.sel.1200, 2, function(x) ASE.cox(x$sel2, type="GCV")),
    AGLASSO = apply(AGLASSO1.sel.1200, 2, function(x) ASE.cox(x$sel2, type="GCV"))
  )
}

######  ASE1.recored
#### 5-fold CV
{
  ASE1.CV.400 = list(
    GLASSO = apply(AGLASSO0.CV.400, 2, function(x) ASE.cox(x$CV1, type="CV")),
    GSCAD = apply(GSCAD.CV.400, 2, function(x)ASE.cox(x$CV1, type="CV")),
    GMCP = apply(GMCP.CV.400, 2, function(x) ASE.cox(x$CV1, type="CV")),
    AGLASSO = apply(AGLASSO1.CV.400, 2, function(x) ASE.cox(x$CV1, type="CV"))
  )
  
  ASE1.CV.800 = list(
    GLASSO = apply(AGLASSO0.CV.800, 2, function(x) ASE.cox(x$CV1, type="CV")),
    GSCAD = apply(GSCAD.CV.800, 2, function(x)ASE.cox(x$CV1, type="CV")),
    GMCP = apply(GMCP.CV.800, 2, function(x) ASE.cox(x$CV1, type="CV")),
    AGLASSO = apply(AGLASSO1.CV.800, 2, function(x) ASE.cox(x$CV1, type="CV"))
  )
  
  ASE1.CV.1200 = list(
    GLASSO = apply(AGLASSO0.CV.1200, 2, function(x) ASE.cox(x$CV1, type="CV")),
    GSCAD = apply(GSCAD.CV.1200, 2, function(x)ASE.cox(x$CV1, type="CV")),
    GMCP = apply(GMCP.CV.1200, 2, function(x) ASE.cox(x$CV1, type="CV")),
    AGLASSO = apply(AGLASSO1.CV.1200, 2, function(x) ASE.cox(x$CV1, type="CV"))
  )
}
#### AIC
{
  ASE1.AIC.400 = list(
    GLASSO = apply(AGLASSO0.sel.400, 2, function(x) ASE.cox(x$sel1, type="AIC")),
    GSCAD = apply(GSCAD.sel.400, 2, function(x)ASE.cox(x$sel1, type="AIC")),
    GMCP = apply(GMCP.sel.400, 2, function(x) ASE.cox(x$sel1, type="AIC")),
    AGLASSO = apply(AGLASSO1.sel.400, 2, function(x) ASE.cox(x$sel1, type="AIC"))
  )
  
  ASE1.AIC.800 = list(
    GLASSO = apply(AGLASSO0.sel.800, 2, function(x) ASE.cox(x$sel1, type="AIC")),
    GSCAD = apply(GSCAD.sel.800, 2, function(x)ASE.cox(x$sel1, type="AIC")),
    GMCP = apply(GMCP.sel.800, 2, function(x) ASE.cox(x$sel1, type="AIC")),
    AGLASSO = apply(AGLASSO1.sel.800, 2, function(x) ASE.cox(x$sel1, type="AIC"))
  )
  
  ASE1.AIC.1200 = list(
    GLASSO = apply(AGLASSO0.sel.1200, 2, function(x) ASE.cox(x$sel1, type="AIC")),
    GSCAD = apply(GSCAD.sel.1200, 2, function(x)ASE.cox(x$sel1, type="AIC")),
    GMCP = apply(GMCP.sel.1200, 2, function(x) ASE.cox(x$sel1, type="AIC")),
    AGLASSO = apply(AGLASSO1.sel.1200, 2, function(x) ASE.cox(x$sel1, type="AIC"))
  )  
}
#### BIC
{
  ASE1.BIC.400 = list(
    GLASSO = apply(AGLASSO0.sel.400, 2, function(x) ASE.cox(x$sel1, type="BIC")),
    GSCAD = apply(GSCAD.sel.400, 2, function(x)ASE.cox(x$sel1, type="BIC")),
    GMCP = apply(GMCP.sel.400, 2, function(x) ASE.cox(x$sel1, type="BIC")),
    AGLASSO = apply(AGLASSO1.sel.400, 2, function(x) ASE.cox(x$sel1, type="BIC"))
  )
  
  ASE1.BIC.800 = list(
    GLASSO = apply(AGLASSO0.sel.800, 2, function(x) ASE.cox(x$sel1, type="BIC")),
    GSCAD = apply(GSCAD.sel.800, 2, function(x)ASE.cox(x$sel1, type="BIC")),
    GMCP = apply(GMCP.sel.800, 2, function(x) ASE.cox(x$sel1, type="BIC")),
    AGLASSO = apply(AGLASSO1.sel.800, 2, function(x) ASE.cox(x$sel1, type="BIC"))
  )
  
  ASE1.BIC.1200 = list(
    GLASSO = apply(AGLASSO0.sel.1200, 2, function(x) ASE.cox(x$sel1, type="BIC")),
    GSCAD = apply(GSCAD.sel.1200, 2, function(x)ASE.cox(x$sel1, type="BIC")),
    GMCP = apply(GMCP.sel.1200, 2, function(x) ASE.cox(x$sel1, type="BIC")),
    AGLASSO = apply(AGLASSO1.sel.1200, 2, function(x) ASE.cox(x$sel1, type="BIC"))
  )  
}
#### GCV
{
  ASE1.GCV.400 = list(
    GLASSO = apply(AGLASSO0.sel.400, 2, function(x) ASE.cox(x$sel1, type="GCV")),
    GSCAD = apply(GSCAD.sel.400, 2, function(x)ASE.cox(x$sel1, type="GCV")),
    GMCP = apply(GMCP.sel.400, 2, function(x) ASE.cox(x$sel1, type="GCV")),
    AGLASSO = apply(AGLASSO1.sel.400, 2, function(x) ASE.cox(x$sel1, type="GCV"))
  )
  
  ASE1.GCV.800 = list(
    GLASSO = apply(AGLASSO0.sel.800, 2, function(x) ASE.cox(x$sel1, type="GCV")),
    GSCAD = apply(GSCAD.sel.800, 2, function(x)ASE.cox(x$sel1, type="GCV")),
    GMCP = apply(GMCP.sel.800, 2, function(x) ASE.cox(x$sel1, type="GCV")),
    AGLASSO = apply(AGLASSO1.sel.800, 2, function(x) ASE.cox(x$sel1, type="GCV"))
  )
  
  ASE1.GCV.1200 = list(
    GLASSO = apply(AGLASSO0.sel.1200, 2, function(x) ASE.cox(x$sel1, type="GCV")),
    GSCAD = apply(GSCAD.sel.1200, 2, function(x)ASE.cox(x$sel1, type="GCV")),
    GMCP = apply(GMCP.sel.1200, 2, function(x) ASE.cox(x$sel1, type="GCV")),
    AGLASSO = apply(AGLASSO1.sel.1200, 2, function(x) ASE.cox(x$sel1, type="GCV"))
  )
}

rm(AGLASSO1.CV.1200)
rm(AGLASSO1.CV.800)
rm(AGLASSO1.CV.400)
rm(AGLASSO1.sel.1200)
rm(AGLASSO1.sel.800)
rm(AGLASSO1.sel.400)
rm(AGLASSO0.CV.1200)
rm(AGLASSO0.CV.800)
rm(AGLASSO0.CV.400)
rm(AGLASSO0.sel.1200)
rm(AGLASSO0.sel.800)
rm(AGLASSO0.sel.400)
rm(GSCAD.CV.1200)
rm(GSCAD.CV.800)
rm(GSCAD.CV.400)
rm(GSCAD.sel.1200)
rm(GSCAD.sel.800)
rm(GSCAD.sel.400)
rm(GMCP.CV.1200)
rm(GMCP.CV.800)
rm(GMCP.CV.400)
rm(GMCP.sel.1200)
rm(GMCP.sel.800)
rm(GMCP.sel.400)
#save.image("E:/Google ���ݵw��/program/yang/group penalized regression 2013/group cox/�T�wweightver2/simulation all/OK/betaresult.cox.RData")
#save.image("E:/Google ���ݵw��/program/yang/group penalized regression 2013/group cox/�T�wweightver2/simulation all/OK/betaresult2.cox.RData")
