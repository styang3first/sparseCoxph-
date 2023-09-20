# path_root = dirname(rstudioapi::getActiveDocumentContext()$path); setwd(path_root)
load("OK/betaresult2.cox.RData")
##############
MSE = function(beta.re, data){
  tmp = sapply(1:100, function(i){
    V = cov(data[,i]$x.m)
    beta = beta.re[,i]
    output = t(beta-beta.true)%*%V%*%(beta-beta.true)
    return(output)
  })
  a = round(mean(tmp),3)
  b = round(sd(tmp),3)
  return( paste0(a, "(", b, ")") )
}
active = function(beta.re){
  A = apply(beta.re, 2, function(x) c(sum(x[-posi]==0), sum(x[posi]==0)) )
  return( apply(A, 1, mean) )
}
table1 = function(beta2, beta1, data){
  output = c(active(beta.re=beta2), 
             MSE(beta.re=beta2, data=data),
             active(beta.re=beta1), 
             MSE(beta.re=beta1, data=data))
  return(output)
}
####table1
r0.400 = c(0, 0, MSE(cox2.400, data400), 0, 0, MSE(cox1.400, data400))
r0.800 = c(0, 0, MSE(cox2.800, data800), 0, 0, MSE(cox1.800, data800))
r0.1200 = c(0, 0, MSE(cox2.1200, data1200), 0, 0, MSE(cox1.1200, data1200))
line = rep(NA, 6)
## CV
tmp1 = table1(beta2.CV.400$GLASSO, beta1.CV.400$GLASSO, data400)
tmp2 = table1(beta2.CV.400$GSCAD, beta1.CV.400$GSCAD, data400)
tmp3 = table1(beta2.CV.400$GMCP, beta1.CV.400$GMCP, data400)
tmp4 = table1(beta2.CV.400$AGLASSO, beta1.CV.400$AGLASSO, data400)
tmp5 = table1(beta2.CV.800$GLASSO, beta1.CV.800$GLASSO, data800)
tmp6 = table1(beta2.CV.800$GSCAD, beta1.CV.800$GSCAD, data800)
tmp7 = table1(beta2.CV.800$GMCP, beta1.CV.800$GMCP, data800)
tmp8 = table1(beta2.CV.800$AGLASSO, beta1.CV.800$AGLASSO, data800)
tmp9 = table1(beta2.CV.1200$GLASSO, beta1.CV.1200$GLASSO, data1200)
tmp10 = table1(beta2.CV.1200$GSCAD, beta1.CV.1200$GSCAD, data1200)
tmp11 = table1(beta2.CV.1200$GMCP, beta1.CV.1200$GMCP, data1200)
tmp12 = table1(beta2.CV.1200$AGLASSO, beta1.CV.1200$AGLASSO, data1200)
table1.CV = rbind(r0.400, tmp1, tmp2, tmp3, tmp4,
                  line,
                  r0.800, tmp5, tmp6, tmp7, tmp8,
                  line,
                  r0.1200, tmp9, tmp10, tmp11, tmp12)
## GCV
tmp1 = table1(beta2.GCV.400$GLASSO, beta1.GCV.400$GLASSO, data400)
tmp2 = table1(beta2.GCV.400$GSCAD, beta1.GCV.400$GSCAD, data400)
tmp3 = table1(beta2.GCV.400$GMCP, beta1.GCV.400$GMCP, data400)
tmp4 = table1(beta2.GCV.400$AGLASSO, beta1.GCV.400$AGLASSO, data400)
tmp5 = table1(beta2.GCV.800$GLASSO, beta1.GCV.800$GLASSO, data800)
tmp6 = table1(beta2.GCV.800$GSCAD, beta1.GCV.800$GSCAD, data800)
tmp7 = table1(beta2.GCV.800$GMCP, beta1.GCV.800$GMCP, data800)
tmp8 = table1(beta2.GCV.800$AGLASSO, beta1.GCV.800$AGLASSO, data800)
tmp9 = table1(beta2.GCV.1200$GLASSO, beta1.GCV.1200$GLASSO, data1200)
tmp10 = table1(beta2.GCV.1200$GSCAD, beta1.GCV.1200$GSCAD, data1200)
tmp11 = table1(beta2.GCV.1200$GMCP, beta1.GCV.1200$GMCP, data1200)
tmp12 = table1(beta2.GCV.1200$AGLASSO, beta1.GCV.1200$AGLASSO, data1200)
table1.GCV = rbind(r0.400, tmp1, tmp2, tmp3, tmp4,
                  line,
                  r0.800, tmp5, tmp6, tmp7, tmp8,
                  line,
                  r0.1200, tmp9, tmp10, tmp11, tmp12)
## BIC
tmp1 = table1(beta2.BIC.400$GLASSO, beta1.BIC.400$GLASSO, data400)
tmp2 = table1(beta2.BIC.400$GSCAD, beta1.BIC.400$GSCAD, data400)
tmp3 = table1(beta2.BIC.400$GMCP, beta1.BIC.400$GMCP, data400)
tmp4 = table1(beta2.BIC.400$AGLASSO, beta1.BIC.400$AGLASSO, data400)
tmp5 = table1(beta2.BIC.800$GLASSO, beta1.BIC.800$GLASSO, data800)
tmp6 = table1(beta2.BIC.800$GSCAD, beta1.BIC.800$GSCAD, data800)
tmp7 = table1(beta2.BIC.800$GMCP, beta1.BIC.800$GMCP, data800)
tmp8 = table1(beta2.BIC.800$AGLASSO, beta1.BIC.800$AGLASSO, data800)
tmp9 = table1(beta2.BIC.1200$GLASSO, beta1.BIC.1200$GLASSO, data1200)
tmp10 = table1(beta2.BIC.1200$GSCAD, beta1.BIC.1200$GSCAD, data1200)
tmp11 = table1(beta2.BIC.1200$GMCP, beta1.BIC.1200$GMCP, data1200)
tmp12 = table1(beta2.BIC.1200$AGLASSO, beta1.BIC.1200$AGLASSO, data1200)
table1.BIC = rbind(r0.400, tmp1, tmp2, tmp3, tmp4,
                  line,
                  r0.800, tmp5, tmp6, tmp7, tmp8,
                  line,
                  r0.1200, tmp9, tmp10, tmp11, tmp12)
## AIC
tmp1 = table1(beta2.AIC.400$GLASSO, beta1.AIC.400$GLASSO, data400)
tmp2 = table1(beta2.AIC.400$GSCAD, beta1.AIC.400$GSCAD, data400)
tmp3 = table1(beta2.AIC.400$GMCP, beta1.AIC.400$GMCP, data400)
tmp4 = table1(beta2.AIC.400$AGLASSO, beta1.AIC.400$AGLASSO, data400)
tmp5 = table1(beta2.AIC.800$GLASSO, beta1.AIC.800$GLASSO, data800)
tmp6 = table1(beta2.AIC.800$GSCAD, beta1.AIC.800$GSCAD, data800)
tmp7 = table1(beta2.AIC.800$GMCP, beta1.AIC.800$GMCP, data800)
tmp8 = table1(beta2.AIC.800$AGLASSO, beta1.AIC.800$AGLASSO, data800)
tmp9 = table1(beta2.AIC.1200$GLASSO, beta1.AIC.1200$GLASSO, data1200)
tmp10 = table1(beta2.AIC.1200$GSCAD, beta1.AIC.1200$GSCAD, data1200)
tmp11 = table1(beta2.AIC.1200$GMCP, beta1.AIC.1200$GMCP, data1200)
tmp12 = table1(beta2.AIC.1200$AGLASSO, beta1.AIC.1200$AGLASSO, data1200)
table1.AIC = rbind(r0.400, tmp1, tmp2, tmp3, tmp4,
                  line,
                  r0.800, tmp5, tmp6, tmp7, tmp8,
                  line,
                  r0.1200, tmp9, tmp10, tmp11, tmp12)
rownames(table1.CV) = rownames(table1.GCV) = rownames(table1.BIC) = rownames(table1.AIC) =  
  c("400&MLE_{full}", "&GLASSO", "&GSCAD", "&GMCP", "&AGLASSO",
    "&", "800&MLE_{full}", "&GLASSO", "&GSCAD", "&GMCP", "&AGLASSO",
    "&", "1200&MLE_{full}", "&GLASSO", "&GSCAD", "&GMCP", "&AGLASSO")

#####
getwd()
# write.csv(table1.CV, "E:/Google ¶³ºÝµwºÐ/program/yang/group penalized regression 2013/group cox/©T©wweightver2/simulation all/table1.CV.csv")
# write.csv(table1.GCV, "E:/Google ¶³ºÝµwºÐ/program/yang/group penalized regression 2013/group cox/©T©wweightver2/simulation all/table1.GCV.csv")
# write.csv(table1.BIC, "E:/Google ¶³ºÝµwºÐ/program/yang/group penalized regression 2013/group cox/©T©wweightver2/simulation all/table1.BIC.csv")
# write.csv(table1.AIC, "E:/Google ¶³ºÝµwºÐ/program/yang/group penalized regression 2013/group cox/©T©wweightver2/simulation all/table1.AIC.csv")
