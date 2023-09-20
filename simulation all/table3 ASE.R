#load("E:/Google ¶³ºÝµwºÐ/program/yang/group penalized regression 2013/group cox/©T©wweightver2/simulation all/OK/betaresult.RData")
table3 = function(sample, asym, posi){
  a=round(apply(sample, 1 ,sd)*1000)
  b=round(apply(asym, 1, mean)*1000)
  c=round(apply(asym, 1, sd)*1000)
  
  tmp = numeric()
  for(j in posi){
    tmp = c(tmp, a[j], paste0(b[j], "(", c[j], ")"))
  }
  return(tmp)
}
r0.400 = table3(cox2.t.400, ASE.cox2.t.400, 1:5)
r0.800 = table3(cox2.t.800, ASE.cox2.t.800, 1:5)
r0.1200 = table3(cox2.t.1200, ASE.cox2.t.1200, 1:5)
line=rep(NA, 10)
## ASE2
## CV
tmp1 = rbind(table3(beta2.CV.400$GLASSO,   ASE2.CV.400$GLASSO, posi),
             table3(beta2.CV.400$GSCAD,    ASE2.CV.400$GSCAD, posi),
             table3(beta2.CV.400$GMCP,     ASE2.CV.400$GMCP, posi),
             table3(beta2.CV.400$AGLASSO,  ASE2.CV.400$AGLASSO, posi))
tmp2 = rbind(table3(beta2.CV.800$GLASSO,   ASE2.CV.800$GLASSO, posi),
             table3(beta2.CV.800$GSCAD,    ASE2.CV.800$GSCAD, posi),
             table3(beta2.CV.800$GMCP,     ASE2.CV.800$GMCP, posi),
             table3(beta2.CV.800$AGLASSO,  ASE2.CV.800$AGLASSO, posi))
tmp3 = rbind(table3(beta2.CV.1200$GLASSO,  ASE2.CV.1200$GLASSO, posi),
             table3(beta2.CV.1200$GSCAD,   ASE2.CV.1200$GSCAD, posi),
             table3(beta2.CV.1200$GMCP,    ASE2.CV.1200$GMCP, posi),
             table3(beta2.CV.1200$AGLASSO, ASE2.CV.1200$AGLASSO, posi))
table3.CV = rbind(r0.400, tmp1,
                  line,
                  r0.800, tmp2,
                  line, 
                  r0.1200, tmp3)
## GCV

tmp1 = rbind(table3(beta2.GCV.400$GLASSO,   ASE2.GCV.400$GLASSO, posi),
             table3(beta2.GCV.400$GSCAD,    ASE2.GCV.400$GSCAD, posi),
             table3(beta2.GCV.400$GMCP,     ASE2.GCV.400$GMCP, posi),
             table3(beta2.GCV.400$AGLASSO,  ASE2.GCV.400$AGLASSO, posi))
tmp2 = rbind(table3(beta2.GCV.800$GLASSO,   ASE2.GCV.800$GLASSO, posi),
             table3(beta2.GCV.800$GSCAD,    ASE2.GCV.800$GSCAD, posi),
             table3(beta2.GCV.800$GMCP,     ASE2.GCV.800$GMCP, posi),
             table3(beta2.GCV.800$AGLASSO,  ASE2.GCV.800$AGLASSO, posi))
tmp3 = rbind(table3(beta2.GCV.1200$GLASSO,  ASE2.GCV.1200$GLASSO, posi),
             table3(beta2.GCV.1200$GSCAD,   ASE2.GCV.1200$GSCAD, posi),
             table3(beta2.GCV.1200$GMCP,    ASE2.GCV.1200$GMCP, posi),
             table3(beta2.GCV.1200$AGLASSO, ASE2.GCV.1200$AGLASSO, posi))
table3.GCV = rbind(r0.400, tmp1,
                   line,
                   r0.800, tmp2,
                   line, 
                   r0.1200, tmp3)
## BIC

tmp1 = rbind(table3(beta2.BIC.400$GLASSO,   ASE2.BIC.400$GLASSO, posi),
             table3(beta2.BIC.400$GSCAD,    ASE2.BIC.400$GSCAD, posi),
             table3(beta2.BIC.400$GMCP,     ASE2.BIC.400$GMCP, posi),
             table3(beta2.BIC.400$AGLASSO,  ASE2.BIC.400$AGLASSO, posi))
tmp2 = rbind(table3(beta2.BIC.800$GLASSO,   ASE2.BIC.800$GLASSO, posi),
             table3(beta2.BIC.800$GSCAD,    ASE2.BIC.800$GSCAD, posi),
             table3(beta2.BIC.800$GMCP,     ASE2.BIC.800$GMCP, posi),
             table3(beta2.BIC.800$AGLASSO,  ASE2.BIC.800$AGLASSO, posi))
tmp3 = rbind(table3(beta2.BIC.1200$GLASSO,  ASE2.BIC.1200$GLASSO, posi),
             table3(beta2.BIC.1200$GSCAD,   ASE2.BIC.1200$GSCAD, posi),
             table3(beta2.BIC.1200$GMCP,    ASE2.BIC.1200$GMCP, posi),
             table3(beta2.BIC.1200$AGLASSO, ASE2.BIC.1200$AGLASSO, posi))
table3.BIC = rbind(r0.400, tmp1,
                   line,
                   r0.800, tmp2,
                   line, 
                   r0.1200, tmp3)
## AIC

tmp1 = rbind(table3(beta2.AIC.400$GLASSO,   ASE2.AIC.400$GLASSO, posi),
             table3(beta2.AIC.400$GSCAD,    ASE2.AIC.400$GSCAD, posi),
             table3(beta2.AIC.400$GMCP,     ASE2.AIC.400$GMCP, posi),
             table3(beta2.AIC.400$AGLASSO,  ASE2.AIC.400$AGLASSO, posi))
tmp2 = rbind(table3(beta2.AIC.800$GLASSO,   ASE2.AIC.800$GLASSO, posi),
             table3(beta2.AIC.800$GSCAD,    ASE2.AIC.800$GSCAD, posi),
             table3(beta2.AIC.800$GMCP,     ASE2.AIC.800$GMCP, posi),
             table3(beta2.AIC.800$AGLASSO,  ASE2.AIC.800$AGLASSO, posi))
tmp3 = rbind(table3(beta2.AIC.1200$GLASSO,  ASE2.AIC.1200$GLASSO, posi),
             table3(beta2.AIC.1200$GSCAD,   ASE2.AIC.1200$GSCAD, posi),
             table3(beta2.AIC.1200$GMCP,    ASE2.AIC.1200$GMCP, posi),
             table3(beta2.AIC.1200$AGLASSO, ASE2.AIC.1200$AGLASSO, posi))
table3.AIC = rbind(r0.400, tmp1,
                   line,
                   r0.800, tmp2,
                   line, 
                   r0.1200, tmp3)
####
# write.csv(table3.CV, "E:/Google ¶³ºÝµwºÐ/program/yang/group penalized regression 2013/group cox/©T©wweightver2/simulation all/table3.cox.CV.csv")
# write.csv(table3.GCV, "E:/Google ¶³ºÝµwºÐ/program/yang/group penalized regression 2013/group cox/©T©wweightver2/simulation all/table3.cox.GCV.csv")
# write.csv(table3.BIC, "E:/Google ¶³ºÝµwºÐ/program/yang/group penalized regression 2013/group cox/©T©wweightver2/simulation all/table3.cox.BIC.csv")
# write.csv(table3.AIC, "E:/Google ¶³ºÝµwºÐ/program/yang/group penalized regression 2013/group cox/©T©wweightver2/simulation all/table3.cox.AIC.csv")
