#load("D:/Google ¶³ºÝµwºÐ/program/yang/group penalized regression 2013/group cox/©T©wweightver2/simulation all/OK/betaresult2.cox.RData")
beta.true2=beta.true[posi]
table4 = function(beta.est, ASE.est){
  ASE.est[which(beta.est==0)]=0
  beta.diff = apply(beta.est, 2, function(x) x-beta.true)
  tmp = (ASE.est*1.96) > abs(beta.diff) 
  output = apply(tmp, 1, sum)[posi]
  for(i in 1:5) output[i] = paste0(output[i], "%")
  return(output)
}
table4.cox = function(beta.est, ASE.est){
  beta.diff = apply(beta.est, 2, function(x) x-beta.true2)
  tmp = (ASE.est*1.96) > abs(beta.diff) 
  output = apply(tmp, 1, sum)
  for(i in 1:5) output[i] = paste0(output[i], "%")
  return(output)
}
line = rep(NA, 5)
##############
tmp1.400 = table4.cox(cox1.t.400, ASE.cox1.t.400)
tmp1.800 = table4.cox(cox1.t.800, ASE.cox1.t.800)
tmp1.1200 = table4.cox(cox1.t.1200, ASE.cox1.t.1200)
tmp2.400 = table4.cox(cox2.t.400, ASE.cox2.t.400)
tmp2.800 = table4.cox(cox2.t.800, ASE.cox2.t.800)
tmp2.1200 = table4.cox(cox2.t.1200, ASE.cox2.t.1200)
#### CV
tmp1 = rbind(table4(beta1.CV.400$GLASSO,   ASE1.CV.400$GLASSO),
             table4(beta1.CV.400$GSCAD,    ASE1.CV.400$GSCAD),
             table4(beta1.CV.400$GMCP,     ASE1.CV.400$GMCP),
             table4(beta1.CV.400$AGLASSO,  ASE1.CV.400$AGLASSO))
tmp2 = rbind(table4(beta1.CV.800$GLASSO,   ASE1.CV.800$GLASSO),
             table4(beta1.CV.800$GSCAD,    ASE1.CV.800$GSCAD),
             table4(beta1.CV.800$GMCP,     ASE1.CV.800$GMCP),
             table4(beta1.CV.800$AGLASSO,  ASE1.CV.800$AGLASSO))
tmp3 = rbind(table4(beta1.CV.1200$GLASSO,  ASE1.CV.1200$GLASSO),
             table4(beta1.CV.1200$GSCAD,   ASE1.CV.1200$GSCAD),
             table4(beta1.CV.1200$GMCP,    ASE1.CV.1200$GMCP),
             table4(beta1.CV.1200$AGLASSO, ASE1.CV.1200$AGLASSO))

tmp4 = rbind(table4(beta2.CV.400$GLASSO,   ASE2.CV.400$GLASSO),
             table4(beta2.CV.400$GSCAD,    ASE2.CV.400$GSCAD),
             table4(beta2.CV.400$GMCP,     ASE2.CV.400$GMCP),
             table4(beta2.CV.400$AGLASSO,  ASE2.CV.400$AGLASSO))
tmp5 = rbind(table4(beta2.CV.800$GLASSO,   ASE2.CV.800$GLASSO),
             table4(beta2.CV.800$GSCAD,    ASE2.CV.800$GSCAD),
             table4(beta2.CV.800$GMCP,     ASE2.CV.800$GMCP),
             table4(beta2.CV.800$AGLASSO,  ASE2.CV.800$AGLASSO))
tmp6 = rbind(table4(beta2.CV.1200$GLASSO,  ASE2.CV.1200$GLASSO),
             table4(beta2.CV.1200$GSCAD,   ASE2.CV.1200$GSCAD),
             table4(beta2.CV.1200$GMCP,    ASE2.CV.1200$GMCP),
             table4(beta2.CV.1200$AGLASSO, ASE2.CV.1200$AGLASSO))
table41.CV = rbind(tmp1.400, tmp1, line, tmp1.800, tmp2, line, tmp1.1200, tmp3)
table42.CV = rbind(tmp2.400, tmp4, line, tmp2.800, tmp5, line, tmp2.1200, tmp6)
table4.CV = cbind(table42.CV, table41.CV)

#### GCV
tmp1 = rbind(table4(beta1.GCV.400$GLASSO,   ASE1.GCV.400$GLASSO),
             table4(beta1.GCV.400$GSCAD,    ASE1.GCV.400$GSCAD),
             table4(beta1.GCV.400$GMCP,     ASE1.GCV.400$GMCP),
             table4(beta1.GCV.400$AGLASSO,  ASE1.GCV.400$AGLASSO))
tmp2 = rbind(table4(beta1.GCV.800$GLASSO,   ASE1.GCV.800$GLASSO),
             table4(beta1.GCV.800$GSCAD,    ASE1.GCV.800$GSCAD),
             table4(beta1.GCV.800$GMCP,     ASE1.GCV.800$GMCP),
             table4(beta1.GCV.800$AGLASSO,  ASE1.GCV.800$AGLASSO))
tmp3 = rbind(table4(beta1.GCV.1200$GLASSO,  ASE1.GCV.1200$GLASSO),
             table4(beta1.GCV.1200$GSCAD,   ASE1.GCV.1200$GSCAD),
             table4(beta1.GCV.1200$GMCP,    ASE1.GCV.1200$GMCP),
             table4(beta1.GCV.1200$AGLASSO, ASE1.GCV.1200$AGLASSO))

tmp4 = rbind(table4(beta2.GCV.400$GLASSO,   ASE2.GCV.400$GLASSO),
             table4(beta2.GCV.400$GSCAD,    ASE2.GCV.400$GSCAD),
             table4(beta2.GCV.400$GMCP,     ASE2.GCV.400$GMCP),
             table4(beta2.GCV.400$AGLASSO,  ASE2.GCV.400$AGLASSO))
tmp5 = rbind(table4(beta2.GCV.800$GLASSO,   ASE2.GCV.800$GLASSO),
             table4(beta2.GCV.800$GSCAD,    ASE2.GCV.800$GSCAD),
             table4(beta2.GCV.800$GMCP,     ASE2.GCV.800$GMCP),
             table4(beta2.GCV.800$AGLASSO,  ASE2.GCV.800$AGLASSO))
tmp6 = rbind(table4(beta2.GCV.1200$GLASSO,  ASE2.GCV.1200$GLASSO),
             table4(beta2.GCV.1200$GSCAD,   ASE2.GCV.1200$GSCAD),
             table4(beta2.GCV.1200$GMCP,    ASE2.GCV.1200$GMCP),
             table4(beta2.GCV.1200$AGLASSO, ASE2.GCV.1200$AGLASSO))
table41.GCV = rbind(tmp1.400, tmp1, line, tmp1.800, tmp2, line, tmp1.1200, tmp3)
table42.GCV = rbind(tmp2.400, tmp4, line, tmp2.800, tmp5, line, tmp2.1200, tmp6)
table4.GCV = cbind(table42.GCV, table41.GCV)

#### BIC
tmp1 = rbind(table4(beta1.BIC.400$GLASSO,   ASE1.BIC.400$GLASSO),
             table4(beta1.BIC.400$GSCAD,    ASE1.BIC.400$GSCAD),
             table4(beta1.BIC.400$GMCP,     ASE1.BIC.400$GMCP),
             table4(beta1.BIC.400$AGLASSO,  ASE1.BIC.400$AGLASSO))
tmp2 = rbind(table4(beta1.BIC.800$GLASSO,   ASE1.BIC.800$GLASSO),
             table4(beta1.BIC.800$GSCAD,    ASE1.BIC.800$GSCAD),
             table4(beta1.BIC.800$GMCP,     ASE1.BIC.800$GMCP),
             table4(beta1.BIC.800$AGLASSO,  ASE1.BIC.800$AGLASSO))
tmp3 = rbind(table4(beta1.BIC.1200$GLASSO,  ASE1.BIC.1200$GLASSO),
             table4(beta1.BIC.1200$GSCAD,   ASE1.BIC.1200$GSCAD),
             table4(beta1.BIC.1200$GMCP,    ASE1.BIC.1200$GMCP),
             table4(beta1.BIC.1200$AGLASSO, ASE1.BIC.1200$AGLASSO))

tmp4 = rbind(table4(beta2.BIC.400$GLASSO,   ASE2.BIC.400$GLASSO),
             table4(beta2.BIC.400$GSCAD,    ASE2.BIC.400$GSCAD),
             table4(beta2.BIC.400$GMCP,     ASE2.BIC.400$GMCP),
             table4(beta2.BIC.400$AGLASSO,  ASE2.BIC.400$AGLASSO))
tmp5 = rbind(table4(beta2.BIC.800$GLASSO,   ASE2.BIC.800$GLASSO),
             table4(beta2.BIC.800$GSCAD,    ASE2.BIC.800$GSCAD),
             table4(beta2.BIC.800$GMCP,     ASE2.BIC.800$GMCP),
             table4(beta2.BIC.800$AGLASSO,  ASE2.BIC.800$AGLASSO))
tmp6 = rbind(table4(beta2.BIC.1200$GLASSO,  ASE2.BIC.1200$GLASSO),
             table4(beta2.BIC.1200$GSCAD,   ASE2.BIC.1200$GSCAD),
             table4(beta2.BIC.1200$GMCP,    ASE2.BIC.1200$GMCP),
             table4(beta2.BIC.1200$AGLASSO, ASE2.BIC.1200$AGLASSO))
table41.BIC = rbind(tmp1.400, tmp1, line, tmp1.800, tmp2, line, tmp1.1200, tmp3)
table42.BIC = rbind(tmp2.400, tmp4, line, tmp2.800, tmp5, line, tmp2.1200, tmp6)
table4.BIC = cbind(table42.BIC, table41.BIC)

#### AIC
tmp1 = rbind(table4(beta1.AIC.400$GLASSO,   ASE1.AIC.400$GLASSO),
             table4(beta1.AIC.400$GSCAD,    ASE1.AIC.400$GSCAD),
             table4(beta1.AIC.400$GMCP,     ASE1.AIC.400$GMCP),
             table4(beta1.AIC.400$AGLASSO,  ASE1.AIC.400$AGLASSO))
tmp2 = rbind(table4(beta1.AIC.800$GLASSO,   ASE1.AIC.800$GLASSO),
             table4(beta1.AIC.800$GSCAD,    ASE1.AIC.800$GSCAD),
             table4(beta1.AIC.800$GMCP,     ASE1.AIC.800$GMCP),
             table4(beta1.AIC.800$AGLASSO,  ASE1.AIC.800$AGLASSO))
tmp3 = rbind(table4(beta1.AIC.1200$GLASSO,  ASE1.AIC.1200$GLASSO),
             table4(beta1.AIC.1200$GSCAD,   ASE1.AIC.1200$GSCAD),
             table4(beta1.AIC.1200$GMCP,    ASE1.AIC.1200$GMCP),
             table4(beta1.AIC.1200$AGLASSO, ASE1.AIC.1200$AGLASSO))

tmp4 = rbind(table4(beta2.AIC.400$GLASSO,   ASE2.AIC.400$GLASSO),
             table4(beta2.AIC.400$GSCAD,    ASE2.AIC.400$GSCAD),
             table4(beta2.AIC.400$GMCP,     ASE2.AIC.400$GMCP),
             table4(beta2.AIC.400$AGLASSO,  ASE2.AIC.400$AGLASSO))
tmp5 = rbind(table4(beta2.AIC.800$GLASSO,   ASE2.AIC.800$GLASSO),
             table4(beta2.AIC.800$GSCAD,    ASE2.AIC.800$GSCAD),
             table4(beta2.AIC.800$GMCP,     ASE2.AIC.800$GMCP),
             table4(beta2.AIC.800$AGLASSO,  ASE2.AIC.800$AGLASSO))
tmp6 = rbind(table4(beta2.AIC.1200$GLASSO,  ASE2.AIC.1200$GLASSO),
             table4(beta2.AIC.1200$GSCAD,   ASE2.AIC.1200$GSCAD),
             table4(beta2.AIC.1200$GMCP,    ASE2.AIC.1200$GMCP),
             table4(beta2.AIC.1200$AGLASSO, ASE2.AIC.1200$AGLASSO))
table41.AIC = rbind(tmp1.400, tmp1, line, tmp1.800, tmp2, line, tmp1.1200, tmp3)
table42.AIC = rbind(tmp2.400, tmp4, line, tmp2.800, tmp5, line, tmp2.1200, tmp6)
table4.AIC = cbind(table42.AIC, table41.AIC)


# ###############
# write.csv(table4.CV, "D:/Google ¶³ºÝµwºÐ/program/yang/group penalized regression 2013/group cox/©T©wweightver2/simulation all/table4.cox.CV.csv")
# write.csv(table4.GCV, "D:/Google ¶³ºÝµwºÐ/program/yang/group penalized regression 2013/group cox/©T©wweightver2/simulation all/table4.cox.GCV.csv")
# write.csv(table4.BIC, "D:/Google ¶³ºÝµwºÐ/program/yang/group penalized regression 2013/group cox/©T©wweightver2/simulation all/table4.cox.BIC.csv")
# write.csv(table4.AIC, "D:/Google ¶³ºÝµwºÐ/program/yang/group penalized regression 2013/group cox/©T©wweightver2/simulation all/table4.cox.AIC.csv")
