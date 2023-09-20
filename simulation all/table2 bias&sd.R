#load("E:/Google ¶³ºÝµwºÐ/program/yang/group penalized regression 2013/group cox/©T©wweightver2/simulation all/OK/betaresult.cox.RData")
table2 = function(beta1, beta2, beta.t, posi){
  tmp = numeric()
  m1 = apply(beta1, 1, mean)-beta.t
  m2 = apply(beta2, 1, mean)-beta.t
  tmp = c(m1[posi], m2[posi])
  tmp = round(tmp, 3)*1000
  return(tmp)
}
beta.true2 = beta.true[posi]


r0.400.full = table2(cox2.400, cox1.400, beta.true, posi = posi)
r0.800.full = table2(cox2.800, cox1.800, beta.true, posi = posi)
r0.1200.full = table2(cox2.1200, cox1.1200, beta.true, posi = posi)
r0.400 = table2(cox2.t.400, cox1.t.400, beta.true2, posi = 1:5)
r0.800 = table2(cox2.t.800, cox1.t.800, beta.true2, posi = 1:5)
r0.1200 = table2(cox2.t.1200, cox1.t.1200, beta.true2, posi = 1:5)

line = rep(NA, 10)
## CV
tmp1 = rbind(r0.400.full,
             table2(beta2.CV.400$GLASSO, beta1.CV.400$GLASSO, beta.true, posi),
             table2(beta2.CV.400$GSCAD, beta1.CV.400$GSCAD, beta.true, posi),
             table2(beta2.CV.400$GMCP, beta1.CV.400$GMCP, beta.true, posi),
             table2(beta2.CV.400$AGLASSO, beta1.CV.400$AGLASSO, beta.true, posi))

tmp2 = rbind(r0.800.full,
             table2(beta2.CV.800$GLASSO, beta1.CV.800$GLASSO, beta.true, posi),
             table2(beta2.CV.800$GSCAD, beta1.CV.800$GSCAD, beta.true, posi),
             table2(beta2.CV.800$GMCP, beta1.CV.800$GMCP, beta.true, posi),
             table2(beta2.CV.800$AGLASSO, beta1.CV.800$AGLASSO, beta.true, posi))

tmp3 = rbind(r0.1200.full,
             table2(beta2.CV.1200$GLASSO, beta1.CV.1200$GLASSO, beta.true, posi),
             table2(beta2.CV.1200$GSCAD, beta1.CV.1200$GSCAD, beta.true, posi),
             table2(beta2.CV.1200$GMCP, beta1.CV.1200$GMCP, beta.true, posi),
             table2(beta2.CV.1200$AGLASSO, beta1.CV.1200$AGLASSO, beta.true, posi))

table2.CV = rbind(tmp1,
                    line,
                    tmp2,
                    line,
                    tmp3)

## GCV
tmp1 = rbind(r0.400.full,
             table2(beta2.GCV.400$GLASSO, beta1.GCV.400$GLASSO, beta.true, posi),
             table2(beta2.GCV.400$GSCAD, beta1.GCV.400$GSCAD, beta.true, posi),
             table2(beta2.GCV.400$GMCP, beta1.GCV.400$GMCP, beta.true, posi),
             table2(beta2.GCV.400$AGLASSO,beta1.GCV.400$AGLASSO, beta.true, posi))

tmp2 = rbind(r0.800.full,
             table2(beta2.GCV.800$GLASSO, beta1.GCV.800$GLASSO, beta.true, posi),
             table2(beta2.GCV.800$GSCAD, beta1.GCV.800$GSCAD, beta.true, posi),
             table2(beta2.GCV.800$GMCP, beta1.GCV.800$GMCP, beta.true, posi),
             table2(beta2.GCV.800$AGLASSO, beta1.GCV.800$AGLASSO, beta.true, posi))

tmp3 = rbind(r0.1200.full,
             table2(beta2.GCV.1200$GLASSO, beta1.GCV.1200$GLASSO, beta.true, posi),
             table2(beta2.GCV.1200$GSCAD, beta1.GCV.1200$GSCAD, beta.true, posi),
             table2(beta2.GCV.1200$GMCP, beta1.GCV.1200$GMCP, beta.true, posi),
             table2(beta2.GCV.1200$AGLASSO, beta1.GCV.1200$AGLASSO, beta.true, posi))

table2.GCV = rbind(tmp1,
                   line,
                   tmp2,
                   line,
                   tmp3)

## BIC
tmp1 = rbind(r0.400.full,
             table2(beta2.BIC.400$GLASSO, beta1.BIC.400$GLASSO, beta.true, posi),
             table2(beta2.BIC.400$GSCAD, beta1.BIC.400$GSCAD, beta.true, posi),
             table2(beta2.BIC.400$GMCP, beta1.BIC.400$GMCP, beta.true, posi),
             table2(beta2.BIC.400$AGLASSO, beta1.BIC.400$AGLASSO, beta.true, posi))

tmp2 = rbind(r0.800.full,
             table2(beta2.BIC.800$GLASSO, beta1.BIC.800$GLASSO, beta.true, posi),
             table2(beta2.BIC.800$GSCAD, beta1.BIC.800$GSCAD, beta.true, posi),
             table2(beta2.BIC.800$GMCP, beta1.BIC.800$GMCP, beta.true, posi),
             table2(beta2.BIC.800$AGLASSO, beta1.BIC.800$AGLASSO, beta.true, posi))

tmp3 = rbind(r0.1200.full,
             table2(beta2.BIC.1200$GLASSO, beta1.BIC.1200$GLASSO, beta.true, posi),
             table2(beta2.BIC.1200$GSCAD, beta1.BIC.1200$GSCAD, beta.true, posi),
             table2(beta2.BIC.1200$GMCP, beta1.BIC.1200$GMCP, beta.true, posi),
             table2(beta2.BIC.1200$AGLASSO, beta1.BIC.1200$AGLASSO, beta.true, posi))

table2.BIC = rbind(tmp1,
                   line,
                   tmp2,
                   line,
                   tmp3)
## AIC
tmp1 = rbind(r0.400.full,
             table2(beta2.AIC.400$GLASSO, beta1.AIC.400$GLASSO, beta.true, posi),
             table2(beta2.AIC.400$GSCAD, beta1.AIC.400$GSCAD, beta.true, posi),
             table2(beta2.AIC.400$GMCP, beta1.AIC.400$GMCP, beta.true, posi),
             table2(beta2.AIC.400$AGLASSO, beta1.AIC.400$AGLASSO, beta.true, posi))

tmp2 = rbind(r0.800.full,
             table2(beta2.AIC.800$GLASSO, beta1.AIC.800$GLASSO, beta.true, posi),
             table2(beta2.AIC.800$GSCAD, beta1.AIC.800$GSCAD, beta.true, posi),
             table2(beta2.AIC.800$GMCP, beta1.AIC.800$GMCP, beta.true, posi),
             table2(beta2.AIC.800$AGLASSO, beta1.AIC.800$AGLASSO, beta.true, posi))

tmp3 = rbind(r0.1200.full,
             table2(beta2.AIC.1200$GLASSO, beta1.AIC.1200$GLASSO, beta.true, posi),
             table2(beta2.AIC.1200$GSCAD, beta1.AIC.1200$GSCAD, beta.true, posi),
             table2(beta2.AIC.1200$GMCP, beta1.AIC.1200$GMCP, beta.true, posi),
             table2(beta2.AIC.1200$AGLASSO, beta1.AIC.1200$AGLASSO, beta.true, posi))

table2.AIC = rbind(tmp1,
                   line,
                   tmp2,
                   line,
                   tmp3)
####
table2.BIC
# write.csv(table2.CV, "E:/Google ¶³ºÝµwºÐ/program/yang/group penalized regression 2013/group cox/©T©wweightver2/simulation all/table2.cox.CV.csv")
# write.csv(table2.GCV, "E:/Google ¶³ºÝµwºÐ/program/yang/group penalized regression 2013/group cox/©T©wweightver2/simulation all/table2.cox.GCV.csv")
# write.csv(table2.BIC, "E:/Google ¶³ºÝµwºÐ/program/yang/group penalized regression 2013/group cox/©T©wweightver2/simulation all/table2.cox.BIC.csv")
# write.csv(table2.AIC, "E:/Google ¶³ºÝµwºÐ/program/yang/group penalized regression 2013/group cox/©T©wweightver2/simulation all/table2.cox.AIC.csv")
