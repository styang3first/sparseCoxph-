load("E:/Google ¶³ºÝµwºÐ/program/yang/group penalized regression 2013/group cox/©T©wweightver2/simulation all/OK/allresult.RData")

beta1.AIC.400 = list(
  GLASSO = apply(aglasso0.sel.400, 2, function(x) x$sel1$AIC$beta),
#   GSCAD = apply(result.gSCAD.sel400, 2, function(x) x$sel1$AIC$beta),
#   GMCP = apply(result.gMCP.sel400, 2, function(x) x$sel1$AIC$beta),
  AGLASSO1 = apply(aglasso1.sel.400, 2, function(x) x$sel1$AIC$beta),
  AGLASSO2 = apply(aglasso2.sel.400, 2, function(x) x$sel1$AIC$beta)
)

beta1.AIC.800 = list(
  GLASSO = apply(result.glasso.sel800, 2, function(x) x$sel1$AIC$beta),
  GSCAD = apply(result.gSCAD.sel800, 2, function(x) x$sel1$AIC$beta),
  GMCP = apply(result.gMCP.sel800, 2, function(x) x$sel1$AIC$beta),
  AGLASSO = apply(result.aglasso.sel800, 2, function(x) x$sel1$AIC$beta)
)
beta1.AIC.1200 = list(
  GLASSO = apply(result.glasso.sel1200, 2, function(x) x$sel1$AIC$beta),
  GSCAD = apply(result.gSCAD.sel1200, 2, function(x) x$sel1$AIC$beta),
  GMCP = apply(result.gMCP.sel1200, 2, function(x) x$sel1$AIC$beta),
  AGLASSO = apply(result.aglasso.sel1200, 2, function(x) x$sel1$AIC$beta)
)