# path_root = dirname(rstudioapi::getActiveDocumentContext()$path); setwd(path_root)
load("datasetall.RData")

library(foreach)
library(doParallel)
library(survival)
cpu = makeCluster(6)
registerDoParallel(cpu)

GMCP.sel.400 = foreach( r=1:100,.packages="survival", .combine=cbind ) %dopar%{
  x1 = x2 = data400[,r]$x.m
  y1 = y2 = data400[,r]$y
  y1[,1] = 0
  group = c(rep(1:5, each = 2), 6:10)
  lambda.seq = from=NULL; to=0;  n.lambda=200; tol=1e-8
  penalty = "MCP"
  true.A = c(1, 2, 5, 6, 11)
  gamma = 0
  save1 = mygrpcox.path(method=1.2, x=x1, y=y1, group = group, gamma=gamma, from=from, to=to, n.lambda=n.lambda, tol = tol, penalty=penalty, true.A = true.A)
  save2 = mygrpcox.path(method=1.2, x=x2, y=y2, group = group, gamma=gamma, from=from, to=to, n.lambda=n.lambda, tol = tol, penalty=penalty, true.A = true.A)
  sel1 = select6(save1)
  sel2 = select6(save2)
  list(sel1=sel1, sel2=sel2) 
}

GMCP.sel.800 = foreach( r=1:100,.packages="survival", .combine=cbind ) %dopar%{
  x1 = x2 = data800[,r]$x.m
  y1 = y2 = data800[,r]$y
  y1[,1] = 0
  group = c(rep(1:5, each = 2), 6:10)
  lambda.seq = from=NULL; to=0;  n.lambda=200; tol=1e-8
  penalty = "MCP"
  true.A = c(1, 2, 5, 6, 11)
  gamma = 0
  save1 = mygrpcox.path(method=1.2, x=x1, y=y1, group = group, gamma=gamma, from=from, to=to, n.lambda=n.lambda, tol = tol, penalty=penalty, true.A = true.A)
  save2 = mygrpcox.path(method=1.2, x=x2, y=y2, group = group, gamma=gamma, from=from, to=to, n.lambda=n.lambda, tol = tol, penalty=penalty, true.A = true.A)
  sel1 = select6(save1)
  sel2 = select6(save2)
  list(sel1=sel1, sel2=sel2) 
}

GMCP.sel.1200 = foreach( r=1:100,.packages="survival", .combine=cbind ) %dopar%{
  x1 = x2 = data1200[,r]$x.m
  y1 = y2 = data1200[,r]$y
  y1[,1] = 0
  group = c(rep(1:5, each = 2), 6:10)
  lambda.seq = from=NULL; to=0;  n.lambda=200; tol=1e-8
  penalty = "MCP"
  true.A = c(1, 2, 5, 6, 11)
  gamma = 0
  save1 = mygrpcox.path(method=1.2, x=x1, y=y1, group = group, gamma=gamma, from=from, to=to, n.lambda=n.lambda, tol = tol, penalty=penalty, true.A = true.A)
  save2 = mygrpcox.path(method=1.2, x=x2, y=y2, group = group, gamma=gamma, from=from, to=to, n.lambda=n.lambda, tol = tol, penalty=penalty, true.A = true.A)
  sel1 = select6(save1)
  sel2 = select6(save2)
  list(sel1=sel1, sel2=sel2) 
}
stopCluster(cl = cpu)
save.image("GMCP.CV.RData")


