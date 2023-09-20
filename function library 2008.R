library(survival)
library(MASS)
data.generate = function(N, Sigma, b0, t1, t2, t3, t4, c.max=10){
  # This function generates left-truncated and right-censored survival data with
  p = 10
  x.o = trunc.o = failure.o = numeric()
  n=5000
  repeat{
    if( length(trunc.o)<N ){
      tmp = mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma, tol = 1e-15)
      a = qnorm(1/3, 0, 1); b = qnorm(2/3, 0, 1);
      posi0 = which(tmp[,1:5] <= a)
      posi1 = which(tmp[,1:5] > b)
      x = tmp; x[,1:5]=2
      x[posi0]=0; x[posi1]=1
      for(j in 6:10){ 
        tmp1 = x[,j]
        x[,j] = sign(tmp1) * ifelse(abs(tmp1)>2, 2, abs(tmp1))
      }
      #beta=c(-1, 1, 0, 0, 1, 1, 0, 0, 1, 1, rep(0, 10))
      tmp =  0.3*(x[,1]!=2) + -0.35*(x[,3]!=2) + 0.2*x[,6]        
      
      tmp2 = b0 + t1*tmp +
        t2*(x[,2]==2) + t2*(x[,4]==2) + t2*(x[,5]==2) +
        -t4*(x[,7]) + t4*(x[,8])  + 
        -t4*(x[,9]) + t4*(x[,10])
      #t3*(x[,6]>0) + t4*(x[,7]>0) + t4*(x[,8]>0)  + t4*(x[,9]>0)  + t4*(x[,10]>0)
      
      
      trunc = (-log(runif(n = n, 0, 1)) * exp(-tmp2))^0.5
      failure = (-log(runif(n = n, 0, 1)) * exp(-tmp))^0.5
      
      index = which(trunc < failure)
      
      trunc.o = c(trunc.o, trunc[index])
      failure.o = c(failure.o, failure[index])
      x.o = rbind(x.o, x[index,])
      tr = 1-length(index)/n
      cat(tr, ' ')
    } else if(length(trunc.o)>N){
      failure.o = failure.o[1:N]
      trunc.o = trunc.o[1:N]
      x.o = x.o[1:N,]
      break;
    }
  }
  if(c.max == Inf) c.o = rep(Inf, N) else c.o = trunc.o + runif(N, 0, c.max)
  status.o = ifelse(failure.o < c.o, 1, 0)
  time.o = ifelse(failure.o < c.o, failure.o, c.o)
  
  x.o = as.data.frame(x.o)
  for(j in 1:5) x.o[,j] = as.factor(x.o[,j])
  for(j in 1:5){ contrasts(x.o[,j]) = contr.treatment(n=3, base = 3)}
  
  x.m = model.matrix(~V1+V2+V3+V4+V5+V6+V7+V8+V9+V10, data = x.o)
  x.m = x.m[,-1]
  return(list(x.m=x.m, y=cbind(trunc.o, time.o, status.o), censored.rate=1-sum(status.o)/N, truncated.rate=tr))
}

data.generate2 = function(N, Sigma, b0, t1, t2, t3, t4, c.max=10){
  # This function generates left-truncated and right-censored survival data with
  # trunc.o = truncation time
  # status.o = idicator for right censored
  p = ncol(Sigma)
  x.o = trunc.o = failure.o = numeric()
  n=5000
  repeat{
    if( length(trunc.o)<N ){
      tmp = mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma, tol = 1e-15)
      a = qnorm(1/3, 0, 1); b = qnorm(2/3, 0, 1);
      posi0 = which(tmp[,1:5] <= a)
      posi1 = which(tmp[,1:5] > b)
      x = tmp; x[,1:5]=2
      x[posi0]=0; x[posi1]=1
      for(j in 6:10){ 
        tmp1 = x[,j]
        x[,j] = sign(tmp1) * ifelse(abs(tmp1)>2, 2, abs(tmp1))
      }
      #beta=c(-1, 1, 0, 0, 1, 1, 0, 0, 1, 1, rep(0, 10))
      tmp =  1*(x[,1]!=2) + 0.7*x[,6]        
      
      tmp2 = b0 + t1*tmp +
        t2*(x[,2]==2) - t2*(x[,3]==2) + t2*(x[,4]==2) - t2*(x[,5]==2) +
        t4*(x[,7]) - t4*(x[,8])  + t4*(x[,9]) - t4*(x[,10])
      trunc = (-log(runif(n = n, 0, 1)) * exp(-tmp2))^0.5
      failure = (-log(runif(n = n, 0, 1)) * exp(-tmp))^0.5
      
      index = which(trunc < failure)
      
      trunc.o = c(trunc.o, trunc[index])
      failure.o = c(failure.o, failure[index])
      x.o = rbind(x.o, x[index,])
      tr = 1-length(index)/n
      cat(tr, ' ')
    } else if(length(trunc.o)>=N){
      failure.o = failure.o[1:N]
      trunc.o = trunc.o[1:N]
      x.o = x.o[1:N,]
      break;
    }
  }
  if(c.max == Inf) c.o = rep(Inf, N) else c.o = trunc.o + runif(N, 0, c.max)
  status.o = ifelse(failure.o < c.o, 1, 0)
  time.o = ifelse(failure.o < c.o, failure.o, c.o)
  
  x.o = as.data.frame(x.o)
  for(j in 1:5) x.o[,j] = as.factor(x.o[,j])
  for(j in 1:5){ contrasts(x.o[,j]) = contr.treatment(n=3, base = 3) }
  
  x.m = model.matrix(~V1+V2+V3+V4+V5+V6+V7+V8+V9+V10, data = x.o)
  x.m = x.m[,-1]
  return(list(x.m=x.m, y=cbind(trunc.o, time.o, status.o), censored.rate=1-sum(status.o)/N, truncated.rate=tr))
}

stand = function(x){
  # This function standardize x-variables
  n = nrow(x); p = ncol(x)
  center = apply(x, 2, mean)
  x.tmp = x - matrix(1, n, 1)%*%center
  scale = sqrt(apply( (x.tmp)^2, 2, mean))
  x.norm = x.tmp%*%diag(1/scale)
  return(list(x.norm, center, scale))
}

RC.list=function(y){
  # This function computes list of censorship w.r.t. ordered survival time
  # because of left-truncation, RC.list does not have nested structure
  Y = as.matrix(y)
  trunc = as.numeric(Y[,1])
  time = as.numeric(Y[,2])
  status = as.numeric(Y[,3])
  od = order(time)
  n = nrow(y)
  R.list = lapply(1:n, function(i) which(trunc<=time[i] & time[i]<=time))
  C.list=list(); tmp = matrix(0, n, n)
  for(i in 1:n) tmp[i,R.list[[i]]] = 1
  for(i in 1:n) C.list[[i]]=which(tmp[,i]==1)
  return(list(R.list=R.list, C.list=C.list, RC.m=tmp))
}

working = function(x, Y, beta, family=c("gaussian","cox"), RC.m){
  # This function computed "working data" for generalized linear model (using "link function")
  family = family[1];
  n = nrow(x);
  if(family == "gaussian"){
    eta = Y
    w = rep(1, n)
    res = Y - x %*% beta
  }else if (family == "cox"){
    status = Y[,3]
    eta = as.numeric(x %*% beta)
    e.eta = exp(eta)
    tmp1 = (RC.m%*%e.eta)[,1]
    tmp2 = I(status==1)/tmp1
    tmp3 = e.eta * (t(RC.m)%*%tmp2)[,1]
    l. = tmp3 - status
    w = tmp3 - e.eta^2 * (t(RC.m)%*%tmp2^2)[,1]
    res = -l./w; res[which(w==0)] = 0
  }
  return(list(eta = eta, w = w, res=res))
}

dp = function(para, lambda, gamma = 3.7, penalty = c("lasso", "SCAD", "MCP")){
  # This function computes (sub)gradients of penalty functions
  if( sum(para<0)!=0 ) stop("beta should be positive")
  penalty = penalty[1];  output = para;
  po = which(para!=0); para.diff = para[po];
  
  output[po] = switch(penalty,
                 "lasso" = rep(lambda, length(para.diff)),
                 "SCAD" = {
                   tmp = gamma*lambda-para.diff;
                   tmp = ifelse(tmp>0, tmp, 0)
                   ifelse(para.diff <= lambda, lambda, tmp/(gamma-1))
                   },
                 "MCP" = ifelse(para.diff <= lambda*gamma, lambda - para.diff/gamma, 0) )
  output[-po]=NA
  return(output)
}

S = function(z, lambda, gamma = 3.7, penalty = c("lasso", "SCAD", "MCP")){
  # This function computes penalty value
  penalty = penalty[1]; 
  
  z.abs = abs(z)
  if(penalty!="SCAD"){
    tmp = ifelse(z.abs > lambda, z.abs-lambda, 0)
    output = switch(penalty,
                    "lasso" = tmp,
                    "MCP" = ifelse(z.abs <= gamma*lambda, tmp/(1-1/gamma), z.abs))
  } else if(penalty == "SCAD"){
    output = z.abs
    a = which(z.abs <=2*lambda)
    b = setdiff(which(z.abs <= gamma*lambda), a)
    
    output[a] = ifelse(z.abs[a] > lambda, z.abs[a]-lambda, 0)
    lambda. = gamma*lambda/(gamma-1)
    output[b] = ifelse(z.abs[b] > lambda., z.abs[b]-lambda., 0)/(1-1/(gamma-1))
  } else {stop("no such penalty")}
  
  return(sign(z) * output)
}

grpcd1 = function(V, VTV, u, n, beta.ini, lambdaj, penalty, group.list, tol=1e-7){
  # This algorithm uses local quadratic approximation (LQA) to approximate 
  # partial likelihood of grouped-variables then optimized them. 
  # Notice that LQA neither gives exactly optimal solution, nor sparse estimation.
  # For non-grouped variables, the coordinate descent algorithm gives exactly optimal solution.
  J = length(lambdaj)
  beta.old = beta.new = beta.ini
  iter = 0
  res = u - V %*% beta.ini
  posi = 1:J
  repeat{
    for(j in posi){
      iter=iter+1
      gj = group.list[[j]]
      pj = length(gj)
      vj = as.matrix(V[, gj])
      betaj. = betaj = beta.new[gj]
      
      resj = res + vj%*%betaj
      vj_resj = (t(vj) %*% resj ) / n        
      vj_vj = VTV[gj, gj] / n
      vj_resj.norm = sqrt(sum(vj_resj^2))
      if(pj==1){ # coordinated descent: exactly optimal solution
        if(penalty == "lasso"){
          betaj. = S(z=vj_resj.norm, lambda=lambdaj[j], penalty=penalty)*sign(vj_resj[1,1])*n/(t(vj)%*%vj)
        } else {
          zj = (t(vj) %*%(res))/n + betaj
          zj.norm = sqrt(sum(zj^2))
          betaj. = S(z=zj.norm, lambda=lambdaj[j], penalty=penalty)*sign(zj)
        }
        betaj. = as.numeric(betaj.)
      } else { # grouped LQA descent: not exactly optimal solution
        if(vj_resj.norm<lambdaj[j]+tol){
          betaj. = rep(0, pj)
        } else {
          if(all(betaj==0)) betaj. = vj_resj
          betaj.old = betaj.
          repeat{
            betaj.norm = sqrt(sum(betaj.old^2))
            lambda.LQA = dp(betaj.norm, lambda=lambdaj[j], penalty=penalty )
            betaj. = solve(vj_vj + diag(lambda.LQA/betaj.norm,  pj))%*%vj_resj
            if(all(abs(betaj.-betaj.old)<tol)) break;
            iter=iter+1
            betaj.old=betaj.
          }
        }
      }                
      beta.diff = betaj. - betaj
      if( any(beta.diff!=0) ){
        beta.new[group.list[[j]]] = betaj.
        res = res - vj%*%beta.diff
      }
    }
    if( all( abs(beta.old-beta.new)<tol ) ) break;
    beta.old = beta.new
  }
  return(list(beta.hat = beta.new, iter=iter))
}

grpcd1.2 = function(V, VTV, u, n, beta.ini, beta.mle, lambdaj, penalty, group.list, tol=1e-7){
  # This algorithm uses local linear approximation (LLA) to approximate 
  # partial likelihood of all variables then optimized them. 
  # Although LLA does not give exactly optimal solution, it gives sparse estimation.
  J = length(lambdaj)
  beta.old = beta.new = beta.ini
  
  iter = 0
  res = u - V %*% beta.ini
  posi = 1:J
  repeat{
    for(j in posi){
      iter=iter+1
      gj = group.list[[j]]
      pj = length(gj)
      vj = as.matrix(V[, gj])
      betaj = beta.new[gj]
      vj_vj = VTV[gj, gj] / n
      bj = max(eigen(vj_vj)$value)+1e-6
      
      Uj = (t(vj) %*% res ) / n        
      Uj.norm = sqrt(sum(Uj^2))
      zj = Uj + bj*betaj
      zj.norm = sum(zj^2)^0.5
      betaj.norm = max(sqrt(sum(betaj^2)), tol)
      lambda.LLA = dp(betaj.norm, lambdaj[j], penalty = penalty)
      betaj. = S(zj.norm, lambda.LLA, penalty="lasso" ) * zj/zj.norm/bj
     
      beta.diff = betaj. - betaj
      if( any(beta.diff!=0) ){
        beta.new[group.list[[j]]] = betaj.
        res = res - vj%*%beta.diff
      }
    }
    if( all( abs(beta.old-beta.new)<tol ) ) break;
    beta.old = beta.new
  }
  return(list(beta.hat = beta.new, iter=iter))
}


grpcd2 = function(x.orthog, Y, beta.ini, family=c("gaussian", cox),
                  RC.m, lambdaj, penalty, group.list, tol=1e-7){
  # This function optimize group penalized partial likelihood when each grouped-x-variables
  # is orthogonalized.
  family = family[1]
  J = length(lambdaj); n = nrow(x.orthog); p = ncol(x.orthog)
  beta.old = beta.new = beta.ini
  iter = 0
  
  posi = 1:J
  repeat{
    work = working(x=x.orthog, Y=Y, beta=beta.new, family=family, RC.m=RC.m)
    w = work$w
    res = work$res
    for(j in posi){
      xj = as.matrix(x.orthog[, group.list[[j]]])
      betaj = beta.new[group.list[[j]]]
      b=1
      zj = (t(xj) %*% (res*w)) / n  + b * betaj
      zj.norm = sqrt(sum(zj^2))
      
      betaj. = S(z = zj.norm, lambda = lambdaj[j], penalty = penalty) * zj/zj.norm/b
      betaj. = as.numeric(betaj.)
      
      beta.diff = betaj. - betaj
      if( any(beta.diff!=0) ){
        beta.new[group.list[[j]]] = betaj.
        res = res - xj%*%beta.diff
      }
    }
    if( all( abs(beta.old-beta.new)<tol ) ) break;
    iter = iter+1
    beta.old = beta.new
  } 
  return(list(beta.hat = beta.new, iter=iter))
}

mybind=function(list1, glist, n, p){
  J = length(list1)
  tmp = matrix(0, n, p)
  for(j in 1:J) tmp[,glist[[j]]] = list1[[j]]
  return(tmp)  
}
mygrpcox.path = function(x, y, group, gamma, true.A=NULL,
                    n.lambda = 100, tol = 1e-10, 
                    from = NULL, to = 1e-2,
                    lambda.seq = NULL, 
                    penalty = c("lasso", "SCAD", "MCP"), method = 1){
  # This function optimizes group penalized partial likelihood (coxph model) 
  # for survival data along with increasing lambda.seq.
  # Penalty function can be 1. lasso, 2. SCAD, 3. MCP
  cd.tol = tol
  penalty = penalty[1]
  family = "cox"
  x = as.matrix(x)
  n = nrow(x); p = ncol(x)
  
  Y=as.matrix(y)
  trunc = Y[,1]; time = Y[,2]; status = Y[,3]
  tmp = RC.list(y)
  RC.m = tmp$RC.m
  rm(tmp)
  
  group.table = table(group)
  J = length(group.table)
  group.type = as.numeric(names(group.table))
  group.list = lapply(1:J, function(j) which(group==group.type[j]))
  
  x.trans = stand(x)
  x.norm = x.trans[[1]]
  center = x.trans[[2]]
  scale = x.trans[[3]]
  
  eigen.list = lapply(1:J, function(j){
    tmp = x.norm[, group.list[[j]]]
    return( eigen( t(tmp) %*% tmp / n  ) )
  })
  x.orthogonal.list = lapply(1:J, function(j){
    tmp = eigen.list[[j]]
    return( x.norm[, group.list[[j]]] %*% tmp$vectors %*% diag( 1/sqrt(tmp$values), group.table[j] ) )
  })
  x.orthog = mybind(x.orthogonal.list, group.list, n, p)  
  colnames(x.orthog) = colnames(x)
  
  cox.save = coxph(Surv(trunc, time, status)~x.orthog)
  beta.mle.stand = cox.save$coef
  beta.mle.origin = coxph(Surv(trunc, time, status)~x)$coef
  if(!is.null(true.A)){
    beta.mle.stand.true = coxph(Surv(trunc, time, status)~x.orthog[,true.A])$coef
    beta.mle.origin.true = coxph(Surv(trunc, time, status)~x[,true.A])$coef
  } else beta.mle.stand.true = beta.mle.origin.true = NULL
  wj = sqrt(group.table) * sqrt( sapply(1:J, function(j) sum(beta.mle.stand[group.list[[j]]]^2)) )^(-gamma)
  
  VTV = solve(cox.save$var) 
  V=chol(VTV)
  u = V%*%beta.mle.stand
  
  if(is.null(lambda.seq)){
    if(is.null(from)){
      lambda.max = sapply(1:J, function(j){
        vj = as.matrix(V[,group.list[[j]]])
        a = t(vj)%*%u
        return( sqrt(sum(a^2))/n/wj[j] )
      })
      from = max( lambda.max )
    }
    to = from*to
    lambda.seq = seq( from = from, to = to, by = -(from-to)/(n.lambda-1))
  }
  n.lambda = length(lambda.seq)
  
  iter = result = lambdaj.seq = numeric(0)
  beta.ini = rep(0, p)
  
  t1 = proc.time()
  for(lambda in lambda.seq){
    lambdaj = lambda * wj
    lambdaj.seq = cbind(lambdaj.seq, lambdaj)
    
    if(method == 1){
      grpcd.info = grpcd1(V=V, VTV=VTV, u=u, n=n, beta.ini=beta.ini,
                         lambdaj=lambdaj, penalty=penalty, 
                         group.list=group.list, tol=cd.tol)
    } else if(method == 1.2){
      grpcd.info = grpcd1.2(V=V, VTV=VTV, u=u, n=n, beta.ini=beta.ini, beta.mle=beta.mle.stand,
                            lambdaj=lambdaj, penalty=penalty, 
                            group.list=group.list, tol=cd.tol)
    }   else if(method == 2){
      grpcd.info = grpcd2(x.orthog = x.orthog, Y=Y, beta.ini=beta.ini,
                          family=family, RC.m=RC.m,
                          lambdaj = lambdaj, penalty = penalty, 
                          group.list=group.list, tol=cd.tol)
    }   
    result = cbind(result, grpcd.info$beta.hat)
    beta.ini = grpcd.info$beta.hat
    iter = c(iter, grpcd.info$iter)
    cat(grpcd.info$iter, " ")
  }
  t2 = proc.time()
  print(paste0("time cost: ", round((t2-t1)[1], 3), " sec" ))
  
  result[abs(result)<tol]=0
  result2 = result
  for(j in 1:J){
    tmp = eigen.list[[j]]
    beta.tmp = tmp$vectors %*% diag( 1/sqrt(tmp$values), group.table[j] ) %*% result[ group.list[[j]], ]
    result2[ group.list[[j]], ] = beta.tmp
  }
  result3 = (diag(1/scale)%*%result2)
  colnames(result) = colnames(result3) = 1:n.lambda
  
  output = list(x=x, x.orthog=x.orthog, y=y, wj=wj, uV=list(u=u, V=V),
                group=group, group.list=group.list,
                beta.mle = cbind(mle.origin = beta.mle.origin, mle.stand = beta.mle.stand),
                beta.mle.true = cbind(mle.origin.true = beta.mle.origin.true, mle.stand.true = beta.mle.stand.true),
                betas = result,
                beta = result3, 
                lambda.seq = lambda.seq, lambdaj.seq=lambdaj.seq,
                iter=iter, penalty=penalty, gamma=gamma, family=family)  
  return(output)  
}  

path.plot = function(save, true=NULL, main=NULL, a=NULL, b=NULL){
  beta.path = save$beta
  J=length(save$w)
  group = save$group
  group.table = table(group)
  group.type = as.numeric(names(group.table))
  group.list = lapply(1:J, function(j) which(group==group.type[j]))
  if(is.null(a) ) a = max((beta.path))
  if(is.null(b) ) b = min((beta.path))
  
  plot(0, type="n", xlim=c(0, max(save$lambda.seq)), ylim=c(b, a), xlab=expression(lambda), ylab="coefficient", main=main)
  
  if(is.null(true)){
    for(j in 1:J){
      for(j2 in group.list[[j]]){
        points(save$lambda.seq, beta.path[j2,], col=j, "l")
      }
    }
  } else {
    for(j in setdiff(1:J, true)){
      for(j2 in group.list[[j]]){
        points(save$lambda.seq, beta.path[j2,], col=1, "l")
      }
    }
    for(j in 1:length(true)){
      for(j2 in group.list[[true[j]]]){
        points(save$lambda.seq, beta.path[j2,], col=(j+1), "l")
      }
    }
  }
  abline(h=0, cex=3)
}

part.likeli = function(x, y, beta, RC.m, lossonly=F){
  x = as.matrix(x)
  trunc = as.numeric(y[,1])
  time = as.numeric(y[,2])
  status = as.numeric(y[,3])
  n = nrow(x); p = ncol(x)
  exb=as.numeric(exp(x%*%beta))
  posi = which(status==1)
  R.list2 = t(RC.m[posi,])
  s0 = (t(R.list2)%*%exb)[,1]
  
  if(lossonly){
    if(p==1){
      lbeta = sum(x[posi]*beta - log(s0))
      l.beta = NULL
      l..beta = NULL
    } else {
      lbeta = sum((x[posi,]%*%beta - log( s0 )))
      l.beta = NULL
      l..beta = NULL
    }
    return( lbeta=lbeta )
  } else {
    s1 = t(x)%*%diag(exb)%*%R.list2
    s2.tmp = array(0, dim=c(p, p, n))
    for(i in 1:n) s2.tmp[,,i] = exb[i]*x[i,]%*%t(x[i,])
    s22.tmp = sapply(1:n, function(i) as.numeric(s2.tmp[,,i]))
    s2 = s22.tmp %*% R.list2
    if(p==1){
      lbeta = sum(x[posi]*beta - log(s0))
      l.beta = sum(x[posi])-sum(s1/s0)
      l..beta = sum(s1^2/s0^2) - s2%*%(1/s0)
    } else {
      lbeta = sum((x[posi,]%*%beta - log( s0 )))
      l.beta = apply(x[posi,], 2, sum) - apply(s1%*%diag(1/s0), 1, sum)  
      l..beta = s1%*%diag(1/s0^2)%*%t(s1) - matrix(s2%*%(1/s0), p, p)
    }
    return( list(lbeta=lbeta, l.beta=l.beta, l..beta=as.matrix(l..beta)) )
  }
}

select3 = function(save, path.plot=F, value.plot=F){
  criterion = c("GCV", "BIC", "AIC")
  x.orthog = as.matrix(save$x.orthog)
  y = as.matrix(save$y)
  RC.m = RC.list(y)$RC.m
  
  betas = save$betas
  lambdaj.seq = save$lambdaj.seq
  lambda.seq = save$lambda.seq
  group = save$group
  penalty = save$penalty
  n = nrow(x.orthog); p = ncol(x.orthog); R = ncol(betas)
  
  group.table = table(group)
  J = length(group.table)
  group.type = as.numeric(names(group.table))
  group.list = lapply(1:J, function(j) which(group==group.type[j]))
  
  
  cri = sapply(1:R, function(r){
    beta = betas[,r]; cat(r)
    eps = 1e-9
    
    A = which(sapply(1:J, function(j) any(beta[group.list[[j]]]!=0) ))
    Ac = setdiff(1:J, A)
    for(j in Ac) beta[group.list[[j]]] = eps*sqrt(1/length(group.list[[j]]))
    
    part = part.likeli(x=x.orthog, y=y, RC.m=RC.m, beta=beta)
    Q = -part$l..beta/n; 
    M = array(0, dim=dim(Q))
    for(j in 1:J){
      gj=group.list[[j]]
      betaj = beta[gj]
      betaj.norm = sqrt(sum(betaj^2))
      dp.lambdaj = dp(para=betaj.norm, lambda=lambdaj.seq[j, r], penalty=penalty)
      a =  dp.lambdaj/betaj.norm * ( diag(1, length(betaj)) - betaj%*%t(betaj)/betaj.norm^2)
      M[gj, gj] = M[gj, gj] + a
    } 
    df = sum(diag(Q%*%solve(Q+M)))
    #cat((t2-t1)[1], " ", (t3-t2)[1], " ", (t4-t3)[1], " \n")
    return(list(loss=part$lbeta, df=df))
  })
  
  value = numeric()
  for(name in criterion){
    value.tmp = switch(name,
                       "GCV" = apply(cri, 2, function(x) -x$loss/(1 - x$df/n)),
                       "AIC" = apply(cri, 2, function(x) -x$loss/n + x$df/n),
                       "BIC" = apply(cri, 2, function(x) -x$loss/n + x$df*log(n)/n) )  
    value = rbind(value, value.tmp)
  }
  rownames(value) = criterion
  
  min = apply(value, 1, which.min)
  if(path.plot){ path.plot(save); abline(v=lambda.seq[min])}
  #if(value.plot){plot(lambda.seq, value[1,]); abline(v=lambda.seq[min])}
  AIC = BIC = GCV = NULL
  for(i in 1:length(criterion)){
    assign(criterion[i], list(lambda=lambda.seq[min[i]], betas=betas[,min[i]], beta=save$beta[,min[i]] ))
  }
  output = list(info=save, cri=cri, value=value, min=min, AIC=AIC, BIC=BIC, GCV=GCV )
  return(output)
}

select4 = function(save, path.plot=F, value.plot=F){
  criterion = c("GCV", "BIC", "AIC")
  x.orthog = as.matrix(save$x.orthog)
  y = as.matrix(save$y)
  RC.m = RC.list(y)$RC.m
  
  betas = save$betas
  lambdaj.seq = save$lambdaj.seq
  lambda.seq = save$lambda.seq
  group = save$group
  penalty = save$penalty
  n = nrow(x.orthog); p = ncol(x.orthog); R = ncol(betas)
  
  group.table = table(group)
  J = length(group.table)
  group.type = as.numeric(names(group.table))
  group.list = lapply(1:J, function(j) which(group==group.type[j]))
  
  cri = sapply(1:R, function(r){
    beta = betas[,r]; cat(r)
    eps = 1e-9
    
    A = which(sapply(1:J, function(j) any(beta[group.list[[j]]]!=0) ))
    Ac = setdiff(1:J, A)
    for(j in Ac) beta[group.list[[j]]] = eps*sqrt(1/length(group.list[[j]]))
    
    part = part.likeli(x=x.orthog, y=y, RC.m=RC.m, beta=beta)
    Q = -part$l..beta/n; 
    M = array(0, dim=dim(Q))
    for(j in 1:J){
      gj=group.list[[j]]
      betaj = beta[gj]
      if(length(gj)==1){
        dp.lambdaj = dp(para=abs(betaj), lambda=lambdaj.seq[j, r], penalty=penalty)
        M[gj, gj] = dp.lambdaj/abs(betaj)
      } else {
        betaj.norm = sqrt(sum(betaj^2))
        dp.lambdaj = dp(para=betaj.norm, lambda=lambdaj.seq[j, r], penalty=penalty)
        a =  dp.lambdaj/betaj.norm * ( diag(1, length(betaj)) - betaj%*%t(betaj)/betaj.norm^2)
        M[gj, gj] = M[gj, gj] + a
      }
    } 
    df = sum(diag(Q%*%solve(Q+M)))
    #cat((t2-t1)[1], " ", (t3-t2)[1], " ", (t4-t3)[1], " \n")
    return(list(loss=-part$lbeta, df=df))
  })
  
  value = numeric()
  for(name in criterion){
    value.tmp = switch(name,
                       "GCV" = apply(cri, 2, function(x) x$loss/(1 - x$df/n)),
                       "AIC" = apply(cri, 2, function(x) x$loss/n + x$df/n),
                       "BIC" = apply(cri, 2, function(x) x$loss/n + x$df*log(n)/n) )  
    value = rbind(value, value.tmp)
  }
  rownames(value) = criterion
  
  min = apply(value, 1, which.min)
  if(path.plot){ path.plot(save); abline(v=lambda.seq[min])}
  #if(value.plot){plot(lambda.seq, value[1,]); abline(v=lambda.seq[min])}
  AIC = BIC = GCV = NULL
  for(i in 1:length(criterion)){
    assign(criterion[i], list(lambda=lambda.seq[min[i]], betas=betas[,min[i]], beta=save$beta[,min[i]] ))
  }
  output = list(info=save, cri=cri, value=value, min=min, AIC=AIC, BIC=BIC, GCV=GCV )
  return(output)
}

select5 = function(save, path.plot=F, value.plot=F){
  criterion = c("GCV", "BIC", "AIC")
  x.orthog = as.matrix(save$x.orthog)
  y = as.matrix(save$y)
  RC.m = RC.list(y)$RC.m
  
  betas = save$betas
  lambdaj.seq = save$lambdaj.seq
  lambda.seq = save$lambda.seq
  group = save$group
  penalty = save$penalty
  n = nrow(x.orthog); p = ncol(x.orthog); R = ncol(betas)
  
  group.table = table(group)
  J = length(group.table)
  group.type = as.numeric(names(group.table))
  group.list = lapply(1:J, function(j) which(group==group.type[j]))
  
  cri = sapply(1:R, function(r){
    beta = betas[,r]; cat(r)
    eps = 1e-9
    
    A = which(sapply(1:J, function(j) any(beta[group.list[[j]]]!=0) ))
    Ac = setdiff(1:J, A)
    for(j in Ac) beta[group.list[[j]]] = eps*sqrt(1/length(group.list[[j]]))
    
    part = part.likeli(x=x.orthog, y=y, RC.m=RC.m, beta=beta)
    work = working(x=x.orthog, Y=y, beta=beta, RC.m=RC.m, family="cox")
    loss1 = -(as.numeric(t(work$res)%*%diag(work$w)%*%work$res))/2
    loss2 = part$lbeta*2
    Q = t(x.orthog)%*%diag(work$w)%*%x.orthog/n #-part$l..beta/n; 
    M = array(0, dim=dim(Q))
    for(j in 1:J){
      gj=group.list[[j]]
      betaj = beta[gj]
      if(length(gj)==1){
        dp.lambdaj = dp(para=abs(betaj), lambda=lambdaj.seq[j, r], penalty=penalty)
        M[gj, gj] = dp.lambdaj/abs(betaj)
      } else {
        betaj.norm = sqrt(sum(betaj^2))
        dp.lambdaj = dp(para=betaj.norm, lambda=lambdaj.seq[j, r], penalty=penalty)
        a =  dp.lambdaj/betaj.norm * ( diag(1, length(betaj)) - betaj%*%t(betaj)/betaj.norm^2)
        M[gj, gj] = M[gj, gj] + a
      }
    } 
    df = sum(diag(Q%*%solve(Q+M)))
    #cat((t2-t1)[1], " ", (t3-t2)[1], " ", (t4-t3)[1], " \n")
    return(list(loss=loss1, loss2=loss2, df=df))
  })
  
  value = numeric()
  for(name in criterion){
    value.tmp = switch(name,
                       "GCV" = apply(cri, 2, function(x) -x$loss/(1 - x$df/n)),
                       "AIC" = apply(cri, 2, function(x) -x$loss/n + x$df/n),
                       "BIC" = apply(cri, 2, function(x) -x$loss/n + x$df*log(n)/n) )  
    value = rbind(value, value.tmp)
  }
  rownames(value) = criterion
  
  min = apply(value, 1, which.min)
  if(path.plot){ path.plot(save); abline(v=lambda.seq[min])}
  AIC = BIC = GCV = NULL
  for(i in 1:length(criterion)){
    assign(criterion[i], list(lambda=lambda.seq[min[i]], betas=betas[,min[i]], beta=save$beta[,min[i]] ))
  }
  output = list(info=save, cri=cri, value=value, min=min, AIC=AIC, BIC=BIC, GCV=GCV )
  return(output)
}

select6 = function(save, path.plot=F, value.plot=F, loss=1, df=1, df.est=NULL, loss.est=NULL ){
  y = save$y
  RC.m = RC.list(y=y)$RC.m
  V = save$uV$V
  u = save$uV$u
  betas = save$betas
  lambdaj.seq = save$lambdaj.seq
  lambda.seq = save$lambda.seq
  group = save$group
  penalty = save$penalty
  x.orthog = save$x.orthog
  n = nrow(x.orthog); p = ncol(V); R = ncol(betas)
  n.V = nrow(V)
  group.table = table(group)
  J = length(group.table)
  group.type = as.numeric(names(group.table))
  group.list = lapply(1:J, function(j) which(group==group.type[j]))
  cox = coxph(Surv(y[,1], y[,2], y[,3])~x.orthog)
  R = ncol(betas)
  
  if(is.null(df.est)){
    df.est = sapply(1:R, function(r){
      beta. = beta = betas[,r]; # cat(r)
      eps = 1e-9
      
      A = which(sapply(1:J, function(j) any(beta[group.list[[j]]]!=0) ))
      Ac = setdiff(1:J, A)
      for(j in Ac) beta[group.list[[j]]] = eps*sqrt(1/length(group.list[[j]]))
      
      Q = t(V)%*%V/n
      M = array(0, dim=dim(Q))
      for(j in 1:J){
        gj = group.list[[j]]
        pj = length(gj)
        betaj = beta[gj]
        betaj. = beta.[gj]
        betaj.norm = sqrt(sum(betaj^2))
        if(pj==1){
          dp.lambdaj = dp(para=abs(betaj), lambda=lambdaj.seq[j, r], penalty=penalty)
          M[gj, gj] = dp.lambdaj/betaj.norm
        } else {
          dp.lambdaj = dp(para=betaj.norm, lambda=lambdaj.seq[j, r], penalty=penalty)
          if(df==1){
            a =  dp.lambdaj/betaj.norm * ( diag(1, pj) - betaj.%*%t(betaj.)/betaj.norm^2)
          } else if(df==2){
            a = diag(dp.lambdaj, pj)/betaj.norm         
          }
          M[gj, gj] = M[gj, gj] + a
        }
      } 
      H = (V%*% solve(Q+M) %*%t(V))/n #projection matrix
      df = sum(diag(H))
      if(df > p) df = p
      return(df=df)
    })
  }
  
  if(is.null(loss.est)){
    loss.est = apply(betas, 2, function(beta){
      if(loss==1){
        loss.est = -part.likeli(x=x.orthog, y=y, RC.m=RC.m, beta=beta, lossonly=T )
      } else if(loss==2){
        loss.est = sum((u-V%*%beta)^2)/2 - cox$loglik[2]
      }
      return(loss.est)
    })
  }
  
  value = rbind("GCV" = loss.est/(1-df.est/n)^2,
                "BIC" = 2*loss.est + df.est*log(n),
                "AIC" = 2*loss.est + 2*df.est)
  
  min = apply(value, 1, which.min)
  if(path.plot){ path.plot(save); abline(v=lambda.seq[min])}
  #if(value.plot){plot(lambda.seq, value[1,]); abline(v=lambda.seq[min])}
  GCV = list(min=min[1], value=value[1,], lambda=lambda.seq[min[1]], betas=betas[,min[1]], beta=save$beta[,min[1]] )
  BIC = list(min=min[2], value=value[2,], lambda=lambda.seq[min[2]], betas=betas[,min[2]], beta=save$beta[,min[2]] )
  AIC = list(min=min[3], value=value[3,], lambda=lambda.seq[min[3]], betas=betas[,min[3]], beta=save$beta[,min[3]] )
  
  output = list(info=save, df=df.est, loss=loss.est, AIC=AIC, BIC=BIC, GCV=GCV )
  return(output)
}

cv.mygrpcox.path = function(x, y, group, gamma, K=5, true.A=NULL,
                         n.lambda = 100, tol = 1e-8, 
                         from = NULL, to = 1e-2,
                         lambda.seq = NULL, 
                         penalty = c("lasso", "SCAD", "MCP"), method){
  # This function perform K-fold cross-validation
  RC.m = RC.list(y)$RC.m
  save.tmp = mygrpcox.path(method=method, x=x, y=y, group = group, gamma=gamma, from=from, to=to, n.lambda=n.lambda, tol = 1, penalty=penalty, true.A = true.A)
  lambda.seq = save.tmp$lambda.seq
  partition = matrix(sample(1:nrow(x), nrow(x), replace=F), nrow = K)
  
  model.err = matrix(0, K, n.lambda)
  #model.err2 = matrix(0, K, n.lambda)
  for(k in 1:K){
    x.train = x[-partition[k,],]
    y.train = y[-partition[k,],]
    RC.m.train=RC.list(y = y.train)$RC.m
    save.train = mygrpcox.path(method=method, x=x.train, y=y.train, group = group, gamma=gamma, 
                            lambda.seq = lambda.seq, tol = tol, penalty = penalty)
    beta.trains = save.train$betas
        
    l1.1 = apply(beta.trains, 2, function(x) part.likeli(x=save.tmp$x.orthog, y=y, beta=x, RC.m=RC.m, lossonly=T) )
    l2.1 = apply(beta.trains, 2, function(x) part.likeli(x=save.train$x.orthog, y=y.train, beta=x, RC.m=RC.m.train, lossonly=T) )
    model.err[k,] = l1.1-l2.1
  }  
  CV = apply(model.err, 2, sum)
  plot(lambda.seq, CV)
  lambda.select = lambda.seq[which.max(CV)]
  output = mygrpcox.path(method=method, x=x, y=y, group = group, gamma=gamma, lambda.seq=lambda.select, tol=tol, penalty=penalty, true.A = true.A)
  return(list(info=output, CV=CV, beta.select=output$beta, lambda.select=lambda.select))
}


kkt.check=function(save, r, roun=5){
  # This function checks whether the solution satisfies KKT conditions
  g = save$group.list
  beta=save$betas[,r]
  lambdaj=save$lambdaj.seq[,r]
  u = save$uV$u
  V = save$uV$V
  beta.mle = save$beta.mle[,2]
  n = nrow(save$x)
  l.beta2 = t(V)%*%V%*%(beta.mle - beta)
  
  value=numeric()
  for(j in 1:length(g)){
    gj = g[[j]]
    z.norm=(sum(beta[gj]^2))^0.5
    if( z.norm==0 ){ 
      value[gj] = lambdaj[j]
    } else {
      value[gj] = beta[gj]/z.norm * dp(para=z.norm, lambdaj[j], penalty=save$penalty)
    }
  }
  value=value*n
  return( round(cbind(l.beta2, value, beta), roun) )
  return(output)
  
}

ASE.est = function(object, type = NULL){
  # This function computes asymptotic standard error
  cat(1, " ")
  info = object$info
  group = info$group
  group.table = table(group)
  J = length(group.table)
  group.type = as.numeric(names(group.table))
  group.list = lapply(1:J, function(j) which(group==group.type[j]))
  x = info$x
  n = nrow(x); p =ncol(x)
  
  x.trans = stand(x)
  x.norm = x.trans[[1]]
  center = x.trans[[2]]
  scale = x.trans[[3]]
  
  eigen.list = lapply(1:J, function(j){
    tmp = x.norm[, group.list[[j]]]
    return( eigen( t(tmp) %*% tmp / n  ) )
  })
  u = info$uV$u
  V = info$uV$V
  penalty = info$penalty
  
  index = switch(type, 
                 "CV" = 1, 
                 "AIC" = object$AIC$min,
                 "GCV" = object$GCV$min,
                 "BIC" = object$BIC$min)
  
  betas = beta = info$betas[,index]
  lambdaj.seq = info$lambdaj.seq[,index]
  
  eps = 1e-9
  A = which(sapply(1:J, function(j) any(beta[group.list[[j]]]!=0) ))
  Ac = setdiff(1:J, A)
  for(j in Ac) beta[group.list[[j]]] = eps*sqrt(1/length(group.list[[j]]))
  
  Q = t(V)%*%V
  D = M = array(0, dim=dim(Q))
  for(j in 1:J){
    gj=group.list[[j]]
    betaj = beta[gj]
    betaj.norm = sqrt(sum(betaj^2))
    dp.lambdaj = dp(para=betaj.norm, lambda=lambdaj.seq[j], penalty=penalty)
    a =  dp.lambdaj/betaj.norm * ( diag(1, length(betaj)) )
    M[gj, gj] = M[gj, gj] + a
    if( all(betas[gj]!=0) ) D[gj, gj] = D[gj, gj] + a
  } 
  tmp1 = solve(Q + n*M)%*%(Q + n*D)
  ASES = tmp1 %*% solve(Q) %*% t(tmp1) 
  
  M=diag(1, p)
  for(j in 1:J){
    gj = group.list[[j]]
    if(length(gj)!=1){
      tmp = eigen.list[[j]]
      tmp2 = (tmp$vectors) %*% diag( 1/sqrt(tmp$values))
      M[ gj,gj ] = tmp2
    }
  }
  ASE = (M) %*% (ASES) %*% t(M)
  ASE = sqrt(diag(ASE))/scale
  return( ASE )
}

ASE.cox = function(object, type = NULL){
  # This function computes asymptotic standard error using coxph after variable selection
  cat(1, " ")
  info = object$info
  index = switch(type, 
                 "CV" = 1, 
                 "AIC" = object$AIC$min,
                 "GCV" = object$GCV$min,
                 "BIC" = object$BIC$min)
  x = info$x
  y = info$y
  beta = info$beta[,index]
  beta. = rep(0, ncol(x))
  A = which(beta!=0)
  if( sum(A) != 0){
    cox = coxph(Surv(y[,1], y[,2], y[,3])~x[,A])
    beta.[A] = cox$coef
  }
  ASE = sqrt(-diag(solve(part.likeli(x = x, y = y, RC.m = RC.list(y)$RC.m, beta=beta.)$l..beta)))
  return(ASE)
}
beta.cox = function(object, type = NULL){
  # This function computes beta estimation using coxph after variable selection
  cat(1, " ")
  info = object$info
  index = switch(type, 
                 "CV" = 1, 
                 "AIC" = object$AIC$min,
                 "GCV" = object$GCV$min,
                 "BIC" = object$BIC$min)
  x = info$x
  y = info$y
  beta = info$beta[,index]
  beta. = rep(0, ncol(x))
  A = which(beta!=0)
  if( sum(A) != 0){
    cox = coxph(Surv(y[,1], y[,2], y[,3])~x[,A])
    beta.[A] = cox$coef
  }
  return(beta.)
}  