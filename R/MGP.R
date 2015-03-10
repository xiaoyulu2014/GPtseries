MGP = function(dat,N,Q,num_gen,lag,max,min,Pc,Pm,x_test,step){ 
 
  library(base)
#   library(foreach)
#   library(doParallel)
#   cl <- makeCluster(4)
#   registerDoParallel(cl)
  
  #step1
  x = matrix(,lag,N)
  for (i in 1:N){
    x[,i] = dat[(i+1):(i+lag)]
  }
  
  x_tilda = cbind(t(x),rep(1,N))
  
  #step-ahead prediction
  y = dat[(1+lag+step):(N+lag+step)]
  
  ##GA
  #form Q strings
  theta=matrix(,Q,3)
  str=list()
  for (i in 1:3) str[[i]]=matrix(rbinom(n=L[i]*Q,size=1,prob=0.5),Q,L[i])
  
  t=1
  while (t < num_gen){
  #decoding
    for (i in 1:3){
      for (j in 1:Q){
        R = strtoi(paste(str[[i]][j,],collapse=""),base=2)
        r = (log(max[i],base=10)-log(min[i],base=10))*R/(2^L[i]-1)+log(min[i],base=10)
        theta[j,i] = 10^r
      }
    }
    #covariance matrix and estimate weights and fitness
    K = array(,c(Q,N,N))
    w = matrix(,Q,lag+1)
    fitness = c()
  
    for (i in 1:Q){
      for (j in 1:N){
        for (k in 1:N){
          K[i,j,k] = theta[i,1]^2*exp(-sum((x[,j]-x[,k])^2)/(2*theta[i,2]^2))         
        }
      }
      
      K[i,,] = K[i,,] + theta[i,3]^2*diag(N)
      cK = chol(K[i,,])
      iK = chol2inv(cK)
      w[i,] = solve(t(x_tilda)%*%iK%*%x_tilda)%*%t(x_tilda)%*%iK%*%y
      fitness[i] = -(0.5*2*sum(log(diag(cK))) + 0.5*t(y-x_tilda%*%w[i,])%*%iK%*%(y-x_tilda%*%w[i,]) + N/2*log(2*pi))
    }
      
    fitness = fitness + abs(min(fitness))
    
    #reproduction
    S = cbind(str[[1]],str[[2]],str[[3]])
    index = sample(1:Q,size=Q,replace=T,prob=fitness/sum(fitness))
    parent = sample(index)[1:2]
    if (runif(1)<Pc){
      crosspt = sample(1:sum(L))[1]
      if (crosspt == sum(L)) crosspt = crosspt-1
      tmp = S[parent[1],(crosspt+1):sum(L)]
      S[parent[1],] = c(S[parent[1],1:crosspt],S[parent[2],(crosspt+1):sum(L)])
      S[parent[2],] = c(S[parent[2],1:crosspt],tmp)
    }
    
    indexmu = sample(1:Q)[1:(Pm*Q)]
    for (j in indexmu){
      mut = sample(1:sum(L))[1]
      S[j,mut] = 1 - S[j,mut]
    }   
    
    str[[1]] = S[,1:L[1]];str[[2]]=S[,(L[1]+1):(L[1]+L[2])];str[[3]]=S[,(L[1]+L[2]+1):sum(L)]
    t = t+1
  }  
    #suboptimal problem
  best = which.max(fitness)
  theta_best = theta[best,]
  w_best = w[best,]
  K_best = K[best,,]
  mean_best = x_tilda%*%w_best
  
  
  #prediction
  mean_pos = cbind(t(x_test),rep(1,N)) %*% w_best
  Sigma_pos = matrix(,ncol(x_test),N)
  for (i in 1:ncol(x_test)){
    for (j in 1:N){
      Sigma_pos[i,j] = theta_best[1]^2*exp(-sum((x_test[,i]-x[,j])^2)/(2*theta_best[2]^2))
    }
  }
  y_hat = mean_pos + Sigma_pos%*%chol2inv(chol(K_best))%*%(y-mean_best)
  return(y_hat)
}















