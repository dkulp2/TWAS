info.matrix <- function(y.w, r.w, r.wAll, weights, weights.all, geno.w, geno.wAll, n.fx, n.tx, n.snps, tx.max, beta.hat, gamma.hat, sige.hat, phi.hat, eta.hat) {
  # first term of eq. 19 is eq 20-29
  # score variables is U in paper

  y.tilde <- y.w - rep(1,length(weights))*beta.hat-r.w*gamma.hat
  genoB.w <- cbind(rep(1,length(weights)),geno.w)
  mu <- exp(genoB.w%*%eta.hat)

  # score1-score5 are eq 8-12

  score1 <- rep(0,n.fx)
  for (i in 1:n.fx){
    tal1 <- 1 + (i-1)*tx.max
    tal2 <- tal1+tx.max - 1
    score1[i] <- sum(weights[tal1:tal2]*(r.w[tal1:tal2]*y.tilde[tal1:tal2]/sige.hat^2))
  }

  score2 <- rep(0,n.fx)
  for (i in 1:n.fx){
    tal1 <- 1 + (i-1)*tx.max
    tal2 <- tal1+tx.max - 1
    score2[i] <- sum(weights[tal1:tal2]*(y.tilde[tal1:tal2]/sige.hat^2))
  }


  score3 <- rep(0,n.fx)
  for (i in 1:n.fx){
    tal1 <- 1 + (i-1)*tx.max
    tal2 <- tal1+tx.max - 1
    score3[i] <- sum(weights[tal1:tal2]*(-1/(2*sige.hat^2) + (y.tilde[tal1:tal2]^2)/(2*sige.hat^4)))
  }

  score4 <- matrix(0,nrow=n.fx,ncol=n.snps+1)
  for (i in 1:n.fx){
    tal1 <- 1 + (i-1)*tx.max
    tal2 <- tal1+tx.max - 1
    score4[i,] <- apply(as.vector(weights[tal1:tal2]*as.vector(r.w[tal1:tal2]-mu[tal1:tal2])/as.vector(1+mu[tal1:tal2]/phi.hat))*genoB.w[tal1:tal2,],2,sum)
  }

  score5 <- rep(0,n.fx)
  for (i in 1:n.fx){
    tal1 <- 1 + (i-1)*tx.max
    tal2 <- tal1+tx.max - 1
    score5[i] <- sum(weights[tal1:tal2]*(digamma(r.w[tal1:tal2]+phi.hat) - digamma(phi.hat) -(r.w[tal1:tal2]+phi.hat)/(mu[tal1:tal2]+phi.hat) - log(mu[tal1:tal2]/phi.hat +1) + 1))
  }

  # 5x5 matrix + eta terms. each element is either a variance or covariance (off diags)
  m1.1 <- outer(score1,score1)
  info1.1 <- sum(weights*(r.w^2/sige.hat^2 - r.w^2*y.tilde^2/sige.hat^4)) - 2*sum(m1.1[upper.tri(m1.1)])


  mat1.2=outer(score1,score2)
  diag(mat1.2)<-0
  info1.2 <- sum(weights*(r.w/sige.hat^2 - r.w*y.tilde^2/sige.hat^4)) - sum(mat1.2)


  mat1.3=outer(score1,score3)
  diag(mat1.3)<-0
  info1.3 <- sum(weights*(r.w*y.tilde/sige.hat^4)*(1+ 1/2 - y.tilde^2/(2*sige.hat^2))) - sum(mat1.3)

  info1.4a <- -apply(as.vector(weights*r.w*y.tilde*as.vector(r.w-mu)/as.vector(sige.hat^2*(1+mu/phi.hat)))*cbind(rep(1,length(weights)),geno.w),2,sum)

  info1.4b <- rep(0,n.snps+1)
  for (i in 1:n.fx){
    for (j in 1:n.fx){
      info1.4b <- info1.4b + score1[i]*score4[j,]*(i!=j)
    }
  }

  info1.4 <- info1.4a - info1.4b

  mat1.5=outer(score1,score5)
  diag(mat1.5)<-0
  info1.5 <- -sum(weights*as.vector(r.w*y.tilde/sige.hat^2)*(digamma(r.w+phi.hat) - digamma(phi.hat) -(r.w+phi.hat)/(mu+phi.hat) - log(mu/phi.hat +1) + 1))  - sum(mat1.5)


  m2.2 <- outer(score2,score2)
  info2.2 <- sum(weights*(1/sige.hat^2 - y.tilde^2/sige.hat^4)) - 2*sum(m2.2[upper.tri(m2.2)])

  mat2.3=outer(score2,score3)
  diag(mat2.3)<-0
  info2.3 <- sum(weights*(y.tilde/sige.hat^4)*(1+ 1/2 - y.tilde^2/(2*sige.hat^2))) - sum(mat2.3)

  info2.4a <- -apply(as.vector(weights*y.tilde*as.vector(r.w-mu)/as.vector(sige.hat^2*(1+mu/phi.hat)))*cbind(rep(1,length(weights)),geno.w),2,sum)

  info2.4b <- rep(0,n.snps+1)
  for (i in 1:n.fx){
    for (j in 1:n.fx){
      info2.4b <- info2.4b + score2[i]*score4[j,]*(i!=j)
    }
  }

  info2.4 <- info2.4a - info2.4b

  mat2.5=outer(score2,score5)
  diag(mat2.5)<-0
  info2.5 <- -sum(weights*as.vector(y.tilde/sige.hat^2)*(digamma(r.w+phi.hat) - digamma(phi.hat) -(r.w+phi.hat)/(mu+phi.hat) - log(mu/phi.hat +1) + 1)) - sum(mat2.5)

  m3.3 <- outer(score3,score3)
  info3.3 <- sum(weights*(y.tilde^2/sige.hat^6 - 1/(2*sige.hat^4) -(-1/(2*sige.hat^2) + (y.tilde^2)/(2*sige.hat^4))^2)) - 2*sum(m3.3[upper.tri(m3.3)])

  info3.4a <- -apply(as.vector(weights*(-1/(2*sige.hat^2) + y.tilde^2/(2*sige.hat^4))*as.vector(r.w-mu)/as.vector(1+mu/phi.hat))*cbind(rep(1,length(weights)),geno.w),2,sum)

  info3.4b <- rep(0,n.snps+1)
  for (i in 1:n.fx){
    for (j in 1:n.fx){
      info3.4b <- info3.4b + score3[i]*score4[j,]*(i!=j)
    }
  }

  info3.4 <- info3.4a - info3.4b

  mat3.5=outer(score3,score5)
  diag(mat3.5)<-0
  info3.5 <- -sum(weights*(-1/(2*sige.hat^2) + as.vector(y.tilde^2/(2*sige.hat^4)))*(digamma(r.w+phi.hat) - digamma(phi.hat) -(r.w+phi.hat)/(mu+phi.hat) - log(mu/phi.hat +1) + 1)) - sum(mat3.5)


  #####################################################
  ### ALL DATA CALCULATION

  genoB.w <- cbind(rep(1,length(weights.all)),geno.wAll)
  mu <- exp(genoB.w%*%eta.hat)

  score4 <- matrix(0,nrow=n.fx+n.tx,ncol=n.snps+1)

  for (i in 1:(n.fx)){
    tal1 <- 1 + (i-1)*tx.max
    tal2 <- tal1+tx.max - 1
    score4[i,] <- apply(as.vector(weights.all[tal1:tal2]*as.vector(r.wAll[tal1:tal2]-mu[tal1:tal2])/as.vector(1+mu[tal1:tal2]/phi.hat))*genoB.w[tal1:tal2,],2,sum)
  }
  for (i in (n.fx+1):(n.fx+n.tx)){
    tal1 <- n.fx*tx.max + i - n.fx
    tal2 <- n.fx*tx.max + i - n.fx
    score4[i,] <- as.vector(weights.all[tal1:tal2]*as.vector(r.wAll[tal1:tal2]-mu[tal1:tal2])/as.vector(1+mu[tal1:tal2]/phi.hat))*genoB.w[tal1:tal2,]
  }
  
  score5 <- rep(0,n.fx+n.tx)
  for (i in 1:(n.fx)){
    tal1 <- 1 + (i-1)*tx.max
    tal2 <- tal1+tx.max - 1
    score5[i] <- sum(weights.all[tal1:tal2]*(digamma(r.wAll[tal1:tal2]+phi.hat) - digamma(phi.hat) -(r.wAll[tal1:tal2]+phi.hat)/(mu[tal1:tal2]+phi.hat) - log(mu[tal1:tal2]/phi.hat +1) + 1))
  }
  for (i in (n.fx+1):(n.fx+n.tx)){
    tal1 <- n.fx*tx.max + i - n.fx
    tal2 <- n.fx*tx.max + i - n.fx
    score5[i] <- sum(weights.all[tal1:tal2]*(digamma(r.wAll[tal1:tal2]+phi.hat) - digamma(phi.hat) -(r.wAll[tal1:tal2]+phi.hat)/(mu[tal1:tal2]+phi.hat) - log(mu[tal1:tal2]/phi.hat +1) + 1))
  }

  info4.4a <- matrix(0,nrow=n.snps+1,ncol=n.snps+1)
  for (i in seq(1,length(weights.all))){
    temp<- weights.all[i]*outer(genoB.w[i,],genoB.w[i,])*(phi.hat*mu[i]*(r.wAll[i]+phi.hat)/((phi.hat+mu[i])^2) - ((r.wAll[i]-mu[i])/(1+mu[i]/phi.hat))^2)
    info4.4a <- info4.4a + temp
  }

  info4.4b <- matrix(0,nrow=n.snps+1,ncol=n.snps+1)
  for (i in seq(1,n.fx+n.tx)){
    for (j in seq(1,n.fx+n.tx)){
      info4.4b <- info4.4b + outer(score4[i,],score4[j,])*(i!=j)
    }
  }
  info4.4 <- info4.4a-info4.4b


  info4.5a <- apply((as.vector(weights.all*(r.wAll-mu)*(-mu/(phi.hat+mu)^2 - (digamma(r.wAll+phi.hat) - digamma(phi.hat) -(r.wAll+phi.hat)/(mu+phi.hat) - log(mu/phi.hat +1) + 1)/(1+mu/phi.hat)))*genoB.w),2,sum)

  info4.5b <- rep(0,n.snps+1)
  for (i in seq(1,n.fx+n.tx)){
    for (j in seq(1,n.fx+n.tx)){
      info4.5b <- info4.5b + score4[j,]*score5[i]*(i!=j)
    }
  }

  info4.5 <- info4.5a - info4.5b

  m5.5 <- outer(score5,score5)
  info5.5 <- sum((weights.all*(-trigamma(r.wAll+phi.hat)+trigamma(phi.hat)-(mu^2+r.wAll*phi.hat)/(phi.hat*(mu+phi.hat)^2) - (digamma(r.wAll+phi.hat) - digamma(phi.hat) -(r.wAll+phi.hat)/(mu+phi.hat) - log(mu/phi.hat +1) + 1)^2))) - 2*sum(m5.5[upper.tri(m5.5)])


  #####################################################


  Info <- matrix(NA,nrow=5+n.snps,ncol=5+n.snps)
  Info[1,1:3] <- c(info1.1,info1.2,info1.3)
  Info[1:3,1] <- c(info1.1,info1.2,info1.3)
  Info[2:2,2:2]<- info2.2
  Info[3,1:3]<-c(info1.3,info2.3,info3.3)
  Info[1:3,3]<-c(info1.3,info2.3,info3.3)

  Info[4:(4+n.snps),4:(4+n.snps)] <- info4.4

  Info[1,4:(5+n.snps)] <- c(info1.4,info1.5)
  Info[4:(5+n.snps),1] <- c(info1.4,info1.5)
  Info[2,4:(5+n.snps)] <- c(info2.4,info2.5)
  Info[4:(5+n.snps),2] <- c(info2.4,info2.5)
  Info[3,4:(5+n.snps)] <- c(info3.4,info3.5)
  Info[4:(5+n.snps),3] <- c(info3.4,info3.5)
  
  Info[(5+n.snps),4:(5+n.snps)] <- c(info4.5,info5.5)
  Info[4:(5+n.snps),(5+n.snps)] <- c(info4.5,info5.5)

  return(Info)
}
