library(mvtnorm)
library(BiocParallel)
source("FunctionsM.R")
source("Matrix_ComputationM.R")



d2opt <- function(data_1, data_2, K = 100000, lambda = 0.01,c = 0.15, weight.vec = rep(1/6, 6), delta.min = 0.12){
        dat <- rbind(data_1, data_2) ##complete data
        
        ##estimators considered
        est.four <- NULL
        for(i in 1:4){
                est.four[i] <- mean(dat$Y[dat$A==1 & dat$Risk<=i/4 & dat$Risk>(i-1)/4])-mean(dat$Y[dat$A==0 & dat$Risk<=i/4 & dat$Risk>(i-1)/4])        
        }
        
        ##sample size matrices
        interval.consid <- matrix(0,2,4)
        interval.consid[1,] <- c(0,0.25,0.5,0.75)
        interval.consid[2,] <- c(0.25,0.5,0.75,1)
        
        samplesize.matrix <- matrix(0,2,4)
        for(i in 1:2){
                for(j in 1:4){
                        samplesize.matrix[i,j] <- sum(dat$Risk<=interval.consid[2,j] & dat$Risk>interval.consid[1,j] & dat$A==i-1)        
                }
        }
        
        ###estimated sigma
        phat <- sum(dat$Y==1)/nrow(dat)
        sigma <- sqrt(phat*(1-phat))
        
        ####generate vcov matrix
        diag.vec <- sigma^2/samplesize.matrix[1,] + sigma^2/samplesize.matrix[2,]
        vcov <- diag(diag.vec)
        
        ##true mean treatment effect under 6 settings considered 
        mat.returned <- matrix.return(c = c, alpha = 0.05, beta = 0.2, delta.min = delta.min, sigma = sigma)
        mean.dat <- t(mat.returned$mean.dat)
        ##calculate posterior weights
        dmv.val <- NULL
        for(i in 1:ncol(mean.dat)){
                dmv.val[i] <-  dmvnorm(est.four,mean=mean.dat[,i],sigma=vcov,log=T) + log(weight.vec[i])       
        }
        
        diff <- dmv.val[1:5]-dmv.val[6]
        ratio.prob <- exp(diff)
        ratio.prob <- c(ratio.prob,1)
        prob.vec <- ratio.prob/sum(ratio.prob)
        
        ##calculate posterior expectation
        sample.ind <- sample(1:6, K , prob = prob.vec, replace = T)
        loss.matrix <- matrix(0, K, ncol = 5)
        integral.matrix <- mat.returned$integral.matrix
        duration <- mat.returned$duration
        power <- mat.returned$power
        
        ###utility function computation is below  "loss" represent the utility here####
        for(i in 1:5){
                for(j in 1:length(sample.ind)){
                        loss.matrix[j,i] <- integral.matrix[sample.ind[j],i] * power[sample.ind[j],i] - lambda*duration[i]       
                }
        }
        loss.expectation <- apply(loss.matrix,2,mean)
        action <- which.max(loss.expectation)
        return(list(action.recom=action)) 
}

d1opt <- function(seed = 314, data_1, mu_0 = 0.125, K1 = 100, K = 1000, lambda = 0.01, c = 0.15, weight.vec = rep(1/6,6), delta.min = 0.2, sigma = 3, n = 64){
        dat <- data_1
        
        stage1.loss.act <- function(act,seed = 314){
                set.seed(act*100 + seed)
                tmp.loss.mat <- NULL
                num.pieces <- round((ideal.pop.prime[2,act]-ideal.pop.prime[1,act])/0.25)
                for(i in seq_along(sample.ind)){
                        
                        risk.score <- matrix(0,n/num.pieces,num.pieces)
                        for(j in 1:num.pieces){
                                risk.score[,j] <- runif(n/num.pieces,(ideal.pop.prime[1,act]+0.25*(j-1)),(ideal.pop.prime[1,act]+0.25*j))
                        }
                        risk.score <- as.vector(risk.score)
                        
                        treatment.ind <- matrix(0,n/num.pieces,num.pieces)
                        for(j in 1:num.pieces){
                                treatment.ind[sample(1:(n/num.pieces),n/(num.pieces*2),replace = F),j] <- 1
                        }
                        treatment.ind <- as.vector(treatment.ind)
                        tmp.model <- sample.ind[i]
                        mu_r <- mu_0 + get(paste0("b",tmp.model))(risk.score)
                        data_2 <- data.frame(Risk=risk.score,mu_0=mu_0,mu_r=mu_r,A=treatment.ind)
                        
                        ##generate response by risk score
                        data_2$Y <- NA
                        data_2$Y[data_2$A==0] <- rbinom(sum(data_2$A==0),1,mu_0)
                        data_2$Y[data_2$A==1] <- rbinom(sum(data_2$A==1),1,data_2$mu_r[data_2$A==1])
                        data_2 <- data_2[,c(1,4,5)]
                        action <- d2opt(data_1 = dat,data_2 = data_2, K = K, lambda = lambda, c = c, weight.vec = weight.vec, 
                                        delta.min = delta.min)$action.recom        
                        tmp.loss.mat[i] <- integral.matrix[sample.ind[i], action] * power[sample.ind[i], action] - lambda * duration[action]  
                        
                }
                return(tmp.loss.mat)
        }
        
        ###determine the posterior weights b,c|X1
        est.four.stage1 <- NULL
        for(i in 1:4){
                est.four.stage1[i] <- mean(dat$Y[dat$A==1 & dat$Risk<=i/4 & dat$Risk>(i-1)/4])-mean(dat$Y[dat$A==0 & dat$Risk<=i/4 & dat$Risk>(i-1)/4]) 
        }
        
        
        ##sample size matrices
        interval.consid <- matrix(0, 2, 4)
        interval.consid[1, ] <- c(0, 0.25, 0.5, 0.75)
        interval.consid[2, ] <- c(0.25, 0.5, 0.75, 1)
        
        samplesize.matrix <- matrix(0, 2, 4)
        for(i in 1:2){
                for(j in 1:4){
                        samplesize.matrix[i,j] <- sum(dat$Risk<=interval.consid[2,j] & dat$Risk>interval.consid[1,j] & dat$A==i-1)        
                }
        }
        
        ###estimated sigma
        phat <- sum(dat$Y==1)/nrow(dat)
        sigma <- sqrt(phat*(1-phat))
        
        ####generate vcov matrix
        diag.vec <- sigma^2/samplesize.matrix[1, ] + sigma^2/samplesize.matrix[2, ]
        vcov <- diag(diag.vec)
        matrix.return.list <- matrix.return(c = c, alpha = 0.05, beta = 0.2, delta.min = delta.min, sigma = sigma)
        mean.dat <- t(matrix.return.list$mean.dat)
        dmv.val <- NULL
        for(i in 1:ncol(mean.dat)){
                dmv.val[i] <-  dmvnorm(est.four.stage1,mean=mean.dat[,i],sigma=vcov,log=T) + log(weight.vec[i])       
        }
        
        diff <- dmv.val[1:5]-dmv.val[6]
        ratio.prob <- exp(diff)
        ratio.prob <- c(ratio.prob,1)
        prob.vec <- ratio.prob/sum(ratio.prob)
        
        sample.ind <- sample(1:6, K1, prob = prob.vec, replace = T)
        integral.matrix <- matrix.return.list$integral.matrix
        power <- matrix.return.list$power
        duration <- matrix.return.list$duration
        ideal.pop.prime <- matrix.return.list$ideal.pop.prime
        ideal.pop.prime <- ideal.pop.prime[,-5]
        tmp.loss <- bplapply(1:4, stage1.loss.act, seed = seed)
        loss <- matrix(unlist(tmp.loss), nrow = K1)
        loss.1.expectation <- apply(loss,2,mean)
        if(diff(range(loss.1.expectation)) == 0 & loss.1.expectation[1] == 0){
                action_1 <- 4
        }else{
                action_1 <- which.max(loss.1.expectation)
        }
        return(action_1)
        
}

choice.population.two.stage <- function(data.simul, seed = 314, n = 64, K1 = 100, K = 1000, mu_0 = 0.125, delta.min = 0.12, lambda = 0.01, c = 0.15, weight.vec = c(rep(0.1,5), 0.5)){
        set.seed(seed)
        
        phat <- sum(data.simul$Y==1)/nrow(data.simul)
        sigma <- sqrt(phat*(1-phat))
        
        data.list <- list()
        data.sampled <- list()
        index.list <- list()
        data1.list <- list()
        data2.list <- list()
        for(i in 1:2){
                for(j in 1:4){
                        data.list[[(i-1)*4+j]] <- data.simul[data.simul$Risk<=j/4 & data.simul$Risk>(j-1)/4 & data.simul$A==(i-1),]
                        M <- nrow(data.list[[(i-1)*4+j]])
                        data.sampled[[(i-1)*4+j]] <- data.list[[(i-1)*4+j]][sample(1:M, 16, replace = T),]
                        index.list[[(i-1)*4+j]] <- sample(1:16,8)
                        data1.list[[(i-1)*4+j]] <- data.sampled[[(i-1)*4+j]][index.list[[(i-1)*4+j]],]
                        data2.list[[(i-1)*4+j]] <- data.sampled[[(i-1)*4+j]][-index.list[[(i-1)*4+j]],]
                }
        }
        dat1 <- do.call(rbind,data1.list)
        dat2 <- do.call(rbind,data2.list)
        
        action_1 <- d1opt(data_1 = dat1, n = n, mu_0 = mu_0, K1 = K1, K = K,
                          lambda = lambda, c = c, weight.vec = weight.vec, seed = seed, delta.min = delta.min)
        
        set.seed(seed*10)
        ideal.pop.prime <- matrix.return(c = c, alpha = 0.05, beta = 0.2, delta.min = delta.min, sigma = sigma)$ideal.pop.prime
        ideal.pop.prime <- ideal.pop.prime[,-5]
        ##continue to enroll action_1 patient populaiton###
        data.list.2 <- list()
        data_2.list <- list()
        num.pieces <- round((ideal.pop.prime[2,action_1]-ideal.pop.prime[1,action_1])/0.25)
        for(i in 1:2){
                for(j in 1:num.pieces){
                        data.list.2[[(i-1)*4+j]] <- dat2[dat2$Risk<=(ideal.pop.prime[1,action_1] + 0.25*j) & 
                                                                 dat2$Risk>(ideal.pop.prime[1,action_1] + 0.25*(j-1))& dat2$A==(i-1),]
                        index.list[[(i-1)*4+j]] <- sample(1:8,ifelse((num.pieces == 2), 16, 8), replace = T)
                        data_2.list[[(i-1)*4+j]] <- data.list.2[[(i-1)*4+j]][index.list[[(i-1)*4+j]],]
                }
        }
        data_2 <- do.call(rbind, data_2.list)
        
        action <- d2opt(data_1 = dat1, data_2 = data_2, K = K, lambda = lambda, c = c, weight.vec = weight.vec, 
                        delta.min = delta.min)$action.recom 
        return(list(action_1=action_1,action_2=action))
        
}



###bootstrap 200 times for second stage and compute the proportion########
###data.mistie is the MISTIE trial data that contains three columns#######
###column 1 Risk --- Pre-randomization ICH volume#########################
###column 2 A --- Treatment Indicator#####################################
###column 3 Y --- The Primary Outcome (binary in the MISTIE)##############
set.seed(321)
Boot.no <- 200
result.dat <- matrix(0, Boot.no, 2)
for(i in 1 : Boot.no){
        seed = i + 321
        result.list <- choice.population.two.stage(data.simul = data.mistie, seed = seed, lambda = 0.01, c = 0.1, weight.vec = c(rep(0.1,5), 0.5))
        result.dat[i, ] <- as.numeric(unlist(result.list))
}

res <- NULL
for(i in 1:5){
        res[i] <- sum(result.dat[,2]==i)/Boot.no        
}
for(i in 6:9){
        res[i] <- sum(result.dat[,1]==(i-5))/Boot.no
}
res
