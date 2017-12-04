library(mvtnorm)
library(BiocParallel)
source("Matrix_Computation.R")
source("Functions.R")


d2opt <- function(data_1, data_2, K = 100000,lambda = 1,c = 1, weight.vec = rep(1/6, 6), delta.min = 0.2, sigma = 3){
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
        contrib.part <- matrix(0, K, 4)
        integral.matrix <- mat.returned$integral.matrix
        duration <- mat.returned$duration
        power <- mat.returned$power
        # health.cost <- mat.returned$health.cost
        # integral.matrixo <- mat.returned$integral.matrixo
        # 
        ###loss function computation is below####
        for(j in 1:length(sample.ind)){
                for(i in 1:5){
                        loss.matrix[j,i] <- integral.matrix[sample.ind[j],i] * power[sample.ind[j],i] - lambda*duration[i]
                }
        }
        loss.expectation <- apply(loss.matrix,2,mean)
        action <- which.max(loss.expectation)  ###loss here refers to utility
        return(list(action.recom=action)) 
}


d1opt <- function(seed = 314, data_1, mu_0 = 3, K1 = 100, K = 1000, lambda = 1, c = 1, n = 264, weight.vec=rep(1/6,6), delta.min = 0.2, sigma = 3){
        dat <- data_1
        
        stage1.loss.act <- function(act){
                #set.seed(act*100 + seed)
                tmp.loss.mat <- NULL
                num.pieces <- round((ideal.pop.prime[2,act]-ideal.pop.prime[1,act])/0.25)
                for(i in seq_along(sample.ind)){
                        #print(i)
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
                        data_2$Y[data_2$A==0] <- rnorm(n/2,mu_0,sigma)
                        data_2$Y[data_2$A==1] <- rnorm(n/2,data_2$mu_r[data_2$A==1],sigma)
                        # risk.score <- runif(n,ideal.pop.prime[1,act],ideal.pop.prime[2,act])
                        # treatment.ind <- rbinom(n,1,0.5)
                        # tmp.model <- sample.ind[i]
                        # mu_r <- mu_0 + get(paste0("b",tmp.model))(risk.score)
                        # data_2 <- data.frame(Risk=risk.score,mu_0=mu_0,mu_r=mu_r,A=treatment.ind)
                        
                        ##generate response by risk score
                        #data_2$Y <- rnorm(n, mean = ifelse(data_2$A == 1, data_2$mu_r, mu_0), sd = 3)
                        data_2 <- data_2[, c(1,4,5)]
                        action <- d2opt(data_1 = dat,data_2 = data_2, K = K, lambda = lambda, c = c, weight.vec = weight.vec, 
                                        delta.min = delta.min, sigma = sigma)$action.recom        
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
        interval.consid <- matrix(0,2,4)
        interval.consid[1,] <- c(0,0.25,0.5,0.75)
        interval.consid[2,] <- c(0.25,0.5,0.75,1)
        
        samplesize.matrix <- matrix(0,2,4)
        for(i in 1:2){
                for(j in 1:4){
                        samplesize.matrix[i,j] <- sum(dat$Risk<=interval.consid[2,j] & dat$Risk>interval.consid[1,j] & dat$A==i-1)        
                }
        }
        
        ####generate vcov matrix
        diag.vec <- sigma^2/samplesize.matrix[1,] + sigma^2/samplesize.matrix[2,]
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
        #tmp.loss <- bplapply(1:4, stage1.loss.act, seed = seed)
        set.seed(seed)
        loss <- matrix(0, nrow = K1, ncol = 4)
        for(act in 1:4) {
                loss[, act] <- stage1.loss.act(act = act)
        }
        loss.1.expectation <- apply(loss,2,mean)
        action_1 <- which.max(loss.1.expectation) # "loss" here means utility
        return(action_1)
        
}


choice.population.two.stage <- function(seed = 314, n = 264, sigma = 3, K1 = 100, K = 1000, mu_0 = 3, lambda = 0.1, c = 0,
                                        weight.vec = rep(1/6,6), true.model = 1, delta.min = 0.3){        
        set.seed(seed)
        Risk <- matrix(0,n/4,4)
        for(j in 1:4){
                Risk[,j] <- runif(n/4,(j-1)/4,j/4)
        }
        Risk <- as.vector(Risk) 
        treatment.ind <- matrix(0,n/4,4)
        for(j in 1:4){
                treatment.ind[sample(1:(n/4),n/8,replace = F),j] <- 1
        }
        treatment.ind <- as.vector(treatment.ind)
        
        mu_r <- mu_0 + get(paste0("b",true.model))(Risk)
        dat1 <- data.frame(Risk=Risk,mu_0=mu_0,mu_r=mu_r,A=treatment.ind)
        
        ##generate response by risk score
        dat1$Y <- NA
        dat1$Y[dat1$A==0] <- rnorm(n/2,mu_0,sigma)
        dat1$Y[dat1$A==1] <- rnorm(n/2,dat1$mu_r[dat1$A==1],sigma)
        dat1 <- dat1[, c(1, 4, 5)]
        
        action_1 <- d1opt(data_1 = dat1, n = n, sigma = sigma, mu_0 = mu_0, K1 = K1, K = K,
                          lambda = lambda, c = c, weight.vec = weight.vec, seed = seed, delta.min = delta.min)
        
        set.seed(seed*10)
        matrix.return.list <- matrix.return(c = c, alpha = 0.05, beta = 0.2, delta.min = delta.min, sigma = sigma)
        ideal.pop.prime <- matrix.return.list$ideal.pop.prime
        ideal.pop.prime <- ideal.pop.prime[, -5]
        ##continue to enroll action_1 patient populaiton###
        num.pieces <- round((ideal.pop.prime[2,action_1]-ideal.pop.prime[1,action_1])/0.25)
        risk.score <- matrix(0,n/num.pieces,num.pieces)
        for(j in 1:num.pieces){
                risk.score[,j] <- runif(n/num.pieces,(ideal.pop.prime[1,action_1]+0.25*(j-1)),(ideal.pop.prime[1,action_1]+0.25*j))
        }
        risk.score <- as.vector(risk.score) 
        treatment.ind <- matrix(0,n/num.pieces,num.pieces)
        for(j in 1:num.pieces){
                treatment.ind[sample(1:(n/num.pieces),n/(num.pieces*2),replace = F),j] <- 1
        }
        treatment.ind <- as.vector(treatment.ind) 
        mu_r <- mu_0 + get(paste0("b",true.model))(risk.score) 
        data_2 <- data.frame(Risk=risk.score,mu_0=mu_0,mu_r=mu_r,A=treatment.ind)
        
        ##generate response by risk score
        data_2$Y <- NA
        data_2$Y[data_2$A==0] <- rnorm(n/2,mu_0,sigma)
        data_2$Y[data_2$A==1] <- rnorm(n/2,data_2$mu_r[data_2$A==1],sigma)
        data_2 <- data_2[, c(1,4,5)]
        action_2 <- d2opt(data_1 = dat1, data_2 = data_2, K = K, lambda = lambda, c = c, weight.vec = weight.vec, 
                          delta.min = delta.min, sigma = sigma)$action.recom 
        return(list(action_1=action_1,action_2=action_2))
        
}


#choice.population.two.stage(seed = 354, true.model = 2, c = 0.32, lambda = 0.01, K1 = 100, K = 1000)
simul.setting <- expand.grid(k = 1:200, model = 1:6)
phase2_simulation <- function(line, lambda = 0, c = 0, weight.vec = c(rep(0.1,5), 0.5)) {
        simul_no <- simul.setting[line, 1]
        true.model <- simul.setting[line, 2]
        return(choice.population.two.stage(seed = true.model * 1000 + simul_no, lambda = lambda, c = c , weight.vec = weight.vec, true.model = true.model)) 
}

result_list <- bplapply(1 : 1200, phase2_simulation, lambda = 0.01, c = 0.34)
result_matrix <- do.call(rbind, result_list)
result.mat <- matrix(0, 6, ncol = 9)
for(i in 1:6){
        Rep.matrix <- result_matrix[((1 + 200 * (i - 1)) : (200 * i)), ]
        for(j in 1:5){
                result.mat[i,j] <- sum(Rep.matrix[,2]==j)/200
        }
        for(j in 6:9){
                result.mat[i,j] <- sum(Rep.matrix[,1]==(j-5))/200
        }
}
result.mat

result_adp_032 <- result.mat
result_fix_032 <- result.mat
result_adp_034 <- result.mat

###bayes risk###
#########compute expected utility#########
Kprime = 15000; seed = 316; n = 264; sigma = 3; K1 = 100; K = 1000; mu_0 = 3; lambda = 0.01; c = 0.32;
weight.vec = c(rep(0.1,5), 0.5); delta.min = 0.3; 

set.seed(seed)
sample.index <- sample(1:6, Kprime, prob = weight.vec, replace = T)
mat.returned <- matrix.return(c = c, alpha = 0.05, beta = 0.2, delta.min = delta.min, sigma = sigma)
integral.matrix <- mat.returned$integral.matrix
duration <- mat.returned$duration
power <- mat.returned$power
integral.matrixo <- mat.returned$integral.matrixo
health.cost <- mat.returned$health.cost
ideal.pop.prime <- mat.returned$ideal.pop.prime
ideal.pop.prime <- ideal.pop.prime[, -5]
util <- NULL

util_value <- function(i, two.stage = T) {
        set.seed(seed * i)
        
        Risk <- matrix(0,n/4,4)
        for(j in 1:4){
                Risk[,j] <- runif(n/4,(j-1)/4,j/4)
        }
        Risk <- as.vector(Risk) 
        treatment.ind <- matrix(0,n/4,4)
        for(j in 1:4){
                treatment.ind[sample(1:(n/4),n/8,replace = F),j] <- 1
        }
        treatment.ind <- as.vector(treatment.ind) 
        true.model <- sample.index[i]
        mu_r <- mu_0 + get(paste0("b",true.model))(Risk)
        dat1 <- data.frame(Risk=Risk,mu_0=mu_0,mu_r=mu_r,A=treatment.ind)
        
        ##generate response by risk score
        dat1$Y <- NA
        dat1$Y[dat1$A==0] <- rnorm(n/2,mu_0,sigma)
        dat1$Y[dat1$A==1] <- rnorm(n/2,dat1$mu_r[dat1$A==1],sigma)
        dat1 <- dat1[, c(1, 4, 5)]
        
        if(two.stage == T) {
                action_1 <- d1opt(data_1 = dat1, n = n, sigma = sigma, mu_0 = mu_0, K1 = K1, K = K,
                                  lambda = lambda, c = c, weight.vec = weight.vec, seed = seed + i, delta.min = delta.min)
        }else {
                action_1 <- 1        
        }
        
        num.pieces <- round((ideal.pop.prime[2,action_1]-ideal.pop.prime[1,action_1])/0.25)
        risk.score <- matrix(0,n/num.pieces,num.pieces)
        for(j in 1:num.pieces){
                risk.score[,j] <- runif(n/num.pieces,(ideal.pop.prime[1,action_1]+0.25*(j-1)),(ideal.pop.prime[1,action_1]+0.25*j))
        }
        risk.score <- as.vector(risk.score) 
        treatment.ind <- matrix(0,n/num.pieces,num.pieces)
        for(j in 1:num.pieces){
                treatment.ind[sample(1:(n/num.pieces),n/(num.pieces*2),replace = F),j] <- 1
        }
        treatment.ind <- as.vector(treatment.ind) 
        mu_r <- mu_0 + get(paste0("b",true.model))(risk.score) 
        data_2 <- data.frame(Risk=risk.score,mu_0=mu_0,mu_r=mu_r,A=treatment.ind)
        
        ##generate response by risk score
        data_2$Y <- NA
        data_2$Y[data_2$A==0] <- rnorm(n/2,mu_0,sigma)
        data_2$Y[data_2$A==1] <- rnorm(n/2,data_2$mu_r[data_2$A==1],sigma)
        data_2 <- data_2[, c(1,4,5)]
        action_2 <- d2opt(data_1 = dat1, data_2 = data_2, K = K, lambda = lambda, c = c, weight.vec = weight.vec, 
                          delta.min = delta.min, sigma = sigma)$action.recom 
        util <- integral.matrix[sample.index[i], action_2] * power[sample.index[i], action_2] - lambda*duration[action_2]
        total_benefit <- integral.matrixo[sample.index[i], action_2]
        health_cost <- health.cost[action_2]
        power_sq <- power[sample.index[i], action_2]
        duration_val <- duration[action_2]
        return(c(util, total_benefit, health_cost, power_sq, duration_val))
}

fix_util <- bplapply(1 : Kprime, util_value, two.stage = F)
adp_util <- bplapply(1 : Kprime, util_value)

# in case bplapply fails
# adp_util <- NULL
# for(i in 1723 : 3000) {
#         temp <- bplapply(((i - 1) * 5 + 1) : (i * 5), util_value)
#         adp_util <- rbind(adp_util, do.call(rbind, temp))
#         print(i)
# }
