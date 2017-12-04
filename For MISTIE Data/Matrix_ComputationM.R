#######calculate the integral matrices and other needed matrices for computation
matrix.return <- function(c, alpha, beta, delta.min, sigma){
        ##function 1#### ideal group (0,1)
        ###Six Functions##########
        b1 <- function(r){
                return(0.5)
        }
        
        ##function 2#### ideal group (0.5,1)
        b2 <- function(r){
                return(1*(r-1/2)*(r>=1/2))
        }
        
        ##function 3#### ideal group (0,0.5)
        b3 <- function(r){
                return(-1*(r-0.5)*(r<=0.5))
        }
        ##function 4#### ideal group (0.25,0.75)
        b4 <- function(r){
                return(3*(r-0.25)*(r>=0.25 & r <= 0.5) + 
                               (-3)*(r-0.75)*(r>0.5 & r<=0.75))
        }
        ##function 5#### ideal group (0,0)
        b5 <- function(r){
                return(-0.125)
        }
        
        ##function 6#### ideal group (0,0)
        b6 <- function(r){
                return(0)
        }
        
        ideal.pop1 <- c(0,1)
        ideal.pop2 <- c(0.5,1)
        ideal.pop3 <- c(0,0.5)
        ideal.pop4 <- c(0.25,0.75)
        ideal.pop5 <- c(0,0)
        ideal.pop6 <- c(0,0)
        
        ###get rid of duplicated selection##########################
        ideal.pop <- NULL
        for(i in 1:6){
                ideal.pop <- cbind(ideal.pop,get(paste0("ideal.pop",i)))
        }
        
        dup.index <- NULL
        for(i in 1:5){
                for(j in (i+1):6){
                        if(identical(ideal.pop[,i],ideal.pop[,j])) dup.index <- c(dup.index,j)        
                }
        }
        dup.index <- dup.index[!duplicated(dup.index)]
        ideal.pop.prime <- ideal.pop[,-dup.index]
        
        ############integral matrix################
        integral.matrix <- matrix(0,6,5)
        for(i in 1:6){
                for(j in 1:5){
                        integral.matrix[i,j] <- integrate(Vectorize(get(paste0("b",i))),ideal.pop.prime[1,j],ideal.pop.prime[2,j])$value
                }
        }
        
        ###########duration##################
        duration <- 1/(ideal.pop.prime[2,] - ideal.pop.prime[1,])
        duration[ideal.pop.prime[1,]==0 & ideal.pop.prime[2,]==0] <- 0
        
        ##########average treatment effect####################
        avg.trteff <- matrix(0,6,5)
        avg.trteff <- sweep(integral.matrix, 2, STATS = (ideal.pop.prime[2,] - ideal.pop.prime[1,]), FUN = "/")
        avg.trteff[, 5] <- 0
        
        ##########integral matrix with health care cost###########
        integral.matrix <-  sweep(integral.matrix, 2, STATS = c * (ideal.pop.prime[2,] - ideal.pop.prime[1,]), FUN = "-")
        
        
        #########average treatment effect for posterior sampling####
        average.treat <- matrix(0, 6, 4)
        interval.consid <- matrix(0, 2, 4)
        interval.consid[1, ] <- c(0, 0.25, 0.5, 0.75)
        interval.consid[2, ] <- c(0.25, 0.5, 0.75, 1)
        for(i in 1:6){
                for(j in 1:4){
                        average.treat[i,j] <- integrate(Vectorize(get(paste0("b",i))),interval.consid[1,j],interval.consid[2,j])$value/(interval.consid[2,j]-interval.consid[1,j])
                }
        }
        
        #########compute probability of phase iii success#######
        n.min <- 2 * (qnorm(1 - alpha) + qnorm(1 - beta))^2/(delta.min^2/sigma^2)
        power <- pnorm(sqrt(n.min) * avg.trteff/sqrt(2 * sigma^2) - qnorm(1 - alpha))  ##type i error default 0.05 for phase iii, beta = 0.2.
        
        
        return(list(integral.matrix=integral.matrix,duration=duration,avg.trteff=avg.trteff,mean.dat=average.treat,ideal.pop.prime=ideal.pop.prime,power=power^2))
}

#matrix.return(c = 0, alpha = 0.1, beta = 0.2, delta.min = 0.2, sigma = 3)

