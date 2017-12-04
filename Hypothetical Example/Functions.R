###Six Functions##########
b1 <- function(r){
        return(1)
}

##function 2#### ideal group (0.5,1)
b2 <- function(r){
        return(2*(r-1/2)*(r>=1/2))
}

##function 3#### ideal group (0,0.5)
b3 <- function(r){
        return(-2*(r-0.5)*(r<=0.5))
}
##function 4#### ideal group (0.25,0.75)
b4 <- function(r){
        return(4*(r-0.25)*(r>=0.25 & r <= 0.5) + 
                       (-4)*(r-0.75)*(r>0.5 & r<=0.75))
}
##function 5#### ideal group (0,0)
b5 <- function(r){
        return(-1)
}

##function 6#### ideal group (0,0)
b6 <- function(r){
        return(0)
}

