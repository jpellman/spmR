# return error covariance constraints (for serially correlated data)
# FORMAT [C] = spm_Ce(v,a)
# v  - (1 x l) v(i) = number of observations for i-th block
# a  - AR coefficient expansion point  (default a = [])
# 
# a  = [] (default) - block diagonal identity matrices specified by v:
#
#   C{i}  = blkdiag( zeros(v(1),v(1)),...,AR(0),...,zeros(v(end),v(end)))
#   AR(0) = eye(v(i),v(i))
#
# otherwise:
#
#   C{i}     = AR(a) - a*dAR(a)/da;
#   C{i + 1} = AR(a) + a*dAR(a)/da;
#
# See also: spm_Q.m

spm_Ce <- function(v, a=numeric(0)) {

    C <- list()
    l <- length(v)
    n <- sum(v)
    k <- 0

    if(l > 1) {
        ### FIXME!!!
        ### not tested, probably not correct!
        c.idx <- 1
        for(i in 1:l) {
            dCda <- spm_Ce(v[i], a)
            for(j in 1:length(dCda)) {
                out <- which(dCda[[j]] != 0, arr.ind=TRUE)
                x <- out[,"row"]; y <- out[,"col"]; q <- dCda[[j]][out]   
                x <- x + k
                y <- y + k
                tmp <- diag(n); tmp[x,y] <- q
                C[[x.idx]] <- tmp
                c.idx <- c.idx + 1
            }
        } 
    } else {
        # dCda
        if(length(a) != 0) {
            Q      <- spm_Q(a,v)
            #dQda   <- spm_diff("spm_Q", a, v, 1)
            dQda   <- spm_Qdiff(a, v)
            C[[1]] <- Q - dQda * a
            C[[2]] <- Q + dQda * a
        } else {
            C[[1]] <- diag(v)
        }
    }

    C
}

spm_Qdiff <- function(a, v) {

    h <- exp(-8)
    Q <- spm_Q(a, v)
    Q.right <- spm_Q(a+h, v)
    
    Q.diff <- (Q.right - Q)/h    
    
    Q.diff
}

