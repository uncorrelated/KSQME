#' gensys
#'
#' Solves multivariate linear rational expectations models
#'
#' @param g0,g1,c0,psi,pi matrices that form the system to be solved
#' @param div Roots as large as or larger than this in absolute value are
#' suppressed in the solution.  The \code{div=-1} default value allows
#' roots with unit absolute value and excludes those larger.
#' 
#' @details
#' System given as
#'     \deqn{g0 y(t) = g1 y(t-1) + c0 + psi z(t) +
#'          pi eta(t)}{g0 \%*\% y(t) = g1 \%*\% y(t-1) + c0 +
#'             psi \%*\% z(t) + pi \%*\% eta(t)},
#' with \code{z} an exogenous variable process and \code{eta} being
#' endogenously determined one-step-ahead expectational errors.
#'
#' @return
#' A list with elements
#' \describe{
#'   \item{\code{G1, C, impact, ywt, fmat, fwt, loose}}{matrices forming the solution
#'      system}
#'   \item{\code{eu}}{\code{eu(1)=1} for existence, \code{eu(2)=1} for
#'          uniqueness.  \code{eu=[-2,-2]} for coincident zeros, implying an
#'          incomplete system.}
#'   \item{gev}{Ratio of second column to first column of this matrix is the
#'      set of generalized eigenvalues of the system.}
#'   }
#' Returned system is
#'       \deqn{y(t) = G1  y(t-1) + C + impact  z(t) +
#'           ywt \sum_{s=1}^\infty fmat^s  fwt z(t+s)
#'              + loose eta(t)}{y(t)=G1 \%*\% y(t-1) + C + impact \%*\% z(t) +
#'           ywt \%*\% inv(I - fmat \%*\% inv(L)) \%*\% fmat \%*\% fwt \%*\% z(t+1)
#'              + loose \%*\% eta }
#' If \code{z(t)} is i.i.d., the term involving \code{ywt, fmat} and \code{fwt}
#' drops out.  If the solution is unique (eu[2]==1) there is no "loose" term.
#' Otherwise loose characterizes the dimensions along which there is
#' non-uniqueness.
#' @export
#' @import QZ
#' 
gensys <- function(g0, g1, c0=matrix(0,dim(g0)[1],1), psi, pi, div=-1) {
    eu <- c(0,0)
    realsmall <- 1e-7
    fixdiv <- (div>0)
    n <- dim(g0)[1]
    nshock <- if (is.matrix(psi)) dim(psi)[2] else if (is.null(psi)) 0 else 1
    g0 <- g0 + 0+0i
    g1 <- g1 + 0+0i
    qzl <- QZ::qz(g0, g1)
    ## The 0+0i terms are necessary to make the arguments of qz complex.  The
    ## QZ package qz() returns the real generalized Schur, with 2x2 real blocks
    ## on the diagonal corresponding to complex roots, when its arguments
    ## are real, but the full complex QZ for complex arguments.
    ## Lines below have also been adjusted to use the qz() from the QZ package
    ## The original cas qz wrapper delivers a list with elements a, b, q, z, gev
    ## and rc.  These correspond to a=S, b=T, q=Q, z=Z, 
    ## gev=matrix(c(ALPHA,BETA), ncol=2), and rc=INFO.
    ##-------- Translation from QZ package qz output to local qz ----------------
    qzl <- with(qzl, list(a=S, b=T, q=Q, z=Z, gev=matrix(c(ALPHA,BETA),ncol=2), rc=INFO))
    ##-----------------------------------------------------------
    zxz <- any((abs(diag(qzl$a))<realsmall) & (abs(diag(qzl$b))<realsmall))
    if (zxz) {
        "Coincident zeros.  Indeterminacy and/or nonexistence.\n"
        eu <- c(-2,-2)
        gev <- qzl$gev
        return(list(eu=eu,gev=gev))
    }
    zeroax <- abs(diag(qzl$a)) < realsmall
    unstabx <- abs(diag(qzl$a)) < (1-realsmall)*abs(diag(qzl$b)) # near unit roots don't count
    unstabx <- (! zeroax) & unstabx
    if (! fixdiv) {
        if (! any(unstabx)){
            div <- 1.01
        } else
            {
    div <- .5*(min(abs(diag(qzl$b)[unstabx]/diag(qzl$a)[unstabx]))+1)
}
    }
    unstabx <- div*abs(diag(qzl$a))<= abs(diag(qzl$b))
    nunstab <- sum(unstabx)
    qzl <- qzdiv(div,qzl)
    qq <- t(Conj(qzl$q))                # to match matlab convention 
    gev <- qzl$gev
    ## note that this means that gev is not simply the diagonals of a nd b.  qzdiv
    ## changes the numbers on the diagonals (though not their ratios), but merely reorders
    ## the original gev.  
    if (nunstab==n){
        six <- NULL
        uix <- 1:n
    } else
        {
    if (nunstab==0){
        uix <- NULL
        six <- 1:n
    } else
        {
    uix <- (n-nunstab+1):n
    six <- 1:(n-nunstab)
}
}
    q1 <- qq[six,,drop=FALSE]
    q2 <- qq[uix,,drop=FALSE]
    z1 <- t(Conj(qzl$z[,six,drop=FALSE]))
    z2 <- t(Conj(qzl$z[,uix,drop=FALSE]))
    a2 <- qzl$a[uix,uix,drop=FALSE]
    b2 <- qzl$b[uix,uix,drop=FALSE]
    ## debug
    ## browser()
    etawt <- q2 %*% pi
    neta <- if (is.matrix(pi)) dim(pi)[2] else if (is.null(pi)) 0 else 1
    ndeta <- min(nunstab,neta)
    if(ndeta==0){
        ueta <- matrix(0,nunstab,0)
        deta <- vector("numeric",0)
        veta <- matrix(0,neta,0)
        bigev <- vector("logical",0)
    } else {
          sd <- svd(etawt)
          ueta <- sd$u; deta <- sd$d; veta <- sd$v
          bigev <- deta>realsmall
          ueta<-ueta[,bigev,drop=FALSE]
          veta<-veta[,bigev,drop=FALSE]
          deta<-deta[bigev]
      }
    eu[1] <- sum(bigev) >= nunstab
    ##----------------------------------------------------
    ## Note that existence and uniqueness are not just matters of comparing
    ## numbers of roots and numbers of endogenous errors.  These counts are
    ## reported below because usually they point to the source of the problem.
    ##------------------------------------------------------
    etawt1 <- q1 %*% pi
    ndeta1 <- min(n-nunstab,neta)
    if(ndeta1==0){
        ueta1 <- matrix(0,n-nunstab,0)
        deta1 <- vector("numeric",0)
        veta1 <- matrix(0,neta,0)
        bigev1 <- vector("logical",0)
    } else {
          sd <- svd(etawt1)
          ueta1<-sd$u
          deta1 <- sd$d
          veta1 <- sd$v
          bigev1 <- deta1 > realsmall
      }
    if (any(bigev1)) { #needed because empty dimensions are dropped after select
        ueta1 <- ueta1[,bigev1,drop=FALSE]
        veta1 <- veta1[,bigev1,drop=FALSE]
        deta1 <- deta1[bigev1]
        loose <- veta1-veta %*% t(Conj(veta)) %*% veta1
        svdl <- svd(loose)
        loose <- sum(abs(svdl$d)>realsmall*n)
        unq <- (loose==0)
    } else {
          ueta1 <- matrix(1,n-nunstab,0)
          veta1 <- matrix(1,neta,0)
          deta1 <- vector("complex",0)
          unq <- TRUE
      }
    if (unq) {
        eu[2] <- 1
    } else
        {
    cat("Indeterminacy.", loose, "loose endog errors.\n")
}
    ## Note: if v is a vector of length n and m is an nxp matrix,
    ## v*m==diag(v)%*%m, m/v==solve(diag(v),m)==diag(v)\m (matlab notation)
    ##
    tmat <- cbind(diag(n-nunstab),
                  -t(Conj((ueta %*% (t(Conj(veta))/deta)) %*% veta1 %*% (deta1 * t(Conj(ueta1)))))  )  
    G0<- rbind( tmat %*% qzl$a, cbind(matrix(0,nunstab,n-nunstab), diag(nunstab)))
    G1<- rbind(tmat %*% qzl$b, matrix(0,nunstab,n))
    ##----------------------
    ## G0 is always non-singular because by construction there are no zeros on
    ## the diagonal of a[1:(n-nunstab),1:(n-nunstab)], which forms G0's ul corner.
    ##-----------------------
    G0I <- solve(G0)
    G1 <- G0I%*%G1
    ##----------- uix can be empty, e.g. in indeterminate systems with no unstable roots ------------
    if(is.null(uix)){
        C <- G0I %*% tmat %*% qq %*% c0
        fmat <- matrix(0,0,0)
        fwt <- matrix(0, 0, nshock)
        impact <- G0I %*% tmat %*% qq %*% psi
    }else{
         C <- G0I %*% rbind(tmat%*% qq %*%c0,solve(qzl$a[uix,uix,drop=FALSE]-qzl$b[uix,uix,drop=FALSE],q2%*%c0) )
         impact <- G0I %*% rbind(tmat %*% qq %*% psi, matrix(0,nunstab, nshock))
         fmat <- solve(qzl$b[uix,uix,drop=FALSE],qzl$a[uix,uix,drop=FALSE])
         fwt <- -solve(qzl$b[uix,uix,drop=FALSE],q2 %*% psi)
     }
    ywt <- G0I[,uix,drop=FALSE]
    ##loose <- G0I %*% etawt1 %*% (diag(neta) - veta %*% t(Conj(veta)))
    loose <- G0I %*% qq %*% pi %*% (diag(neta) - veta %*% t(Conj(veta)))
    ## loose <- G0I %*% rbind(loose,matrix(0,nunstab,neta))  #(think the above is a mistaken remnant)
    ##-------------------- above are output for system in terms of z'y -------
    G1 <- (qzl$z %*% G1 %*% t(Conj(qzl$z)))
    C <- (qzl$z %*% C)
    impact <- (qzl$z %*% impact)
    ywt <- qzl$z %*% ywt
    loose <- (qzl$z %*% loose)
    vn <- dimnames(g0)[[2]]
    dimnames(G1) <- list(vn,vn)
    dimnames(C) <- list(vn,NULL)
    dimnames(impact)[[1]] <- vn
    dimnames(ywt)[[1]] <- vn
    dimnames(loose)[[1]] <- vn
    return(list(G1=G1,C=C,impact=impact,fmat=fmat,fwt=fwt,ywt=ywt,gev=gev,eu=eu,loose=loose))
}

