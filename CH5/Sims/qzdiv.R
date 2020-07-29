#' qzdiv
#'
#' Sort the roots in a qz decomposition around a value
#' 
#' Takes a list containing matrices \code{a, b}, orthonormal
#' matrices \code{q,z}, and rearranges them so that all cases of
#' \code{abs(b[i,i] / a[i,i]) > stake} are in lower right corner, while
#' preserving U.T. and orthonormal properties and \code{q \%*\% a \%*\% t(z)} and
#' \code{q \%*\% b \%*\% t(z)}. 
#'
#' @param stake the real number about which the roots are sorted
#' @param qzlist the qz transform list.  The \code{qz} from the QZ package
#'   uses different names in its output list, which have to be translated
#'   as at the start of \code{gensys} before being passed to this program.
#' @param flip If \code{flip} is TRUE, then cases with
#'   \code{abs(b[i,i]/a[i,i]) < stake} are in lower right.
#'
#' @return The input list, with all matrices re-ordered.
#' @seealso \link{gensys}, \link{QZ::qz}, \link{qzswitch}
#' @export 
qzdiv <- function (stake,qzlist, flip=FALSE)  {
  ##
  a <- qzlist$a
  b <- qzlist$b
  q <- qzlist$q
  z <- qzlist$z
  n  <- dim(a)[1]
  root <- abs(cbind(diag(a), diag(b)))
  if (flip) {
    root <-  root[ , 2:1]
    stake <- 1/stake
  }
  root[,1] <- root[,1]-(root[,1]<1.e-13)*(root[,1]+root[,2])
  root[,2] <- root[,2]/root[,1]
  for (i in  n:1) {
    m <- 0
    for (j in i:1) {
      if (root[j,2] > stake || root[j,2] < -.1) {
        m <- j
        break                           #found an unstable root.  Now check for stable ones below.
      }
    }
    if (m == 0) {
      break                             #quit sort because only stable roots left
    } else {
      if (m<i) {
        for (k in m:(i-1)) {
          qzlist <- qzswitch(k,qzlist)
          tmp <- root[k,2]
          root[k,2] <- root[k+1,2]
          root[k+1,2] <- tmp
        }
      }         
    }
  }
  return(qzlist)
}
