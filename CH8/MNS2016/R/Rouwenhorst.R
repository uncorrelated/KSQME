# Creates a discrete approximation to a first order autoregressive process
# with serial correlation p + q - 1.  See Rouwenhorst (1995) (pps. 325-329)
# in Cooley <i>Frontiers of Business Cycle Research</i> Princeton.

Rouwenhorst <- function(rho, sigmas, znum=2){ # default is 2x2 transition matrix

	p = (rho + 1)/2;
	q = p;

	hlag = matrix(1, 1, 1);

	for(i in seq(2, znum, 1)){
		n <- nrow(hlag)
		m <- i

		h1 <- h2 <- h3 <- h4 <- matrix(0, m, m) # m×mのゼロ行列をつくる

		idx1 <- 1:n
		idx2 <- (m - n + 1):m
		h1[idx1, idx1] <- hlag # 左上にhlagをコピー
		h2[idx1, idx2] <- hlag # 右上にhlagをコピー
		h3[idx2, idx1] <- hlag # 左下にhlagをコピー
		h4[idx2, idx2] <- hlag # 右下にhlagをコピー
		h <- p*h1 + (1-p)*h2 + (1-q)*h3 + q*h4

		if(2 < i){ # MatlabやJuliaは行列の添字2:jで2>jだとnull指定になるが、Rはj:2指定になる
			idx3 <- i - 1
			h[idx3, ] <- h[idx3, ] / 2
		}

		hlag <- h
	}

	PI = h;

	# symmetrically and evenly spaced between [-epsilon, epsilon] with h elements.  
	# When p = q, then then variance of shock is epsilon^2/(h-1).  
	zvar = (sigmas^2)/( 1 - rho^2);
	epsilon = sqrt((znum - 1)*zvar);

	Z = seq(-epsilon, epsilon, length.out=znum)
	list(p=PI, state_values=Z)
}

testRouwenhorst <- function(){
	r <- Rouwenhorst(0.96566, 0.01695^(1/2), 3)

	p <- matrix(c(0.9659548089, 0.03375038220000003, 0.00029480890000000066, 0.016875191100000016, 0.9662496178, 0.016875191100000016, 0.00029480890000000066, 0.03375038220000003, 0.9659548089), 3, 3, byrow=TRUE)

	if(any(1e-6 < abs(r$p - p))){
		print("coding error: p is not same as original.")
	} 

	sv <- seq(-0.7086723748665338, 0.7086723748665338, 0.7086723748665338)

	if(any(1e-6 < abs(r$state_values - sv))){
		print("coding error: state_values is not same as original.")
	}

}

