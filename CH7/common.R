#
# 共通関数
#

# nleqslvパッケージがインストールされていない場合、先にインストールするように促して止まる
if(!any(suppressWarnings(library(quietly=TRUE, verbose=FALSE)$results[,"Package"] == "nleqslv"))){
	stop("Do install.packages(\"nleqslv\") before runnning this script.")
}
library(nleqslv)

calcOLG <- function(dparam){
	with(dparam, {
		meaJ <- rep(1/Nj, Nj) # AGE DISTRIBUTION (no population growth) 
		theta <- numeric(Nj) 
		theta[1:NjW] <- 1.0
		
		L <- sum(meaJ * theta) # 1:NjWは要るのか？
		
		ss0 <- numeric(Nj)
		ss0[(NjW + 1):Nj] <- rho # MatlabとRでオペラントの優先順序が違うかも知れず
		
		tau <- rho*(Nj-NjW)/NjW;
		
		# INITIAL GUESS FOR ASSET DECISIONS
		# a(2) to a(Nj) (NOTE a(1)=a(Nj+1)=0)
		stepsize <- 0.01
		a_g <- seq(0.01, stepsize*(Nj - 1), stepsize) # 0.01刻みで0.01～0.60
		
		policyA <- function(aX){
			a <- c(0, aX, 0) # aのサイズはNj+1, 死ぬ瞬間は資産ゼロ
			K <- sum(a[1:Nj]*meaJ)
			r <- alpha*(K/L)^(alpha - 1) - delta # INTEREST RATE
			w <- (1-alpha)*(K/L)^alpha # WAGE
			ss <- ss0 * w
			# EULER EQUATION (u(c)=log(c))
			ff <- beta*(1+r)*( (1+r)*a[1:(Nj-1)] + w*theta[1:(Nj-1)]*(1-tau) + ss[1:(Nj-1)] - a[2:Nj])
			f <- (1+r)*a[2:Nj] + w*theta[2:Nj]*(1-tau) + ss[2:Nj] - a[3:(Nj + 1)] - ff
		}
		
		# FIND SOLUTIONS : aX(1:Nj-1) 
		r_nleqslv <- nleqslv(a_g, policyA) # 初期値はa_g; 非線形連立方程式を解くnleqslvは、Matlabのfsolve相当になる
		if(1 < r_nleqslv$termcd){
			stop(r_nleqslv$message) # 収束しなかったら、止まる
		}
		aX <- r_nleqslv$x
		
		a <- numeric(Nj + 1)
		a[2:Nj] <- aX
		K <- sum(a[1:Nj]*meaJ)
		r <- alpha*(K/L)^(alpha - 1) - delta # INTEREST RATE
		w <- (1-alpha)*(K/L)^alpha # WAGE
		ss <- ss0*w
		c <- w*theta[1:Nj]*(1-tau) + (1+r)*a[1:Nj] + ss[1:Nj] - a[2:(Nj + 1)]
		list(a=a, K=K, c=c, r=r, w=w)
	})
}

