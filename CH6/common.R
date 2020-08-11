#==========================================================================
# SET PARAMETER VALUES
#==========================================================================

mu     = 3;               # risk aversion (=3 baseline)             
beta   = 0.96;            # subjective discount factor 
delta  = 0.08;            # depreciation
alpha  = 0.36;            # capital's share of income
b      = 3;               # borrowing limit

#==========================================================================
# COMPUTE TRANSITION MATRIX OF LABOR PRODUCTIVITY
#==========================================================================

# ROUTINE tauchen.m TO COMPUTE TRANSITION MATRIX, GRID OF AN AR(1) AND
# STATIONARY DISTRIBUTION
# approximate labor endowment shocks with seven states Markov chain
# log(s_{t}) = rho*log(s_{t-1})+e_{t} 
# e_{t}~ N(0,sig^2)

Nl       = 7;             # number of discretized states
rho      = 0.6;           # first-order autoregressive coefficient
sig      = 0.4;           # intermediate value to calculate sigma (=0.4 BASE)

M=2; 

source("tauchen.R")

#
# σから遷移行列などを計算する
# args: sigma
# return: 
# 	prob   : transition matrix of the Markov chain
# 	logs   : the discretized states of log labor earnings
# 	invdist: the invariant distribution of Markov chain
#	labor  : 定常状態の労働供給量？
#
calcTM <- function(sigma){
	r_tauchen <- tauchen(Nl,0,rho,sigma,M) # 第2引数は定常状態の値でAiyagari (1994)の設定ではゼロ
	invdist <- with(r_tauchen, {
		# 固有値問題を解くより単純にループさせた方が速い
		tol <- 1e-7
		# 定常状態の値を計算する
		x0 <- matrix(1, 1, ncol(TransitionMatrix)) / ncol(TransitionMatrix)
		for(i in 1:1000){
			x1 <- x0 %*% TransitionMatrix
			if(tol > sum(abs(x1 - x0))){
				break
			}
			x0 <- x1
		}
		matrix(x1, ncol(TransitionMatrix), 1)
	})
	prob <- r_tauchen$TransitionMatrix;
	logs <- r_tauchen$Grid;
	s <- matrix(exp(logs), 1, length(logs))
	labor <- c(s %*% invdist)
	return(list(
		prob = prob,
		logs = logs,
		invdist = invdist,
		s = s,
		labor = labor	
	))
}

#
# 実験コードを読み込む（機能はしない）
# source("experiment.R")
#

#
# フラグチェック関数
# args: fpos;フラグのビット位置
# return: true/false
#
chkflag <- function(n, fpos){
	if(fpos < 1 || !is.numeric(fpos)){
		stop("Wrong number as fpos!")
	}
	m <- 2^(fpos - 1)
	1 == (as.integer(n/m)%%2)
}

#
# パッケージがインストールされていなければ止まる
#
chkPkg <- function(pkgnames){
	for(pkgname in pkgnames){
		if(!any(suppressWarnings(library(quietly=TRUE, verbose=FALSE)$results[,"Package"] == pkgname))){
			stop(paste("Do install.packages(\"", pkgname, "\") before runnning this script.", sep=""))
		}
	}
}

#
# Cのモジュールを呼び出す関数 
#
loadCModule <- function(){
	dll <- dyn.load(paste(vfi_type, .Platform$dynlib.ext, sep = ""))
	function(r, alpha, beta, delta, mu, b, s, prob){
		.Call("VFI", r, alpha, beta, delta, mu, b, s, prob)
	}
}

