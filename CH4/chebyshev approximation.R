#
# 経済セミナー連載「定量的マクロ経済学と数値計算」
# 第4回「オイラー方程式と多項式近似」のソースコードのRへの野良移植: チェビシェフ多項式
#

#################################################################
### Warning: Package 'chebpol' の存在に気づいてはいけません。 ###
#################################################################

#
# y = f(x)をチェビシェフ多項式で近似するための関数群
#

#
# パラメーターθを作成するためのxの評価点を取得する
# args: xmin:xの定義域の最小値, xmax:xの定義域の最大値, n:求める評価点の数
#
ChebyshevCollocationPoints <- function(xmin, xmax, n, type="extrema"){
	if("extrema" == type){
		i <- 0:(n - 1)
		s <- -cos((pi/(n - 1)) * i)
	} else if("zeros" == type){
		i <- 1:(n - 1)
		s <- c(0, -cos((pi/2/(n-1))*(2*i - 1)))
	} else {
		stop("unknown type")
	}
	(xmax - xmin)*(s + 1)/2 + xmin
}

# 
# 基底関数の行列を求める
# args: xmin:xの定義域の最小値, xmax:xの定義域の最大値, d:次数, x:近似する関数に実際に引数として与えられたベクトル
# 
ChebyshevBasis <- function(xmin, xmax, d, x){
	n <- d + 1
	s <- (2/(xmax - xmin))*(x - xmin) - 1
	T <- matrix(0, length(x), n)
	T[, 1] <- 1
	T[, 2] <- s
	for(i in 3:n){
		T[, i] <- 2*s*T[, i-1] - T[, i-2]
	}
	T
}

#
# Collocation pointsになるcp_xとcp_yの対応から関数の近似のためのパラメーターΘを作成し、xの定義域とセットにしてリストにまとめる
#
ChebyshevMapApproximation <- function(cp_x, cp_y, xmin, xmax, d){
	T <- ChebyshevBasis(xmin, xmax, d, cp_x)
	theta <- solve(T, cp_y) # = solve(T) %*% cp_y
	list(theta=theta, xmin=xmin, xmax=xmax, cp_x=cp_x, cp_y=cp_y)
}

#
# 関数の近似のためのθを求める
# args: f:近似する対象の関数, xmin:xの定義域の最小, xmax:xの定義域の最大, d:次数
#
ChebyshevApproximation <- function(f, xmin, xmax, d){
	cp_x <- ChebyshevCollocationPoints(xmin, xmax, d + 1)
	cp_y <- f(cp_x)
	ChebyshevMapApproximation(cp_x, cp_y, xmin, xmax, d)
}

#
# θから予測値をつくる
# args: xmin:xの定義域の最小, xmax:xの定義域の最大, theta:θ, x:予測する点
#
ChebyshevPredict <- function(caparam, x){
	with(caparam, {
		T <- ChebyshevBasis(xmin, xmax, length(theta) - 1, x)
		predict_y <- T %*% theta
	})
}

