### 相対的リスク回避度一定効用関数（Constant Relative Risk Aversion Utility Function） ###
CRRA <- function(cons, gamma){
       if(1.0==gamma){
		return(log(cons))
       } else {
		return(cons^(1-gamma)/(1-gamma))
       }
}

# 消費の限界効用
dCRRA <- function(cons, gamma) cons^(-gamma)

### ベルマン方程式を定義 ###
BellmanEq <- function(kprime, k, vnext, alpha, beta, gamma){
	if(kprime <= 0){
		return(-1e+6) # トリック(1): k'は正の値しか取らないので、ペナルティを与えてその値が選ばれないようにする
	}
	wealth = k^alpha;
	cons = wealth - kprime;
	util = ifelse(0<cons, CRRA(cons, gamma), -1e+6); # 消費が負値の場合、ペナルティを与えてその値が選ばれないようにする
	util + beta*vnext(kprime)
}

# Discretized Dynamic Programming.Rのオイラー方程式の誤差を保存するファイル名
file_err_ddp <- "err_ddp.csv"

