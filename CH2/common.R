#
# 複数の手法で共通する設定と関数
#

# グリッドのパラメーター
# number:グリッドの数, min:最小値, max:最大値
gparam_w <- list(number = 10, min=0.1, max=1.0)
gparam_a <- list(number = 40, min=0.025, max=1.0)

# グリッド生成関数
init_grid <- function(gparam){
       with(gparam,{
               seq(min, max, length.out=number)
       })
}

#
# ディープ・パラメーター
# beta:割引因子, gamma:相対的危険回避度, nir: 純利子率
# 
dparam <- list(beta=0.985^30, gamma=2.0, nir=1.025^30-1.0)

#
# 相対的リスク回避度一定効用関数（Constant Relative Risk Aversion Utility Function）
# args: 消費, 相対的危険回避度
# return: 効用水準
# 
CRRA <- function(cons, gamma=dparam[["gamma"]]){
       if(1.0==gamma){
               log(cons)
       } else {
               cons^(1-gamma)/(1-gamma)
       }
}

# 消費の限界効用
dCRRA <- function(cons, gamma=dparam[["gamma"]]) cons^(-gamma)

# オイラー方程式の変形
euler_eq <- function(a, w){
	cons_0 <- w - a # 0期の消費
	if(any(cons_0 <= 0)){
		return(-1e+10) # -∞の代わりを返す
	}
	with(dparam, {
		cons_1 <- (1.0 + nir)*a # 1期の消費
		beta*(1.0 + nir)*dCRRA(cons_1)/dCRRA(cons_0)
	})
}

