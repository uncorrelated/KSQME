dparam <- list(
	r_star = 0.75,	# 状態H→Lの確率:p_H=0のときの、定常状態での名目金利の値
	phi = 5.0,	# テイラールールにおける名目金利の期待インフレ率に対する反応係数φ(注: 小さいとpi_L=0にならない)
	p_L = 0.75,	# 状態L（危機？）の継続確率
	s_L=-1.5625,	# 状態Lでの自然利子率の値
	kappa=0.0091	# フィリップス曲線の傾き
)
dparam["s_H"] <- dparam$r_star # 状態Hでの自然利子率の値
dparam["beta"] = 1/(1+dparam$r_star/100); # 割引率(オイラー方程式の定常状態より)

# source("find parameters.R") # ディープパラメーターの調整に用いたコードだが、計算結果は初期設定に反映しているため、実行不要

# 計算結果をまとめる行列をつくる
makePolicyMatrix <- function(y_H, y_L, pi_H, pi_L, r_H, r_L=0){
	policy <- matrix(c(y_H, y_L, pi_H, pi_L, r_H, r_L), 2, 3)
	colnames(policy) <- c("y","pi","r") # 列が変数
	rownames(policy) <- c("H", "L") # 行が状態
	policy
}

