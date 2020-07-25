#
# 経済セミナー連載「定量的マクロ経済学と数値計算」
# 第5回「ニューケインジアン・モデルの新展開」のソースコードのRへの野良移植: 2状態モデル（解析的解）
#

# 共通の設定を読み込む
source("common.R")
#
# 解析的解を求める関数
# args: p_H:状態H→Lへの転落確率
#
closed_form_solution <- function(p_H){
	with(dparam, {
		# x := t(y_H, y_L, pi_H, pi_L)として、A %*% x = bが雑誌本文の代数的解の式と同値になるようにAとbを定める
		# 本文の式のうちr_L=0とr_Hは代入して消され、方程式は4本の減っている
		A <- matrix(c(
			-1+(1-p_H), p_H, -(phi-1)*(1-p_H), -(phi-1)*p_H,
			kapper, 0, -1+beta*(1-p_H), beta*p_H,
			(1 - p_L), -1+p_L, 1-p_L, p_L,
			0, kapper, beta*(1-p_L), -1+beta*p_L
		),4,4, byrow=TRUE)
		b <- matrix(c(r_star - s_H, 0, -s_L, 0), 4, 1)

		x <- solve(A, b) # = solve(A) %*% b
		makePolicyMatrix(x[1], x[2], x[3], x[4], r_star + phi*((1 - p_H)*x[3] + p_H * x[4]))
	})
}

closed_form_solution(0)
closed_form_solution(0.025)

