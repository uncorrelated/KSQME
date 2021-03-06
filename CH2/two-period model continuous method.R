#
# 経済セミナー連載「定量的マクロ経済学と数値計算」
# 第2回「2期間モデルと数値計算の概観」のソースコードのRへの野良移植: 求める値は連続的な方法, 図2, 代数的解, 最適化アルゴリズム
#


# 他の手法と共通する設定などを読み込む
source("common.R")

# 計算用グリッドを初期化
grid_w <- init_grid(gparam_w) # 所得wを示すグリッド

# 生涯効用関数
# args: a:貯蓄, w:所得
luf <- function(a, w){
	cons_0 <- w - a # 0期の消費
	if(cons_0 <= 0){
		return(-1e+10)
	}
	cons_1 <- (1.0 + dparam$nir)*a # 1期の消費
	utility <- CRRA(cons_0) + dparam$beta*CRRA(cons_1)
	utility
}

# 所得wごとの貯蓄aを求めてデータフレームにまとめる
df_fig2 <- data.frame(
	 wage = grid_w,
	 
	 policy_nlm = sapply(grid_w, function(w){
		 # 最適化アルゴリズムで効用関数の符号を反転させた目的関数の最小値を探す
		 # なおoptimとoptimize関数はwが小さいときに計算に失敗する ← 引数の問題？
		 r_nlm <- nlm(function(a){
			 -1*luf(a, w)
		 }, gparam_a$min)
		 # 生涯効用関数を最大化する貯蓄a
		 r_nlm$estimate
	 })
)

# 図2をプロットする
par(mar=c(4, 4.5, 2, 1), bg="white")
plot(policy_nlm ~ wage, xlab=expression(paste("Income: ",w)), ylab=expression(paste("Saving: ",a)), data=df_fig2, type="o", pch=21, bg='white', ylim=c(0, 0.5))

# プロットした図をepsで保存する
dev.copy2eps(file="fig 2.eps", width=6, height=4)


