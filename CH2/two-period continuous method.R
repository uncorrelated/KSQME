#
# 経済セミナー連載「定量的マクロ経済学と数値計算」
# 第2回「2期間モデルと数値計算の概観」のソースコードのRへの野良移植: 求める値は連続的な方法, 図2, 代数的解, 最適化アルゴリズム, 求根法
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

df_others <- data.frame(
	wage = grid_w,

	# 代数的解を求める
	policy_closed_form =  with(dparam, {
		1.0/(1.0+(beta*(1+nir))^(-1/gamma)*(1+nir))*grid_w
	}),
	 
	policy_uniroot = sapply(grid_w, function(w){
		# 求根法で0期と1期の消費の限界代替率が1になる貯蓄aを探す
		r_uniroot <- uniroot(function(a){
			euler_eq(a, w) - 1
		}, upper=gparam_a$max, lower=gparam_a$min)
		# 生涯効用関数を最大化を満たす貯蓄a
		r_uniroot$root
	})
)

# 図2をプロットする
par(mar=c(4, 4.5, 2, 1), bg="white")
plot(policy_nlm ~ wage, xlab=expression(paste("Income: ",w)), ylab=expression(paste("Saving: ",a)), data=df_fig2, type="o", pch=21, bg='white', ylim=c(0, 0.5))

# プロットした図をepsで保存する
dev.copy2eps(file="fig 2.eps", width=6, height=4)


