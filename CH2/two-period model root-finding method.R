#
# 経済セミナー連載「定量的マクロ経済学と数値計算」
# 第2回「2期間モデルと数値計算の概観」のソースコードのRへの野良移植: 求根法
#

# 他の手法と共通する設定などを読み込む
source("common.R")

# 計算用グリッドを初期化
grid_w <- init_grid(gparam_w) # 所得wを示すグリッド

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

# ソースを読み込んでも無反応なので、結果を表示するようにする
print(df_others)

