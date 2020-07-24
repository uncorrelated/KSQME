# ディープパラメーター
dparam <- list(
	beta  = 0.96, # 割引因子
	gamma = 1.0, # 相対的危険回避度(異時点間の代替の弾力性の逆数)
	alpha = 0.40, # 資本分配率
	delta = 1.00 # 固定資本減耗(delta=1.0のときは解析解が存在)
)

# グリッド生成関数
init_grid <- function(gparam){
       with(gparam,{
               seq(min, max, length.out=number)
       })
}

#
# 最適状態であれば戻り値が0なるオイラー方程式
# args: c1:来期消費, k1:来期資本ストック, c0:今期消費
#
euler_eq <- function(c1, k1, c0){
	with(dparam, {
		1/c0 - beta*(1/c1)*(alpha*k1^(alpha-1) + (1-delta))
	})
}

maxit <- 100 # 最大ループ回数
tol <- 1e-6 # 最小更新幅

