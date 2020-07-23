#
# 経済セミナー連載「定量的マクロ経済学と数値計算」
# 第4回「オイラー方程式と多項式近似」のソースコードのRへの野良移植: 分権経済と時間反復法
#

# ディープパラメーター
dparam <- list(
	beta  = 0.96, # 割引因子
	gamma = 1.0, # 相対的危険回避度(異時点間の代替の弾力性の逆数)
	alpha = 0.40, # 資本分配率
	delta = 1.00 # 固定資本減耗(delta=1.0のときは解析解が存在)
)

gparam_k <- list(number = 21, min=0.05, max=0.5)

# グリッド生成関数
init_grid <- function(gparam){
       with(gparam,{
               seq(min, max, length.out=number)
       })
}

grid_k <- init_grid(gparam_k) # 資本ストックを表すグリッド

# 計算目的である資本ストックに対応した消費を定める政策関数を表すグリッド
grid_policy <- grid_k # 初期値は資本ストックに等しくする

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
msteps <- numeric(maxit) # 更新幅を保存するベクトル

for(i in 1:maxit){
	# 資本ストックごとの生産量（f(k)） + 残資本ストック = 富を保存
	grid_wealth <- with(dparam, { grid_k^alpha + (1 - delta)*grid_k })

	# 政策関数をスプライン関数で近似
	apf <- smooth.spline(grid_k, grid_policy)

	# 富グリッドのセルごとに最適な消費量を求める
	grid_c <- sapply(grid_wealth, function(w){
		# Rのビルトイン関数の求根法を用いる（ので初期値がいらない）
		r_uniroot <- uniroot(function(c0){
			# 今期の資源制約式と今期の消費、来期の資本ストックを計算（負の値をとらないようにする）			
			k1 <- max(w - c0, gparam_k$min)
			# 政策関数の近似から、来期の消費を計算
			c1 <- predict(apf, k1)$y
			euler_eq(c1, k1, c0)
		}, upper=w, lower=0) # 消費は0〜富の間になる
		r_uniroot$root
	})

	# 繰り返し計算誤差と言うか、grid_policyとgrid_cの各セルの差の最大値を求める
	msteps[i] <- max(abs(grid_c - grid_policy))

	# 計算された最適消費で政策関数を更新
	grid_policy <- grid_c

	# 更新幅がtol未満であれば終了
	if(tol > msteps[i]){
		break;
	}
}

# 計算された政策関数から、貯蓄関数を表すグリッドを計算
grid_saving = with(dparam, { grid_k^alpha + (1 - delta)*grid_k - grid_policy})

# 貯蓄関数の代数的解
grid_saving_by_closed_form = with(dparam, { beta*alpha*(grid_k^alpha) })

# 資本ストックのグリッドを細かくしてオイラー方程式で計算誤差を求める
err <- with(dparam, {
	# 精度評価用グリッドのパラメーター（postfixのatはaccuracy test）
	gparam_k_at <- with(gparam_k,{ list(number = 20*(number - 1) + 1, min=min, max=max) }) 
	grid_k0_at <- init_grid(gparam_k_at)

	# 政策関数をスプライン関数で近似はループ中ですでに計算している
	# apf <- smooth.spline(grid_k, grid_policy)
	grid_c0_at <- predict(apf, grid_k0_at)$y

	LHS <- grid_c0_at^gamma # 消費の限界効用を計算

	grid_k1_at <- grid_k0_at^alpha + (1 - delta)*grid_k0_at - grid_c0_at
	grid_c1_at <-predict(apf, grid_k1_at)$y

	# 資本レンタルコスト
	r <- alpha*grid_k1_at^(alpha - 1) - delta
	RHS <- beta*(1 + r)*grid_c1_at^gamma

	RHS / LHS - 1
})

# 図1(a)をプロットする
par(mar=c(4, 5, 1, 1), bg="white")
col <- c("blue", "red", "black")
lty <- c(1, 2, 3)
lwd <- c(2, 2, 1)
plot(grid_k, grid_saving, col="blue", xlab=expression(paste("capital stock for the this term: ", k[t])), ylab=expression(paste("capital stock for the next term: ", k=k[t+1])), xlim=c(gparam_k$min, 0.5), ylim=c(0,0.5), type="l", lty=lty[1], lwd=lwd[2])
lines(grid_k, grid_saving_by_closed_form, col="red", lty=lty[2], lwd=lwd[2])
lines(grid_k, grid_k, col="black", lty=lty[3], lwd=lwd[3])
legend("topleft", c("approximate solution", "closed-form solution", "45-degree line"), col=col, lty=lty, lwd=lwd, bg='white', bty="n", seg.len=3)
# プロットした図をepsで保存する
dev.copy2eps(file="fig 1-a.eps", width=6, height=4)

# Dynamic Programmingとの比較図はDPのコードを書いていないので省略

