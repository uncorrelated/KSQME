#
# 経済セミナー連載「定量的マクロ経済学と数値計算」
# 第4回「オイラー方程式と多項式近似」のソースコードのRへの野良移植: 分権経済と時間反復法
#

source("common.R") # 共通するパラメーターや関数などを読み込む

gparam_k <- list(number = 21, min=0.05, max=0.5)
grid_k <- init_grid(gparam_k) # 資本ストックを表すグリッド

# 計算目的である資本ストックに対応した消費を定める政策関数を表すグリッド
grid_policy <- grid_k # 初期値は資本ストックに等しくする

msteps <- numeric(maxit) # 更新幅を保存するベクトル

for(i in 1:maxit){
	# 資本ストックごとの生産量（f(k)） + 残資本ストック = 富を保存
	grid_wealth <- with(dparam, { grid_k^alpha + (1 - delta)*grid_k })

	# 政策関数をスプライン関数で近似（splinefunを使う方がよいかも）
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
err_ti <- with(dparam, {
	# 精度評価用グリッドのパラメーター（postfixのatはaccuracy test）
	gparam_k_at <- with(gparam_k,{ list(number = 20*(number - 1) + 1, min=min, max=max) }) 
	grid_k0_at <- init_grid(gparam_k_at)

	# 政策関数をスプライン関数で近似はループ中ですでに計算している
	grid_c0_at <- predict(apf, grid_k0_at)$y
	
	LHS <- grid_c0_at^(-gamma) # 消費の限界効用を計算

	grid_k1_at <- grid_k0_at^alpha + (1 - delta)*grid_k0_at - grid_c0_at
	grid_c1_at <- predict(apf, grid_k1_at)$y

	# 資本レンタルコスト
	r <- alpha*grid_k1_at^(alpha - 1) - delta
	RHS <- beta*(1 + r)*grid_c1_at^(-gamma)

	data.frame(
		kgrid = grid_k0_at,
		error = RHS / LHS - 1
	)
})

# Dynamic Programmingの誤差データのファイルがあればロードする
if(exists("err_pdp")) rm(err_pdp)
path_err_pdp <- "../CH3/err_pdp.csv"
if(file.exists(path_err_pdp)){
	err_pdp <- read.table(path_err_pdp, sep=",", header=TRUE)

	# 2つを同時にプロットできるように準備
	dev.new(width=10, height=5)
	par(oma=c(0, 0, 0, 0), mfrow=c(1, 2), mar=c(4.5, 4.5, 1, 1), bg="white")
} else {
	# DPの誤差を保存したファイルが無ければ、プロット1つにあわせる
	par(mar=c(4, 5, 1, 1), bg="white")
}

# 図1(a)をプロットする
col <- c("blue", "red", "black")
lty <- c(1, 2, 3)
lwd <- c(2, 2, 1)
plot(grid_k, grid_saving, col="blue", xlab=expression(paste("capital stock for the this term: ", k[t])), ylab=expression(paste("capital stock for the next term: ", k=k[t+1])), xlim=c(gparam_k$min, 0.5), ylim=c(0,0.5), type="l", lty=lty[1], lwd=lwd[2])
lines(grid_k, grid_saving_by_closed_form, col="red", lty=lty[2], lwd=lwd[2])
lines(grid_k, grid_k, col="black", lty=lty[3], lwd=lwd[3])
legend("topleft", c("approximate solution", "closed-form solution", "45-degree line"), col=col, lty=lty, lwd=lwd, bg='white', bty="n", seg.len=3)

# Dynamic Programmingとの比較図
if(exists("err_pdp")){
	# 図1(b)をプロットする
	lty <- c(1, 2)
	lwd <- c(2, 2)
	col <- c("black", "red")
	legends <- c("Time Iteration", "Value Function Iteration")
	plot(err_ti$kgrid, abs(err_ti$error), main="", xlab=expression(paste("Current Priod Assets: ", k[t])), ylab="absolute value of error of the Euler equation", type="l", lty=lty[1], lwd=lwd[1], col=col[1])
	lines(err_pdp$kgrid, abs(err_pdp$error), lty=lty[2], lwd=lwd[2], col=col[2])
	legend("topright", lty=lty, lwd=lwd, col=col, legend=legends, bty="n", y.intersp=1.1, seg.len=2.5, inset=0.02)

} else {
	# DPの誤差を保存したファイルが無ければプロットを一つepsで保存する
	dev.copy2eps(file="fig 1-a.eps", width=6, height=4)
}

