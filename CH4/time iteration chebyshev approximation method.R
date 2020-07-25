#
# 経済セミナー連載「定量的マクロ経済学と数値計算」
# 第4回「オイラー方程式と多項式近似」のソースコードのRへの野良移植: 時間反復法への多項式近似の応用
#

source("common.R") # 共通するパラメーターや関数などを読み込む
source("chebyshev approximation.R")

ngrid <- c(3, 5, 9) # ベンチマークをとるグリッドのサイズ
ngrid <- c(ngrid, 17, 33) # 微妙過ぎて違いが分からないので足す

# ベンチマーク結果を格納するデータフレーム
r_benckmark <- data.frame(
	no.of.grid = ngrid, # グリッドのサイズ
	log10_error.mean = numeric(length(ngrid)), # 平均誤差
	log10_error.max = numeric(length(ngrid)), # 最大誤差
	elapsed.time.absolute = numeric(length(ngrid)), # 経過時間
	elapsed.time.relative = numeric(length(ngrid)) # 相対経過時間
)

for(j in 1:length(ngrid)){
	gc(); gc(); # ベンチマーク中にガーベッジコレクターが動かないように準備
	r_benckmark$elapsed.time.absolute[j] <- system.time({

		n <- ngrid[j]
		gparam_k <- list(number = n, min=0.05, max=0.5) # グリッドのサイズはループごとに変わる
		grid_k <- with(gparam_k, { ChebyshevCollocationPoints(min, max, number) }) # 資本ストックを表すグリッド

		# 計算目的である資本ストックに対応した消費を定める政策関数を表すグリッド
		grid_policy <- grid_k # 初期値は資本ストックに等しくする

		msteps <- numeric(maxit) # 更新幅を保存するベクトル

		for(i in 1:maxit){
			# 資本ストックごとの生産量（f(k)） + 残資本ストック = 富を保存
			grid_wealth <- with(dparam, { grid_k^alpha + (1 - delta)*grid_k })

			# 政策関数をチェビシェフ多項式で近似
			r_ca <- ChebyshevMapApproximation(grid_k, grid_policy, gparam_k$min, gparam_k$max, n - 1)

			# 富グリッドのセルごとに最適な消費量を求める
			grid_c <- sapply(grid_wealth, function(w){
				# Rのビルトイン関数の求根法を用いる（ので初期値がいらない）
				r_uniroot <- uniroot(function(c0){
					# 今期の資源制約式と今期の消費、来期の資本ストックを計算（負の値をとらないようにする）			
					k1 <- max(w - c0, gparam_k$min)
					# 政策関数の近似から、来期の消費を計算
					c1 <- ChebyshevPredict(r_ca, k1)
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

			# 政策関数をチェビシェフ多項式で近似するのはループ中ですでに終えている
			grid_c0_at <- ChebyshevPredict(r_ca, grid_k0_at)

			LHS <- grid_c0_at^gamma # 消費の限界効用を計算

			grid_k1_at <- grid_k0_at^alpha + (1 - delta)*grid_k0_at - grid_c0_at
			grid_c1_at <- ChebyshevPredict(r_ca, grid_k1_at)

			# 資本レンタルコスト
			r <- alpha*grid_k1_at^(alpha - 1) - delta
			RHS <- beta*(1 + r)*grid_c1_at^gamma

			RHS / LHS - 1
		})

		# 精度を保存
		r_benckmark$log10_error.mean[j] <- log10(mean(abs(err)))
		r_benckmark$log10_error.max[j] <- log10(max(abs(err)))

	})[3] # system.time()の終わり
}

# 相対時間を計算
r_benckmark$elapsed.time.relative <- r_benckmark$elapsed.time.absolute /  r_benckmark$elapsed.time.absolute[1]

# 画面に表示
print(r_benckmark)

