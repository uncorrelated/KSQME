#
# 経済セミナー連載「定量的マクロ経済学と数値計算」
# 第5回「ニューケインジアン・モデルの新展開」のソースコードのRへの野良移植: 2状態モデル（時間反復法）, 図1
#

# 共通の設定を読み込む
source("common.R")

#
# 時間反復法で計算
# args: p_H:状態H→Lへの転落確率
#
time_iteration_method <- function(p_H){
	with(dparam, {
		# 状態の種類を表すグリッド
		grid_states <- matrix(c(s_H, s_L), 2, 1)
		names(grid_states) <- c("H", "L")

		# マルコフ連鎖に使う遷移行列
		transition <- matrix(c(1-p_H, p_H, 1-p_L, p_L), 2, 2, byrow=TRUE)
		colnames(transition) <- c("(1-p)", "p")
		rownames(transition) <- c("H", "L")

		# 政策関数をあらわすグリッド
		policy_new <- policy_old <- makePolicyMatrix(0, 0, 0, 0, 0)

		tol <- 1e-06 # 最小更新幅
		maxit <- 5000 # 最大ループ回数
		msteps <- numeric(maxit) # 更新幅を保存するベクトル

		for(i in 1:maxit){
			# 古い政策関数から期待値を計算
			E <- transition %*% policy_old

			# 期待値を所与として最適化
			r_s = pmax(r_star + phi*E[,2], 0)
			y_s = E[,1] - (r_s - E[,2] - grid_states);
			pi_s = kappa*y_s + beta*E[,2]

			# 新しい政策関数を保存
			policy_new[,] <- c(y_s, pi_s, r_s)

			# 繰り返し計算誤差と言うか、policy_newとpolicy_oldの各セルの差の最大値を求める
			msteps[i] <- max(abs(policy_new - policy_old))

			# 政策関数を更新
			policy_old <- policy_new

			# 更新幅がtol未満であれば終了
			if(tol > msteps[i]){
				break;
			}
		}

		# 収束しなかった場合は止める
		if(maxit <= i){
			stop("iteration limit exceeded.")
		}

		policy_new
	})
}

r_tim_000 <- time_iteration_method(0)
r_tim_025 <- time_iteration_method(0.025)

plot_policy <- function(i){
	par(mar=c(4.5, 4.5, 1, 1), bg="white")
	ylab <- c("output gap", "inflation rate", "policy interest rate")[i]
	x <- c(0, 1)
	ylim <- c(floor(min(r_tim_000[, i], r_tim_025[, i])), ceiling(max(r_tim_000[, i],r_tim_025[, i])))
	plot(x, r_tim_000[, i], axes=FALSE, xlab="State", ylab=ylab, type="l", ylim=ylim)
	points(x, r_tim_000[, i], pch=19, col="black")
	lines(x, r_tim_025[, i], lty=2)
	points(x, r_tim_025[, i], pch=20, col="black", bg="white")
	axis(1, at=x, labels=c("H", "L"))
	y_at <- seq(ylim[1], ylim[2], 1)
	axis(2, at=y_at, labels=y_at)
}

plot_policy(1)
dev.copy2eps(file="fig 1-y.eps", width=6, height=4)

plot_policy(2)
dev.copy2eps(file="fig 1-pi.eps", width=6, height=4)

plot_policy(3)
dev.copy2eps(file="fig 1-r.eps", width=6, height=4)

