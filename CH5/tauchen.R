#
# tauchen (1986) のマルコフ連鎖に使うAR(1)な遷移行列を作成する関数
# xₜ₊₁=c+ρxₜ+ϵₜ （ϵₜ∈𝒩(0,σ²)）
# args:
# 	N: グリッドポイントの数
# 	μ: 定常状態
# 	ρ: 慣性（persistence）
# 	σ: 標準偏差
# 	m: グリッドポイントが代理する区間（標準偏差単位）
# return:
#	G: グリッド
#	P: 遷移行列（transition matrix）
#
tauchen <- function(N, mu, rho, sigma, m){

	P <- matrix(0, N, N) # 遷移確率の行列
	c <- (1 - rho)*mu # 定数項

	# 最大値と最小値
	gmax <- m*sqrt(sigma^2/(1-rho^2)) + mu
	gmin <- -gmax + mu

	# 等間隔のグリッドを定める
	G <- seq(gmin, gmax, length.out=N)

	# グリッド間の間隔
	w <- (gmax - gmin)/(N - 1)

	# 遷移行列を作成
        P[, 1] <- pnorm(G[1] - c - rho*G[1:N] + w/2, sd=sigma)
        P[, 2:(N-1)] <- sapply(2:(N-1), function(i){
		pnorm(G[i] - c - rho*G[1:N] + w/2, sd=sigma) - pnorm(G[i] - c - rho*G[1:N] - w/2, sd=sigma)
        })
	P[, N] <- 1 - pnorm(G[N] - c - rho*G[1:N] - w/2, sd=sigma)

	list(Grid=G, TransitionMatrix=P)
}

