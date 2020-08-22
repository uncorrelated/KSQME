#
# 経済セミナー連載「定量的マクロ経済学と数値計算」
# 第5回「ニューケインジアン・モデルの新展開」のソースコードのRへの野良移植: 準線形モデル（時間反復法）, 図3
#
# Replication of Adam and Billi (2007,JME)
# January 2018, Takeki Sunakawa
#
# ported to R by @uncorrelated
#

# ディープパラメーター
dparam <- list(
	r_star = 3.5/4,
	sigma = 6.0,
	alpha = 0.66,
	theta = 7.66,
	omega = 0.47,
	lambda = 0.048/16
)

dparam$beta = 1/(1 + dparam$r_star/100)
dparam$kappa = with(dparam, (1-alpha)*(1-alpha*beta)/alpha*(1/sigma+omega)/(1+omega*theta) )

# 需要者{の状態/へのショック}のパラメーター
gparam <- list(
	rho = 0.8,
	sigma = 1.524,
	N = 31
)

# 供給者{の状態/へのショック}のパラメーター
uparam <- list(
	rho = 0.0,
	sigma = 1e-5,
	N = 31
)

# {状態/ショック}グリッド生成関数を読み込む
source("tauchen.R")

# 需要者{の状態/へのショック}を代理できるグリッドを作成
g_r <- with(gparam, tauchen(N, dparam$sigma*dparam$r_star, rho, sigma, 3.0))
# 供給者{の状態/へのショック}を代理できるグリッドを作成
u_r <- with(uparam, tauchen(N, 0, rho, sigma, 3.0)) 

# 2つの状態遷移行列のクロネッカー積をとり、2つの{状態/ショック}をまとめて扱えるようにする
TM <- g_r$TransitionMatrix %x% u_r$TransitionMatrix
Grid <- matrix(c(rep(g_r$Grid, each=uparam$N), rep(u_r$Grid, gparam$N)), ncol(TM), 2)

policy <- with(dparam, {

	# 政策関数をあらわすグリッド
	policy_new <- policy_old <- matrix(0, ncol(TM), 3)
	colnames(policy_new) <- colnames(policy_old) <- c("y","pi","r")

	tol <- 1e-06 # 最小更新幅
	maxit <- 1000 # 最大ループ回数
	msteps <- numeric(maxit) # 更新幅を保存するベクトル

	for(i in 1:maxit){
		# 古い政策関数から期待値を計算
		E <- TM %*% policy_old

		# 期待値を所与として最適化
        	pi_s <- (beta*E[,2] + Grid[, 2]) / (1 + kappa^2 / lambda)
		y_s <- (-kappa / lambda) * pi_s
		r_s <- (1 / sigma)*(E[, 1] - y_s + Grid[, 1]) + E[, 2] # 記事ではσ=1が暗に仮定されていた

		# マイナス金利のときは、ゼロ金利として再計算
		f <- 0 > r_s # fには条件文（と言うか論理式）の結果になるTRUE/FALSEのベクトルが入る
		y_s[f] <- (E[, 1] - sigma*(0 - E[, 2]) + Grid[, 1])[f] # 記事ではσ=1が暗に仮定されていた
		pi_s[f] <- (kappa*y_s + beta*E[, 2] + Grid[, 2])[f]
		r_s[f] <- 0

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

# クロネッカー積でまとまったショックへの政策関数のグリッドなので、分解しないと把握できない
# 行が需要者{の状態/へのショック}、列が供給者の{の状態/へのショック}に対応する
policy_y <- matrix(policy[, 1], gparam$N, uparam$N)
policy_pi <- matrix(policy[, 2], gparam$N, uparam$N)
policy_r <- matrix(policy[, 3], gparam$N, uparam$N)

# 雑誌記事の図3に相当する図をプロット
# replicate Figures 4-5 in the paperとのこと
i <- ceiling(length(g_r$Grid) / 2);
g <- g_r$Grid/dparam$sigma
xlab <- "The natural rate of interest"
par(oma=c(0, 0, 0, 0), mfrow=c(3,1), mar=c(4.5, 4.5, 1, 1), bg="white") 
plot(g, policy_y[i, ] / 4, xlab=xlab, ylab="output gap: y", type="l") # 四半期データにするための /4 なので、本文の文脈では無くてよい
plot(g, policy_pi[i, ], xlab=xlab, ylab=expression(paste("inflation rate: " , pi)), type="l")
plot(g, policy_r[i, ], xlab=xlab, ylab=expression(paste("policy interest rate: " , r[n])), type="l")

# 図を保存する
dev.copy2eps(file="fig 3.eps", width=6, height=4)

