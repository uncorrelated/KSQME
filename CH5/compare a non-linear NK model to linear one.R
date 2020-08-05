#
# 経済セミナー連載「定量的マクロ経済学と数値計算」
# 第5回「ニューケインジアン・モデルの新展開」のソースコードのRへの野良移植: Githubにあげられている非線形モデルと線形モデルの比較プロット
# 

# 非線形モデルを解く
source("non-linear NK model index-function approach.R")

# 線形モデルを解く
source("linear NK model to compare.R")

attach(dparam)

# 横軸になる前期金利; 線形モデル、非線形モデルの両方の計算に利用
n <- 101
Rpast <- with(bound,seq(Rmin, Rmax, length.out=n))

#
# 線形モデル、非線形モデルともに、推定パラメーターからそれぞれの予測値(から定常状態の値で割って対数化して100倍した値)を計算
#

linear <- with(bound, {
	c <- rep(NA, n)
	y <- numeric(n)
	pi <- numeric(n)
	R <- numeric(n) # Juliaのコードではrnに対応

	for(i in 1:n){
		gnow <- 0
		znow <- zmin
		rnow <- 0

		X0 = numeric(8)
		X0[3] = log(Rpast[i]/ss$R)
		X0[5] = gnow
		X0[6] = znow    

		Z0 = numeric(3) # shocks
		X1 = Re(r_gensys$G1 %*% X0 + r_gensys$impact %*% Z0)

		R[i] = X1[3] * 100
		pi[i] = X1[2] * 100
		y[i] = X1[1] * 100
	}

	data.frame(Rpast=Rpast, c=c, y=y, pi=pi, R=R)
})

non_linear <- with(bound, {
	c <- numeric(n)
	y <- numeric(n)
	pi <- numeric(n)
	R <- numeric(n) # Juliaのコードではrnに対応

	for(i in 1:n){
		gnow <- 0
		znow <- zmin
		rnow <- 0

		y_star <- calcY_star(gnow)

		gp  <- rho_g*gnow
		zp  <- rho_z*znow
		rp  <- rep(0, length(rnow))

		X <- matrix(c(Rpast[i], gp, zp, rp), 4, byrow=TRUE) * slopecon[, 1] + slopecon[, 2]
		rownames(X) <- c("R", "gp", "zp", "rp")
		p2s <- poly2s(t(X))
		P <- p2s %*% r_solve$coef

		fc1 <- P[, "c0n"] # coeffcn
		fp1 <- P[, "pi0n"] # coeffpn

		c[i] <- calcC(fc1)
		pi[i] <- calcPi(c[i], fp1)
		y[i] <- calcY(c[i], pi[i], gp)
	#	print(sprintf("y:%f, pi:%f, Rpast:%f, rnow:%f, ystar:%f", y[i], pi[i], Rpast[i], rnow, y_star))
		R[i] <- calcR(y[i], pi[i], Rpast[i], rnow, y_star)

		if (is.na(R[i]) || R[i]<1.0){
			fc1 <- P[, "c0b"] # coeffcb
			fp1 <- P[, "pi0b"] # coeffpb

			c[i] <- calcC(fc1)
			pi[i] <- calcPi(c[i], fp1)
			y[i] <- calcY(c[i], pi[i], gp)
			R[i] <- 1
		}
	}
	data.frame(Rpast=Rpast, c=log(c/ss$c) * 100, y=log(y/ss$y) * 100, pi=log(pi/ss$pi) * 100, R=log(R/ss$R) * 100)
})

detach(dparam)

# プロットする
par(oma=c(0, 0, 0, 0), mfrow=c(3,1), mar=c(4.5, 4.5, 1, 1), bg="white") # mfrow=c(3,1)は、3つのプロットを横3分割して表示するためのオプション

drawPlot <- function(ylab, y1, y2, lpos){
	col <- c("red", "blue")
	lwd <- c(2, 2)
	lty <- c(1, 2)
	labels <- c("non-linear NK model", "linear NK model")
	x <- log(Rpast/ss$R)*100
	ylim <- c(floor(min(y1, y2)), ceiling(max(y1, y2)))
	plot(x, y1, xlab="", ylab=ylab, type="l", ylim=ylim, col=col[1], lty=lty[1], lwd=lwd[1])
	lines(x, y2, col=col[2], lty=lty[2], lwd=lwd[2])
	legend(lpos, lty=lty, lwd=lwd, col=col, legend=labels, bty="n", y.intersp=1.5, seg.len=3, inset=0.05)
}

drawPlot("Policy Rate", non_linear$R, linear$R, "bottomright")
drawPlot("Inflation", non_linear$pi, linear$pi, "topright")
drawPlot("Output Gap", non_linear$y, linear$y, "topright")

# 図を保存する
dev.copy2eps(file="fig X1-comparing.eps", width=6, height=9)
# dev.copy(png, "fig X1-comparing.png", width=600, height=900, type="cairo", bg="white"); dev.off()

