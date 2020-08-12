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

	gnow <- 0
	znow <- zmin
	rnow <- 0

	X0 <- matrix(0, 8, n)
	X0[3, ] = log(Rpast/ss$R)
	X0[5, ] = gnow
	X0[6, ] = znow    

	Z0  <- matrix(0, 3, n) # shocks
	X1 = Re(r_gensys$G1 %*% X0 + r_gensys$impact %*% Z0) # 実数部のみを扱う

	data.frame(Rpast=Rpast, c=c, y=X1[1, ] * 100, pi=X1[2, ] * 100, R=X1[3, ] * 100)
})

non_linear <- with(bound, {
	c <- numeric(n)
	y <- numeric(n)
	pi <- numeric(n)
	R <- numeric(n) # Juliaのコードではrnに対応

	gnow <- rep(0, n)
	znow <- rep(zmin, n)
	rnow <- rep(0, n)

	y_star <- calcY_star(gnow)

	gp  <- rho_g*gnow
	zp  <- rho_z*znow
	rp  <- rep(0, length(rnow))

	# non-ZLBとZLBのfc1とfp1を同時に計算
	X <- matrix(c(Rpast, gp, zp, rp), 4, byrow=TRUE) * slopecon[, 1] + slopecon[, 2]
	rownames(X) <- c("R", "gp", "zp", "rp")
	p2s <- poly2s(t(X))
	P <- p2s %*% r_solve$coef

	fc1 <- P[, "c0n"] # coeffcn
	fp1 <- P[, "pi0n"] # coeffpn

	# next period's c and pi (obtained by next period's fc and fp)
	c <- calcC(fc1)
	pi <- calcPi(c, fp1)
	y <- calcY(c, pi, gp)
	R <- calcR(y, pi, Rpast, rnow, y_star)

	flag <- is.na(R) | R<1.0 # ZLBに引っかかったフラグ

	fc1 <- P[, "c0b"] # coeffcb
	fp1 <- P[, "pi0b"] # coeffpb

	c[flag] <- calcC(fc1[flag])
	pi[flag] <- calcPi(c[flag], fp1[flag])
	y[flag] <- calcY(c[flag], pi[flag], gp[flag])
	R[flag] <- 1

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

