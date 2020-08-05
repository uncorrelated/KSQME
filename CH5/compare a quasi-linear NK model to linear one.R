#
# 経済セミナー連載「定量的マクロ経済学と数値計算」
# 第5回「ニューケインジアン・モデルの新展開」のソースコードのRへの野良移植: Githubにあげられている準線形モデルと線形モデルの比較プロット
# 

# 準線形モデルを解く
source("quasi-linear NK model index-function approach.R")

# 横軸になる前期金利; 準線形モデル、線形モデルの両方の計算に利用
n <- 101
Rpast <- with(bound,seq(Rmin, Rmax, length.out=n))

attach(dparam)
quasi_linear <- with(bound, {
	c <- numeric(n)
	y <- numeric(n)
	pi <- numeric(n)
	R <- numeric(n) # Juliaのコードではrnに対応

	for(i in 1:n){
		gnow <- 0
		znow <- zmin
		rnow <- 0

		# non-ZLBとZLBのfc1とfp1を同時に計算
		X <- matrix(c(Rpast[i], gnow, znow, rnow), 4, byrow=TRUE) * slopecon[, 1] + slopecon[, 2]
		rownames(X) <- c("R", "gp", "zp", "rp")
		p2s <- poly2s(t(X))
		P <- p2s %*% r_solve$coef

		# first assume the ZLB is not binding, and use coeffcn and coeffpn
		fc1 <- P[, "c0n"] # coeffcn
		fp1 <- P[, "pi0n"] # coeffpn

		# next period's c and pi (obtained by next period's fc and fp)
		c1 <- calcC(fc1, fp1, Rpast[i], znow, rnow)
		pi1 <- calcPi(c1, fp1)
		y1 <- calcY(c1, gnow)
		R1 <- calcR(Rpast[i], c1, pi1, rnow)

		if (is.na(R1) || R1 < -1*ss$R){
			c1 <- calcC(fc1, fp1, Rpast[i], znow, rnow, TRUE)
			pi1 <- calcPi(c1, fp1)
			y1 <- calcY(c1, gnow)
			R1 <- -1*ss$R
		}

		R[i] = R1
		pi[i] = pi1
		y[i] = y1
	}
	data.frame(Rpast=Rpast, c=c, y=y, pi=pi, R=R)
})

detach(dparam)

# 線形モデルを解く
source("linear NK model to compare.R")

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
		X0[3] = Rpast[i]
		X0[5] = gnow
		X0[6] = znow    

		Z0 = numeric(3) # shocks
		X1 = Re(r_gensys$G1 %*% X0 + r_gensys$impact %*% Z0)

		R[i] = X1[3]
		pi[i] = X1[2] 
		y[i] = X1[1] 
	}

	data.frame(Rpast=Rpast, c=c, y=y, pi=pi, R=R)
})



# プロットする
par(oma=c(0, 0, 0, 0), mfrow=c(3,1), mar=c(4.5, 4.5, 1, 1), bg="white")  # mfrow=c(3,1)は、3つのプロットを横3分割して表示するためのオプション

drawPlot <- function(ylab, y1, y2, lpos){
	col <- c("red", "blue")
	lwd <- c(2, 2)
	lty <- c(1, 2)
	labels <- c("quasi-linear NK model", "linear NK model")
	x <- Rpast
	ylim <- c(floor(min(y1, y2)), ceiling(max(y1, y2)))
	plot(x, y1, xlab="", ylab=ylab, type="l", ylim=ylim, col=col[1], lty=lty[1], lwd=lwd[1])
	lines(x, y2, col=col[2], lty=lty[2], lwd=lwd[2])
	legend(lpos, lty=lty, lwd=lwd, col=col, legend=labels, bty="n", y.intersp=1.5, seg.len=3, inset=0.05)
}

drawPlot("Policy Rate", quasi_linear$R, linear$R, "bottomright")
drawPlot("Inflation", quasi_linear$pi, linear$pi, "topright")
drawPlot("Output Gap", quasi_linear$y, linear$y, "topright")

# 図を保存する
dev.copy2eps(file="fig X2-comparing.eps", width=6, height=9)
# dev.copy(png, "fig X2-comparing.png", width=600, height=900, type="cairo", bg="white"); dev.off()

