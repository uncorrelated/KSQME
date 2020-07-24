#
# 経済セミナー連載「定量的マクロ経済学と数値計算」
# 第4回「オイラー方程式と多項式近似」のソースコードのRへの野良移植: 図4
#

f <- function(x) 1/(1 + 25*x^2) # 近似対象のルンゲ関数
fn <- expression(frac(1,1+25*x^2)) # プロットの表記

# 引数xの定義域
xmin <- -1
xmax <- 1

# xを生成
x <- seq(xmin, xmax, length.out=1001)

# xに対応するyを計算
y <- f(x)


###### スプライン関数による近似 ######

# 等間隔に評価点を次元数とる
cp_x <- seq(xmin, xmax, length.out=11)
cp_y <- f(cp_x)

# スプライン関数で予測を行う
r_ss <- smooth.spline(cp_x, cp_y)
yp <- predict(r_ss, x)$y

# プロットする
par(mar=c(4.5, 4.5, 1, 1), bg="white")
lwd <- c(1, 2, 3)
lty <- c(1, 2, -1)
col <- c("black", "red", "red")
pch <- c(-1, -1, 21)
cex <- c(-1, -1, 2)
plot(y ~ x, type="l", col=col[1], lwd=lwd[1], lty=lty[1], ylim=c(-0.25, 1.5))
lines(yp ~ x, col=col[2], lwd=lwd[2], lty=lty[2])
points(cp_y ~ cp_x, col=col[3], pch=pch[3], lwd=lwd[3], cex=cex[3])
legend("topright", c(fn, "Spline", "Collocation points"), col=col, lwd=lwd, lty=lty, pch=pch, pt.cex=cex, bg='white', bty="n", seg.len=3)

dev.copy2eps(file="spline approximation.eps", width=6, height=4)


###### 通常の多項式近似 ######

d <- 10 # 10次

# 等間隔に評価点を次元数とる
cp_x <- seq(xmin, xmax, length.out=d + 1)
cp_y <- f(cp_x)

X <- sapply(0:d, function(i){
	cp_x^i
})
b <- solve(X, cp_y)

# xに対する予測値をつくる
X <- sapply(0:d, function(i){
	x^i
})
yp <- X %*% b

# プロットする
par(mar=c(4.5, 4.5, 1, 1), bg="white")
lwd <- c(1, 2)
lty <- c(1, 2)
col <- c("black", "red")
plot(y ~ x, type="l", col=col[1], lwd=lwd[1], lty=lty[1], ylim=c(-0.25, 2))
lines(yp ~ x, col=col[2], lwd=lwd[2], lty=lty[2])
legend("top", c(fn, "Ordinary polynomial"), col=col, lwd=lwd, lty=lty, bg='white', bty="n", seg.len=3)

dev.copy2eps(file="fig 2-b.eps", width=6, height=4)


###### チェビシェフの多項式近似 ######

# チェビシェフ多項式近似のための関数群を読み込む
source("chebyshev approximation.R")

#
# 引数で指定した関数とそのチェビシェフ多項式近似をプロットする
# args: d:字数, f:関数, fn:関数名（凡例での表記）, x:引数
#
plotChebyshevApproximation <- function(d, f, fn, x, y=f(x)){
	# チェビシェフ多項式近似
	r_ca <- ChebyshevApproximation(f, xmin, xmax, d)
	# パラメーターΘの推定と、それからの予測の関数が分かれているので、推定パラメーターを使いまわす事ができるが、今回は使い捨て
	yp <- ChebyshevPredict(r_ca, x)
	
	par(mar=c(4.5, 4.5, 1, 1), bg="white")
	lwd <- c(1, 2, 2)
	pch <- c(-1, -1, 4)
	col <- c("black", "blue", "blue")
	cex <- c(-1, -1, 2)
	lty <- c(1, 1, -1)
	plot(y ~ x, type="l")
	lines(yp ~ x, col=col[2], lwd=lwd[2])
	points(r_ca$cp_y ~ r_ca$cp_x, col=col[3], pch=pch[3], lwd=lwd[3], cex=cex[3])
	legend("topright", c(fn, "Chebyshev polynomial", "Collocation points"), col=col, lty=lty, lwd=lwd, pch=pch, pt.cex=cex, bg='white', bty="n", seg.len=3)
}

# ルンゲ関数を10次で近似
plotChebyshevApproximation(10, f, fn, x, y)

# プロットした図をepsで保存する
dev.copy2eps(file="fig 4-a.eps", width=6, height=4)

# ルンゲ関数を20次で近似
plotChebyshevApproximation(20, f, fn, x, y)

# プロットした図をepsで保存する
dev.copy2eps(file="fig 4-b.eps", width=6, height=4)

