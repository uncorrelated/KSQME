#
# 経済セミナー連載「定量的マクロ経済学と数値計算」
# 第2回「2期間モデルと数値計算の概観」のソースコードのRへの野良移植: 射影法, 図3
#

# 他の手法と共通する設定などを読み込む
source("common.R")

# 計算用グリッドを初期化
grid_w <- init_grid(gparam_w) # 所得wを示すグリッド
grid_a <- init_grid(gparam_a) # 貯蓄a計算用

#
# 多項式近似関数（Polynomial Approximation Function）
# args: 所得w, パラメーターtheta
# return: 貯蓄a
#
paf <- function(w, theta){
	dim <- length(theta) - 1
	p <- 0:dim
	sapply(w, function(w){
		sum(theta*w^p)
	})
}

#
# オイラー条件との乖離を表す残差を計算
#
residuals <- function(a, w){
	euler_eq(a, w) - 1
}

#
# 距離関数Ρ; pafとresidualsとrhoをまとめた関数を作る方が速い気がするが、見通しを優先
#
rho <- function(theta){
	a <- paf(grid_w, theta)
	sum(residuals(a, grid_w)^2)
}


# 基底ΨがないせいかOLSで推定した方が早くて確実感があるものの、Nelder-Mead法で計算する
# なお、レーベンバーグ・マルカート法はminpack.lmパッケージで使える
init <- c(0.01, 0.2) # 係数の初期値を当て推量(initial guess)と説明があった初期値. {ここに並べたパラメーターの数-1}次の多項式になる
r_optim <- optim(init, rho, method="Nelder-Mead")

# 所得wごとの貯蓄aを求めてデータフレームにまとめる
df_projection <- data.frame(
	wage = grid_w,
	policy_projection = paf(grid_w, r_optim$par)
)

# データフレームに所得wを定めたときの貯蓄と残差の関係をまとめる
df_fig3 <- data.frame(a=seq(0.1, 0.5 - 1e-06, length.out=50)) # 図3
for(w in c(0.5, 0.8, 1.0)){
	cname <- sprintf("w=%.1f", w)
	df_fig3[cname] <- residuals(df_fig3$a, w)
}

# 図3をプロットする
pch <- c(21, 24, 22)
col <- c("blue", "darkgreen", "red")
lty <- c(1, 1, 1)
par(mar=c(4, 5, 1, 1), bg="white")
plot(df_fig3$a, df_fig3[,colnames(df_fig3)[2]], xlab=expression(paste("Saving: ",a)), ylab=expression(paste("Residuals: ",R(w,theta))), ylim=c(-1, 1), type="o", lty=lty[1], pch=pch[1], col=col[1], bg='white')
lines(df_fig3$a, df_fig3[,colnames(df_fig3)[3]], type="o", lty=lty[2], pch=pch[2], col=col[2], bg='white')
lines(df_fig3$a, df_fig3[,colnames(df_fig3)[4]], type="o", lty=lty[3], pch=pch[3], col=col[3], bg='white')
legend("topright", colnames(df_fig3)[2:4], col=col, lty=lty, lwd=1, pch=pch, pt.bg='white', bg='white', bty="n")

x_margin <-  (par()$usr[2]-par()$usr[1])/100
abline(h=0.0, lty=3)
abline(v=0.17, lty=3)
text(0.17 - x_margin, -1, expression(a=0.17), adj=c(1,0))

# プロットした図をepsで保存する
dev.copy2eps(file="fig 3.eps", width=6, height=4)

