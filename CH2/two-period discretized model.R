#
# 経済セミナー連載「定量的マクロ経済学と数値計算」
# 第2回「2期間モデルと数値計算の概観」のソースコードのRへの野良移植
#

#
# 計算用グリッドを初期化
#

# グリッドのパラメーター
# number:グリッドの数, min:最小値, max:最大値
gparam_w <- list(number = 10, min=0.1, max=1.0)
gparam_a <- list(number = 40, min=0.025, max=1.0)

# グリッド生成関数
init_grid <- function(gparam){
	with(gparam,{
		seq(min, max, length.out=number)
	})
}

grid_w <- init_grid(gparam_w) # 所得w計算用
grid_a <- init_grid(gparam_a) # 貯蓄a計算用

#
# ディープ・パラメーター
# beta:割引因子, gamma:相対的危険回避度, nir: 純利子率
# 
dparam <- list(beta=0.985^30, gamma=2.0, nir=1.025^30-1.0)

#
# 相対的リスク回避度一定効用関数（Constant Relative Risk Aversion Utility Function）
# args: 消費, 相対的危険回避度
# return: 効用水準
# 
CRRA <- function(cons, gamma=dparam[["gamma"]]){
	if(1.0==gamma){
		log(cons)
	} else {
		cons^(1-gamma)/(1-gamma)
	}
}

# 所得wと貯蓄aの組み合わせに応じた消費と効用を計算
cons_0 <- rep(grid_w, each=gparam_a$number) - rep(grid_a, gparam_w$number) # 0期の消費
cons_1 <- (1.0 + dparam$nir)*rep(grid_a, gparam_w$number) # 1期の消費
utility <- CRRA(cons_0) + dparam$beta*CRRA(cons_1)
utility[0>=cons_0] <- ncu <- -1e+10 # 0期の消費が負のときの効用は-1e+10

# 目的関数を代理する行列
# 列ごとに所得w、行ごとに貯蓄aが異なる
mobjf <- matrix(utility, gparam_a$number, gparam_w$number)

# 列（所得w）ごとに目的関数が最大になる行（貯蓄a）の位置を探す
gp_a <- apply(mobjf, 2, which.max) 

# 列（所得w）ごとに目的関数が最大になる貯蓄aを求める
optimal_a <- sapply(1:gparam_w$number, function(gp_w){ mobjf[gp_a[gp_w], gp_w] })

# データフレームに結果をまとめる
df_fig1_a <- data.frame(saving=grid_a) # 図1(a)
for(w in c(0.5, 0.8, 1.0)){
	gp_w <- which(w == grid_w)
	cname <- sprintf("w=%.1f", w)
	df_fig1_a[cname] <- mobjf[, gp_w]
	df_fig1_a[ncu==mobjf[, gp_w], cname] <- NA # 消費が負のとき
}

df_fig1_b <- data.frame(wage=grid_w, policy=optimal_a) # 図1(b)

# 図1(a)をプロットする
pch <- c(21, 24, 22)
col <- c("blue", "darkgreen", "red")
lty <- c(1, 1, 1)
par(mar=c(4, 5, 1, 1), bg="white") # フォントを指定しないとエラー
plot(df_fig1_a$saving, df_fig1_a[,colnames(df_fig1_a)[2]], xlab=expression(paste("Saving: ",a)), ylab=expression(paste("Lifetime Utility: ",U(c[1],c[2]))), ylim=c(-10, 0), type="o", lty=lty[1], pch=pch[1], col=col[1], bg='white')
lines(df_fig1_a$saving, df_fig1_a[,colnames(df_fig1_a)[3]], type="o", lty=lty[2], pch=pch[2], col=col[2], bg='white')
lines(df_fig1_a$saving, df_fig1_a[,colnames(df_fig1_a)[4]], type="o", lty=lty[3], pch=pch[3], col=col[3], bg='white')
legend("topright", colnames(df_fig1_a)[2:4], col=col, lty=lty, lwd=1, pch=pch, pt.bg='white', bg='white', bty="n")

# プロットした図をepsで保存する
dev.copy2eps(file="fig 1-a.eps", width=6, height=4)

# 図1(b)をプロットする
par(mar=c(4, 4, 2, 1), bg="white")
plot(policy ~ wage, main="Policy Function", xlab=expression(paste("Income: ",w)), ylab=expression(paste("Saving: ",a)), data=df_fig1_b, type="o", pch=21, bg='white')

# プロットした図をepsで保存する
dev.copy2eps(file="fig 1-b.eps", width=6, height=4)

