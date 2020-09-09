#
# 経済セミナー連載「定量的マクロ経済学と数値計算」
# 第3回「動的計画法」のソースコードのRへの野良移植: 第2節ロビンソン・クルーソーとベルマン方程式, 図2
#

#
# 最適化(optimization)と内挿法(interpolation)をつかってロビンソン・クルーソー経済を解く.
#

source("common.R")

# パラメーター
param <- list(
	### カリブレーション ###
	beta  = 0.96, # 割引因子
	gamma = 1.0,  # 相対的危険回避度(異時点間の代替の弾力性の逆数)
	alpha = 0.4,  # 資本分配率

	### 離散化用のパラメータ ###
	nk   = 11,    # グリッドの数
	kmax =  1.0,  # 資本グリッドの最大値
	kmin =  0.05, # 資本グリッドの最小値

	# 無人島に滞在する期間
	TT   = 10
)

# 数値演算関数
calcRCModel <- function(param){
	with(param, {

		### グリッドポイントを計算 ###
		kgrid = seq(kmin, kmax, length.out=nk)

		### 変数を定義 ###
		vfcn = matrix(0, nk, TT); # 価値関数
		pfcn = matrix(0, nk, TT); # 政策関数
		cfcn = matrix(0, nk, TT); # 消費関数

		### 最終期(全てを消費) ###
		pfcn[, TT] = 0; # 全て消費するので貯蓄はゼロ
		cfcn[, TT] = kgrid^alpha; # 生産=消費
		vfcn[, TT] = CRRA(cfcn[, TT], gamma); # 消費から得られる効用

		### メインループ ###
		for(t in seq(TT-1, 1, by=-1)){
			cat(sprintf("period %d: \n", t))

			# 次期の価値関数をスプライン関数で近似する
			vnext <- splinefun(kgrid, vfcn[,t + 1])

			for(i in 1:nk){
				# 価値関数が最大値となる来季資本（kprime）を初期値初期値0.01から探す
				r_nlm <- nlm(function(p, k, vnext, alpha, beta, gamma){ -1 * BellmanEq(p, k, vnext, alpha, beta, gamma) }, 0.01, k=kgrid[i], vnext=vnext, alpha=alpha, beta=beta, gamma=gamma)
				pfcn[i, t] <- r_nlm$estimate
				vfcn[i, t] <- r_nlm$minimum * -1
			}

			# 消費関数を計算(必ずしも計算する必要はない)
			cfcn[, t] = kgrid^alpha - pfcn[, t];
		}

		# 代数的解
		p_true <- sapply(1:TT, function(t){
			i <- 1:nk
			alpha*beta*( (1-(alpha*beta)^(TT-t)) / (1-(alpha*beta)^(TT-t+1)) )*(kgrid[i]^alpha)
		})

		list(kgrid=kgrid, TT=TT, p_true=p_true, vfcn=vfcn, pfcn=pfcn, cfcn=cfcn)
	})
}

plotRCModel <- function(r_RCModel){
	with(r_RCModel, {
		plotNorth <- function(title, xlab, ylab, y){
			T <- c(1, 7, 8, 9, 10)
			lty <- c(1, 3, 2, 5, 6)
			lwd <- c(3, 2, 2, 2, 2)
			col <- c("black", "red", "blue", "purple", "black")
			ylim <- c(min(y), max(y)*1.2)
			plot(kgrid, y[, T[1]], main=title, xlab=xlab, ylab=ylab, ylim=ylim, type="l", lty=lty[1], lwd=lwd[1], col=col[1])
			for(i in 2:length(T)){
				lines(kgrid, y[, T[i]], lty=lty[i], lwd=lwd[i], col=col[i])
			}
			legend("bottomright", lty=lty, lwd=lwd, col=col, legend=parse(text=sprintf("t==%d", T)), bty="o", y.intersp=1.1, seg.len=2.5, bg="white", inset=0.02)
		}


		# グラフを描画
		dev.new(width=9, height=9)
		par(oma=c(0, 0, 3, 0), mfrow=c(2, 2), mar=c(4.5, 4.5, 2, 1), bg="white")

		plotNorth("(a) Value Function", expression(paste("Current Priod Assets: ", k[t])), expression(paste("Value Function: ", V[t](k[t]))), vfcn)
		plotNorth("(b) Policy Function", expression(paste("Current Priod Assets: ", k[t])), expression(paste("Next Priod Assets: ", k[t+1])), pfcn)
		plotNorth("(c) Consumption Function", expression(paste("Current Priod Assets: ", k[t])), expression(paste("Consumption: ", c[t])), cfcn)

		with(new.env(), {
			T <- c(1, 1, 9, 9)
			lty <- c(3, 1, 1, 2)
			lwd <- c(3, 2, 2, 3)
			col <- c("black", "red", "black", "red")
			method <- c("approximation", "polynominal", "approximation", "polynominal")
			ylim <- c(min(pfcn, p_true), max(pfcn, p_true)*1.2)
			plot(kgrid, pfcn[, T[1]], main="(d) Precision of the approximation", xlab=expression(paste("Current Priod Assets: ", k[t])), ylab=expression(paste("Next Priod Assets: ", k[t+1])), ylim=ylim, type="l", lty=lty[1], lwd=lwd[1], col=col[1])
			for(i in 2:length(T)){
				if(1 == i %% 2){
					lines(kgrid, pfcn[, T[i]], lty=lty[i], lwd=lwd[i], col=col[i])
				} else{
					lines(kgrid, p_true[, T[i]], lty=lty[i], lwd=lwd[i], col=col[i])
				}
			}
			legend("topleft", lty=lty, lwd=lwd, col=col, legend=sapply(1:length(T), function(i){
				eval( parse(text=sprintf("expression(paste(\"%s: \", t==%d))", method[i], T[i])) )
			}), bty="o", y.intersp=1.1, seg.len=2.5, bg="white", inset=0.01)
		})

		mtext(side=3, line=1, outer=T, text="Fig.2 Robinson Cruesoe Model", cex=1.5)
	})
}

r_RCModel <- calcRCModel(param) # 計算する
plotRCModel(r_RCModel) # プロットする
dev.copy2eps(file="fig 2.eps", width=8, height=8) # 図を保存する
# dev.copy(png, "fig 2.png", width=800, height=800, type="cairo", bg="white"); dev.off()

