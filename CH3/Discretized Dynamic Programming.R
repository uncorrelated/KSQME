#
# 経済セミナー連載「定量的マクロ経済学と数値計算」
# 第3回「動的計画法」のソースコードのRへの野良移植: 第5節無限期間モデルの解法, 図3
#

# ここは本家にもR版のコードがあるが、Matlab坂と挙動が異なり、本文図4などの作成できないので修正版をつくった

source("common.R")

param <- list(
	### カリブレーション ###
	beta  = 0.96, # 割引因子
	gamma = 1.0,  # 相対的危険回避度(異時点間の代替の弾力性の逆数)
	alpha = 0.40, # 資本分配率
	delta = 1.00, # 固定資本減耗

	### 離散化用のパラメータ ###
	nk   = 1001, # グリッドの数 ← Matlabのコードでは10001だが、設定するとメモリー的に計算不能かも
	kmax = 0.5,  # 資本グリッドの最大値
	kmin = 0.05  # 資本グリッドの最小値 (0にすると生産が出来なくなる)
)

# 固定資本減耗を変えるときは以下2行をアンコメント
# param$delta <- 0.08
# param$kmax <- 10.0

### 収束の基準 ###
maxit = 1000;    # 繰り返し計算の最大値
tol  = 1.0e-005; # 許容誤差(STEP 2)
dif1 = 1.0;      # 価値関数の繰り返し誤差
dif2 = 1.0;      # 政策関数の繰り返し誤差

# パラメーターを変えて比較するかも知れないので、数値計算部分を関数化する
calcDDP <- function(param){
	with(param, {
		cat("\n-+- Solve a neoclassical growth model -+-\n");
		cat("\n-+- PARAMETER VALUES -+-\n");
		cat(sprintf("beta=%5.2f, gamma=%5.2f, alpha=%5.2f, delta=%5.2f\n", beta, gamma, alpha, delta))
		cat(sprintf("kmin=%5.2f, kmax=%5.2f, #grid=%d\n", kmin, kmax, nk))

		## STEP 1(a): グリッド生成
		kgrid = seq(kmin, kmax, length.out=nk);

		# STEP 1(b):価値関数・政策関数の初期値を設定
		vfcn  <- numeric(nk)
		pfcn  <- numeric(nk)
		Tvfcn <- numeric(nk)
		Tpfcn <- numeric(nk)
		vkp   <- matrix(data = 0.0, nrow = nk, ncol = nk)

		# STEP 3: 効用関数の組み合わせ

		# 効用関数の初期値 (消費が0以下になる組み合わせにはペナルティ)
		wealth <- t(replicate(nk, kgrid^alpha + (1.0-delta)*kgrid)) # 第2引数の縦ベクトルをnk回繰り替えした行列
		cons <- wealth - replicate(nk, kgrid) # 第2項のkgridは自動的に行列に拡大され計算される
		util <- matrix(-1e+6, nrow = nk, ncol = nk) # 初期値は大きな負の値
		util[0<cons] <- CRRA(cons[0<cons], gamma) # 消費が正値になる(k,k')の組み合わせについて効用を計算

		# STEP 4: 価値関数を繰り返し計算
		for (it in 1:maxit) { # itがループのカウンターになる

			# ベルマン方程式: V(k;k')
			vkp = util + beta*vfcn

			# 最適化: 各kについてV(k;k')を最大にするk'を探す
			for(i in 1:nk){
				j <- which.max(vkp[, i])
				Tvfcn[i] <- vkp[j, i]
				Tpfcn[i] <- kgrid[j]
			}

			# 繰り返し計算誤差を確認
			dif1 <- max(abs((Tvfcn-vfcn)/vfcn))
			dif2 <- max(abs((Tpfcn-pfcn)/pfcn))

			# 価値関数・効用関数をアップデート
			vfcn <- Tvfcn
			pfcn <- Tpfcn

			cat(" iteration index:", it, ", iter dif of value:", dif1, ", iter dif of policy:", dif2, "\n")

			# 繰り返し計算誤差がtol未満になったら、収束したとして終了
			if (dif1 < tol){
				cat(sprintf("The iteration finished successfully: dif1=%f < %f", dif1, tol))
				break
			}
		}

		# 繰り返し上限に達していたらエラー
		if(maxit <= it){
			stop("iteration limit exceeded.")
		}

		# 最終的な政策関数が得られてから消費関数を計算
		cfcn = kgrid^alpha + (1-delta)*kgrid - pfcn

		## 代数的解の計算
		AA = (1.0-beta)^(-1) * (log(1.0-alpha*beta) + ((alpha*beta)/(1.0-alpha*beta))*log(alpha*beta));
		BB = alpha/(1.0-alpha*beta);
		v_true = AA + BB*log(kgrid);
		p_true = beta*alpha*(kgrid^alpha);

		# オイラー方程式から誤差を測定
		pfcn0 <- pfcn[seq(1, nk, length.out=21)] # pfcnの一部の要素を抜き出す
		kgrid2 <- seq(kmin, kmax, length.out=length(pfcn0))

		cons = kgrid2^alpha + (1-delta)*kgrid2 - pfcn0
		LHS  = dCRRA(cons, gamma);
		kp   = pfcn0
		interp1 <- splinefun(kgrid2, pfcn0)
		kpp  = interp1(kp);
		cons = kp^alpha + (1-delta)*kp - kpp;
		rent = alpha*kp^(alpha-1.0) - delta;
		RHS  = beta*(1+rent)*dCRRA(cons, gamma);
		err  = RHS/LHS-1.0;

		# 計算結果はリストにまとめておく
		list(kgrid=kgrid, v_true=v_true, p_true=p_true, vfcn=vfcn, pfcn=pfcn, cfcn=cfcn, dif=dif1, err=err, kgrid2=kgrid2)
	})
}

# 実行時間の計測と表示
print(system.time({
	r_ddp <- calcDDP(param)
}))

# 他で比較するために、誤差を保存
write.table(data.frame(
	kgrid = r_ddp$kgrid2, 
	error = r_ddp$err
), file=file_err_ddp, col.names=TRUE, row.names=FALSE, sep=",")

# 図を描く
with(r_ddp, {

	dev.new(width=8, height=4)
	par(oma=c(0, 0, 3, 0), mfrow=c(1, 2), mar=c(4.5, 4.5, 2, 1), bg="white")

	lty <- c(1, 2)
	lwd <- c(2, 2)
	col <- c("black", "red")
	legends <- c("approximation", "polynominal")
	plot(kgrid, vfcn,
		main="(a) Value Function",
		xlab=expression(paste("Current Priod Assets: ", k[t])),
		ylab=expression(paste("Value Function: ", V[t](k[t]))),		
		type="l", lty=lty[1], lwd=lwd[1], col=col[1])
	lines(kgrid, v_true, lty=lty[2], lwd=lwd[2], col=col[2])
	legend("bottomright", lty=lty, lwd=lwd, col=col, legend=legends, bty="n", y.intersp=1.1, seg.len=2.5, inset=0.02)

	lty <- c(lty, 3)
	lwd <- c(lwd, 1)
	col <- c(col, "gray")
	legends <- c(legends, "45-degree line")
	y <- c(pfcn, p_true)
	ylim <- c(0, max(y))
	plot(kgrid, pfcn,
		main="(b) Policy Function",
		xlab=expression(paste("Current Priod Assets: ", k[t])),
		ylab=expression(paste("Next Priod Assets: ", k[t+1])),		
		type="l", lty=lty[1], lwd=lwd[1], col=col[1], ylim=ylim)
	lines(kgrid, p_true, lty=lty[2], lwd=lwd[2], col=col[2])
	lines(kgrid, kgrid, lty=lty[3], lwd=lwd[3], col=col[3])
	legend("bottomright", lty=lty, lwd=lwd, col=col, legend=legends, bty="n", y.intersp=1.1, seg.len=2.5, inset=0.02)

	mtext(side=3, line=1, outer=T, text="Fig.3 Discretized Dynamic Programming", cex=1.5)

	dev.copy2eps(file="fig 3.eps", width=8, height=4) # 図を保存する
	# dev.copy(png, "fig 3.png", width=800, height=400, type="cairo", bg="white"); dev.off()
})

