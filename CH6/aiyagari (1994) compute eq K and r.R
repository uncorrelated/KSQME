#
# 経済セミナー連載「定量的マクロ経済学と数値計算」
# 第6回「異質な個人を組み込んだマクロモデル」のソースコードのRへの野良移植: 均衡点の資本と金利の計算と図3のプロット
#
#==========================================================================
# Model of Aiyagari (1994)
# Written for Keizai Seminar #6
# By Sagiri Kitao
# 
# ported to R by uncorrelated
#==========================================================================
#
# 共通コードを読み込む
source("common.R")

# 価値関数のタイプ（高精度/低精度）を定める
# (1) VALUE FUNCTION (USE THE SAME GRID FOR STATE AND CONTROL)
vfi_type <- "aiyagari_vfi2"

# (2) VALUE FUNCTION (USE FINER GRID FOR CONTROL)
# vfi_type <- "aiyagari_vfi2" # 時間がかかるが精度がアップ

#
# （1と2の処理は別ファイルに移した）
# =3 compute eq K and r : method 1 : search over r-grid from the bottom 
# =4 compute eq K and r : method 2 : update new guess of r in two ways
#
# =5 compute eq K and r : method X : bisection method
#
indE = 5;

# 演算ルーチンを指定する（3を指定しても連続実行はしない）
# 1:R, 2:C
r_c_mp_switch <- 2

# 価値関数を読み込む
if(chkflag(r_c_mp_switch,1)){
	print("R Code:")
	source(paste(vfi_type, "R", sep=".")) 
}else if(chkflag(r_c_mp_switch,2)){
	print("C Code:")
	vfi <- loadCModule()
}else{
	stop("r_c_mp_switch is invald.")
}


#==========================================================================
# COMPUTE K and r in EQ
#==========================================================================

print(system.time({
	if(3==indE){

		rate0=0.02;    # initial guess (START WITH A VALUE LESS THAN EQ VALUE)
		adj  =0.001;

		tm_baseline <- calcTM(sig)    

		ind=0;
		while(0==ind){
			K0=tm_baseline$labor*(alpha/(rate0+delta))^(1/(1-alpha));

			r_vfi <- with(tm_baseline, {
				vfi(rate0, alpha, beta, delta, mu, b, s, prob)
			})

			K1 <- r_vfi$a
			if(K0<K1){
			    ind=1;
			}
			rate0=rate0+adj;
			#        [ind rate0 K0 K1 K0-K1] ← 謎
		}

		# INTEREST RATE AND CAPITAL IN EQUILIBRIUM (SOLUTIONS)
		EQ_rate=rate0
		EQ_K=K0
		# 明示されていないがvif関数内で計算したgridkとkfunもプロットで使う   

	}else if(4==indE){

		rate0= 0.025; # INITIAL GUESS

		tm_baseline <- calcTM(sig)    

		err=1;
		errTol=0.001;
		maxiter=200;
		iter=1;

		adj=0.001;

		for(i in 1:maxiter){

			K0=tm_baseline$labor*(alpha/(rate0+delta))^(1/(1-alpha));

			r_vfi <- with(tm_baseline, {
				vfi(rate0, alpha, beta, delta, mu, b, s, prob)
			})

			K1 <- r_vfi$a

			err = abs(K0-K1);
				
			# (1) UPDATE GUESS AS (r0+r(K1))/2
			#rtemp=alpha*((K1/labor)^(alpha-1))-delta;
			#rate0=(rtemp+rate0)/2;

			# (2) UPDATE GUESS AS r0+adj*(K0-K1)
			rate0=rate0+adj*(K0-K1);

			if(errTol >= err){
				break;
			}      
		}    

		if(maxiter==i){
			stop("iteration limit exceeded.")
		}
		
		# INTEREST RATE AND CAPITAL IN EQUILIBRIUM (SOLUTIONS)
		EQ_rate=rate0
		EQ_K=K0
		# 明示されていないがvif関数内で計算したgridkとkfunもプロットで使う   
	}else if(5==indE){
		#
		# Matlabコードには無かったが、3==indEが総当たりアルゴリズムで、4==indEがなぜか収束しないので二分探索法での探索を用意する
		#
		tm_baseline <- calcTM(sig)    
		calcK <- function(rate){
			tm_baseline$labor*(alpha/(rate+delta))^(1/(1-alpha))
		}

		rate_l = 0.0;
		rate_h = 1.0;

		calcVFI <- function(rate0){
			r_vfi <- with(tm_baseline, {
				vfi(rate0, alpha, beta, delta, mu, b, s, prob)
			})
			r_vfi$K0 =  calcK(rate0);
			r_vfi$K1 = r_vfi$a
			r_vfi$flag = r_vfi$K0 > r_vfi$K1
			r_vfi$err = abs(r_vfi$K0 - r_vfi$K1)
			r_vfi$rate = rate0
			r_vfi
		}

		r_vfi_l <- calcVFI(rate_l)
		r_vfi_h <- calcVFI(rate_h)

		if(r_vfi_l$flag == r_vfi_h$flag){
			stop("ill-condtion for the bisection method.")
		}

		r_vfi_m <- calcVFI((rate_l + rate_h)/2)

		errTol = 0.001;
		steptol = 1e-4
		step <- 1e+4

		for(i in 1:100){
			if(r_vfi_m$err < errTol || step < steptol){
				break;
			}
			err0 <- r_vfi_m$err
			if(r_vfi_m$flag == r_vfi_l$flag){
				r_vfi_l <- r_vfi_m
				rate_l <- r_vfi_m$rate
			} else {
				r_vfi_h <- r_vfi_m
				rate_h <- r_vfi_m$rate
			}
			r_vfi_m <- calcVFI((rate_l + rate_h)/2)
			step <- abs(r_vfi_m$err - err0)
#			print(r_vfi_m$rate)
		}

		r_vfi <- r_vfi_m # 原理的に解付近に到達するので

		# INTEREST RATE AND CAPITAL IN EQUILIBRIUM (SOLUTIONS)
		EQ_rate=r_vfi$rate
		EQ_K=r_vfi$K0
		# 明示されていないがvif関数内で計算したgridkとkfunもプロットで使う   

	}else{
		stop("indE is invalid.")
	}
}))

print(sprintf("EQ_rate=%f EQ_K=%f", EQ_rate, EQ_K)) # 計算した均衡点を表示しておく

# r_vfiをバイナリ保存（グラフを描き直すときにreadRDS("eq_vfi.obj")で読み込めば、上の時間のかかる計算を省略できる） 
saveRDS(r_vfi, file="eq_vfi.obj")

with(r_vfi, {
	par(bg="white")
	col <- c("black", "red", "blue", NA, NA)
	lty <- c(4, 2, 1, NA, NA)
	lwd <- c(2, 2, 2, NA, NA)
	expr1 <- parse(text = sprintf("hat(r)==%.3f", EQ_rate)) # 反則気味だが、均衡値を凡例に書いておく
	expr2 <- parse(text = sprintf("hat(K)==%.3f", EQ_K)) 
	plot(gridk, kfun[1,], main="Policy Function", xlab=expression(a[t]), ylab=expression(a[t+1]==g(a[t],l)), col=col[1], lty=lty[1], lwd=lwd[1], type="l")
	lines(gridk, kfun[4,], col=col[2], lty=lty[2], lwd=lwd[2])
	lines(gridk, kfun[7,], col=col[3], lty=lty[3], lwd=lwd[3])
	legend("bottomright", lty=lty, lwd=lwd, col=col, legend=c(expression(l[low]), expression(l[middle]), expression(l[high]), expr1, expr2), bty="n", y.intersp=1.5, seg.len=3, inset=0.05)
})

# 図を保存する
dev.copy2eps(file="fig 3.eps", width=6, height=6)
# dev.copy(png, "fig 3.png", width=480, height=480, type="cairo", bg="white"); dev.off()

