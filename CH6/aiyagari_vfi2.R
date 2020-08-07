#
# Aiyagari (1994) の需要供給曲線の計算とプロットで用いるVFの計算の高精度版
# (2) aiyagari_vfi2(r) : grid search over a finer grid for the control (grida2)
#
# 資本グリッドが細かくなり、途中の計算で使う一時グリッドのサイズをさらに大きくしている
#
vfi <- function(r, alpha, beta, delta, mu, b, s, prob){

	# write wage as a function of interest rate 
	wage = (1-alpha)*((alpha/(r+delta))^alpha)^(1/(1-alpha));

	# borrowing limit
	phi <- wage*s[1]/r # s := exp(GRID)
	phi[r<=0 || phi>b] <- b

	# -phi is borrowing limit, b is adhoc
	# the second term is natural limit

	# capital grid (need define in each iteration since it depends on r/phi)
	Nk   = 150;                     # grid size for state
	maxK = 20;                      # maximum value of capital grid  
	minK = -phi;                    # borrowing constraint
	intK = (maxK - minK)/(Nk - 1)
	gridk= seq(minK, maxK, length.out=Nk)      # state of assets 

	Nk2   = 500;                    # grid size for CONTROL
	gridk2 = seq(minK, maxK, length.out=Nk2) # vtemp用のグリッド

	#  initialize some variables
	kfunG  = matrix(0, Nl, Nk);
	v      = matrix(0, Nl, Nk); 
	tv     = matrix(0, Nl, Nk);

	kfun = kfunG_old = kfunG;

	err     = 10;
	maxiter = 2000;

	j <- 0

	for(iter in 1:maxiter){

		#  tabulate the utility function such that for zero or negative
		#  consumption utility remains a large negative number so that
		#  such values will never be chosen as utility maximizing

		for(kc in 1:Nk){
			for(lc in 1:Nl){            

				vtemp = -1000000*matrix(1, Nk2, 1); # グリッドサイズを大きめにとる
				kccmax <- Nk2 # Matlabコードでは初期化されておらず、ループ中の前の処理結果の値が入る可能性がある

				for(kcc in 1:Nk2){
					cons = s[lc]*wage + (1+r)*gridk[kc] - gridk2[kcc];
					if(cons<=0){
						kccmax=kcc-1;
						break
					}
					util = (cons^(1-mu))/(1-mu);                                

					# gridをgrid2に重み付き平均で対応させるために、2点kcc1とkcc2を選んで、ウェイトpr1, pr2を計算
					if(kcc==Nk2){
						# kccが一時グリッドの最大値をさしているときは、2点をとれない
						kcc1=Nk;
						kcc2=Nk;
						pr1=1;
						pr2=0;
					}else{
						temp=(phi+gridk2[kcc])/intK;
						kcc1=floor(temp)+1;
						kcc2=kcc1+1;
						pr2=(gridk2[kcc]-gridk[kcc1])/intK;
						pr1=1-pr2;
					}

					vpr=0;
					for(lcc in 1:Nl){
						vpr = vpr + prob[lc,lcc]*(pr1*v[lcc,kcc1] + pr2*v[lcc,kcc2]); 

					}                
					vtemp[kcc] = util + beta*vpr;
				}

				t2 <- which.max(vtemp[1:kccmax])
				t1 <- vtemp[t2]
				tv[lc, kc] = t1;
				kfunG[lc,kc] = t2;
				kfun[lc,kc] = gridk2[t2];            
			}
		}
	    
		v=tv;

		err=max(abs(kfunG-kfunG_old));
#		print(sprintf("err: %d", err))

		kfunG_old=kfunG;
		if(err <=0){
			break
		}
	}

#	print(sprintf("iter: %d", iter))


	mea0 <- matrix(1, Nl, Nk) / (Nl*Nk);
	mea1 <- matrix(0, Nl, Nk);
	errTol = 0.00001;
	maxiter = 2000;

	for(iter in 1:maxiter){
    		for(kc in 1:Nk){
			for(lc in 1:Nl){
				kcc = kfunG[lc,kc];
				
				# gridをgrid2に重み付き平均で対応させるために、2点kcc1とkcc2を選んで、ウェイトpr1, pr2を計算
				if(kcc==Nk2){
					# kccが一時グリッドの最大値をさしているときは、2点をとれない
					kcc1=Nk;
					kcc2=Nk;
					pr1=1;
					pr2=0;
				}else{
					temp=(phi+gridk2[kcc])/intK;
					kcc1=floor(temp)+1;
					kcc2=kcc1+1;
					pr2=(gridk2[kcc]-gridk[kcc1])/intK;
					pr1=1-pr2;
				}
		
				for(lcc in 1:Nl){
					mea1[lcc, kcc1] = mea1[lcc, kcc1] + prob[lc,lcc]*pr1*mea0[lc,kc];
					mea1[lcc, kcc2] = mea1[lcc, kcc2] + prob[lc,lcc]*pr2*mea0[lc,kc];
				}
			}        
		}

		err = max(abs(mea1-mea0))
		mea0 = mea1;
		mea1 = matrix(0, Nl, Nk);

		if(err <= errTol){
			break
		}
	}

#	print(sprintf("iter: %d\n", iter))
	
	list(a=sum(mea0*kfun), gridk=gridk, kfun=kfun)
}

