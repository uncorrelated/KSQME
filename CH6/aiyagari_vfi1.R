#
# Aiyagari (1994) の需要供給曲線の計算とプロットで用いるVFの計算の低精度版
# (1) aiyagari_vfi1(r) : grid search over the same grid as the state (grida)
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
	Nk   = 100;                     # grid size for state and control
	maxK = 20;                      # maximum value of capital grid  
	minK = -phi;                    # borrowing constraint
	gridk= seq(minK, maxK, length.out=Nk)      # state of assets 

	#  initialize some variables
	kfunG  = matrix(0, Nl, Nk);
	v      = matrix(0, Nl, Nk); 
	tv     = matrix(0, Nl, Nk);

	kfun = kfunG_old=kfunG;

	err     = 10;
	maxiter = 2000;

	j <- 0

	for(iter in 1:maxiter){

		#  tabulate the utility function such that for zero or negative
		#  consumption utility remains a large negative number so that
		#  such values will never be chosen as utility maximizing

		for(kc in 1:Nk){
			for(lc in 1:Nl){            
				vtemp = -1000000*matrix(1, Nk, 1);
				kccmax <- Nk # Matlabコードでは初期化されておらず、ループ中の前の処理結果の値が入る可能性がある
				for(kcc in 1:Nk){
					cons = s[lc]*wage + (1+r)*gridk[kc] - gridk[kcc];
					if(cons<=0){
						kccmax=kcc-1;
						break
					}
					util = (cons^(1-mu))/(1-mu);                                
					vpr=0;
					for(lcc in 1:Nl){
						vpr = vpr + prob[lc,lcc]*v[lcc,kcc]; 
					}                
					vtemp[kcc] = util + beta*vpr;                
				}

				t2 <- which.max(vtemp[1:kccmax])
				t1 <- vtemp[t2]

				tv[lc, kc] = t1;
				kfunG[lc,kc] = t2;
				kfun[lc,kc] = gridk[t2];            
			}
		}
	    
		v=tv;       

		err=max(abs(kfunG-kfunG_old));

		kfunG_old=kfunG;
		if(err <=0){
			break
		}
	}

	mea0 <- matrix(1, Nl, Nk) / (Nl*Nk);
	mea1 <- matrix(0, Nl, Nk);
	errTol = 0.00001;
	maxiter = 2000;

	for(iter in 1:maxiter){
		if(FALSE && j<3){
			print(c(mea0));
			j = j + 1
		}

    		for(kc in 1:Nk){
			for(lc in 1:Nl){
				kcc = kfunG[lc,kc];
				for(lcc in 1:Nl){
					mea1[lcc, kcc] = mea1[lcc, kcc] + prob[lc,lcc]*mea0[lc,kc];
					if(FALSE && j<100){
						print(sprintf("j: %d", j))
						print(sprintf("prob[%d,%d]: %f", lc-1, lcc-1, prob[lc,lcc]))
						print(sprintf("mea0[p(%d, %d)]: %f", lc-1, kc-1, mea0[lc, kc]))
						print(sprintf("mea1[p(%d, %d)]: %f", lcc-1, kcc-1, mea1[lcc, kc]))
						j = j + 1
					}

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

	list(a=sum(mea0*kfun), gridk=gridk, kfun=kfun)
}


