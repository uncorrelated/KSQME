#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include <float.h>

#define Nk 150 /* grid size for state */
#define Nk2 500 /* grid size for control */
#define maxK 20 /* maximum value of capital grid */
#define p(r,c) (r+Nl*c)
#define maxiter 2000
#define errTol 0.00001

SEXP VFI(SEXP r, SEXP alpha, SEXP beta, SEXP delta, SEXP mu, SEXP b, SEXP s, SEXP prob)
{
/* C++と異なりC言語は関数の冒頭に変数宣言をまとめる必要がある */
	SEXP  rv1, rv2, rv3, lst, lstnames, dim;
	double wage, phi, *gridk, *gridk2, minK, intK, intK2, *v, *tv, *kfun, *vtemp, cons, vpr, t1, err=0, diff, util, *mea0, *mea1, result, pr1, pr2;
	int *kfunG, *kfunG_old;
	int nor, noc, Nl, diff_int, err_int;
	 /* ループ用の変数群 */
	unsigned i, iter, kc, lc, lcc, kcc, kcc1, kcc2;
	int kccmax, t2, t3; /* プログラムの構造上、マイナスが代入されるのでint */

/* Rが提供する関数を用いて、第1から第6引数が実数か確認 */
	if(!isNumeric(r) || !isNumeric(alpha) || !isNumeric(beta) || !isNumeric(delta) || !isNumeric(mu) || !isNumeric(b)){
		error("Numbers are required for the arguments between 1st and 5th.");
	}

/* 第7引数がベクターか確認 */
	if(!isVector(s)){
		error("A vector is required for the 6th argument: 's'.");
	}

/* 第8引数が行列か確認 */
	if(!isMatrix(prob)){
		error("A matrix is required for the 7th argument: 'prob'.");
	}

/* 変数Nlはsの長さとする; lengthはRが提供する関数 */
	Nl = length(s);

/* 行列probの行数、列数を得る */
	dim = getAttrib(prob, R_DimSymbol);
	nor = INTEGER(dim)[0]; /* 行数 */
	if(Nl>nor){
		error("the number of row is smaller than the length of s.");
	}
	noc = INTEGER(dim)[1]; /* 列数 */
	if(Nl>noc){
		error("the number of column is smaller than the length of s.");
	}

/* write wage as a function of interest rate */
/* 
	REAL(Rの変数名)で浮動小数点型へのポインターが得られるので、*をつけてアクセスする
	C言語にはa^bのような便利な記法は無いので、Rmath.hをインクルードして、Rが提供するR_pow関数を使う
*/
	wage = (1-*REAL(alpha))*R_pow(R_pow((*REAL(alpha)/(*REAL(r)+*REAL(delta))), *REAL(alpha)), 1/(1-*REAL(alpha)));
	if(*REAL(r) < 0){
		phi = *REAL(b);
	} else {
		phi = wage*REAL(s)[1] / *REAL(r); /* ベクターの場合は[]でポインター記号*と同じ意味になる */
		if(phi > *REAL(b))
			phi = *REAL(b);
	}

/*
	capital grid (need define in each iteration since it depends on r/phi)
*/
	minK = -phi; /* borrowing constraint */
	gridk = Calloc(Nk, double);  /* state of assets; CallocはRが提供する関数 */
	intK = (maxK-minK)/(Nk-1);
	for(i=0; i<Nk; i++){
		gridk[i] = minK + i*intK;
	}

	gridk2 = Calloc(Nk2, double);  /* control of assets; vtemp用のグリッド */
	intK2 = (maxK-minK)/(Nk2-1);
	for(i=0; i<Nk2; i++){
		gridk2[i] = minK + i*intK2;
	}

/*
	initialize some variables
*/

	/* 一次元配列をファイル先頭で定義したマクロp(x,y)で二次元配列として扱っているので注意 */
	kfunG = Calloc(Nl*Nk, int); // Nk行Nk列の行列
	kfunG_old = Calloc(Nl*Nk, int);
	kfun = Calloc(Nl*Nk, double);
	tv = Calloc(Nl*Nk, double);
	v = Calloc(Nl*Nk, double);
	vtemp = Calloc(Nk2, double);
	mea0 = Calloc(Nl*Nk, double);
	mea1 = Calloc(Nl*Nk, double);

	/*
		Matlabのコードでは最初のループで（インデックスの最小値）1が得られる（ようだ）が、
		C言語は最初のループで（インデックスの最小値）0が得られるので、
		初期値0の羅列をkfunG_oldの初期値にしておくと、いきなりループが終了する
	*/
	for(i=0;i<Nk*Nl;i++){
		kfunG_old[i] = -1;
	}

	for(iter=1; maxiter>=iter; iter++){
/*
	tabulate the utility function such that for zero or negative
	consumption utility remains a large negative number so that
	such values will never be chosen as utility maximizing
*/
		for(kc=0; Nk>kc; kc++){ /* kcが0からはじまるのに注意 */ 
			for(lc=0; Nl>lc; lc++){ /* lcが0からはじまるのに注意 */
				for(i=0; i<Nk2; i++){
					vtemp[i] = -1000000;// -DBL_MAX; /* double型の最小値 */
				}

				/* Matlabのコードではkccmaxがここで初期化されないため、cons<=0のグリッドがない場合は、誤差が入る */
				for(kcc=0, kccmax=Nk2-1; Nk2>kcc; kcc++){
					cons = REAL(s)[lc]*wage + (1 + *REAL(r))*gridk[kc] - gridk2[kcc];
					if(cons<=0){
						kccmax = kcc - 1;
						break;
					}
					util = R_pow(cons, (1 - *REAL(mu)))/(1 - *REAL(mu));

					/* gridをgrid2に重み付き平均で対応させるために、2点kcc1とkcc2を選んで、ウェイトpr1, pr2を計算 */
					kcc1 = floor((gridk2[kcc]-minK)/intK); /* Matlabと異なりCのインデックスは0からはじまるので +1 は削除; phiは-minKに変更 */
					kcc2 = kcc1 + 1;
					if(Nk <= kcc2){
						/* kcc2がグリッドの最大値を越えるときの処理 */
						kcc1=Nk - 1;
						kcc2=Nk - 1;
						pr1=1;
						pr2=0;
					} else {
						pr2=(gridk2[kcc]-gridk[kcc1])/intK;
						pr1=1-pr2;
					}

					for(vpr=0, lcc=0; Nl>lcc; lcc++){
						vpr += REAL(prob)[lc + nor*lcc]*(pr1*v[p(lcc, kcc1)] + pr2*v[p(lcc, kcc2)]); /* probは他と縦横サイズが異なる */
					}

					vtemp[kcc] = util + (*REAL(beta))*vpr;
				}

				for(t3=1,t2=0,t1=vtemp[0]; kccmax>=t3; t3++){
					if(vtemp[t3]>t1){
						t1 = vtemp[t3];
						t2 = t3;
					}
				}

				tv[p(lc, kc)] = t1;
				kfunG[p(lc, kc)] = t2;
				kfun[p(lc, kc)] = gridk2[t2];
			}
		}

		/* memcpyを使うべきかも */ 
		for(i=0;i<Nk*Nl;i++){
			v[i] = tv[i];
		}

		/* グリッド・ポジションを指すインデックスの差だからdiff_int/err_intは整数 */
		for(i=0,err_int=0;i<Nk*Nl;i++){
			diff_int = abs(kfunG[i] - kfunG_old[i]);
			if(diff_int > err_int)
				err_int = diff_int;
			kfunG_old[i] = kfunG[i]; /* MatlabのkfunG_old=kfunGの部分も同じループ内で処理 */
		}

		if(err_int <= 0){ /* グリッド・ポジションを指すインデックスの差だから0で終了 */
			break;
		}
	}

	for(i=0;i<Nl*Nk;i++){
		mea0[i] = 1.0F/(Nl*Nk); /* 1.0Fと表記しておかないと浮動小数点にならないかも */
	}

	for(iter=1; maxiter>=iter; iter++){

		for(kc=0; kc<Nk; kc++){
			for(lc=0; lc<Nl; lc++){
				kcc = kfunG[p(lc,kc)];

				/* gridをgrid2に重み付き平均で対応させるために、2点kcc1とkcc2を選んで、ウェイトpr1, pr2を計算 */
				kcc1 = floor((gridk2[kcc]-minK)/intK); /* Matlabと異なりCのインデックスは0からはじまるので +1 は削除; phiは-minKに変更 */
				kcc2 = kcc1 + 1;
				if(Nk <= kcc2){
					/* kcc2がグリッドの最大値を越えるときの処理 */
					kcc1=Nk - 1;
					kcc2=Nk - 1;
					pr1=1;
					pr2=0;
				} else {
					pr2=(gridk2[kcc]-gridk[kcc1])/intK;
					pr1=1-pr2;
				}

				for(lcc=0;lcc<Nl;lcc++){
					mea1[p(lcc, kcc1)] = mea1[p(lcc, kcc1)] + REAL(prob)[lc + nor*lcc]*pr1*mea0[p(lc, kc)];
					mea1[p(lcc, kcc2)] = mea1[p(lcc, kcc2)] + REAL(prob)[lc + nor*lcc]*pr2*mea0[p(lc, kc)];
				}            
			}        
		}



		for(i=0,err=0,diff=0;i<Nk*Nl;i++){
			diff = fabs(mea1[i] - mea0[i]);
			if(diff > err)
				err = diff;
			mea0[i] = mea1[i];
			mea1[i] = 0;
		}

		if(err <= errTol){
			break;
		}
	}

	rv1 = PROTECT(allocVector(REALSXP, 1)); /* 計算結果を保持するベクターを作成（注意： 値が1つでもベクター扱いになるし、C言語からみると配列） */
	rv2 = PROTECT(allocVector(REALSXP, Nk)); /* gridkを戻すためのベクターを作成 */
	rv3 = PROTECT(allocMatrix(REALSXP, Nl, Nk)); /* kfunを戻すための行列を作成 */

	for(i=0; i<Nk; i++){
		REAL(rv2)[i] = gridk[i];
	}

	for(i=0,result=0; i<Nl*Nk; i++){
		result += mea0[i]*kfun[i]; /* 貯蓄E[a]を計算しているハズ */
		REAL(rv3)[i] = kfun[i]; /* 行列と言ってもただの1次元ベクトル */
	}

	REAL(rv1)[0] = result; /* C言語は配列の添字は0からなのに注意; 確保したサイズを越えたところに代入してもエラーメッセージは出ず処理続行 */

	lst = PROTECT(allocVector(VECSXP, 3)); /* 戻り値になるリストを作成 */

	/* リスト内に計算結果ベクトルなどを入れる */
	SET_VECTOR_ELT(lst, 0, rv1); 
	SET_VECTOR_ELT(lst, 1, rv2);
	SET_VECTOR_ELT(lst, 2, rv3);

	/* リストの要素の名前をリストの設定 */ 
	PROTECT(lstnames = allocVector(STRSXP, 3));
	SET_STRING_ELT(lstnames, 0,  mkChar("a"));
	SET_STRING_ELT(lstnames, 1,  mkChar("gridk"));
	SET_STRING_ELT(lstnames, 2,  mkChar("kfun"));
	setAttrib(lst, R_NamesSymbol, lstnames); 
	UNPROTECT(5); /* PROTECTした回数だけUNPROTECT */

/* 確保した領域を解放; FreeはRが提供する関数 */ 
	Free(mea1);
	Free(mea0);
	Free(vtemp);
	Free(v);
	Free(tv);
	Free(kfun);
	Free(kfunG_old);
	Free(kfunG);
	Free(gridk2);
	Free(gridk);

	return(lst);
}

