#include <Rcpp.h>
#include "pchip.h"
#include <cmath>

using namespace Rcpp;

void RprintMatrix(NumericMatrix m, int s, int n){

	for(int i=s;i<s+n && i<m.nrow(); i++){
		for(int j=0;j<m.ncol();j++){
			Rprintf("%f\t", m(i, j));
		}
		Rprintf("\n");
	}
}


double update(NumericMatrix m1, NumericMatrix m2){
	double max = -1;
	for(int j=0;j<m1.ncol();j++){
		for(int i=0;i<m1.nrow();i++){
			double d = fabs(m1(i,j) - m2(i,j));
			if(max < d){
				max = d;
			}
			m2(i,j) = m1(i,j); // m1とm2を入れ替えたら、不要になるはずのコード
		}
	}
	return max;
}

// [[Rcpp::export]]
List HH_opt_EGM(double beta, List param, NumericVector grid_a, NumericVector grid_z, NumericMatrix prob_z, double wToday, double RToday, double tauToday, double dToday, NumericMatrix pf_c_init, NumericMatrix pf_n_init, NumericMatrix pf_sav_init) {

	NumericMatrix pf_c(clone(pf_c_init)), pf_c_new(clone(pf_c_init)); // コンストラクタにそのまま放り込むと参照先が一緒になり、片方を更新したらもう片方も更新されてしまうので、cloneをしておく
	NumericMatrix pf_n(clone(pf_n_init)), pf_n_new(clone(pf_n_init));
	NumericMatrix pf_sav(clone(pf_sav_init)), pf_sav_new(clone(pf_sav_init));

	// リストで渡される値は展開して型変換しておく
	int max_egm_iter = as<int>(param["max_egm_iter"]);
	int n_a = as<int>(param["n_a"]);
	int n_z = as<int>(param["n_z"]);
	double gamma = as<double>(param["gamma"]);
	double psi = as<double>(param["psi"]);
	NumericVector tau_bar = as<NumericVector>(param["tau_bar"]);
	double bl = as<double>(param["bl"]);
	double egm_err_tol = as<double>(param["egm_err_tol"]);
	int max_egm_const_iter = as<int>(param["max_egm_const_iter"]);
	double egm_const_err_tol = as<double>(param["egm_const_err_tol"]);

	// WHILE LOOP: iterate until the policy functions converge.
	int EGM_iter;
	for(EGM_iter = 0; EGM_iter < max_egm_iter; EGM_iter++){

		NumericMatrix pf_c_int(n_a + 1, n_z);
		NumericMatrix endog_grid_a(n_a + 1, n_z);

		// Assuming that the borrowing constraint is not binding, solve for the policy functions.
		for(int i_z = 0; i_z < n_z; i_z++){
			double zToday      = grid_z[i_z];
			NumericVector cond_prob_z = prob_z(i_z, _); // i_z行をベクトルとして取り出す

			for(int i_a = 0; i_a < n_a; i_a++){
				// EGM: construct a grid for a', rather than a grid for a.
				//     solve for the value of a that would have led to the choice a'.

				double aTomorrow = grid_a[i_a];

				// Compute the expected marginal utility of consumption using a', and derive consumption today from the Euler equation.
				// Note that here we skip a root-finding algorithm, thereby making the algorithm much more efficient than the usual policy function iteration!
				double EMUc = 0.0;
				for(int j_z=0; j_z < n_z; j_z++)
					EMUc += beta*RToday*(cond_prob_z[j_z]*std::pow(pf_c(i_a, j_z),-gamma));

				double cToday = std::pow(EMUc, (-1.0/gamma));

				// Derive labor supply today from the labor supply equation.
				double nToday = std::pow(wToday*zToday*std::pow(cToday,-gamma), 1.0/psi);

				// Derive today's asset from the budget constraint. This current asset level is called the endogenous gridpoints.
				double aToday = cToday + aTomorrow/RToday - wToday*zToday*nToday + tauToday*tau_bar[i_z] - dToday;

				// Store the values
				pf_c_int(i_a, i_z)	= cToday;
				endog_grid_a(i_a, i_z)	= aToday;
			}
		}

		// The following is used for evaluating the new policy functions at the exogenous gridpoints using the pchip interpolation if the exogenous gridpoints are beyond the maximum endogenous gridpoints. Here, pf_c_int(N_A+1,:) is linearly extrapolated.
		for(int i_z = 0; i_z < n_z; i_z++){
			endog_grid_a(n_a + 1 - 1, i_z) = 1e+5; // Note(20190825): Nakata-san set this as 1e8, rather than 1d5.

			pf_c_int(n_a + 1 - 1, i_z) = pf_c_int(n_a - 1, i_z) 
				+ ((pf_c_int(n_a - 1, i_z) - pf_c_int(n_a - 1 - 1, i_z)) / (endog_grid_a(n_a - 1, i_z) - endog_grid_a(n_a - 1 - 1, i_z))) 
				* (endog_grid_a(n_a + 1 - 1, i_z) - endog_grid_a(n_a - 1, i_z)); // C++ではn_aはインデックスの最大値ではなく、外側になるので-1している

			// The coefficient of the second term is just a slope.
		}

		// Evaluate the new policy functions at the exogenous gridpoints (or the original gridpoints).
		for(int i_z = 0; i_z < n_z; i_z++){
			double		zToday      = grid_z[i_z];
//			NumericVector	cond_prob_z = prob_z(i_z, _); // i_z行をベクトルとして取り出す

			for(int i_a = 0; i_a < n_a; i_a++){

				double aToday = grid_a[i_a], cToday, nToday, aTomorrow;

				// In this case, endog_grid_a(1,:) is the value of bond holdings that induces the borrowing constraint to bind next period.
				// This is because the far left gridpoint in this program is set to the borrowing limit.
				if (aToday > endog_grid_a(0, i_z)){ // Matlabとインデックス開始位置が異なる（Matlab:1, C++:0)

					// The borrowing constraint does not bind.
					// shape-preserving spline!!!

					// Matlab互換Hermite補間
					cToday = pchip(endog_grid_a(_, i_z), pf_c_int(_, i_z), aToday); // 別ファイルに用意した関数
					nToday = std::pow(wToday*zToday*std::pow(cToday, -gamma), 1.0/psi);
					aTomorrow = (aToday + wToday*zToday*nToday - tauToday*tau_bar[i_z] + dToday - cToday)*RToday;

				} else {

					// The borrowing costraint binds. Use the subroutine 'EGMConstrained' to compute cToday and nToday when the borrowing constraint is binding.
					aTomorrow = bl;

					// Call 'EGMConstrained' to obtain the values of cToday and nToday when the borrowing constraint is binding.
					// [cToday nToday] = EGMconstrained(p,aToday,zToday,RToday,wToday,tauToday,dToday,i_z);

					// ここでしか使っていない上に、戻り値が二つなので、サブルーチンをインライン展開する
					// Initial guess of the value of labor supply and the associated consumption value.
					nToday = 0.6;
					cToday = aToday + wToday*zToday*nToday - tauToday*tau_bar[i_z] + dToday - bl/RToday;

					double labor_eq_diff = std::pow(cToday,-gamma)*wToday*zToday - std::pow(nToday, psi); // labor_eq_diff denotes the difference of the labor supply equation.

					int EGM_const_iter;

					// WHILE LOOP: iterate until we find a pair of labor supply and consumption that satisfies the labor supply equation.
					//            Here we use the Newton-Raphson method.
					for(EGM_const_iter=0; EGM_const_iter < max_egm_const_iter; EGM_const_iter++){

						double labor_eq_adj  = -gamma*std::pow(cToday, -gamma - 1.0)*std::pow(wToday*zToday, 2.0) - psi*std::pow(nToday, psi-1.0);

						nToday        = nToday - labor_eq_diff/labor_eq_adj;

						cToday        = aToday + zToday*wToday*nToday - tauToday*tau_bar[i_z] + dToday - bl/RToday;

						labor_eq_diff = std::pow(cToday, -gamma)*wToday*zToday - std::pow(nToday, psi);

						double EGM_const_err = fabs(labor_eq_diff);

						/*
						%     ! write(*,*) "------------------------------------------------------"
						%     ! write(*,*) "AT ITERATION   = ", EGM_const_iter
						%     ! write(*,*) "DIFFERENCE = ", EGM_const_err
						%     ! write(*,*) "------------------------------------------------------"
						%     disp([EGM_const_iter EGM_const_err labor_eq_adj nToday cToday]);
						*/
						if (EGM_const_err < egm_const_err_tol)
							break;
					}

					if(EGM_const_iter >= max_egm_const_iter){
						stop("iteration limit exceeded: EGM_const_iter >= max_egm_const_iter");
					}
				}

				/*
				% Obtain the new policy functions for consumption, labor supply, and savings as well as the associated value function.
				% Prepation for the value function
				%             Ev_int = 0d0
				%             do j_z = 1, p%N_Z
				%                 call interp1_pchip(vTomorrow, p%N_A, grid_a, aTomorrow, values(:,j_z))
				%                 Ev_int = Ev_int + cond_prob_z(j_z)*vTomorrow
				%             end do
				%             vToday = cToday**(1d0 - p%GAMMA)/(1d0 - p%GAMMA) - nToday**(1d0 + p%PSI)/(1d0 + p%PSI) + beta*Ev_int
				%             disp([i_a i_z grid_a(i_a) cToday nToday aTomorrow]);
				%             pause
				*/

				pf_c_new(i_a, i_z)   = cToday;
				pf_n_new(i_a, i_z)   = nToday;
				pf_sav_new(i_a, i_z) = aTomorrow;
				// values_new(i_a,i_z) = vToday
			}
		}

		// Evaluate convergence.
		NumericVector EGM_err(3);
		EGM_err[0] = update(pf_c_new, pf_c); 
		EGM_err[1] = update(pf_n_new, pf_n);
		EGM_err[2] = update(pf_sav_new, pf_sav);
		double EGM_maximum_error = max(EGM_err);

//		RprintMatrix(pf_c, 100, 1);

		if(EGM_maximum_error < egm_err_tol)
			break;

		// ! Update the policy functions for consumption, labor supply, and savings as well as the associated value function

/* 値をコピーすると遅いので、参照を入れ替えようかと思ったが、見通しが悪い気がして保留している
		NumericMatrix tmp;
		tmp = pf_c; pf_c = pf_c_new; pf_c_new = tmp;
		tmp = pf_n; pf_n = pf_n_new; pf_n_new = tmp;
		tmp = pf_sav; pf_sav = pf_sav_new; pf_sav_new = tmp;
*/

		//     values = values_new
	}

	if(EGM_iter >= max_egm_iter){
		stop("iteration limit exceeded: EGM_iter >= max_egm_iter");
	}

	return List::create(Named("pf_c")=pf_c_new, Named("pf_n")=pf_n_new, Named("pf_sav")=pf_sav_new);
}

