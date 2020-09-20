#include <Rcpp.h>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <stdio.h>

#include "pchip.h"

using namespace Rcpp;
using namespace Eigen;
using namespace Spectra;

int	which(const NumericVector array, int size, double a);

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
List HH_dist(double beta, List param, NumericVector grid_a, NumericVector grid_a_NS, NumericVector grid_z, NumericMatrix prob_z, double wToday, double RToday, double tauToday, double dToday, NumericMatrix pf_c, NumericMatrix pf_n, NumericMatrix pf_sav) {

	int n_a_ns = as<int>(param["n_a_ns"]);
	int n_z = as<int>(param["n_z"]);
	double min_a_ns = as<double>(param["min_a_ns"]);
	double max_a_ns = as<double>(param["max_a_ns"]);
	double bl = as<double>(param["bl"]);

	SparseMatrix<double> AA(n_a_ns*n_z, n_a_ns*n_z);

	// Non-stochastic simulations for period t.
	for(int i_z=0; i_z < n_z; i_z++){ 

		double zToday = grid_z[i_z];

		for(int i_a_NS = 0; i_a_NS < n_a_ns; i_a_NS++){ // Fix i_a_NS ← 謎コメント

			double aToday = grid_a_NS[i_a_NS]; 
			if(aToday > max_a_ns)
				aToday = max_a_ns;

			// Compute savings of households whose state variables are a_{i_a_NS} and z_{i_z}, using the policy function for savings.
			double aTomorrow = pchip(grid_a, pf_sav(_, i_z), aToday);

			int indexToday = n_a_ns*i_z + i_a_NS; // Matlabとインデックス開始位置が異なる（Matlab:1, C++:0)

			// Find the index j such that grid_a_NS(j) < aTomorrow <= grid_a_NS(j+1).
			int j_a_NS = which(grid_a_NS, n_a_ns - 1, aTomorrow);

			/*
			% Redistribute the current mass, Hist(i_a_NS,i_z), to the points (grid_a_NS(j), z(1)), (grid_a_NS(j), z(2)), (grid_a_NS(j), z(3)),
			% (grid_a_NS(j+1), z(1)), (grid_a_NS(j+1), z(2)), and (grid_a_NS(j+1), z(3)) according to the weights, weight_NS*prob(i_z, 1),  weight_NS*prob(i_z, 2),  weight_NS*prob(i_z, 3),
			% (1d0 - weight_NS)*prob(i_z, 1),  (1d0 - weight_NS)*prob(i_z, 2), and (1d0 - weight_NS)*prob(i_z, 3).
			% Note that 'Hist_up' denotes the end-of-period t distribution, while 'Hist' denotes the beginning-of-period t distribution.
			*/

			if ( j_a_NS == 0 ){

				for(int j_z =0; j_z < n_z; j_z++){
					int indexTomorrow = n_a_ns*j_z; // Matlabとインデックス開始位置が異なる（Matlab:1, C++:0)
					AA.insert(indexTomorrow, indexToday) = prob_z(i_z, j_z);
				}

			} else if ( j_a_NS == n_a_ns ){

				for(int j_z =0; j_z < n_z; j_z++){
					int indexTomorrow = n_a_ns*j_z + n_a_ns - 1; // Matlabとインデックス開始位置が異なる（Matlab:1, C++:0)
					AA.insert(indexTomorrow, indexToday) = prob_z(i_z, j_z);
				}
			} else {

				double weight_NS = 1.0 - ((aTomorrow - grid_a_NS[j_a_NS])/(grid_a_NS[j_a_NS + 1] - grid_a_NS[j_a_NS]));

				for(int j_z=0; j_z < n_z; j_z++){
					int indexTomorrow = n_a_ns*j_z + j_a_NS - 1; // Matlabとインデックス開始位置が異なる（Matlab:1, C++:0)
					AA.insert(indexTomorrow, indexToday) = prob_z(i_z, j_z)*weight_NS;
					AA.insert(indexTomorrow + 1, indexToday) = prob_z(i_z, j_z)*(1.0 - weight_NS);
				}
			}
		}
	}

	// 以下、Spectraのリファレンス通り
	
	// Construct matrix operation object using the wrapper class SparseGenMatProd
	SparseGenMatProd<double> op(AA); // insert時にindexTomorrowとindexTodayを入れ替え.transpose()済み

	// Construct eigen solver object, requesting the largest three eigenvalues
	GenEigsSolver<double, LARGEST_MAGN, SparseGenMatProd<double>> eigs(&op, 1, 6);

	// Initialize and compute
	eigs.init();

	int nconv = eigs.compute();

	// Retrieve results
	if(eigs.info() != SUCCESSFUL){
		stop("eigs.info() != SUCCESSFUL");
	}

//	VectorXcd evalues = eigs.eigenvalues(); // 固有値を取得
	MatrixXcd evectors = eigs.eigenvectors(1); // 固有ベクトル1つの行列を取得（大きい順に並んでるので、最大になる）
	VectorXd evector = evectors.col(0).real(); // ベクトルの実数部を取り出す
	VectorXd elements = evector/evector.sum(); // 合計で割り算する

	// Eigenの疎行列をRcppの行列に変換するのには捻りがいる

	std::vector<double> dv; // doubleのvectorを作る
	dv.resize(evector.size()); // サイズをあらかじめ拡張
	Eigen::Map<VectorXd>(&dv[0], elements.size()) = elements; // doubleのvectorにMapを使って代入
	NumericMatrix dist(n_a_ns, n_z, dv.begin()); // doubleのvectorをRcpp形式の行列に変換

	double meanC = 0.0;
	double meanN = 0.0;
	double meanA = 0.0;
//	meanV = 0.0;

	for(int i_z=0;i_z<n_z;i_z++){

		double zToday = grid_z[i_z];

		for(int i_a_NS = 0; i_a_NS<n_a_ns; i_a_NS++){

			double aToday = grid_a_NS[i_a_NS];
			if(max_a_ns < aToday)
				aToday = max_a_ns;

			double cToday = pchip(grid_a, pf_c(_, i_z), aToday);
			double nToday = pchip(grid_a, pf_n(_, i_z), aToday);
			double aTomorrow = pchip(grid_a, pf_sav(_, i_z), aToday);

			if (aTomorrow < bl)
			    aTomorrow = bl;

			/*
			%         if (i_a_NS == 1)
			%             disp([cToday nToday aTomorrow]); 
			%             pause;
			%         end
			*/

			meanC += dist(i_a_NS, i_z)*cToday;
			meanN += dist(i_a_NS, i_z)*zToday*nToday; // Note: meanN is calculated as the mean value of the 'effective' labor supply, labor supply multiplied by idiosyncratic labor productivity.
			meanA += dist(i_a_NS, i_z)*aTomorrow; // Note: The original code uses aToday instead of aTomorrow, which was wrong.
			// meanV = meanV + (dist(i_a_NS,i_z)/sum(dist))*vToday
		}
	}

	return List::create(Named("dist")=dist, Named("meanC")=meanC, Named("meanN")=meanN, Named("meanA")=meanA);
}

/*
	a ∈ [array[0], array[size-1]]でarray[m] < a ≤ array[m+1]となるmを返す
	a < array[0]のときは0、a > array[size-1]のときはsize-1
*/
int	which(const NumericVector array, int size, double a)
{
	int	l, r, m;

	l = 0;
	r = size - 1;
	/*
		if(a >= array[size-1]) return size-1; をループ前につけて、
		while(l + 1 < r){ にして、if(a < array[m + 1]){ ... }を無くして、
		l = m + 1;をl = m;にした方がたぶん速い。
	*/
	while(l < r){
		m = (l + r)/2;

		if(a > array[m]){
			if(a <= array[m + 1]){ /* 端数切り捨てによりmは必ずm<size-1 */
				return m;
			}
			l = m + 1;/* int型変数を/2をすると切り捨てになるので、下端を動かすときは+1しないと収束しない */
		} else {
			r = m;
		}
	}
	return (l + r)/2;
}

