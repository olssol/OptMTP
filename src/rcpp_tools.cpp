#include <Rcpp.h>
#include <Rmath.h>

using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]

// [[Rcpp::init]]
void my_package_init(DllInfo *dll) {
  // initialization code here
  R_useDynamicSymbols(dll, TRUE);
}


//' Test Rcpp function
//'
//'
//' @param test test parameter
//'
//' @export
// [[Rcpp::export]]
double c_test(double test) {
  return pow(test,2);
}


//' Multiple testing following the graph
//'
//' @param alphas    vector of the original alphas
//' @param mat_g     transition matrix G
//' @param p_values  vector of p-values for the elementary hypothesis
//' @param log       TRUE: print log at each step; FALSE: silent
//'
//' @return Hypothesis rejection status indicator vector
//'
//' @export
// [[Rcpp::export]]
IntegerVector c_mtp(NumericVector p_values, NumericVector alphas,
                    NumericMatrix mat_g, bool log = false) {

		int m = p_values.size();
    IntegerVector h_ind(m, 1);
		IntegerVector rst(m, 0);

		int n_h = m, fstop = 0;
		int i, j, l, k;
		double j_pal, pal;

		//testing algorithm
		while (n_h > 0 & 0 == fstop) {
				// print log
				if (log) {
						Rcout << "-------------------------- \n";
						Rcout << "----Rejection status: \n";
						Rf_PrintValue(rst);
						Rcout << "----Alphas:   \n";
						Rf_PrintValue(alphas);
						Rcout << "----p-values: \n";
						Rf_PrintValue(p_values);
						Rcout << "----G-Matrix: \n";
						Rf_PrintValue(mat_g);
				}

        // arg_min pval / alpha
				j      = 0;
				j_pal  = R_PosInf;
				for (i = 0; i < m; i ++) {
						if (0 == h_ind[i])
								continue;
						pal = p_values[i] / alphas[i];
						if (pal < j_pal) {
								j_pal = pal;
								j     = i;
						}
				}

				// test hypothesis j
				if (p_values[j] > alphas[j]) {
						fstop = 1;
						continue;
				}

				// update tests
				rst[j]   = 1;
				h_ind[j] = 0;
				n_h--;

				if (0 == n_h)
						continue;

				// update alpha
				for (i = 0; i < m; i ++) {
						if (0 == h_ind[i])
								continue;

						alphas[i] += alphas[j] * mat_g(j, i);
				}
				alphas[j] = 0;

				// update G
				for (l = 0; l < m; l++) {
						for (k = 0; k < m; k++) {
								if (0 == h_ind[l] |
										0 == h_ind[k] |
										l == k)
										continue;

								mat_g(l, k) += mat_g(j, k) * mat_g(l, j);
								mat_g(l, k) /= 1 - mat_g(j, l) * mat_g(l, j);
						}
				}

				for (l = 0; l < m; l++) {
						for (k = 0; k < m; k++) {
								if (0 == h_ind[l] |
										0 == h_ind[k] |
										l == k)
										mat_g(l, k) = 0.0;
						}
				}
		}

		// return
		return(rst);
}
