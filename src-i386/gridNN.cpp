#include <Rcpp.h>
using namespace Rcpp;

//' Nearest Neighbor Gridding
//'
//' @export
// [[Rcpp::export]]
NumericVector gridNN(NumericVector gx, NumericVector gy, NumericVector x, NumericVector y, NumericVector z,
                      double p, double xscale, double yscale, double uncertainty) {
  int n = gx.size();
  int nn = x.size();

  NumericVector out(n);
  double temp;
  double temp2;
  double w;

  for(int i = 0; i < n; ++i) {
    temp = pow(abs(x[0] - gx[i]), p) + pow(abs(y[0] - gy[i]), p);
    w = 0;


    for(int j = 1; j < nn; j++) {
      temp2 = pow(abs(x[j] - gx[i]), p) + pow(abs(y[j] - gy[i]), p);

      if (temp2 < temp) {
        w = j;
        temp = temp2;
      }
    }

    out[i] = z[w];
  }
  return out;
}
