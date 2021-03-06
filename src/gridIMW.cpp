#include <Rcpp.h>
using namespace Rcpp;

//' Inverse Metropolis-distance Gridding
//'
//' @export
// [[Rcpp::export]]
NumericVector gridIMW(NumericVector gx, NumericVector gy, NumericVector x, NumericVector y, NumericVector z,
                      double p, double xscale, double yscale, double uncertainty) {
  int n = gx.size();
  int nn = x.size();

  NumericVector out(n);
  double temp;
  double wsum;
  double w;

  double deltamin = (powf(fabs(xscale/2.0), p) + powf(fabs(yscale/2.0), p)) * uncertainty;

  for(int i = 0; i < n; ++i) {
    temp = 0;
    w = 0;
    wsum = 0;

    for(int j = 0; j < nn; j++) {
      w = 1.0 / (fabs(x[j] - gx[i]) + fabs(y[j] - gy[i]) + deltamin);
      temp += w * z[j];
      wsum += w;
    }

    out[i] = temp / wsum;
  }
  return out;
}
