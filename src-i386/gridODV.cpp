#include <Rcpp.h>
using namespace Rcpp;
// @export
// [[Rcpp::export]]
NumericVector gridODV(NumericVector gx, NumericVector gy, NumericVector x, NumericVector y, NumericVector z,
                      double p, double xscale, double yscale, double uncertainty) {
  int n = gx.size();
  int nn = x.size();

  NumericVector out(n);
  double temp;
  double wsum;
  double w;

  double deltamin = (pow(abs(xscale/2), p) + pow(abs(yscale/2), p)) * uncertainty;

  for(int i = 0; i < n; ++i) {
    temp = 0;
    w = 0;
    wsum = 0;

    for(int j = 0; j < nn; j++) {
      w = 1 / exp(-1*pow(abs(x[j] - gx[i]), 2) - pow(abs(y[j] - gy[i]), 2) - deltamin);
      temp += w * z[j];
      wsum += w;
    }

    out[i] = temp / wsum;
  }
  return out;
}
