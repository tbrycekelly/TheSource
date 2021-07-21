#include <Rcpp.h>
using namespace Rcpp;

//' Inverse Distance Gridding (as in ODV)
//'
//' @export
// [[Rcpp::export]]
NumericVector gridODV(NumericVector gx, NumericVector gy, NumericVector x, NumericVector y, NumericVector z,
                      double p, double xscale, double yscale, double uncertainty) {
  int n = gx.size();
  int nn = x.size();

  NumericVector out(n);
  double temp;
  double wsum;
  double w;

  double deltamin = (powf(fabs(xscale/2.0)/xscale, p) + powf(fabs(yscale/2.0)/yscale, p)) * uncertainty;

  for(int i = 0; i < n; ++i) {
    temp = 0;
    w = 0;
    wsum = 0;

    for(int j = 0; j < nn; j++) {
      w = expf(-1*(powf(abs(x[j] - gx[i])/xscale, 2) + powf(abs(y[j] - gy[i])/yscale, 2) + deltamin));
      temp += w * z[j];
      wsum += w;
    }

    out[i] = temp / wsum;
  }
  return out;
}
