#include <itpp/itstat.h>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <ios>

using std::cout;
using std::endl;
using std::fixed;
using std::setprecision;

using namespace itpp;

int main()
{

  bool print_progress = false;

  //
  // first, let's generate some synthetic data

  int N = 100000;  // number of vectors
  int D = 3;       // number of dimensions
  int K = 5;       // number of Gaussians

  Array<vec> X(N);
  for (int n = 0;n < N;n++) { X(n).set_size(D); X(n) = 0.0; }

  // the means

  Array<vec> mu(K);
  mu(0) = "-6, -4, -2";
  mu(1) = "-4, -2,  0";
  mu(2) = "-2,  0,  2";
  mu(3) = " 0, +2, +4";
  mu(4) = "+2, +4, +6";


  // the diagonal variances

  Array<vec> var(K);
  var(0) = "0.1, 0.2, 0.3";
  var(1) = "0.2, 0.3, 0.1";
  var(2) = "0.3, 0.1, 0.2";
  var(3) = "0.1, 0.2, 0.3";
  var(4) = "0.2, 0.3, 0.1";

  cout << fixed << setprecision(3);
  cout << "user configured means and variances:" << endl;
  cout << "mu = " << mu << endl;
  cout << "var = " << var << endl;

  // randomise the order of Gaussians "generating" the vectors
  I_Uniform_RNG rnd_uniform(0, K - 1);
  ivec gaus_id = rnd_uniform(N);

  ivec gaus_count(K);
  gaus_count = 0;
  Array<vec> mu_test(K);
  for (int k = 0;k < K;k++) { mu_test(k).set_size(D); mu_test(k) = 0.0; }
  Array<vec> var_test(K);
  for (int k = 0;k < K;k++) { var_test(k).set_size(D); var_test(k) = 0.0; }

  Normal_RNG rnd_normal;
  for (int n = 0;n < N;n++) {

    int k = gaus_id(n);
    gaus_count(k)++;

    for (int d = 0;d < D;d++) {
      rnd_normal.setup(mu(k)(d), var(k)(d));
      double tmp = rnd_normal();
      X(n)(d) = tmp;
      mu_test(k)(d) += tmp;
    }
  }

  //
  // find the stats for the generated data

  for (int k = 0;k < K;k++) mu_test(k) /= gaus_count(k);

  for (int n = 0;n < N;n++) {
    int k = gaus_id(n);

    for (int d = 0;d < D;d++) {
      double tmp = X(n)(d) - mu_test(k)(d);
      var_test(k)(d) += tmp * tmp;
    }
  }

  for (int k = 0;k < K;k++)  var_test(k) /= (gaus_count(k) - 1.0);

  cout << endl << endl;
  cout << fixed << setprecision(3);
  cout << "stats for X:" << endl;

  for (int k = 0;k < K;k++) {
    cout << "k = " << k << "  count = " << gaus_count(k) << "  weight = " << gaus_count(k) / double(N) << endl;
    for (int d = 0;d < D;d++) cout << "  d = " << d << "  mu_test = " << mu_test(k)(d) << "  var_test = " << var_test(k)(d) << endl;
    cout << endl;
  }


  // make a model with initial values (zero mean and unit variance)
  // the number of gaussians and dimensions of the model is specified here

  MOG_diag mog(K, D);

  cout << endl;
  cout << fixed << setprecision(3);
  cout << "mog.avg_log_lhood(X) = " << mog.avg_log_lhood(X) << endl;

  //
  // find initial parameters via k-means (which are then used as seeds for EM based optimisation)

  cout << endl << endl;
  cout << "running kmeans optimiser" << endl << endl;

  MOG_diag_kmeans(mog, X, 10, 0.5, true, print_progress);

  cout << fixed << setprecision(3);
  cout << "mog.get_means() = " << endl << mog.get_means() << endl;
  cout << "mog.get_diag_covs() = " << endl << mog.get_diag_covs() << endl;
  cout << "mog.get_weights() = " << endl << mog.get_weights() << endl;

  cout << endl;
  cout << "mog.avg_log_lhood(X) = " << mog.avg_log_lhood(X) << endl;


  //
  // EM ML based optimisation

  cout << endl << endl;
  cout << "running ML optimiser" << endl << endl;

  MOG_diag_ML(mog, X, 10, 0.0, 0.0, print_progress);

  cout << fixed << setprecision(3);
  cout << "mog.get_means() = " << endl << mog.get_means() << endl;
  cout << "mog.get_diag_covs() = " << endl << mog.get_diag_covs() << endl;
  cout << "mog.get_weights() = " << endl << mog.get_weights() << endl;

  cout << endl;
  cout << "mog.avg_log_lhood(X) = " << mog.avg_log_lhood(X) << endl;

  return 0;
}
