#include <itpp/itbase.h>

using std::cout;
using std::endl;
using namespace itpp;


double f(const double x)
{
  return x*log(x);
}

int main()
{
  cout << "=========================================" << endl;
  cout << "  Test of numerical integration routines" << endl;
  cout << "=========================================" << endl;

  double res = itpp::quad(f, 1.5, 3.5);
  double res2 = itpp::quadl(f, 1.5, 3.5);

  cout << "Integration of f(x)=x*log(x) over [1.5,3.5]" << endl;
  cout << "quad = " << res << endl;
  cout << "quadl = " << res2 << endl;

  return 0;
}
