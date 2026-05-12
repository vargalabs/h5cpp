#include <itpp/itbase.h>

using std::cout;
using std::endl;
using namespace itpp;

int main(void)
{
  cout << "====================================================" << endl;
  cout << "              Test of fastmath" << endl;
  cout << "====================================================" << endl;

  mat m0("1 2 3;4 5 6;7 8 9"), mv0("2;3;1");
  vec v0("2 3 1");

  cout << "sub_v_vT_m: the slow and fast way" << endl;
  cout << (m0 - mv0*transpose(mv0)*m0) << endl;
  sub_v_vT_m(m0, v0);
  cout << m0 << endl;

  m0 = "1 2 3;4 5 6;7 8 9";

  cout << endl << "sub_m_v_vT: the slow and fast way" << endl;
  cout << (m0 - m0*mv0*transpose(mv0)) << endl;
  sub_m_v_vT(m0, v0);
  cout << m0 << endl;

  return 0;

}
