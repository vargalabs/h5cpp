#include <itpp/itbase.h>

using std::cout;
using std::endl;
using namespace itpp;

int main()
{
//********* Testing Sparse_Vec *********
  cout << endl << "******* Testing Sparse_Vec *******" << endl;

  Sparse_Vec<double> v1(5), v2(5);
  Sparse_Vec<double> v_out(5);
  ivec index_vec;
  vec v;

  v1.set(1, 3);
  v2.set(1, 1);
  v2.set(2, 2);
  cout << "v1=" << v1.full() << endl;
  cout << "v2=" << v2.full() << endl;
  v_out = v1 + v2;
  cout << "v1+v2=" << endl << v_out.full() << endl;
  cout << "v1*v2=" << v1*v2 << endl;
  cout << "density(v1)=" << v1.density() << endl;
  v1.zeros();
  cout << "v1=" << v1.full() << endl;

  index_vec = "0 2 2 3";
  v = "1 2 4 3";
  v1.set(index_vec, v);
  cout << "Set " << v << " to the indices " << index_vec << endl;
  cout << "v1=" << v1.full() << endl;
  v1.set_new(index_vec, v);
  cout << "Set new " << v << " to the indices " << index_vec << endl;
  cout << "Unnoticed error in set_new() if the same index is used several times" << endl;
  cout << "v1=" << v1.full() << endl;

  v1.zeros();
  v1 -= v2;
  cout << "v1(for v1-=v2) = " << v1.full() << endl;
  v1 /= 4;
  cout << "v1(for v1/=4) = " << v1.full() << endl;
  v1 *= 2;
  cout << "v1(for v1*=2) = " << v1.full() << endl;
  v2 /= 2;
  v1 += v2;
  cout << "v1(for v1+=v2/2) = " << v1.full() << endl;
  cout << "number of non-zero elements in v1: nnz = " << v1.nnz() << endl;

  index_vec = "0 2 2 1";
  v = "-1 5 -2 4";
  cout << "Add " << v << " to the indices " << index_vec << endl;
  v1.add(index_vec, v);
  cout << "v1=" << v1.full() << endl;

  v1.clear_elem(2);
  cout << "v1(element 2 cleared) = " << v1.full() << endl;


  Sparse_Vec<double> v3;
  v3.set_size(2);
  v3.set(0, 3);
  cout << "v3=" << endl << v3.full() << endl;
  v3.set_size(5);
  cout << "v3=" << v3.full() << endl;

  Sparse_Vec<double> v4(5), v5(5);
  v4.set(1, 1);
  v5.set(1, 1);
  v4.set(2, 2);
  v5.set(2, 4);
  cout << "v4=" << v4.full() << endl;
  cout << "v5=" << v5.full() << endl;
  if (v4 == v5)
    cout << "v4 equal to v5" << endl;
  else
    cout << "v4 not equal to v5" << endl;


//********* Testing Sparse_Mat *********
  cout << endl << "******* Testing Sparse_Mat *******" << endl;

  Sparse_Mat<double> m1(3, 3), m2(3, 3);

  m1.set(1, 1, 3);
  m2.set(1, 2, 1);
  m2.set(2, 1, 2);
  cout << "m1=" << endl << full(m1) << endl;
  cout << "m2=" << endl << full(m2) << endl;
  cout << "m1+m2=" << endl << full(m1 + m2) << endl;
  cout << "m1*m2=" << endl << full(m1*m2) << endl;
  cout << "density(m1)=" << m1.density() << endl;
  cout << "transpose(m2)=" << endl << full(transpose(m2)) << endl;

  m1.zeros();
  cout << "m1.zeros()=" << endl << full(m1) << endl;
  m1 -= m2;
  cout << "m1(-=m2)=" << endl << full(m1) << endl;
  m1 /= 4;
  cout << "m1(/=4)=" << endl << full(m1) << endl;
  m1 *= 2;
  cout << "m1(*=2)=" << endl << full(m1) << endl;

  m1.add_elem(0, 2, 4);
  cout << "m1.add_elem(0,2,4)=" << endl << full(m1) << endl;
  m1.clear_elem(1, 2);
  cout << "m1.clear_elem(1,2)=" << endl << full(m1) << endl;

  Sparse_Mat<double> m3(2, 3), m4(3, 2);

  m3.set(0, 0, 3);
  m4.set(0, 1, 1);
  m4.set(2, 0, 2);

  cout << "m3=" << endl << full(m3) << endl;
  cout << "m4=" << endl << full(m4) << endl;
  cout << "m3*m4=" << endl << full(m3*m4) << endl;
  cout << "transpose(m4)=" << endl << full(transpose(m4)) << endl;

  cout << "trans_mult(m4,m4)=" << full(trans_mult(m4, m4)) << endl;
  cout << "mult_trans(m4,m4)=" << full(mult_trans(m4, m4)) << endl;
  cout << "mult_trans(m4,m4)=" << full(mult_trans(m4, m4)) << endl;

  v = "1 2 3";
  cout << "v = " << v << endl;
  cout << "m3 * v = " << m3*v << endl;


  Sparse_Mat<double> A(3, 5), B(5, 3), C;
  vec x(3), y(5), z1, z2;

  A = randn(3, 5);
  B = randn(5, 3);
  x = randn(3);
  y = randn(5);

  cout << "A = " << full(A) << endl;
  cout << "B = " << full(B) << endl;
  cout << "x = " << x << endl;
  cout << "y = " << y << endl;

  C = A * B;
  z1 = A * y;
  z2 = x * A;

  cout << "C = " << full(C) << endl;
  cout << "z1 = " << z1 << endl;
  cout << "z2 = " << z2 << endl;

  return 0;
}
