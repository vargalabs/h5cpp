#include <itpp/itbase.h>

using namespace itpp;

//These lines are needed for use of cout and endl
using std::cout;
using std::endl;

int main()
{
  //Declare vectors and matricies:
  vec a, b, c;
  mat A, B;

  //Use the function linspace to define a vector:
  a = linspace(1.0, 2.0, 10);

  //Use a string of values to define a vector:
  b = "0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0";

  //Add two vectors:
  c = a + b;

  //Print results:
  cout << "a = " << a << endl;
  cout << "b = " << b << endl;
  cout << "c = " << c << endl;

  //Use a string to define a matrix:
  A = "1.0 2.0;3.0 4.0";

  //Calculate the inverse of matrix A:
  B = inv(A);

  //Print results:
  cout << "A = " << A << endl;
  cout << "B = " << B << endl;

  //Exit program:
  return 0;

}
