#include <itpp/itbase.h>

using namespace itpp;

//These lines are needed for use of cout and endl
using std::cout;
using std::endl;

int main()
{
  //Declare the scalars used:
  long i, sum, N;

  //Declare tt as an instance of the timer class:
  Real_Timer tt;

  //Initiate the variables:
  N = 1000000;
  sum = 0;

  //Start and reset the timer:
  tt.tic();

  //Do some processing
  for (i = 0; i < N; i++) {
    sum += i;
  }

  // Print the elapsed time
  tt.toc_print();

  //Print the result of the processing:
  cout << "The sum of all integers from 0 to " << N - 1 << " equals " << sum << endl;

  //Exit program:
  return 0;

}
