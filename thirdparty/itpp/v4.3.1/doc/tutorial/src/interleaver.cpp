#include <itpp/itcomm.h>

using namespace itpp;

//These lines are needed for use of cout and endl
using std::cout;
using std::endl;

int main()
{
  //Declare scalars and vectors:
  int rows, cols;
  ivec input, output, deinterleaved;

  //Declare the interleaver. The interleaver classes are templated, and therefore we must specify
  //the type of the data elements. In this example we are using integers:
  Block_Interleaver<int> my_interleaver;

  //Initialize the interleaver class. Note that this can be done already in the declaration by writing
  //Block_Interleaver<int> my_interleaver(rows,cols);
  rows = 4;
  cols = 5;
  my_interleaver.set_rows(rows);
  my_interleaver.set_cols(cols);

  //Define the input to the interleaver:
  input = "1:20";

  //Do the interleaving:
  output = my_interleaver.interleave(input);

  //Do the de-interleaving:
  deinterleaved = my_interleaver.deinterleave(output);

  //Print the results:
  cout << "input = " << input << endl;
  cout << "output = " << output << endl;
  cout << "deinterleaved = " << deinterleaved << endl;

  //Exit program:
  return 0;

}
