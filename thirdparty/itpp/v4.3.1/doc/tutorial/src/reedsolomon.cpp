#include <itpp/itcomm.h>

using namespace itpp;

//These lines are needed for use of cout and endl
using std::cout;
using std::endl;

int main()
{

  //Scalars and vectors:
  int m, t, n, k, q, NumBits, NumCodeWords;
  double p;
  bvec uncoded_bits, coded_bits, received_bits, decoded_bits;

  //Set parameters:
  NumCodeWords = 1000;  //Number of Reed-Solomon code-words to simulate
  p = 0.01;             //BSC Error probability
  m = 3;                //Reed-Solomon parameter m
  t = 2;                //Reed-Solomon parameter t

  cout << "Number of Reed-Solomon code-words to simulate: " <<  NumCodeWords << endl;
  cout << "BSC Error probability : " << p << endl;
  cout << "RS m: " << m << endl;
  cout << "RS t: " << t << endl;

  //Classes:
  Reed_Solomon reed_solomon(m, t);
  BSC bsc(p);
  BERC berc;

  RNG_randomize();

  //Calculate parameters for the Reed-Solomon Code:
  n = round_i(pow(2.0, m) - 1);
  k = round_i(pow(2.0, m)) - 1 - 2 * t;
  q = round_i(pow(2.0, m));

  cout << "Simulating an Reed-Solomon code with the following parameters:" << endl;
  cout << "n = " << n << endl;
  cout << "k = " << k << endl;
  cout << "q = " << q << endl;

  NumBits = m * k * NumCodeWords;
  uncoded_bits = randb(NumBits);
  coded_bits = reed_solomon.encode(uncoded_bits);
  received_bits = bsc(coded_bits);
  decoded_bits = reed_solomon.decode(received_bits);

  berc.count(uncoded_bits, decoded_bits);
  cout << "The bit error probability after decoding is " << berc.get_errorrate() << endl;

  //Exit program:
  return 0;

}
