#include <itpp/itcomm.h>

using namespace itpp;

//These lines are needed for use of cout and endl
using std::cout;
using std::endl;

int main()
{
  //Scalars
  int N;
  double N0;

  //Vectors
  bvec bits, dec_bits;
  vec symbols, rec;

  //Classes
  BPSK bpsk;  //The BPSK modulator/debodulator class
  BERC berc;  //The Bit Error Rate Counter class

  //Init
  N = 500000; //The number of bits to simulate
  N0 = 1;     //0 dB SNR

  //Randomize the random number generator
  RNG_randomize();

  //Generate the bits:
  bits = randb(N);

  //Do the BPSK modulation
  bpsk.modulate_bits(bits, symbols);

  //Add the AWGN
  rec = symbols + sqrt(N0 / 2) * randn(N);

  //Decode the received bits
  bpsk.demodulate_bits(rec, dec_bits);

  //Count the number of errors
  berc.count(bits, dec_bits);

  //Print the results
  cout << "There were " << berc.get_errors() << " received bits in error." << endl;
  cout << "There were " << berc.get_corrects() << " correctly received bits." << endl;
  cout << "The error probability was " << berc.get_errorrate() << endl;
  cout << "The theoretical error probability is " << 0.5*erfc(1.0) << endl;

  //Exit program:
  return 0;

}
