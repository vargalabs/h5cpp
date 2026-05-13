#include <itpp/itcomm.h>

using namespace itpp;

//These lines are needed for use of cout and endl
using std::cout;
using std::endl;

int main()
{
  //Scalars
  int constraint_length, MaxNrofErrors, Nobits, MaxIterations, p, i;
  double Ec, Eb;

  //Vectors
  ivec generators;
  vec EbN0dB, EbN0, N0, ber, trans_symbols, rec_symbols;
  bvec uncoded_bits, coded_bits, decoded_bits;

  //Classes
  BPSK bpsk;
  BERC berc;
  Convolutional_Code conv_code;
  AWGN_Channel channel;

  /*
  Set up the convolutional encoder/decoder class:
  The generators are given in octal form by adding a zero in front of the numbers.
  In this example we will simulate a rate 1/3 code that is listed in J. G. Proakis,
  "Digital communications". The encoder has constraint length 7.
  */
  generators.set_size(3, false);
  generators(0) = 0133;
  generators(1) = 0145;
  generators(2) = 0175;
  constraint_length = 7;
  conv_code.set_generator_polynomials(generators, constraint_length);

  //Init: Calculate some simulation specific parameters:
  Ec = 1.0;
  EbN0dB = linspace(-2, 6, 5);
  EbN0 = inv_dB(EbN0dB);
  Eb = Ec / conv_code.get_rate();
  N0 = Eb * pow(EbN0, -1);
  MaxNrofErrors = 100;
  Nobits = 10000;
  MaxIterations = 10;
  ber.set_size(EbN0dB.length(), false);
  ber.clear();

  //Randomize the random number generators.
  RNG_randomize();

  for (p = 0; p < EbN0dB.length(); p++) {

    cout << "Now simulating point " << p + 1 << " out of " << EbN0dB.length() << endl;
    berc.clear();                 //Clear the bit error rate counter.
    channel.set_noise(N0(p) / 2.0); //Set the noise value of the AWGN channel.

    for (i = 0; i < MaxIterations; i++) {

      uncoded_bits = randb(Nobits);                   //The uncoded bits.
      coded_bits = conv_code.encode(uncoded_bits);    //The convolutional encoder function.
      bpsk.modulate_bits(coded_bits, trans_symbols); //The BPSK modulator.
      rec_symbols = channel(trans_symbols);           //The AWGN channel.
      decoded_bits = conv_code.decode(rec_symbols);   //The Viterbi decoder function.
      berc.count(uncoded_bits, decoded_bits);         //Count the errors.
      ber(p) = berc.get_errorrate();

      //Break the simulation on this point if sufficient number of bit errors were observed:
      if (berc.get_errors() > MaxNrofErrors) {
        cout << "Breaking on point " << p + 1 << " with " << berc.get_errors() << " errors." << endl;
        break;
      }

    }
  }

  //Print the results:
  cout << "BER    = " << ber  << endl;
  cout << "EbN0dB = " << EbN0dB << endl;

  //Exit program:
  return 0;

}
