#include <itpp/itcomm.h>

using namespace itpp;

//These lines are needed for use of cout and endl
using std::cout;
using std::endl;

int main()
{
  //Scalars:
  int SF, Ncode, sc, i, j, NumOfBits, MaxIterations, MaxNrOfErrors, MinNrOfErrors;
  double Eb;

  //Vectors and matrixes:
  vec EbN0dB, EbN0, N0, ber;
  smat all_codes, spread_codesI, spread_codesQ;
  bvec transmited_bits, received_bits;
  cvec transmited_symbols, received_symbols, transmited_chips, received_chips;

  //Classes:
  Multicode_Spread_2d mc_spread;
  AWGN_Channel channel;
  QPSK qpsk;
  BERC berc;

  //Initialize the spreading:
  SF = 4;                            //The spreading factor is 4
  Ncode = 4;                         //Use all four codes in the multi-code spreader
  all_codes = to_smat(hadamard(SF)); //Calculate the spreading codes

  //Set the spreading codes:
  spread_codesI.set_size(Ncode, SF, false);
  spread_codesQ.set_size(Ncode, SF, false);
  for (sc = 0; sc < Ncode; sc++) {
    spread_codesI.set_row(sc, all_codes.get_row(sc));
    spread_codesQ.set_row(sc, all_codes.get_row(sc));
  }
  mc_spread.set_codes(to_mat(spread_codesI), to_mat(spread_codesQ));

  //Specify simulation specific variables:
  EbN0dB = linspace(-2, 14, 17);
  EbN0 = pow(10, EbN0dB / 10);
  Eb = 1.0;
  N0 = Eb * pow(EbN0, -1.0);
  NumOfBits = 50000;
  MaxIterations = 10;
  MaxNrOfErrors = 200;
  MinNrOfErrors = 5;
  ber.set_size(EbN0dB.size(), false);
  ber.clear();

  //Randomize the random number generators:
  RNG_randomize();

  for (i = 0; i < EbN0dB.length(); i++) {

    cout << endl << "Simulating point nr " << i + 1 << endl;
    berc.clear();
    channel.set_noise(N0(i) / 2.0);

    for (j = 0; j < MaxIterations; j++) {

      transmited_bits = randb(NumOfBits);
      transmited_symbols  = qpsk.modulate_bits(transmited_bits);

      //This is where we do the multi-code spreading:
      transmited_chips =  mc_spread.spread(transmited_symbols);

      received_chips = channel(transmited_chips);

      //The multi-code despreading:
      //The second argument tells the despreader that the offset is zero chips.
      //This offset is usefull on channels with delay.
      received_symbols = mc_spread.despread(received_chips, 0);

      received_bits = qpsk.demodulate_bits(received_symbols);

      berc.count(transmited_bits, received_bits);
      ber(i) = berc.get_errorrate();

      cout << "   Itteration " << j + 1 << ": ber = " << berc.get_errorrate() << endl;
      if (berc.get_errors() > MaxNrOfErrors) {
        cout << "Breaking on point " << i + 1 << " with " << berc.get_errors() << " errors." << endl;
        break;
      }

    }

    if (berc.get_errors() < MinNrOfErrors) {
      cout << "Exiting Simulation on point " << i + 1 << endl;
      break;
    }

  }

  //Print results:
  cout << endl;
  cout << "EbN0dB = " << EbN0dB << endl;
  cout << "ber = " << ber << endl;

  //Exit program:
  return 0;

}

