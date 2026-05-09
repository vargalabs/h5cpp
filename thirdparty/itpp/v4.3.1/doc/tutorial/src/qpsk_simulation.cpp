#include <itpp/itcomm.h>

using namespace itpp;

//These lines are needed for use of cout and endl
using std::cout;
using std::endl;

int main()
{
  //Declarations of scalars and vectors:
  int i, Number_of_bits;
  double Ec, Eb;
  vec EbN0dB, EbN0, N0, noise_variance, bit_error_rate; //vec is a vector containing double
  bvec transmitted_bits, received_bits;                 //bvec is a vector containing bits
  cvec transmitted_symbols, received_symbols;           //cvec is a vector containing double_complex

  //Declarations of classes:
  QPSK qpsk;                     //The QPSK modulator class
  AWGN_Channel awgn_channel;     //The AWGN channel class
  it_file ff;                    //For saving the results to file
  BERC berc;                     //Used to count the bit errors
  Real_Timer tt;                 //The timer used to measure the execution time

  //Reset and start the timer:
  tt.tic();

  //Init:
  Ec = 1.0;                      //The transmitted energy per QPSK symbol is 1.
  Eb = Ec / 2.0;                 //The transmitted energy per bit is 0.5.
  EbN0dB = linspace(0.0, 9.0, 10); //Simulate for 10 Eb/N0 values from 0 to 9 dB.
  EbN0 = inv_dB(EbN0dB);         //Calculate Eb/N0 in a linear scale instead of dB.
  N0 = Eb * pow(EbN0, -1.0);     //N0 is the variance of the (complex valued) noise.
  Number_of_bits = 100000;       //One hundred thousand bits is transmitted for each Eb/N0 value

  //Allocate storage space for the result vector.
  //The "false" argument means "Do not copy the old content of the vector to the new storage area."
  bit_error_rate.set_size(EbN0dB.length(), false);

  //Randomize the random number generators in it++:
  RNG_randomize();

  //Iterate over all EbN0dB values:
  for (i = 0; i < EbN0dB.length(); i++) {

    //Show how the simulation progresses:
    cout << "Now simulating Eb/N0 value number " << i + 1 << " of " << EbN0dB.length() << endl;

    //Generate a vector of random bits to transmit:
    transmitted_bits = randb(Number_of_bits);

    //Modulate the bits to QPSK symbols:
    transmitted_symbols = qpsk.modulate_bits(transmitted_bits);

    //Set the noise variance of the AWGN channel:
    awgn_channel.set_noise(N0(i));

    //Run the transmited symbols through the channel using the () operator:
    received_symbols = awgn_channel(transmitted_symbols);

    //Demodulate the received QPSK symbols into received bits:
    received_bits = qpsk.demodulate_bits(received_symbols);

    //Calculate the bit error rate:
    berc.clear();                               //Clear the bit error rate counter
    berc.count(transmitted_bits, received_bits); //Count the bit errors
    bit_error_rate(i) = berc.get_errorrate();   //Save the estimated BER in the result vector

  }

  tt.toc();

  //Print the results:
  cout << endl;
  cout << "EbN0dB = " << EbN0dB << " [dB]" << endl;
  cout << "BER = " << bit_error_rate << endl;
  cout << "Saving results to ./qpsk_result_file.it" << endl;
  cout << endl;

  //Save the results to file:
  ff.open("qpsk_result_file.it");
  ff << Name("EbN0dB") << EbN0dB;
  ff << Name("ber") << bit_error_rate;
  ff.close();

  //Exit program:
  return 0;

}

