#include <itpp/itcomm.h>

using namespace itpp;
using namespace std;


int main(int argc, char **argv)
{
  // -- modulation and channel parameters (taken from command line input) --
  int nC;                    // type of constellation  (1=QPSK, 2=16-QAM, 3=64-QAM)
  int nRx;                   // number of receive antennas
  int nTx;                   // number of transmit antennas
  int Tc;                    // coherence time (number of channel vectors with same H)

  if (argc != 5) {
    cout << "Usage: cm nTx nRx nC Tc" << endl << "Example: cm 2 2 1 100000 (2x2 QPSK MIMO on slow fading channel)" << endl;
    exit(1);
  }
  else {
    sscanf(argv[1], "%i", &nTx);
    sscanf(argv[2], "%i", &nRx);
    sscanf(argv[3], "%i", &nC);
    sscanf(argv[4], "%i", &Tc);
  }

  cout << "Initializing.. " << nTx << " TX antennas, " << nRx << " RX antennas, "
       << (1 << nC) << "-PAM per dimension, coherence time " << Tc << endl;

  // -- simulation control parameters --
  const vec EbN0db = "-5:0.5:50";        // SNR range
  const int Nmethods = 2;                 // number of demodulators to try
  const int Nbitsmax = 50000000;  // maximum number of bits to ever simulate per SNR point
  const int Nu = 1000;                   // length of data packet (before applying channel coding)

  int Nbers, Nfers;              // target number of bit/frame errors per SNR point
  double BERmin, FERmin;         // BER/FER at which to terminate simulation
  if (Tc == 1) {         // Fast fading channel, BER is of primary interest
    BERmin = 0.001;      // stop simulating a given method if BER<this value
    FERmin = 1.0e-10;    // stop simulating a given method if FER<this value
    Nbers = 1000;        // move to next SNR point after counting 1000 bit errors
    Nfers = 200;         // do not stop on this condition
  }
  else {               // Slow fading channel, FER is of primary interest here
    BERmin = 1.0e-15;    // stop simulating a given method if BER<this value
    FERmin = 0.01;       // stop simulating a given method if FER<this value
    Nbers = -1;          // do not stop on this condition
    Nfers = 200;         // move to next SNR point after counting 200 frame errors
  }

  // -- Channel code parameters --
  Convolutional_Code code;
  ivec generator(3);
  generator(0) = 0133;  // use rate 1/3 code
  generator(1) = 0165;
  generator(2) = 0171;
  double rate = 1.0 / 3.0;
  code.set_generator_polynomials(generator, 7);
  bvec dummy;
  code.encode_tail(randb(Nu), dummy);
  const int Nc = length(dummy);      // find out how long the coded blocks are

  // ============= Initialize ====================================

  const int Nctx = (int)(2 * nC * nTx * ceil(double(Nc) / double(2 * nC * nTx)));   // Total number of bits to transmit
  const int Nvec = Nctx / (2 * nC * nTx);    // Number of channel vectors to transmit
  const int Nbitspvec = 2 * nC * nTx;        // Number of bits per channel vector

  // initialize MIMO channel with uniform QAM per complex dimension and Gray coding
  ND_UQAM chan;
  chan.set_M(nTx, 1 << (2*nC));
  cout << chan << endl;

  // initialize interleaver
  Sequence_Interleaver<bin> sequence_interleaver_b(Nctx);
  Sequence_Interleaver<int> sequence_interleaver_i(Nctx);
  sequence_interleaver_b.randomize_interleaver_sequence();
  sequence_interleaver_i.set_interleaver_sequence(sequence_interleaver_b.get_interleaver_sequence());

  //  RNG_randomize();

  Array<cvec> Y(Nvec);        // received data
  Array<cmat> H(Nvec / Tc + 1);   // channel matrix (new matrix for each coherence interval)

  ivec Contflag = ones_i(Nmethods);   // flag to determine whether to run a given demodulator
  if (pow(2.0, nC*2.0*nTx) > 256) {   // ML decoder too complex..
    Contflag(1) = 0;
  }
  if (nTx > nRx) {
    Contflag(0) = 0;                  // ZF not for underdetermined systems
  }
  cout << "Running methods: " << Contflag << endl;

  cout.setf(ios::fixed, ios::floatfield);
  cout.setf(ios::showpoint);
  cout.precision(5);

  // ================== Run simulation =======================
  for (int nsnr = 0; nsnr < length(EbN0db); nsnr++) {
    const double Eb = 1.0; // transmitted energy per information bit
    const double N0 = inv_dB(-EbN0db(nsnr));
    const double sigma2 = N0; // Variance of each scalar complex noise sample
    const double Es = rate * 2 * nC * Eb; // Energy per complex scalar symbol
    // (Each transmitted scalar complex symbol contains rate*2*nC
    // information bits.)
    const double Ess = sqrt(Es);

    Array<BERC> berc(Nmethods);  // counter for coded BER
    Array<BERC> bercu(Nmethods); // counter for uncoded BER
    Array<BLERC> ferc(Nmethods); // counter for coded FER

    for (int i = 0; i < Nmethods; i++) {
      ferc(i).set_blocksize(Nu);
    }

    long int nbits = 0;
    while (nbits < Nbitsmax) {
      nbits += Nu;

      // generate and encode random data
      bvec inputbits = randb(Nu);
      bvec txbits;
      code.encode_tail(inputbits, txbits);
      // coded block length is not always a multiple of the number of
      // bits per channel vector
      txbits = concat(txbits, randb(Nctx - Nc));
      txbits = sequence_interleaver_b.interleave(txbits);

      // -- generate channel and data ----
      for (int k = 0; k < Nvec; k++) {
        /* A complex valued channel matrix is used here. An
           alternative (with equivalent result) would be to use a
           real-valued (structured) channel matrix of twice the
           dimension.
        */
        if (k % Tc == 0) {   // generate a new channel realization every Tc intervals
          H(k / Tc) = Ess * randn_c(nRx, nTx);
        }

        // modulate and transmit bits
        bvec bitstmp = txbits(k * 2 * nTx * nC, (k + 1) * 2 * nTx * nC - 1);
        cvec x = chan.modulate_bits(bitstmp);
        cvec e = sqrt(sigma2) * randn_c(nRx);
        Y(k) = H(k / Tc) * x + e;
      }

      // -- demodulate --
      Array<QLLRvec> LLRin(Nmethods);
      for (int i = 0; i < Nmethods; i++) {
        LLRin(i) = zeros_i(Nctx);
      }

      QLLRvec llr_apr = zeros_i(nC * 2 * nTx);  // no a priori input to demodulator
      QLLRvec llr_apost = zeros_i(nC * 2 * nTx);
      for (int k = 0; k < Nvec; k++) {
        // zero forcing demodulation
        if (Contflag(0)) {
          chan.demodulate_soft_bits(Y(k), H(k / Tc), sigma2, llr_apr, llr_apost,
                                    ND_UQAM::ZF_LOGMAP);
          LLRin(0).set_subvector(k*Nbitspvec, llr_apost);
        }

        // ML demodulation
        if (Contflag(1)) {
          chan.demodulate_soft_bits(Y(k), H(k / Tc), sigma2, llr_apr, llr_apost);
          LLRin(1).set_subvector(k*Nbitspvec, llr_apost);
        }
      }

      // -- decode and count errors --
      for (int i = 0; i < Nmethods; i++) {
        bvec decoded_bits;
        if (Contflag(i)) {
          bercu(i).count(txbits(0, Nc - 1), LLRin(i)(0, Nc - 1) < 0);  // uncoded BER
          LLRin(i) = sequence_interleaver_i.deinterleave(LLRin(i), 0);
          // QLLR values must be converted to real numbers since the convolutional decoder wants this
          vec llr = chan.get_llrcalc().to_double(LLRin(i).left(Nc));
          //   llr=-llr; // UNCOMMENT THIS LINE IF COMPILING WITH 3.10.5 OR EARLIER (BEFORE HARMONIZING LLR CONVENTIONS)
          code.decode_tail(llr, decoded_bits);
          berc(i).count(inputbits(0, Nu - 1), decoded_bits(0, Nu - 1));  // coded BER
          ferc(i).count(inputbits(0, Nu - 1), decoded_bits(0, Nu - 1));  // coded FER
        }
      }

      /* Check whether it is time to terminate the simulation.
       Terminate when all demodulators that are still running have
       counted at least Nbers or Nfers bit/frame errors. */
      int minber = 1000000;
      int minfer = 1000000;
      for (int i = 0; i < Nmethods; i++) {
        if (Contflag(i)) {
          minber = min(minber, round_i(berc(i).get_errors()));
          minfer = min(minfer, round_i(ferc(i).get_errors()));
        }
      }
      if (Nbers > 0 && minber > Nbers) { break;}
      if (Nfers > 0 && minfer > Nfers) { break;}
    }

    cout << "-----------------------------------------------------" << endl;
    cout << "Eb/N0: " << EbN0db(nsnr) << " dB. Simulated " << nbits << " bits." << endl;
    cout << " Uncoded BER: " << bercu(0).get_errorrate() << " (ZF);     " << bercu(1).get_errorrate() << " (ML)" << endl;
    cout << " Coded BER:   " << berc(0).get_errorrate()  << " (ZF);     " << berc(1).get_errorrate()  << " (ML)" << endl;
    cout << " Coded FER:   " << ferc(0).get_errorrate()  << " (ZF);     " << ferc(1).get_errorrate()  << " (ML)" << endl;
    cout.flush();

    /* Check wheter it is time to terminate simulation. Stop when all
    methods have reached the min BER/FER of interest. */
    int contflag = 0;
    for (int i = 0; i < Nmethods; i++) {
      if (Contflag(i)) {
        if (berc(i).get_errorrate() > BERmin)  {  contflag = 1;  }
        else { Contflag(i) = 0; }
        if (ferc(i).get_errorrate() > FERmin)  {  contflag = 1;  }
        else { Contflag(i) = 0; }
      }
    }
    if (contflag) { continue; }
    else {break; }
  }

  return 0;
}
