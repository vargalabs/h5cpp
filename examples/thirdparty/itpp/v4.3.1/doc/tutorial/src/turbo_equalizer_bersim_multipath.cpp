#define USE_PRECODER

#include "itpp/itcomm.h"

using namespace itpp;
using std::cout;
using std::endl;
using std::string;

int main(void)
{
  //general parameters
  double threshold_value = 50;
  string map_metric = "maxlogMAP";
  ivec gen = "07 05";//octal notation
  int constraint_length = 3;
  int ch_nb_taps = 4;//number of channel multipaths
  int nb_errors_lim = 3000;
  int nb_bits_lim = int(1e6);
  int perm_len = pow2i(14);//permutation length
  int nb_iter = 10;//number of iterations in the turbo decoder
  vec EbN0_dB = "0:0.5:10";
  double R = 1.0 / 2.0;//coding rate of FEC
  double Ec = 1.0;//coded bit energy
#ifdef USE_PRECODER
  ivec prec_gen = "03 02";//octal notation
  int prec_gen_length = 2;
#endif

  //other parameters
  int nb_bits_tail = perm_len / gen.length();
  int nb_bits = nb_bits_tail - (constraint_length - 1);//number of bits in a block (without tail)
  vec sigma2 = (0.5 * Ec / R) * pow(inv_dB(EbN0_dB), -1.0);//N0/2
  int nb_blocks;//number of blocks
  int nb_errors;
  bvec bits(nb_bits);//data bits
  bvec nsc_coded_bits(perm_len);//tail is added
  bvec em_bits(perm_len);
  bmat parity_bits;
  ivec perm(perm_len);
  ivec inv_perm(perm_len);
  vec rec(perm_len);
  //SISO equalizer
  vec eq_apriori_data(perm_len);
  vec eq_extrinsic_data;
  //SISO NSC
  vec nsc_intrinsic_coded(perm_len);
  vec nsc_apriori_data(nb_bits_tail);
  nsc_apriori_data.zeros();//always zero
  vec nsc_extrinsic_coded;
  vec nsc_extrinsic_data;
  //decision
  bvec rec_bits(nb_bits_tail);
  int snr_len = EbN0_dB.length();
  mat ber(nb_iter, snr_len);
  ber.zeros();
  register int en, n;

  //CCs
  Convolutional_Code nsc;
  nsc.set_generator_polynomials(gen, constraint_length);
#ifdef USE_PRECODER
  Rec_Syst_Conv_Code prec;
  prec.set_generator_polynomials(prec_gen, prec_gen_length);
#endif

  //BPSK
  BPSK bpsk;

  //AWGN
  AWGN_Channel awgn;

  //multipath channel impulse response (Rayleigh fading) with real coefficients
  vec ch_imp_response(ch_nb_taps);
  vec ini_state = ones(ch_nb_taps);//initial state is zero
  MA_Filter<double, double, double> multipath_channel;

  //SISO blocks
  SISO siso;
  siso.set_generators(gen, constraint_length);
  siso.set_map_metric(map_metric);
#ifdef USE_PRECODER
  siso.set_precoder_generator(prec_gen(0), prec_gen_length);
#endif

  //BER
  BERC berc;

  //Randomize generators
  RNG_randomize();

  //main loop
  for (en = 0;en < snr_len;en++) {
    cout << "EbN0_dB = " << EbN0_dB(en) << endl;
    awgn.set_noise(sigma2(en));
    siso.set_noise(sigma2(en));
    nb_errors = 0;
    nb_blocks = 0;
    while ((nb_errors < nb_errors_lim) && (nb_blocks*nb_bits < nb_bits_lim))//if at the last iteration the nb. of errors is inferior to lim, then process another block
    {
      //permutation
      perm = sort_index(randu(perm_len));
      //inverse permutation
      inv_perm = sort_index(perm);

      //bits generation
      bits = randb(nb_bits);

      //convolutional code
      nsc.encode_tail(bits, nsc_coded_bits);//tail is added here to information bits to close the trellis

      //permutation
      em_bits = nsc_coded_bits(perm);

#ifdef USE_PRECODER
      //precoder
      prec.encode(em_bits, parity_bits);
      em_bits = parity_bits.get_col(0);
#endif

      //BPSK modulation (1->-1,0->+1) + multipath channel
      ch_imp_response = randray(ch_nb_taps);
      ch_imp_response /= sqrt(sum_sqr(ch_imp_response));//normalized power profile
      multipath_channel.set_coeffs(ch_imp_response);
      multipath_channel.set_state(ini_state);//inital state is zero
      rec = awgn(multipath_channel(bpsk.modulate_bits(em_bits)));

      //turbo equalizer
      eq_apriori_data.zeros();//a priori information of emitted symbols
      siso.set_impulse_response(ch_imp_response);
      for (n = 0;n < nb_iter;n++)
      {
        //first decoder
        siso.equalizer(eq_extrinsic_data, rec, eq_apriori_data, false);//no tail
        //deinterleave+threshold
        nsc_intrinsic_coded = SISO::threshold(eq_extrinsic_data(inv_perm), threshold_value);
        //second decoder
        siso.nsc(nsc_extrinsic_coded, nsc_extrinsic_data, nsc_intrinsic_coded, nsc_apriori_data, true);//tail

        //decision
        rec_bits = bpsk.demodulate_bits(-nsc_extrinsic_data);//assume that a priori info is zero
        //count errors
        berc.clear();
        berc.count(bits, rec_bits.left(nb_bits));
        ber(n, en) += berc.get_errorrate();

        //interleave
        eq_apriori_data = nsc_extrinsic_coded(perm);
      }//end iterations
      nb_errors += int(berc.get_errors());//get number of errors at the last iteration
      nb_blocks++;
    }//end blocks (while loop)

    //compute BER over all tx blocks
    ber.set_col(en, ber.get_col(en) / nb_blocks);
  }

  //save results to file
  it_file ff("turbo_equalizer_bersim_multipath.it");
  ff << Name("BER") << ber;
  ff << Name("EbN0_dB") << EbN0_dB;
  ff.close();

  return 0;
}

