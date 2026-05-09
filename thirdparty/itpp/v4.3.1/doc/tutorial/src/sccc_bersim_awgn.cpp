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
  ivec gen = "07 05";//octal form
  int constraint_length = 3;
  int nb_errors_lim = 1500;
  int nb_bits_lim = int(1e6);
  int perm_len = pow2i(14);//permutation length
  int nb_iter = 10;//number of iterations in the turbo decoder
  vec EbN0_dB = "0:0.1:5";
  double R = 1.0 / 4.0;//coding rate (non punctured SCCC)
  double Ec = 1.0;//coded bit energy

  //other parameters
  int nb_bits_tail = perm_len / gen.length();
  int nb_bits = nb_bits_tail - (constraint_length - 1);//number of bits in a block (without tail)
  vec sigma2 = (0.5 * Ec / R) * pow(inv_dB(EbN0_dB), -1.0);//N0/2
  double Lc;//scaling factor for intrinsic information
  int nb_blocks;//number of blocks
  int nb_errors;
  bvec bits(nb_bits);//data bits
  bvec nsc_coded_bits;//tail is added
  bmat rsc_parity_bits;
  ivec perm(perm_len);
  ivec inv_perm(perm_len);
  int rec_len = gen.length() * perm_len;
  bvec coded_bits(rec_len);
  vec rec(rec_len);
  //SISO RSC
  vec rsc_intrinsic_coded(rec_len);
  vec rsc_apriori_data(perm_len);
  vec rsc_extrinsic_coded;
  vec rsc_extrinsic_data;
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

  //Non recursive non Systematic Convolutional Code
  Convolutional_Code nsc;
  nsc.set_generator_polynomials(gen, constraint_length);

  //Recursive Systematic Convolutional Code
  Rec_Syst_Conv_Code rsc;
  rsc.set_generator_polynomials(gen, constraint_length);//initial state should be the zero state

  //BPSK modulator
  BPSK bpsk;

  //AWGN channel
  AWGN_Channel channel;

  //SISO blocks
  SISO siso;
  siso.set_generators(gen, constraint_length);
  siso.set_map_metric(map_metric);

  //BER
  BERC berc;

  //Randomize generators
  RNG_randomize();

  //main loop
  for (en = 0;en < snr_len;en++) {
    cout << "EbN0_dB = " << EbN0_dB(en) << endl;
    channel.set_noise(sigma2(en));
    Lc = -2.0 / sigma2(en);//take into account the BPSK mapping
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

      //serial concatenated convolutional code
      nsc.encode_tail(bits, nsc_coded_bits);//tail is added here to information bits to close the trellis
      nsc_coded_bits = nsc_coded_bits(perm);//interleave
      rsc.encode(nsc_coded_bits, rsc_parity_bits);//no tail added
      for (n = 0;n < perm_len;n++)
      {
        coded_bits(2*n) = nsc_coded_bits(n);//systematic output
        coded_bits(2*n + 1) = rsc_parity_bits(n, 0);//parity output
      }

      //BPSK modulation (1->-1,0->+1) + channel
      rec = channel(bpsk.modulate_bits(coded_bits));

      //turbo decoder
      rsc_intrinsic_coded = Lc * rec;//intrinsic information of coded bits
      rsc_apriori_data.zeros();//a priori LLR for information bits
      for (n = 0;n < nb_iter;n++)
      {
        //first decoder
        siso.rsc(rsc_extrinsic_coded, rsc_extrinsic_data, rsc_intrinsic_coded, rsc_apriori_data, false);

        //deinterleave+threshold
        nsc_intrinsic_coded = SISO::threshold(rsc_extrinsic_data(inv_perm), threshold_value);

        //second decoder
        siso.nsc(nsc_extrinsic_coded, nsc_extrinsic_data, nsc_intrinsic_coded, nsc_apriori_data, true);

        //decision
        rec_bits = bpsk.demodulate_bits(-nsc_extrinsic_data);//suppose that a priori info is zero
        //count errors
        berc.clear();
        berc.count(bits, rec_bits.left(nb_bits));
        ber(n, en) += berc.get_errorrate();

        //interleave
        rsc_apriori_data = nsc_extrinsic_coded(perm);
      }//end iterations
      nb_errors += int(berc.get_errors());//get number of errors at the last iteration
      nb_blocks++;
    }//end blocks (while loop)

    //compute BER over all tx blocks
    ber.set_col(en, ber.get_col(en) / nb_blocks);
  }

  //save results to file
  it_file ff("sccc_bersim_awgn.it");
  ff << Name("BER") << ber;
  ff << Name("EbN0_dB") << EbN0_dB;
  ff.close();

  return 0;
}
