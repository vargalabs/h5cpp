#include <itpp/itcomm.h>

using namespace itpp;

int main()
{
  // Declare my_channel variable as an instance of the Rayleigh_Channel
  // class
  TDL_Channel my_channel;

  // The normalized Doppler frequency is set to 0.1
  double norm_dopp = 0.1;
  my_channel.set_norm_doppler(norm_dopp);

  // Generate nrof_samples of the fading process and store them in ch_coeffs
  // matrix
  int nrof_samples = 10000;
  cmat ch_coeffs;
  my_channel.generate(nrof_samples, ch_coeffs);

  // Open an output file "rayleigh_test.it"
  it_file ff("rayleigh_test.it");

  // Save channel coefficients to the output file
  ff << Name("ch_coeffs") << ch_coeffs;

  // Close the output file
  ff.close();

  // Exit program
  return 0;
}
