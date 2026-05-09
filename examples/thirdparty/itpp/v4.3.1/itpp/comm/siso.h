/*!
 * \file
 * \brief Definitions for Soft Input Soft Output (SISO) modules class
 * \author Bogdan Cristea
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2010  (see AUTHORS file for a list of contributors)
 *
 * This file is part of IT++ - a C++ library of mathematical, signal
 * processing, speech processing, and communications classes and functions.
 *
 * IT++ is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * IT++ is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along
 * with IT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 * -------------------------------------------------------------------------
 */

#ifndef SISO_H
#define SISO_H

#include <itpp/itbase.h> //IT++ base module
#include <itpp/itexports.h>

namespace itpp
{

/*!
  \ingroup fec
  \brief Soft Input Soft Output (%SISO) modules

  The following SISO modules are implemented:
  - decoder for an \f$1/2\f$ Recursive Systematic Convolutional (RSC) code
  - decoder for an \f$1/r\f$ Non-recursive non-Systematic Convolutional (NSC) code
  - equalizer (without and with precoding)
  - descrambler used in Interleave Division Multiple Access (IDMA) systems
  - Multi User Detectors (MUDs) for IDMA systems
  - demappers for Bit Interleaved Coded Modulation (BICM) systems
  - demappers for Space Time (ST) BICM systems

   \note BPSK mapping is realized as follows: 0 -> +1 and 1 -> -1. Thus the xor truth
   table is preserved when multiplying BPSK symbols.

   \note There is some overlap in functionality between the SISO class and the
   Modulator_ND class as follows:
   - When used for ST-BICM systems, both SISO and Modulator_ND classes can be used
   for iterative (turbo reception), but the for the SISO class the emitter part
   is implemented by the STC class, while for the Modulator_ND class the emitter
   is implemented by the same class.
   - When used for reception, the SISO class is more generic than the Modulator_ND
   class since it allows the use of several ST codes following Hassibi's model
   (see STC class) and several reception algorithms.
   - The SISO demapper for ST-BICM systems when using V-BLAST as ST code can be
   replaced by the ND demodulator (see Modulator_ND class). The best performance
   could be achieved with the ND demodulator using FULL_ENUM_LOGMAP algorithm
   followed by the SISO demapper with Hassibi_maxlogMAP algorithm (less complex
   than FULL_ENUM_LOGMAP). There is no configuration in which the ND demodulator
   and the SISO demapper can be considered completely equivalent (from an
   implementation point of view).
 */
class ITPP_EXPORT SISO
{
public:
    //! %SISO class constructor
    /*! Internal variables are initialized with default values:
     * - the trellis used in MAP algorithm is not terminated
     * - maxlogMAP metric is used
     * - simplified GCD is selected
     * - no scrambler is used
     * - GA demapper is selected, etc.
     */
    SISO();
    //! Sets the metric for MAP algorithm (convolutional codes and multipath channels)
    /*! Possible input values are:
     * - logMAP
     * - maxlogMAP
     * - SOVA
     * - Viterbi
     *
     * \note Soft Output Viterbi Algorithm (SOVA) is equivalent to the MAP algorithm
     * only when the a priori information is zero.
     */
    void set_map_metric(const std::string &in_MAP_metric);
    //! Sets the precoder generator polynomial for turbo equalizer
    /*! The precoder is used in turbo equalization.
     * The generator polynomial describes the feedback connections of the precoder.
     */
    void set_precoder_generator(const itpp::bvec &in_prec_gen);
    void set_precoder_generator(const int &in_prec_gen, const int &constraint_length);
    //! Sets convolutional code generator polynomials
    /*!
     * The generator polynomials are specified as rows of the input binary matrix.
     */
    void set_generators(const itpp::bmat &in_gen);
    void set_generators(const itpp::ivec &in_gen, const int &constraint_length);
    //! Signals whether the trellis used in the MAP algorithm is terminated or not (only for convolutional codes and multipath channels)
    /*!
     * If the input value is true, the trellis is terminated. In order to terminate
     * the trellis an ending bit tail is added to the input stream of the encoder.
     */
    void set_tail(const bool &in_tail);
    //! Sets the length of the trellis used by the SOVA
    void set_viterbi_win_len(const int &win_len);
    //! Sets the scaling factor used to normalize the reliability value computed by the SOVA
    void set_sova_scaling_factor(const double &scaling_factor);
    //! Sets the threshold value used to limit the reliability value computed by SOVA
    void set_sova_threshold(const double &threshold);
    //! Sets the Viterbi algorithm scaling factors
    void set_viterbi_scaling_factors(const double &matching_scaling_factor, //!< scaling factor for matching bits
                                     const double &nonmatching_scaling_factor //!< scaling factor for non matching bits
                                    );
    //! Sets the Viterbi algorithm hard output flag (true when only the hard output is needed)
    void set_viterbi_hard_output_flag(const bool &flag);
    //channel setup functions
    //! Sets Additive White Gaussian Noise variance for each dimension
    void set_noise(const double &in_sigma2);
    //! Sets channel impulse response for equalizer
    /*! The input is a real channel impulse response.
     */
    void set_impulse_response(const itpp::vec &h);
    //! Sets channel impulse response for Multi-User Detectors
    /*! The input is a real matrix with each row represented by the channel impulse response of one user.
     */
    void set_impulse_response(const itpp::mat &H);
    //! Sets the channel attenuations for demappers (when only a modulator is used)
    /*! The input is the vector of channel attenuations for each received symbol
     */
    void set_impulse_response(const itpp::cvec &h);
    //! Sets channel attenuations for demappers (when Space-Time codes are used)
    /*! The input is a complex matrix of dimension \f$MN\times tx\_duration\f$,
     * where \f$M\f$ and \f$N\f$ is the number of emission and reception antennas,
     * respectively and \f$tx\_duration\f$ is the transmission duration expressed
     * in symbol durations.
     *
     * This input matrix is formed as follows:
     * - the starting point is an \f$M\times N\f$ complex matrix (channel matrix)
     * with one row represented by the atenuations seen by each of \f$N\f$ reception
     * antennas during one symbol duration, when the signal is emitted by a given
     * emission antenna. The row number represents the number of emission antenna.
     * - the \f$M\times N\f$ channel matrix is then transformed into a vector of
     * length \f$MN\f$, with the first \f$M\f$ elements the first column of the channel matrix
     * - the vector of length \f$MN\f$ represents one column of the input matrix
     * - in the input matrix, the vector is repeated \f$\tau_c\f$ times and
     * \f$tx\_duration/\tau_c\f$ different vectors are used. Thus, the channel
     * is supposed constant over \f$\tau_c\f$ symbol durations (channel coherence time)
     * and \f$tx\_duration/\tau_c\f$ different channel realisations are used.
     * - in our implementation \f$\tau_c\f$ must be and integer multiple of \f$T\f$,
     * where \f$T\f$ is the ST block code duration expressed in symbol durations.
     * This means that the channel matrix must be constant over at least \f$T\f$ symbol durations.
     */
    void set_impulse_response(const itpp::cmat &cH);
    //! Sets scrambler pattern
    /*! The scrambler pattern must be a sequence of \f$\pm1\f$ and is used by
     * the %SISO NSC module in IDMA systems reception. At emission side the bits
     * are first encoded by an NSC code, BPSK modulated and then scrambled with the given pattern.
     */
    void set_scrambler_pattern(const itpp::vec &phi);
    void set_scrambler_pattern(const itpp::bvec &phi);
    //! Sets Multi-User Detector method
    /*! Possible input values are:
     * - maxlogMAP
     * - GCD
     * - sGCD (simplified GCD)
     */
    void set_mud_method(const std::string &method);
    //demodulator and MIMO demapper setup
    //! Sets symbol constellation
    /*! The correspondence between each symbol and its binary representation is given.
     * Note that if the symbols are normalized to the square root of the number of emission antenna
     * you should use the normalized constellation as input for this method.
     */
    void set_constellation(const int &in_nb_bits_symb, //!< the number of symbols
                           const itpp::cvec &in_constellation, //!< all possible symbols as a complex vector
                           const itpp::bmat &in_bin_constellation //!< binary representations of symbols as a binary matrix (each row corresponds to one symbol)
                          );
    void set_constellation(const int &in_nb_bits_symb, //!< the number of symbols
                           const itpp::cvec &in_constellation, //!< all possible symbols as a complex vector
                           const itpp::ivec &in_int_constellation //!< integer representations of symbols as a vector
                          );
    //! Sets Space-Time block code parameters
    /*! ST block codes are generated using Hassibi's model.
     */
    void set_st_block_code(const int &Q, //!< the number of symbols per block
                           const itpp::cmat &A, //!< generator matrices
                           const itpp::cmat &B, //!< of the ST block code
                           const int &N //!< the number of reception antennas
                          );
    //! Sets demapper method
    /*! Possible input values are:
     * - Hassibi_maxlogMAP (maxlogMAP algorithm applied for ST block codes represented
     * using Hassibi's model)
     * - GA
     * - sGA (simplified GA)
     * - mmsePIC
     * - zfPIC (simplified mmsePIC)
     * - Alamouti_maxlogMAP (maxlogMAP algorithm applied to Alamouti code using matched-filter
     * reception method)
     */
    void set_demapper_method(const std::string &method);
    //! %SISO decoder for RSC codes
    void rsc(itpp::vec &extrinsic_coded, //!< extrinsic information of coded bits
             itpp::vec &extrinsic_data, //!< extrinsic information of data bits
             const itpp::vec &intrinsic_coded, //!< intrinsic information of coded bits
             const itpp::vec &apriori_data //!< a priori information of data bits
            );
    //! %SISO decoder for RSC codes (tail is set through input)
    void rsc(itpp::vec &extrinsic_coded, //!< extrinsic information of coded bits
             itpp::vec &extrinsic_data, //!< extrinsic information of data bits
             const itpp::vec &intrinsic_coded, //!< intrinsic information of coded bits
             const itpp::vec &apriori_data, //!< a priori information of data bits
             const bool &tail //!< if true the trellis is terminated
            );
    //! %SISO decoder for NSC codes
    void nsc(itpp::vec &extrinsic_coded, //!< extrinsic information of coded bits
             itpp::vec &extrinsic_data, //!< extrinsic information of data bits
             const itpp::vec &intrinsic_coded, //!< intrinsic information of coded bits
             const itpp::vec &apriori_data //!< a priori information of data bits
            );
    //! %SISO decoder for NSC codes (tail is set through input)
    void nsc(itpp::vec &extrinsic_coded, //!< extrinsic information of coded bits
             itpp::vec &extrinsic_data, //!< extrinsic information of data bits
             const itpp::vec &intrinsic_coded, //!< intrinsic information of coded bits
             const itpp::vec &apriori_data, //!< a priori information of data bits
             const bool &tail //!< if true the trellis is terminated
            );
    //! %SISO equalizer
    /*! Channel trellis is generated so that BPSK mapping is assumed: 0->+1 and 1->-1
     * (xor truth table is preserved)
     */
    void equalizer(itpp::vec &extrinsic_data, //!< extrinsic informations of input symbols
                   const itpp::vec &rec_sig, //!< received signal
                   const itpp::vec &apriori_data //!< a priori informations of input symbols
                  );
    //! %SISO equalizer (tail is set through input)
    /*! Channel trellis is generated so that BPSK mapping is assumed: 0->+1 and 1->-1
     * (xor truth table is preserved)
     */
    void equalizer(itpp::vec &extrinsic_data, //!< extrinsic informations of input symbols
                   const itpp::vec &rec_sig, //!< received signal
                   const itpp::vec &apriori_data,  //!< a priori informations of input symbols
                   const bool &tail //!< if true the trellis is terminated
                  );
    //! %SISO descrambler
    void descrambler(itpp::vec &extrinsic_coded, //!< extrinsic information of scrambled bits
                     itpp::vec &extrinsic_data, //!< extrinsic information of informational bits
                     const itpp::vec &intrinsic_coded, //!< intrinsic information of scrambled bits
                     const itpp::vec &apriori_data //!< a priori information of informational bits
                    );
    //! %SISO Multi-User Detector
    void mud(itpp::mat &extrinsic_data, //!< extrinsic informations of emitted chips from all users
             const itpp::vec &rec_sig, //!< received signal
             const itpp::mat &apriori_data //!< a priori informations of emitted chips from all users
            );
    //! %SISO demapper (when only a modulator is used)
    void demapper(itpp::vec &extrinsic_data, //!< extrinsic informations of emitted bits
                  const itpp::cvec &rec_sig, //!< received signal
                  const itpp::vec &apriori_data //!< a priori informations of emitted bits
                 );
    //! %SISO demapper (when Space-Time codes are used)
    void demapper(itpp::vec &extrinsic_data, //!< extrinsic informations of emitted bits
                  const itpp::cmat &rec_sig, //!< received signal
                  const itpp::vec &apriori_data //!< a priori informations of emitted bits
                 );
    //! Functions used to limit values at a given +- threshold
    static double threshold(const double &x, const double &value);
    static itpp::vec threshold(const itpp::vec &in, const double &value);
    static itpp::mat threshold(const itpp::mat &in, const double &value);
private:
    //! SISO::rsc using logMAP algorithm
    void rsc_logMAP(itpp::vec &extrinsic_coded, itpp::vec &extrinsic_data,
                    const itpp::vec &intrinsic_coded, const itpp::vec &apriori_data);
    //! SISO::rsc using maxlogMAP algorithm
    void rsc_maxlogMAP(itpp::vec &extrinsic_coded, itpp::vec &extrinsic_data,
                       const itpp::vec &intrinsic_coded, const itpp::vec &apriori_data);
    //! SISO::rsc using %SOVA
    void rsc_sova(itpp::vec &extrinsic_data, //!< extrinsic information of data bits
                  const itpp::vec &intrinsic_coded, //!< intrinsic information of coded bits
                  const itpp::vec &apriori_data, //!< a priori information of data bits
                  const int &win_len //!< window length used to represent the trellis
                 );
    //! SISO::rsc using Viterbi algorithm
    void rsc_viterbi(itpp::vec &extrinsic_coded, //!< extrinsic information of coded bits
                     itpp::vec &extrinsic_data, //!< extrinsic information of data bits
                     const itpp::vec &intrinsic_coded, //!< intrinsic information of coded bits
                     const itpp::vec &apriori_data, //!< a priori information of data bits
                     const int &win_len //!< window length used to represent the trellis
                    );
    //! SISO::nsc using logMAP algorithm
    void nsc_logMAP(itpp::vec &extrinsic_coded, itpp::vec &extrinsic_data,
                    const itpp::vec &intrinsic_coded, const itpp::vec &apriori_data);
    //! SISO::nsc using maxlogMAP algorithm
    void nsc_maxlogMAP(itpp::vec &extrinsic_coded, itpp::vec &extrinsic_data,
                       const itpp::vec &intrinsic_coded, const itpp::vec &apriori_data);
    //! SISO::equalizer using logMAP algorithm
    void equalizer_logMAP(itpp::vec &extrinsic_data, const itpp::vec &rec_sig,
                          const itpp::vec &apriori_data);
    //! SISO::equalizer using maxlogMAP algorithm
    void equalizer_maxlogMAP(itpp::vec &extrinsic_data, const itpp::vec &rec_sig,
                             const itpp::vec &apriori_data);
    //! SISO::mud using maxlogMAP algorithm
    void mud_maxlogMAP(itpp::mat &extrinsic_data, const itpp::vec &rec_sig,
                       const itpp::mat &apriori_data);
    //! SISO::mud using maxlogMAP algorithm based on T-BCJR
    void mud_maxlogTMAP(itpp::mat &extrinsic_data, const itpp::vec &rec_sig,
                        const itpp::mat &apriori_data, const double &threshold=-5);
    //! SISO::mud using Gaussian Chip Detector (GCD)
    void GCD(itpp::mat &extrinsic_data, const itpp::vec &rec_sig,
             const itpp::mat &apriori_data);
    //! SISO::mud using simplified Gaussian Chip Detector (sGCD)
    void sGCD(itpp::mat &extrinsic_data, const itpp::vec &rec_sig,
              const itpp::mat &apriori_data);
    //! SISO::demapper using maxlogMAP algorithm for ST block codes described using Hassibi's model
    void Hassibi_maxlogMAP(itpp::vec &extrinsic_data, const itpp::cmat &rec_sig,
                           const itpp::vec &apriori_data);
    //! SISO::demapper using Gaussian Approximation (GA) algorithm
    void GA(itpp::vec &extrinsic_data, const itpp::cmat &rec_sig,
            const itpp::vec &apriori_data);
    //! SISO::demapper using simplified Gaussian Approximation (sGA) algorithm
    void sGA(itpp::vec &extrinsic_data, const itpp::cmat &rec_sig,
             const itpp::vec &apriori_data);
    //! SISO::demapper using MMSE Parallel Interference Canceller (PIC)
    void mmsePIC(itpp::vec &extrinsic_data, const itpp::cmat &rec_sig,
                 const itpp::vec &apriori_data);
    //! SISO::demapper using ZF Parallel Interference Canceller (PIC)
    void zfPIC(itpp::vec &extrinsic_data, const itpp::cmat &rec_sig,
               const itpp::vec &apriori_data);
    //! SISO::demapper using maxlogMAP algorithm and matched filter receiver for Alamouti ST code
    void Alamouti_maxlogMAP(itpp::vec &extrinsic_data, const itpp::cmat &rec_sig,
                            const itpp::vec &apriori_data);
    //! SISO::demapper using logMAP algorithm for complex modulators
    void demodulator_logMAP(itpp::vec &extrinsic_data, const itpp::cvec &rec_sig,
                            const itpp::vec &apriori_data);
    //! SISO::demapper using maxlogMAP algorithm for complex modulators
    void demodulator_maxlogMAP(itpp::vec &extrinsic_data, const itpp::cvec &rec_sig,
                               const itpp::vec &apriori_data);
    //! Prints an error message to standard output
    /*! If the %SISO class is used in a mex file, this function ensures that
     * the proper function is used for displaying the error message
     */
    void print_err_msg(const std::string &msg) const;

    // MAP algorithm variables
    //! MAP algorithm metric
    struct ITPP_EXPORT MAP_Metrics
    {
      enum Type {Unknown, logMAP, maxlogMAP, SOVA, Viterbi};
      MAP_Metrics() : _t(Unknown) {}
      MAP_Metrics(Type t) : _t(t) {}
      operator Type () const {return _t;}
    private:
      Type _t;
      template<typename T> operator T () const;
    };
    MAP_Metrics MAP_metric;
    //! Generator polynomials for convolutional codes (CC)
    itpp::bmat gen;
    //! Precoder generator polynomial
    itpp::bvec prec_gen;
    //! True if trellis of CC is terminated
    bool tail;
    // SOVA & Viterbi variables
    //! Viterbi trellis window length
    int Viterbi_win_len;
    //! SOVA scaling factor used to multiply the reliabiliy information
    double SOVA_scaling_factor;
    //! SOVA threshold used to limit the reliability information
    double SOVA_threshold;
    //! Viterbi scaling factors
    double Viterbi_scaling_factor[2];
    //! Viterbi hard output flag (true when only hard output is needed)
    bool Viterbi_hard_output_flag;
    //channel variables
    //! AWGN noise variance
    double sigma2;
    //! Real channel impulse response
    itpp::mat impulse_response;
    //! Complex channel impulse response
    itpp::cmat c_impulse_response;
    //! Scrambler pattern
    itpp::bvec scrambler_pattern;
    //! MUD method
    struct ITPP_EXPORT MUD_Methods
    {
      enum Type {Unknown, sGCD, maxlogMAP, GCD};
      MUD_Methods() : _t(Unknown) {}
      MUD_Methods(Type t) : _t(t) {}
      operator Type () const {return _t;}
    private:
      Type _t;
      template<typename T> operator T () const;
    };
    MUD_Methods MUD_method;
    //constellation variables
    //! Number of bits/symbol
    int nb_bits_symb;
    //! Complex constellation
    itpp::cvec constellation;
    //! Binary constellation
    itpp::bmat bin_constellation;
    //Space Time block code variables
    //! Number of symbols/block
    int symbols_block;
    //! Number of emission antennas
    int nb_em_ant;
    //! Number of reception antennas
    int nb_rec_ant;
    //! ST code block duration
    int block_duration;
    //! ST generator matrix 1
    itpp::cmat ST_gen1;
    //! ST generator matrix 2
    itpp::cmat ST_gen2;
    //! Demapper method
    struct ITPP_EXPORT Demapper_Methods
    {
      enum Type {Unknown, GA, Hassibi_maxlogMAP, sGA, mmsePIC, zfPIC, Alamouti_maxlogMAP};
      Demapper_Methods() : _t(Unknown) {}
      Demapper_Methods(Type t) : _t(t) {}
      operator Type () const {return _t;}
    private:
      Type _t;
      template<typename T> operator T () const;
    };
    Demapper_Methods demapper_method;

    //internal variables and functions
    //! FIR filter for a zero padded signals (\f$L\f$ zeros are added at the end of the signal, where \f$L\f$ is the order of the filter)
    void zpFIRfilter(itpp::vec& filt, //!< filtered signal
                     const itpp::vec &h, //!< filter impulse response
                     const itpp::vec &sig //!< signal to filter
                    );
    //! Generates (precoded) channel trellis
    void gen_chtrellis(void);
    //! Generates (precoded) hyper channel trellis
    void gen_hyperTrellis(void);
    //! (Hyper) Channel trellis
    struct
    {
        int numInputSymbols;//!< number of input symbols
        int stateNb;//!< number of states
        int* prevState;//!< previous states
        int* nextState;//!< next states
        double* output;//!< output
        int* input;//!< input
    } chtrellis;
    //! Generates Recursive and Systematic Convolutional (RSC) code trellis
    void gen_rsctrellis(void);
    //! RSC code trellis
    struct
    {
        int numStates;//!< number of states
        int* prevStates;//!< previous states
        int* nextStates;//!< next states
        double* PARout;//!< parity output bit
        itpp::bin* fm;//! feedback memory
    } rsctrellis;
    //! Generates Non recursive and non Systematic Convolutional (NSC) code trellis
    void gen_nsctrellis(void);
    //! NSC code trellis
    struct
    {
        int stateNb;//!< number of states
        int* prevState;//!< previous states
        int* nextState;//!< next states
        double* output;//!< output
        int* input;//!< input
    } nsctrellis;
    //! Finds half constellations
    void find_half_const(int &select_half, itpp::vec &re_part,
                         itpp::bmat &re_bin_part, itpp::vec &im_part, itpp::bmat &im_bin_part);
    //! Finds equivalent received signal with real coefficients
    void EquivRecSig(itpp::vec &x_eq, const itpp::cmat &rec_sig);
    //! Finds equivalent channel with real coefficients
    void EquivCh(itpp::mat &H_eq, const itpp::cvec &H);
    //! Computes equivalent symbols statistics (mean and variance of the real and imaginary part)
    void compute_symb_stats(itpp::vec &Es, itpp::vec &Vs,
	 		    int ns, int select_half, const itpp::vec &apriori_data,
		   	    const itpp::vec &re_part, const itpp::vec &im_part,
			    const itpp::bmat &re_bin_part, const itpp::bmat &im_bin_part);
    static MAP_Metrics map_metric_from_string(const std::string &in_MAP_metric);
    static MUD_Methods mud_method_from_string(const std::string &in_mud_method);
    static Demapper_Methods demapper_method_from_string(const std::string &in_dem_method);
};

inline SISO::SISO()
{
    tail = false;
    MAP_metric = MAP_Metrics::maxlogMAP;
    MUD_method = MUD_Methods::sGCD;
    scrambler_pattern = "0";//corresponds to +1 using BPSK mapping
    prec_gen = "1";
    demapper_method = Demapper_Methods::GA;
    Viterbi_win_len = 20;//should be set according to the generator polynomials
    SOVA_scaling_factor = 0.8;//set according to Wang [2003]
    SOVA_threshold = 10;//according to Wang [2003] an adaptive value should be used
    Viterbi_scaling_factor[0] = 1.4;//according to Kerner [2009]
    Viterbi_scaling_factor[1] = 0.4;
    Viterbi_hard_output_flag = false;
}

inline SISO::MAP_Metrics SISO::map_metric_from_string(const std::string &in_MAP_metric)
{
    if (in_MAP_metric=="logMAP")
    {
        return MAP_Metrics::logMAP;
    } else if (in_MAP_metric=="maxlogMAP")
    {
        return MAP_Metrics::maxlogMAP;
    } else if (in_MAP_metric=="SOVA")
    {
        return MAP_Metrics::SOVA;
    } else if (in_MAP_metric=="Viterbi")
    {
        return MAP_Metrics::Viterbi;
    } else
    {
        return MAP_Metrics::Unknown;
    }
}

inline SISO::MUD_Methods SISO::mud_method_from_string(const std::string &in_mud_method)
{
    if (in_mud_method=="maxlogMAP")
    {
        return MUD_Methods::maxlogMAP;
    } else if (in_mud_method=="sGCD")
    {
        return MUD_Methods::sGCD;
    } else if (in_mud_method=="GCD")
    {
        return MUD_Methods::GCD;
    } else
    {
        return MUD_Methods::Unknown;
    }
}

inline SISO::Demapper_Methods SISO::demapper_method_from_string(const std::string &in_dem_method)
{
    if (in_dem_method=="Hassibi_maxlogMAP")
    {
        return Demapper_Methods::Hassibi_maxlogMAP;
    } else if (in_dem_method=="Alamouti_maxlogMAP")
    {
        return Demapper_Methods::Alamouti_maxlogMAP;
    } else if (in_dem_method=="GA")
    {
        return Demapper_Methods::GA;
    } else if (in_dem_method=="sGA")
    {
        return Demapper_Methods::sGA;
    } else if (in_dem_method=="mmsePIC")
    {
        return Demapper_Methods::mmsePIC;
    } else if (in_dem_method=="zfPIC")
    {
        return Demapper_Methods::zfPIC;
    } else
    {
        return Demapper_Methods::Unknown;
    }
}

inline void SISO::set_map_metric(const std::string &in_MAP_metric)
{
    MAP_metric = map_metric_from_string(in_MAP_metric);
}

inline void SISO::set_precoder_generator(const itpp::bvec &in_prec_gen)//set precoder polynomial
{
    prec_gen = in_prec_gen;
}

inline void SISO::set_precoder_generator(const int &in_prec_gen,
        const int &constraint_length)//set precoder polynomial
{
    prec_gen = itpp::dec2bin(constraint_length, in_prec_gen);
}

inline void SISO::set_generators(const itpp::bmat &in_gen)
{
    gen = in_gen;
}

inline void SISO::set_generators(const itpp::ivec &in_gen,
                                 const int &constraint_length)
{
    int nb_outputs = in_gen.length();
    gen.set_size(nb_outputs, constraint_length);
    for (int n=0; n<nb_outputs; n++)
        gen.set_row(n, itpp::dec2bin(constraint_length, in_gen(n)));
}

inline void SISO::set_tail(const bool &in_tail)
{
    tail = in_tail;
}

inline void SISO::set_viterbi_win_len(const int &win_len)
{
    Viterbi_win_len = win_len;
}

inline void SISO::set_sova_scaling_factor(const double &scaling_factor)
{
    SOVA_scaling_factor = scaling_factor;
}

inline void SISO::set_sova_threshold(const double &threshold)
{
    SOVA_threshold = threshold;
}

inline void SISO::set_viterbi_scaling_factors(const double &matching_scaling_factor,
        const double &nonmatching_scaling_factor)
{
    Viterbi_scaling_factor[0] = matching_scaling_factor;
    Viterbi_scaling_factor[1] = nonmatching_scaling_factor;
}

inline void SISO::set_viterbi_hard_output_flag(const bool &flag)
{
    Viterbi_hard_output_flag = flag;
}

inline void SISO::set_noise(const double &in_sigma2)
{
    sigma2 = in_sigma2;
}

inline void SISO::set_impulse_response(const itpp::vec &h)
{
    impulse_response.set_size(1, h.length());
    impulse_response.set_row(0, h);
}

inline void SISO::set_impulse_response(const itpp::mat &H)
{
    impulse_response = H;
}

inline void SISO::set_impulse_response(const itpp::cvec &h)
{
    c_impulse_response.set_size(1, h.length());
    c_impulse_response.set_row(0, h);
}

inline void SISO::set_impulse_response(const itpp::cmat &cH)
{
    c_impulse_response = cH;
}

inline void SISO::set_scrambler_pattern(const itpp::vec &phi)
{
    int phi_len = phi.length();
    scrambler_pattern.set_size(phi_len);
    //scrambler_pattern = to_bvec((1-phi)/2);//BPSK mapping: 0->+1 and 1->-1
    register int n;
    for (n=0; n<phi_len; n++)
        scrambler_pattern(n) = itpp::bin((1-int(phi(n)))/2);//BPSK mapping: 0->+1 and 1->-1
}

inline void SISO::set_scrambler_pattern(const itpp::bvec &phi)
{
    scrambler_pattern = phi;
}

inline void SISO::set_mud_method(const std::string &method)
{
    MUD_method = mud_method_from_string(method);
}

inline void SISO::set_constellation(const int &in_nb_bits_symb,
                                    const itpp::cvec &in_constellation, const itpp::bmat &in_bin_constellation)
{
    nb_bits_symb = in_nb_bits_symb;
    constellation = in_constellation;
    bin_constellation = in_bin_constellation;
}

inline void SISO::set_constellation(const int &in_nb_bits_symb,
                                    const itpp::cvec &in_constellation, const itpp::ivec &in_int_constellation)
{
    nb_bits_symb = in_nb_bits_symb;
    int nb_symb = in_constellation.length();
    constellation.set_size(nb_symb);
    bin_constellation.set_size(nb_symb, nb_bits_symb);
    for (int n=0; n<nb_symb; n++)
    {
        constellation(n) = in_constellation(in_int_constellation(n));
        bin_constellation.set_row(n, itpp::dec2bin(nb_bits_symb, n));
    }
}

inline void SISO::set_st_block_code(const int &Q, const itpp::cmat &A,
                                    const itpp::cmat &B, const int &N)
{
    symbols_block = Q;
    nb_em_ant = A.cols();
    nb_rec_ant = N;
    block_duration = A.rows()/Q;
    ST_gen1 = A;
    ST_gen2 = B;
}

inline void SISO::set_demapper_method(const std::string &method)
{
    demapper_method = demapper_method_from_string(method);
}

inline void SISO::rsc(itpp::vec &extrinsic_coded, itpp::vec &extrinsic_data,
                      const itpp::vec &intrinsic_coded, const itpp::vec &apriori_data, const bool &tail)
{
    set_tail(tail);
    rsc(extrinsic_coded, extrinsic_data, intrinsic_coded, apriori_data);
}

inline void SISO::rsc(itpp::vec &extrinsic_coded, itpp::vec &extrinsic_data,
                      const itpp::vec &intrinsic_coded, const itpp::vec &apriori_data)
{
    if (gen.size()==0)
    {
        print_err_msg("SISO::rsc: generator polynomials not initialized");
        return;
    }

    if (MAP_metric==MAP_Metrics::logMAP)
    {
        rsc_logMAP(extrinsic_coded, extrinsic_data, intrinsic_coded, apriori_data);
    } else if (MAP_metric==MAP_Metrics::maxlogMAP)
    {
        rsc_maxlogMAP(extrinsic_coded, extrinsic_data, intrinsic_coded, apriori_data);
    } else if (MAP_metric==MAP_Metrics::SOVA)
    {
        //no extrinsic information for coded bits is provided
        rsc_sova(extrinsic_data, intrinsic_coded, apriori_data, Viterbi_win_len);
    } else if (MAP_metric==MAP_Metrics::Viterbi)
    {
        rsc_viterbi(extrinsic_coded, extrinsic_data, intrinsic_coded, apriori_data, Viterbi_win_len);
    } else
    {
        print_err_msg("SISO::rsc: unknown MAP metric. The MAP metric should be either logMAP or maxlogMAP or SOVA or Viterbi");
    }
}

inline void SISO::nsc(itpp::vec &extrinsic_coded, itpp::vec &extrinsic_data,
                      const itpp::vec &intrinsic_coded, const itpp::vec &apriori_data, const bool &tail)
{
    set_tail(tail);
    nsc(extrinsic_coded, extrinsic_data, intrinsic_coded, apriori_data);
}

inline void SISO::nsc(itpp::vec &extrinsic_coded, itpp::vec &extrinsic_data,
                      const itpp::vec &intrinsic_coded, const itpp::vec &apriori_data)
{
    if (gen.size()==0)
    {
        print_err_msg("SISO::nsc: generator polynomials not initialized");
        return;
    }

    if (MAP_metric==MAP_Metrics::logMAP)
        nsc_logMAP(extrinsic_coded, extrinsic_data, intrinsic_coded, apriori_data);
    else if (MAP_metric==MAP_Metrics::maxlogMAP)
        nsc_maxlogMAP(extrinsic_coded, extrinsic_data, intrinsic_coded, apriori_data);
    else
        print_err_msg("SISO::nsc: unknown MAP metric. The MAP metric should be either logMAP or maxlogMAP");
}

inline void SISO::equalizer(itpp::vec &extrinsic_data, //!< extrinsic informations of input symbols
                            const itpp::vec &rec_sig, //!< received signal
                            const itpp::vec &apriori_data,  //!< a priori informations of input symbols
                            const bool &tail //!< if true the trellis is terminated
                           )
{
    set_tail(tail);
    equalizer(extrinsic_data, rec_sig, apriori_data);
}

inline void SISO::equalizer(itpp::vec &extrinsic_data, const itpp::vec &rec_sig,
                            const itpp::vec &apriori_data)
{
    if (impulse_response.size()==0)
    {
        print_err_msg("SISO::equalizer: channel impulse response not initialized");
        return;
    }
    if ((impulse_response.size()==1)&&(prec_gen.length()==1))
    {
        print_err_msg("SISO::equalizer: flat fading channel and no precoder. Use the soft output of the channel (no need for a priori information)");
        return;
    }

    if (MAP_metric==MAP_Metrics::logMAP)
        equalizer_logMAP(extrinsic_data, rec_sig, apriori_data);
    else if (MAP_metric==MAP_Metrics::maxlogMAP)
        equalizer_maxlogMAP(extrinsic_data, rec_sig, apriori_data);
    else
        print_err_msg("SISO::equalizer: unknown MAP metric. The MAP metric should be either logMAP or maxlogMAP");
}

inline void SISO::mud(itpp::mat &extrinsic_data, const itpp::vec &rec_sig,
                      const itpp::mat &apriori_data)
{
    if (impulse_response.size()==0)
    {
        print_err_msg("SISO::mud: channel impulse response not initialized");
        return;
    }
    if (impulse_response.rows()!=apriori_data.rows())
    {
        print_err_msg("SISO::mud: channel impulse response must have the same number of rows as a priori info.");
        return;
    }

    if (MUD_method==MUD_Methods::maxlogMAP)
        mud_maxlogMAP(extrinsic_data, rec_sig, apriori_data);
    else if (MUD_method==MUD_Methods::GCD)
        GCD(extrinsic_data, rec_sig, apriori_data);
    else if (MUD_method==MUD_Methods::sGCD)
        sGCD(extrinsic_data, rec_sig, apriori_data);
    else
        print_err_msg("SISO::mud: unknown MUD method. The MUD method should be either maxlogMAP, GCD or sGCD");
}

inline void SISO::demapper(itpp::vec &extrinsic_data, const itpp::cvec &rec_sig,
                           const itpp::vec &apriori_data)
{
    if (c_impulse_response.size()==0)
    {
        print_err_msg("SISO::demapper: channel impulse response not initialized");
        return;
    }
    if ((constellation.size()==0) || (bin_constellation.size()==0))
    {
        print_err_msg("SISO::demapper: constellation not initialized");
        return;
    }
    if (MAP_metric==MAP_Metrics::logMAP)
        demodulator_logMAP(extrinsic_data, rec_sig, apriori_data);
    else if (MAP_metric==MAP_Metrics::maxlogMAP)
        demodulator_maxlogMAP(extrinsic_data, rec_sig, apriori_data);
    else
        print_err_msg("SISO::demapper: unknown MAP metric. The MAP metric should be either logMAP or maxlogMAP");
}

inline void SISO::demapper(itpp::vec &extrinsic_data, const itpp::cmat &rec_sig,
                           const itpp::vec &apriori_data)
{
    if (c_impulse_response.size()==0)
    {
        print_err_msg("SISO::demapper: channel impulse response not initialized");
        return;
    }
    if ((ST_gen1.size()==0) || (ST_gen2.size()==0))
    {
        print_err_msg("SISO::demapper: Space-Time generator polynomials not initialized");
        return;
    }
    if ((constellation.size()==0) || (bin_constellation.size()==0))
    {
        print_err_msg("SISO::demapper: constellation not initialized");
        return;
    }

    if (demapper_method==Demapper_Methods::Hassibi_maxlogMAP)
        Hassibi_maxlogMAP(extrinsic_data, rec_sig, apriori_data);
    else if (demapper_method==Demapper_Methods::GA)
        GA(extrinsic_data, rec_sig, apriori_data);
    else if (demapper_method==Demapper_Methods::sGA)
        sGA(extrinsic_data, rec_sig, apriori_data);
    else if (demapper_method==Demapper_Methods::mmsePIC)
        mmsePIC(extrinsic_data, rec_sig, apriori_data);
    else if (demapper_method==Demapper_Methods::zfPIC)
        zfPIC(extrinsic_data, rec_sig, apriori_data);
    else if (demapper_method==Demapper_Methods::Alamouti_maxlogMAP)
        Alamouti_maxlogMAP(extrinsic_data, rec_sig, apriori_data);
    else
        print_err_msg("SISO::demapper: unknown demapper method. The demapper method should be either Hassibi_maxlogMAP, GA, sGA, mmsePIC, zfPIC or Alamouti_maxlogMAP");
}

inline void SISO::print_err_msg(const std::string &msg) const
{
#ifdef mex_h
    mexErrMsgTxt(msg.c_str());
#else
    std::cout << msg << std::endl;
#endif
}

inline double SISO::threshold(const double &x, const double &value)
{
    if ((x>value)||(x<-value))
        return (x>0?value:-value);
    return x;
}

inline itpp::vec SISO::threshold(const itpp::vec &in, const double &value)
{
    itpp::vec out(in.length());
    register int n;
    for (n=0; n<in.length(); n++)
        out(n) = threshold(in(n), value);
    return out;
}

inline itpp::mat SISO::threshold(const itpp::mat &in, const double &value)
{
    itpp::mat out(in.rows(),in.cols());
    register int n;
    for (n=0; n<in.rows(); n++)
        out.set_row(n, threshold(in.get_row(n), value));
    return out;
}

}

#endif /*SISO_H_*/
