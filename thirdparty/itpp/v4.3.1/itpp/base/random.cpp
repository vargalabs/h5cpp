/*!
 * \file
 * \brief Implementation of classes for random number generators
 * \author Tony Ottosson and Adam Piatyszek
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

#include <itpp/base/random.h>
#include <itpp/base/itcompat.h>
#include <itpp/base/math/elem_math.h>


namespace itpp
{

namespace random_details
{

//Thread-local context for thread-safe RNGs
static ActiveDSFMT::Context thread_local_context;
#pragma omp threadprivate(thread_local_context)
//Thread-local context initialization flag
static bool is_thread_local_context_initialized = false;
#pragma omp threadprivate(is_thread_local_context_initialized)

ActiveDSFMT::Context& lc_get()
{
  return thread_local_context;
}

bool lc_is_initialized()
{
  return is_thread_local_context_initialized;
}

void lc_mark_initialized()
{
  is_thread_local_context_initialized = true;
}

/*!
* \brief Get an unsigned int from time variables t and c.
*
* Better than uint(x) in case x is floating point in [0,1]
* Based on code by Lawrence Kirby (fred@genesis.demon.co.uk)
*/
static unsigned int hash_time_to_seed(time_t t, clock_t c)
{
  static unsigned int differ = 0; // guarantee time-based seeds will change

  unsigned int h1 = 0;
  unsigned char *p = (unsigned char *) &t;
  for(size_t i = 0; i < sizeof(t); ++i) {
    h1 *= std::numeric_limits<unsigned char>::max() + 2U;
    h1 += p[i];
  }
  unsigned int h2 = 0;
  p = (unsigned char *) &c;
  for(size_t j = 0; j < sizeof(c); ++j) {
    h2 *= std::numeric_limits<unsigned char>::max() + 2U;
    h2 += p[j];
  }
  return (h1 + differ++) ^ h2;
}



// ----------------------------------------------------------------------
// ActiveDSFMT (DSFMT_19937_RNG)
// ----------------------------------------------------------------------
template <>
const bool ActiveDSFMT::bigendian = is_bigendian();

#if defined(__SSE2__)
template <>
const __m128i ActiveDSFMT::sse2_param_mask = _mm_set_epi32(ActiveDSFMT::MSK32_3, ActiveDSFMT::MSK32_4, ActiveDSFMT::MSK32_1, ActiveDSFMT::MSK32_2);
#endif // __SSE2__
}

/*
 *Global Seed Provider class definition.
 *
 * Provides unique seeds for thread-safe generators running in each thread.
 */
class GlobalSeedProvider
{
  static const unsigned int default_first_seed = 4257U;
  typedef random_details::ActiveDSFMT DSFMT;
public:
  //constructor
  GlobalSeedProvider(): _dsfmt(_c), _first_seed_given(false) {
    _dsfmt.init_gen_rand(default_first_seed);
  }
  //set new seed
  void set_seed(unsigned int s) {_dsfmt.init_gen_rand(s);}
  //get preset seed
  unsigned int get_seed() {return _c.last_seed;}
  //reset provider state to previously set one
  void reset() { if(_first_seed_given) _dsfmt.init_gen_rand(get_seed());}
  //set previously saved state from ivec
  void set_state(const ivec& st) {
    int size = (DSFMT::N + 1) * 4;
    it_assert(st.size() == size + 1, "GlobalSeedProvider::state(): "
              "Invalid state initialization vector");
    uint32_t *psfmt = &_c.status[0].u32[0];
    for(int i = 0; i < size; ++i) {
      psfmt[i] = static_cast<uint32_t>(st(i));
    }
    _c.idx = st(size);
    _first_seed_given = true;
  }
  //get current provider state
  ivec get_state() {
    int size = (DSFMT::N + 1) * 4;
    uint32_t *psfmt = &_c.status[0].u32[0];
    ivec state(size + 1); // size + 1 to save idx variable in the same vec
    for(int i = 0; i < size; ++i) {
      state(i) = static_cast<int>(psfmt[i]);
    }
    state(size) = _c.idx;
    return state;
  }
  //randomize current provider state with system time
  void randomize() {_dsfmt.init_gen_rand(random_details::hash_time_to_seed(time(0), clock())); _first_seed_given = true;}
  //generate new seed for random number generators
  unsigned int generate() {
    if(_first_seed_given)
      return _dsfmt.genrand_uint32();
    else {
      //return default seed on first request.
      //it is done in order not to breake the old-style itpp tests.
      //Some tests rely on the default state equal to 4257U
      _first_seed_given = true;
      return default_first_seed;
    }
  }
private:
  //
  DSFMT _dsfmt;
  //
  DSFMT::Context _c;
  //
  bool _first_seed_given;
};


//Global seed provider instance
GlobalSeedProvider& global_seed_provider()
{
  static GlobalSeedProvider global_seed_provider_instance;
  return global_seed_provider_instance;
}


void GlobalRNG_reset(unsigned int seed)
{
  #pragma omp critical
  {
    global_seed_provider().set_seed(seed);
  }
}


void GlobalRNG_reset()
{
  #pragma omp critical
  {
    global_seed_provider().reset();
  }
}

unsigned int GlobalRNG_get_local_seed()
{
  unsigned int s;
  #pragma omp critical
  {
    s = global_seed_provider().generate();
  }
  return s;
}


void GlobalRNG_randomize()
{
  #pragma omp critical
  {
    global_seed_provider().randomize();
  }
}


void GlobalRNG_get_state(ivec &state)
{
  #pragma omp critical
  {
    state = global_seed_provider().get_state();
  }
}


void GlobalRNG_set_state(const ivec &state)
{
  #pragma omp critical
  {
    global_seed_provider().set_state(state);
  }
}

void RNG_reset(unsigned int seed)
{
  random_details::ActiveDSFMT dsfmt(random_details::lc_get());
  dsfmt.init_gen_rand(seed);
  random_details::lc_mark_initialized();
}


void RNG_reset()
{
  random_details::ActiveDSFMT dsfmt(random_details::lc_get());

  if(random_details::lc_is_initialized()) {
    //already initialized. Reinit with last set seed;
    dsfmt.init_gen_rand(random_details::lc_get().last_seed);
  }
  else {
    //query global seed provider for new seed and init with it
    dsfmt.init_gen_rand(GlobalRNG_get_local_seed());
    random_details::lc_mark_initialized();
  }
}


void RNG_randomize()
{
  random_details::ActiveDSFMT dsfmt(random_details::lc_get());
  dsfmt.init_gen_rand(random_details::hash_time_to_seed(time(0), clock()));
  random_details::lc_mark_initialized();
}


void RNG_get_state(ivec &state)
{
  int size = (random_details::ActiveDSFMT::N + 1) * 4;
  uint32_t *psfmt = &random_details::lc_get().status[0].u32[0];
  state.set_size(size + 1); // size + 1 to save idx variable in the same vec
  for(int i = 0; i < size; ++i) {
    state(i) = static_cast<int>(psfmt[i]);
  }
  state(size) = random_details::lc_get().idx;
}


void RNG_set_state(const ivec &state)
{
  int size = (random_details::ActiveDSFMT::N + 1) * 4;
  it_assert(state.size() == size + 1, "RNG_set_state: "
            "Invalid state initialization vector");
  uint32_t *psfmt = &random_details::lc_get().status[0].u32[0];
  for(int i = 0; i < size; ++i) {
    psfmt[i] = static_cast<uint32_t>(state(i));
  }
  random_details::lc_get().idx = state(size);
}

///////////////////////////////////////////////
// I_Uniform_RNG
///////////////////////////////////////////////

I_Uniform_RNG::I_Uniform_RNG(int min, int max)
{
  setup(min, max);
}

void I_Uniform_RNG::setup(int min, int max)
{
  if(min <= max) {
    lo = min;
    hi = max;
  }
  else {
    lo = max;
    hi = min;
  }
}

void I_Uniform_RNG::get_setup(int &min, int &max) const
{
  min = lo;
  max = hi;
}

ivec I_Uniform_RNG::operator()(int n)
{
  ivec vv(n);

  for(int i = 0; i < n; i++)
    vv(i) = sample();

  return vv;
}

imat I_Uniform_RNG::operator()(int h, int w)
{
  imat mm(h, w);
  int i, j;

  for(i = 0; i < h; i++)
    for(j = 0; j < w; j++)
      mm(i, j) = sample();

  return mm;
}

///////////////////////////////////////////////
// Uniform_RNG
///////////////////////////////////////////////

Uniform_RNG::Uniform_RNG(double min, double max)
{
  setup(min, max);
}

void Uniform_RNG::setup(double min, double max)
{
  if(min <= max) {
    lo_bound = min;
    hi_bound = max;
  }
  else {
    lo_bound = max;
    hi_bound = min;
  }
}

void Uniform_RNG::get_setup(double &min, double &max) const
{
  min = lo_bound;
  max = hi_bound;
}

///////////////////////////////////////////////
// Exp_RNG
///////////////////////////////////////////////

Exponential_RNG::Exponential_RNG(double lambda)
{
  setup(lambda);
}

vec Exponential_RNG::operator()(int n)
{
  vec vv(n);

  for(int i = 0; i < n; i++)
    vv(i) = sample();

  return vv;
}

mat Exponential_RNG::operator()(int h, int w)
{
  mat mm(h, w);
  int i, j;

  for(i = 0; i < h; i++)
    for(j = 0; j < w; j++)
      mm(i, j) = sample();

  return mm;
}

///////////////////////////////////////////////
// Gamma_RNG
///////////////////////////////////////////////
void Gamma_RNG::init_state()
{
  const static double sqrt32 = 5.656854;

  const static double q1 = 0.04166669;
  const static double q2 = 0.02083148;
  const static double q3 = 0.00801191;
  const static double q4 = 0.00144121;
  const static double q5 = -7.388e-5;
  const static double q6 = 2.4511e-4;
  const static double q7 = 2.424e-4;

  double r = 1.0 / alpha;
  _scale = 1.0 / beta;

  it_error_if(!std::isfinite(alpha) || !std::isfinite(_scale) || (alpha < 0.0)
              || (_scale <= 0.0), "Gamma_RNG::init_state() - wrong parameters");

  _s2 = alpha - 0.5;
  _s = std::sqrt(_s2);
  _d = sqrt32 - _s * 12.0;

  _q0 = ((((((q7 * r + q6) * r + q5) * r + q4) * r + q3) * r
          + q2) * r + q1) * r;

  /* Approximation depending on size of parameter alpha */
  /* The constants in the expressions for _b, _si and _c */
  /* were established by numerical experiments */
  if(alpha <= 3.686) {
    _b = 0.463 + _s + 0.178 * _s2;
    _si = 1.235;
    _c = 0.195 / _s - 0.079 + 0.16 * _s;
  }
  else if(alpha <= 13.022) {
    _b = 1.654 + 0.0076 * _s2;
    _si = 1.68 / _s + 0.275;
    _c = 0.062 / _s + 0.024;
  }
  else {
    _b = 1.77;
    _si = 0.75;
    _c = 0.1515 / _s;
  }

}
vec Gamma_RNG::operator()(int n)
{
  vec vv(n);
  for(int i = 0; i < n; i++)
    vv(i) = sample();
  return vv;
}

mat Gamma_RNG::operator()(int r, int c)
{
  mat mm(r, c);
  for(int i = 0; i < r * c; i++)
    mm(i) = sample();
  return mm;
}

double Gamma_RNG::sample()
{
  // A copy of rgamma code from the R package, adapted to IT++ by Vasek
  // Smidl

  /* Constants : */
  const static double exp_m1 = 0.36787944117144232159;/* exp(-1) = 1/e */

  const static double a1 = 0.3333333;
  const static double a2 = -0.250003;
  const static double a3 = 0.2000062;
  const static double a4 = -0.1662921;
  const static double a5 = 0.1423657;
  const static double a6 = -0.1367177;
  const static double a7 = 0.1233795;

  double e, p, q, t, u, v, w, x, ret_val;
  double a = alpha;
  double scale = _scale;

  if(a < 1.) {  /* GS algorithm for parameters a < 1 */
    if(a == 0)
      return 0.;
    e = 1.0 + exp_m1 * a;
    for(;;) {  //VS repeat
      p = e * RNG.genrand_open_open();
      if(p >= 1.0) {
        x = -std::log((e - p) / a);
        if(-std::log(RNG.genrand_open_close()) >= (1.0 - a) * std::log(x))
          break;
      }
      else {
        x = std::exp(std::log(p) / a);
        if(-std::log(RNG.genrand_open_close()) >= x)
          break;
      }
    }
    return scale * x;
  }

  /* --- a >= 1 : GD algorithm --- */

  /* Step 1: t = standard normal deviate, x = (s,1/2) -normal deviate. */
  /* immediate acceptance (i) */
  t = NRNG.sample();
  x = _s + 0.5 * t;
  ret_val = x * x;
  if(t >= 0.0)
    return scale * ret_val;

  /* Step 2: u = 0,1 - uniform sample. squeeze acceptance (s) */
  u = RNG.genrand_close_open();
  if((_d * u) <= (t * t * t))
    return scale * ret_val;


  /* Step 3: no quotient test if x not positive */
  if(x > 0.0) {
    /* Step 4: calculation of v and quotient q */
    v = t / (_s + _s);
    if(std::fabs(v) <= 0.25)
      q = _q0 + 0.5 * t * t * ((((((a7 * v + a6) * v + a5) * v + a4) * v
                                 + a3) * v + a2) * v + a1) * v;
    else
      q = _q0 - _s * t + 0.25 * t * t + (_s2 + _s2) * log(1.0 + v);

    /* Step 5: quotient acceptance (q) */
    if(log(1.0 - u) <= q)
      return scale * ret_val;
  }

  for(;;) {  //VS repeat
    /* Step 6: e = standard exponential deviate
     *         u =  0,1 -uniform deviate
     *         t = (b,si)-double exponential (laplace) sample */
    e = -std::log(RNG.genrand_open_close()); //see Exponential_RNG
    u = RNG.genrand_open_close();
    u = u + u - 1.0;
    if(u < 0.0)
      t = _b - _si * e;
    else
      t = _b + _si * e;
    /* Step 7: rejection if t < tau(1) = -0.71874483771719 */
    if(t >= -0.71874483771719) {
      /* Step 8:  calculation of v and quotient q */
      v = t / (_s + _s);
      if(std::fabs(v) <= 0.25)
        q = _q0 + 0.5 * t * t *
            ((((((a7 * v + a6) * v + a5) * v + a4) * v + a3) * v
              + a2) * v + a1) * v;
      else
        q = _q0 - _s * t + 0.25 * t * t + (_s2 + _s2) * log(1.0 + v);
      /* Step 9: hat acceptance (h) */
      /* (if q not positive go to step 6) */
      if(q > 0.0) {
        // Try to use w = expm1(q); (Not supported on w32)
        w = expm1(q);
        /*  ^^^^^ original code had approximation with rel.err < 2e-7 */
        /* if t is rejected sample again at step 6 */
        if((_c * std::fabs(u)) <= (w * std::exp(e - 0.5 * t * t)))
          break;
      }
    }
  } /* repeat .. until `t' is accepted */
  x = _s + 0.5 * t;
  return scale * x * x;
}

///////////////////////////////////////////////
// Normal_RNG
///////////////////////////////////////////////

void Normal_RNG::get_setup(double &meanval, double &variance) const
{
  meanval = mean;
  variance = sigma * sigma;
}

// (Ziggurat) tabulated values for the heigt of the Ziggurat levels
const double Normal_RNG::ytab[128] = {
  1, 0.963598623011, 0.936280813353, 0.913041104253,
  0.892278506696, 0.873239356919, 0.855496407634, 0.838778928349,
  0.822902083699, 0.807732738234, 0.793171045519, 0.779139726505,
  0.765577436082, 0.752434456248, 0.739669787677, 0.727249120285,
  0.715143377413, 0.703327646455, 0.691780377035, 0.68048276891,
  0.669418297233, 0.65857233912, 0.647931876189, 0.637485254896,
  0.62722199145, 0.617132611532, 0.607208517467, 0.597441877296,
  0.587825531465, 0.578352913803, 0.569017984198, 0.559815170911,
  0.550739320877, 0.541785656682, 0.532949739145, 0.524227434628,
  0.515614886373, 0.507108489253, 0.498704867478, 0.490400854812,
  0.482193476986, 0.47407993601, 0.466057596125, 0.458123971214,
  0.450276713467, 0.442513603171, 0.434832539473, 0.427231532022,
  0.419708693379, 0.41226223212, 0.404890446548, 0.397591718955,
  0.390364510382, 0.383207355816, 0.376118859788, 0.369097692334,
  0.362142585282, 0.355252328834, 0.348425768415, 0.341661801776,
  0.334959376311, 0.328317486588, 0.321735172063, 0.31521151497,
  0.308745638367, 0.302336704338, 0.29598391232, 0.289686497571,
  0.283443729739, 0.27725491156, 0.271119377649, 0.265036493387,
  0.259005653912, 0.253026283183, 0.247097833139, 0.241219782932,
  0.235391638239, 0.229612930649, 0.223883217122, 0.218202079518,
  0.212569124201, 0.206983981709, 0.201446306496, 0.195955776745,
  0.190512094256, 0.185114984406, 0.179764196185, 0.174459502324,
  0.169200699492, 0.1639876086, 0.158820075195, 0.153697969964,
  0.148621189348, 0.143589656295, 0.138603321143, 0.133662162669,
  0.128766189309, 0.123915440582, 0.119109988745, 0.114349940703,
  0.10963544023, 0.104966670533, 0.100343857232, 0.0957672718266,
  0.0912372357329, 0.0867541250127, 0.082318375932, 0.0779304915295,
  0.0735910494266, 0.0693007111742, 0.065060233529, 0.0608704821745,
  0.056732448584, 0.05264727098, 0.0486162607163, 0.0446409359769,
  0.0407230655415, 0.0368647267386, 0.0330683839378, 0.0293369977411,
  0.0256741818288, 0.0220844372634, 0.0185735200577, 0.0151490552854,
  0.0118216532614, 0.00860719483079, 0.00553245272614, 0.00265435214565
};

/*
 * (Ziggurat) tabulated values for 2^24 times x[i]/x[i+1], used to accept
 * for U*x[i+1]<=x[i] without any floating point operations
 */
const unsigned int Normal_RNG::ktab[128] = {
  0, 12590644, 14272653, 14988939,
  15384584, 15635009, 15807561, 15933577,
  16029594, 16105155, 16166147, 16216399,
  16258508, 16294295, 16325078, 16351831,
  16375291, 16396026, 16414479, 16431002,
  16445880, 16459343, 16471578, 16482744,
  16492970, 16502368, 16511031, 16519039,
  16526459, 16533352, 16539769, 16545755,
  16551348, 16556584, 16561493, 16566101,
  16570433, 16574511, 16578353, 16581977,
  16585398, 16588629, 16591685, 16594575,
  16597311, 16599901, 16602354, 16604679,
  16606881, 16608968, 16610945, 16612818,
  16614592, 16616272, 16617861, 16619363,
  16620782, 16622121, 16623383, 16624570,
  16625685, 16626730, 16627708, 16628619,
  16629465, 16630248, 16630969, 16631628,
  16632228, 16632768, 16633248, 16633671,
  16634034, 16634340, 16634586, 16634774,
  16634903, 16634972, 16634980, 16634926,
  16634810, 16634628, 16634381, 16634066,
  16633680, 16633222, 16632688, 16632075,
  16631380, 16630598, 16629726, 16628757,
  16627686, 16626507, 16625212, 16623794,
  16622243, 16620548, 16618698, 16616679,
  16614476, 16612071, 16609444, 16606571,
  16603425, 16599973, 16596178, 16591995,
  16587369, 16582237, 16576520, 16570120,
  16562917, 16554758, 16545450, 16534739,
  16522287, 16507638, 16490152, 16468907,
  16442518, 16408804, 16364095, 16301683,
  16207738, 16047994, 15704248, 15472926
};

// (Ziggurat) tabulated values of 2^{-24}*x[i]
const double Normal_RNG::wtab[128] = {
  1.62318314817e-08, 2.16291505214e-08, 2.54246305087e-08, 2.84579525938e-08,
  3.10340022482e-08, 3.33011726243e-08, 3.53439060345e-08, 3.72152672658e-08,
  3.8950989572e-08, 4.05763964764e-08, 4.21101548915e-08, 4.35664624904e-08,
  4.49563968336e-08, 4.62887864029e-08, 4.75707945735e-08, 4.88083237257e-08,
  5.00063025384e-08, 5.11688950428e-08, 5.22996558616e-08, 5.34016475624e-08,
  5.44775307871e-08, 5.55296344581e-08, 5.65600111659e-08, 5.75704813695e-08,
  5.85626690412e-08, 5.95380306862e-08, 6.04978791776e-08, 6.14434034901e-08,
  6.23756851626e-08, 6.32957121259e-08, 6.42043903937e-08, 6.51025540077e-08,
  6.59909735447e-08, 6.68703634341e-08, 6.77413882848e-08, 6.8604668381e-08,
  6.94607844804e-08, 7.03102820203e-08, 7.11536748229e-08, 7.1991448372e-08,
  7.2824062723e-08, 7.36519550992e-08, 7.44755422158e-08, 7.52952223703e-08,
  7.61113773308e-08, 7.69243740467e-08, 7.77345662086e-08, 7.85422956743e-08,
  7.93478937793e-08, 8.01516825471e-08, 8.09539758128e-08, 8.17550802699e-08,
  8.25552964535e-08, 8.33549196661e-08, 8.41542408569e-08, 8.49535474601e-08,
  8.57531242006e-08, 8.65532538723e-08, 8.73542180955e-08, 8.8156298059e-08,
  8.89597752521e-08, 8.97649321908e-08, 9.05720531451e-08, 9.138142487e-08,
  9.21933373471e-08, 9.30080845407e-08, 9.38259651738e-08, 9.46472835298e-08,
  9.54723502847e-08, 9.63014833769e-08, 9.71350089201e-08, 9.79732621669e-08,
  9.88165885297e-08, 9.96653446693e-08, 1.00519899658e-07, 1.0138063623e-07,
  1.02247952126e-07, 1.03122261554e-07, 1.04003996769e-07, 1.04893609795e-07,
  1.05791574313e-07, 1.06698387725e-07, 1.07614573423e-07, 1.08540683296e-07,
  1.09477300508e-07, 1.1042504257e-07, 1.11384564771e-07, 1.12356564007e-07,
  1.13341783071e-07, 1.14341015475e-07, 1.15355110887e-07, 1.16384981291e-07,
  1.17431607977e-07, 1.18496049514e-07, 1.19579450872e-07, 1.20683053909e-07,
  1.21808209468e-07, 1.2295639141e-07, 1.24129212952e-07, 1.25328445797e-07,
  1.26556042658e-07, 1.27814163916e-07, 1.29105209375e-07, 1.30431856341e-07,
  1.31797105598e-07, 1.3320433736e-07, 1.34657379914e-07, 1.36160594606e-07,
  1.37718982103e-07, 1.39338316679e-07, 1.41025317971e-07, 1.42787873535e-07,
  1.44635331499e-07, 1.4657889173e-07, 1.48632138436e-07, 1.50811780719e-07,
  1.53138707402e-07, 1.55639532047e-07, 1.58348931426e-07, 1.61313325908e-07,
  1.64596952856e-07, 1.68292495203e-07, 1.72541128694e-07, 1.77574279496e-07,
  1.83813550477e-07, 1.92166040885e-07, 2.05295471952e-07, 2.22600839893e-07
};

// (Ziggurat) position of right-most step
const double Normal_RNG::PARAM_R = 3.44428647676;

// Get a Normal distributed (0,1) sample
double Normal_RNG::sample()
{
  uint32_t u, sign, i, j;
  double x, y;

  while(true) {
    u = RNG.genrand_uint32();
    sign = u & 0x80;            // 1 bit for the sign
    i = u & 0x7f;               // 7 bits to choose the step
    j = u >> 8;                 // 24 bits for the x-value

    x = j * wtab[i];

    if(j < ktab[i])
      break;

    if(i < 127) {
      y = ytab[i + 1] + (ytab[i] - ytab[i + 1]) * RNG.genrand_close_open();
    }
    else {
      x = PARAM_R - std::log(1.0 - RNG.genrand_close_open()) / PARAM_R;
      y = std::exp(-PARAM_R * (x - 0.5 * PARAM_R)) * RNG.genrand_close_open();
    }

    if(y < std::exp(-0.5 * x * x))
      break;
  }
  return sign ? x : -x;
}


///////////////////////////////////////////////
// Laplace_RNG
///////////////////////////////////////////////

Laplace_RNG::Laplace_RNG(double meanval, double variance)
{
  setup(meanval, variance);
}

void Laplace_RNG::setup(double meanval, double variance)
{
  mean = meanval;
  var = variance;
  sqrt_12var = std::sqrt(variance / 2.0);
}

void Laplace_RNG::get_setup(double &meanval, double &variance) const
{
  meanval = mean;
  variance = var;
}



vec Laplace_RNG::operator()(int n)
{
  vec vv(n);

  for(int i = 0; i < n; i++)
    vv(i) = sample();

  return vv;
}

mat Laplace_RNG::operator()(int h, int w)
{
  mat mm(h, w);
  int i, j;

  for(i = 0; i < h; i++)
    for(j = 0; j < w; j++)
      mm(i, j) = sample();

  return mm;
}

///////////////////////////////////////////////
// AR1_Normal_RNG
///////////////////////////////////////////////

AR1_Normal_RNG::AR1_Normal_RNG(double meanval, double variance, double rho)
{
  setup(meanval, variance, rho);
}

void AR1_Normal_RNG::setup(double meanval, double variance, double rho)
{
  mean = meanval;
  var = variance;
  r = rho;
  factr = -2.0 * var * (1.0 - rho * rho);
  mem = 0.0;
  odd = true;
}

void AR1_Normal_RNG::get_setup(double &meanval, double &variance,
                               double &rho) const
{
  meanval = mean;
  variance = var;
  rho = r;
}

vec AR1_Normal_RNG::operator()(int n)
{
  vec vv(n);

  for(int i = 0; i < n; i++)
    vv(i) = sample();

  return vv;
}

mat AR1_Normal_RNG::operator()(int h, int w)
{
  mat mm(h, w);
  int i, j;

  for(i = 0; i < h; i++)
    for(j = 0; j < w; j++)
      mm(i, j) = sample();

  return mm;
}

void AR1_Normal_RNG::reset()
{
  mem = 0.0;
}

///////////////////////////////////////////////
// Weibull_RNG
///////////////////////////////////////////////

Weibull_RNG::Weibull_RNG(double lambda, double beta)
{
  setup(lambda, beta);
}

void Weibull_RNG::setup(double lambda, double beta)
{
  l = lambda;
  b = beta;
  mean = tgamma(1.0 + 1.0 / b) / l;
  var = tgamma(1.0 + 2.0 / b) / (l * l) - mean;
}


vec Weibull_RNG::operator()(int n)
{
  vec vv(n);

  for(int i = 0; i < n; i++)
    vv(i) = sample();

  return vv;
}

mat Weibull_RNG::operator()(int h, int w)
{
  mat mm(h, w);
  int i, j;

  for(i = 0; i < h; i++)
    for(j = 0; j < w; j++)
      mm(i, j) = sample();

  return mm;
}

///////////////////////////////////////////////
// Rayleigh_RNG
///////////////////////////////////////////////

Rayleigh_RNG::Rayleigh_RNG(double sigma)
{
  setup(sigma);
}

vec Rayleigh_RNG::operator()(int n)
{
  vec vv(n);

  for(int i = 0; i < n; i++)
    vv(i) = sample();

  return vv;
}

mat Rayleigh_RNG::operator()(int h, int w)
{
  mat mm(h, w);
  int i, j;

  for(i = 0; i < h; i++)
    for(j = 0; j < w; j++)
      mm(i, j) = sample();

  return mm;
}

///////////////////////////////////////////////
// Rice_RNG
///////////////////////////////////////////////

Rice_RNG::Rice_RNG(double lambda, double beta)
{
  setup(lambda, beta);
}

vec Rice_RNG::operator()(int n)
{
  vec vv(n);

  for(int i = 0; i < n; i++)
    vv(i) = sample();

  return vv;
}

mat Rice_RNG::operator()(int h, int w)
{
  mat mm(h, w);
  int i, j;

  for(i = 0; i < h; i++)
    for(j = 0; j < w; j++)
      mm(i, j) = sample();

  return mm;
}

} // namespace itpp
