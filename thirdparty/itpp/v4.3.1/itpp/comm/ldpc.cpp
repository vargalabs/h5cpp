/*!
 * \file
 * \brief Implementation of Low-Density Parity Check (LDPC) codes
 * \author Erik G. Larsson, Mattias Andersson and Adam Piatyszek
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

#include <itpp/comm/ldpc.h>
#include <iomanip>
#include <sstream>

namespace itpp
{

/*!
 * \brief Version of the binary file with generator and decoder data
 *
 * This has to be global since it is used in LDPC_Generator and LDPC_Code
 * classes
 */
static const int LDPC_binary_file_version = 2;

// ---------------------------------------------------------------------------
// LDPC_Parity
// ---------------------------------------------------------------------------

// public methods

LDPC_Parity::LDPC_Parity(int nc, int nv): init_flag(false)
{
  initialize(nc, nv);
}

LDPC_Parity::LDPC_Parity(const std::string& filename,
                         const std::string& format): init_flag(false)
{
  if (format == "alist") {
    load_alist(filename);
  }
  else {
    it_error("LDPC_Parity::LDPC_Parity(): Only 'alist' format is supported");
  }
}

LDPC_Parity::LDPC_Parity(const GF2mat_sparse_alist &alist):
    init_flag(false)
{
  import_alist(alist);
}

void LDPC_Parity::initialize(int nc, int nv)
{
  ncheck = nc;
  nvar = nv;
  H = GF2mat_sparse(ncheck, nvar);
  Ht = GF2mat_sparse(nvar, ncheck);
  sumX1 = zeros_i(nvar);
  sumX2 = zeros_i(ncheck);
  init_flag = true;
}

void LDPC_Parity::set(int i, int j, bin x)
{
  it_assert(init_flag, "LDPC_Parity::set(): Object not initialized");
  it_assert_debug((i >= 0) && (i < ncheck),
                  "LDPC_Parity::set(): Wrong index i");
  it_assert_debug((j >= 0) && (j < nvar),
                  "LDPC_Parity::set(): Wrong index j");
  it_assert_debug(H(i, j) == Ht(j, i), "LDPC_Parity:set(): Internal error");

  int diff = static_cast<int>(x) - static_cast<int>(H(i, j));
  sumX1(j) += diff;
  sumX2(i) += diff;

  if (x == 1) {
    H.set(i, j, 1);
    Ht.set(j, i, 1);
  }
  else {
    H.clear_elem(i, j);
    Ht.clear_elem(j, i);
  }

  it_assert_debug(H(i, j) == x, "LDPC_Parity::set(): Internal error");
  it_assert_debug(Ht(j, i) == x, "LDPC_Parity::set(): Internal error");
}

void LDPC_Parity::display_stats() const
{
  it_assert(init_flag,
            "LDPC_Parity::display_stats(): Object not initialized");
  int cmax = max(sumX1);
  int vmax = max(sumX2);
  vec vdeg = zeros(cmax + 1); // number of variable nodes with n neighbours
  vec cdeg = zeros(vmax + 1); // number of check nodes with n neighbours
  for (int col = 0; col < nvar; col++) {
    vdeg(length(get_col(col).get_nz_indices()))++;
    it_assert(sumX1(col) == length(get_col(col).get_nz_indices()),
              "LDPC_Parity::display_stats(): Internal error");
  }
  for (int row = 0; row < ncheck; row++) {
    cdeg(length(get_row(row).get_nz_indices()))++;
    it_assert(sumX2(row) == length(get_row(row).get_nz_indices()),
              "LDPC_Parity::display_stats(): Internal error");
  }

  // from edge perspective
  // number of edges connected to vnodes of degree n
  vec vdegedge = elem_mult(vdeg, linspace(0, vdeg.length() - 1,
                                          vdeg.length()));
  // number of edges connected to cnodes of degree n
  vec cdegedge = elem_mult(cdeg, linspace(0, cdeg.length() - 1,
                                          cdeg.length()));

  int edges = sum(elem_mult(to_ivec(linspace(0, vdeg.length() - 1,
                                    vdeg.length())),
                            to_ivec(vdeg)));

  it_info("--- LDPC parity check matrix ---");
  it_info("Dimension [ncheck x nvar]: " << ncheck << " x " << nvar);
  it_info("Variable node degree distribution from node perspective:\n"
          << vdeg / nvar);
  it_info("Check node degree distribution from node perspective:\n"
          << cdeg / ncheck);
  it_info("Variable node degree distribution from edge perspective:\n"
          << vdegedge / edges);
  it_info("Check node degree distribution from edge perspective:\n"
          << cdegedge / edges);
  it_info("Rate: " << get_rate());
  it_info("--------------------------------");
}


void LDPC_Parity::load_alist(const std::string& alist_file)
{
  import_alist(GF2mat_sparse_alist(alist_file));
}

void LDPC_Parity::save_alist(const std::string& alist_file) const
{
  GF2mat_sparse_alist alist = export_alist();
  alist.write(alist_file);
}


void LDPC_Parity::import_alist(const GF2mat_sparse_alist& alist)
{
  GF2mat_sparse X = alist.to_sparse();

  initialize(X.rows(), X.cols());
  // brute force copy from X to this
  for (int i = 0; i < ncheck; i++) {
    for (int j = 0; j < nvar; j++) {
      if (X(i, j)) {
        set(i, j, 1);
      }
    }
  }
}

GF2mat_sparse_alist LDPC_Parity::export_alist() const
{
  it_assert(init_flag,
            "LDPC_Parity::export_alist(): Object not initialized");
  GF2mat_sparse_alist alist;
  alist.from_sparse(H);
  return alist;
}


int LDPC_Parity::check_connectivity(int from_i, int from_j, int to_i,
                                    int to_j, int godir, int L) const
{
  it_assert(init_flag,
            "LDPC_Parity::check_connectivity(): Object not initialized");
  int i, j, result;

  if (L < 0) {         // unable to reach coordinate with given L
    return (-3);
  }

  // check if reached destination
  if ((from_i == to_i) && (from_j == to_j) && (godir != 0)) {
    return L;
  }

  if (get(from_i, from_j) == 0) {  // meaningless search
    return (-2);
  }

  if (L == 2) {    // Treat this case separately for efficiency
    if (godir == 2) { // go horizontally
      if (get(from_i, to_j) == 1) { return 0; }
    }
    if (godir == 1) { // go vertically
      if (get(to_i, from_j) == 1) { return 0; }
    }
    return (-3);
  }

  if ((godir == 1) || (godir == 0)) {   // go vertically
    ivec cj = get_col(from_j).get_nz_indices();
    for (i = 0; i < length(cj); i++) {
      if (cj(i) != from_i) {
        result = check_connectivity(cj(i), from_j, to_i, to_j, 2, L - 1);
        if (result >= 0) {
          return (result);
        }
      }
    }
  }

  if (godir == 2) { // go horizontally
    ivec ri = get_row(from_i).get_nz_indices();
    for (j = 0; j < length(ri); j++) {
      if (ri(j) != from_j) {
        result = check_connectivity(from_i, ri(j), to_i, to_j, 1, L - 1);
        if (result >= 0) {
          return (result);
        }
      }
    }
  }

  return (-1);
};

int LDPC_Parity::check_for_cycles(int L) const
{
  it_assert(init_flag,
            "LDPC_Parity::check_for_cycles(): Object not initialized");
  // looking for odd length cycles does not make sense
  if ((L&1) == 1) { return (-1); }
  if (L == 0) { return (-4); }

  int cycles = 0;
  for (int i = 0; i < nvar; i++) {
    ivec ri = get_col(i).get_nz_indices();
    for (int j = 0; j < length(ri); j++) {
      if (check_connectivity(ri(j), i, ri(j), i, 0, L) >= 0) {
        cycles++;
      }
    }
  }
  return cycles;
};

//   ivec LDPC_Parity::get_rowdegree()  const
//   {
//     ivec rdeg = zeros_i(Nmax);
//     for (int i=0; i<ncheck; i++)     {
//       rdeg(sumX2(i))++;
//     }
//     return rdeg;
//   };

//   ivec LDPC_Parity::get_coldegree()  const
//   {
//     ivec cdeg = zeros_i(Nmax);
//     for (int j=0; j<nvar; j++)     {
//       cdeg(sumX1(j))++;
//     }
//     return cdeg;
//   };


// ----------------------------------------------------------------------
// LDPC_Parity_Unstructured
// ----------------------------------------------------------------------

int LDPC_Parity_Unstructured::cycle_removal_MGW(int Maxcyc)
{
  it_assert(init_flag,
            "LDPC_Parity::cycle_removal_MGW(): Object not initialized");
  typedef Sparse_Mat<short> Ssmat;
  typedef Sparse_Vec<short> Ssvec;

  Maxcyc -= 2;

  // Construct the adjacency matrix of the graph
  Ssmat G(ncheck + nvar, ncheck + nvar, 5);
  for (int j = 0; j < nvar; j++) {
    GF2vec_sparse col = get_col(j);
    for (int i = 0; i < col.nnz(); i++) {
      if (get(col.get_nz_index(i), j) == 1) {
        G.set(col.get_nz_index(i), j + ncheck, 1);
        G.set(ncheck + j, col.get_nz_index(i), 1);
      }
    }
  }

  Array<Ssmat> Gpow(Maxcyc);
  Gpow(0).set_size(ncheck + nvar, ncheck + nvar, 1);
  Gpow(0).clear();
  for (int i = 0; i < ncheck + nvar; i++) {
    Gpow(0).set(i, i, 1);
  }
  Gpow(1) = G;

  /* Main cycle elimination loop starts here. Note that G and all
     powers of G are symmetric matrices. This fact is exploited in
     the code.
  */
  int r;
  int cycles_found = 0;
  int scl = Maxcyc;
  for (r = 4; r <= Maxcyc; r += 2) {
    // compute the next power of the adjacency matrix
    Gpow(r / 2) = Gpow(r / 2 - 1) * G;
    bool traverse_again;
    do {
      traverse_again = false;
      cycles_found = 0;
      it_info_debug("Starting new pass of cycle elimination, target girth "
                    << (r + 2) << "...");
      int pdone = 0;
      for (int j = 0; j < ncheck + nvar; j++) { // loop over elements of G
        for (int i = 0; i < ncheck + nvar; i++) {
          int ptemp = floor_i(100.0 * (i + j * (ncheck + nvar)) /
                              ((nvar + ncheck) * (nvar + ncheck)));
          if (ptemp > pdone + 10) {
            it_info_debug(ptemp << "% done.");
            pdone = ptemp;
          }

          if (((Gpow(r / 2))(i, j) >= 2)  && ((Gpow(r / 2 - 2))(i, j) == 0)) {
            // Found a cycle.
            cycles_found++;

            // choose k
            ivec tmpi = (elem_mult(Gpow(r / 2 - 1).get_col(i),
                                   G.get_col(j))).get_nz_indices();
            //        int k = tmpi(rand()%length(tmpi));
            int k = tmpi(randi(0, length(tmpi) - 1));
            it_assert_debug(G(j, k) == 1 && G(k, j) == 1,
                            "LDPC_Parity_Unstructured::cycle_removal_MGW(): "
                            "Internal error");

            // determine candidate edges for an edge swap
            Ssvec rowjk = Gpow(r / 2) * (Gpow(r / 2 - 1).get_col(j)
                                         + Gpow(r / 2 - 1).get_col(k));
            int p, l;
            ivec Ce_ind = sort_index(randu(nvar + ncheck)); // random order

            for (int s = 0; s < nvar + ncheck; s++)  {
              l = Ce_ind(s);
              if (rowjk(l) != 0) { continue; }
              ivec colcandi = G.get_col(l).get_nz_indices();
              if (length(colcandi) > 0)  {
                // select a node p which is a member of Ce
                for (int u = 0; u < length(colcandi); u++) {
                  p = colcandi(u);
                  if (p != l) {
                    if (rowjk(p) == 0) {
                      goto found_candidate_vector;
                    }
                  }
                }
              }
            }
            continue; // go to the next entry (i,j)

          found_candidate_vector:
            // swap edges

            if (p >= ncheck) { int z = l; l = p; p = z; }
            if (j >= ncheck) { int z = k; k = j; j = z; }

            // Swap endpoints of edges (p,l) and (j,k)
            // cout << "(" << j << "," << k << ")<->("
            // << p << "," << l << ") " ;
            // cout << ".";
            // cout.flush();

            // Update the matrix
            it_assert_debug((get(j, k - ncheck) == 1) && (get(p, l - ncheck) == 1),
                            "LDPC_Parity_Unstructured::cycle_removal_MGW(): "
                            "Internal error");
            set(j, k - ncheck, 0);
            set(p, l - ncheck, 0);
            it_assert_debug((get(j, l - ncheck) == 0) && (get(p, k - ncheck) == 0),
                            "LDPC_Parity_Unstructured::cycle_removal_MGW(): "
                            "Internal error");
            set(j, l - ncheck, 1);
            set(p, k - ncheck, 1);

            // Update adjacency matrix
            it_assert_debug(G(p, l) == 1 && G(l, p) == 1 && G(j, k) == 1
                            && G(k, j) == 1, "G");
            it_assert_debug(G(j, l) == 0 && G(l, j) == 0 && G(p, k) == 0
                            && G(k, p) == 0, "G");

            // Delta is the update term to G
            Ssmat Delta(ncheck + nvar, ncheck + nvar, 2);
            Delta.set(j, k, -1);
            Delta.set(k, j, -1);
            Delta.set(p, l, -1);
            Delta.set(l, p, -1);
            Delta.set(j, l, 1);
            Delta.set(l, j, 1);
            Delta.set(p, k, 1);
            Delta.set(k, p, 1);

            // update G and its powers
            G = G + Delta;
            it_assert_debug(G(p, l) == 0 && G(l, p) == 0 && G(j, k) == 0
                            && G(k, j) == 0, "G");
            it_assert_debug(G(j, l) == 1 && G(l, j) == 1 && G(p, k) == 1
                            && G(k, p) == 1, "G");

            Gpow(1) = G;
            Gpow(2) = G * G;
            for (int z = 3; z <= r / 2; z++) {
              Gpow(z) = Gpow(z - 1) * G;
            }

            traverse_again = true;
          } // if G()...
        } // loop over i
      }  // loop over j
      if ((!traverse_again) && (cycles_found > 0)) { // no point continue
        scl = r - 2;
        goto finished;
      }
    }
    while (cycles_found != 0);
    scl = r;  // there were no cycles of length r; move on to next r
    it_info_debug("Achieved girth " << (scl + 2)
                  << ". Proceeding to next level.");
  } // loop over r

finished:
  int girth = scl + 2;  // scl=length of smallest cycle
  it_info_debug("Cycle removal (MGW algoritm) finished. Graph girth: "
                << girth << ". Cycles remaining on next girth level: "
                << cycles_found);

  return girth;
}

void LDPC_Parity_Unstructured::generate_random_H(const ivec& C,
    const ivec& R,
    const ivec& cycopt)
{
  // Method based on random permutation. Attempts to avoid placing new
  // edges so that cycles are created. More aggressive cycle avoidance
  // for degree-2 nodes. EGL January 2007.

  initialize(sum(R), sum(C));
  // C, R: Target number of columns/rows with certain number of ones

  // compute number of edges
  int Ne = 0;
  for (int i = 0;i < C.length();i++) {
    for (int j = 0; j < C(i); j++) {
      for (int m = 0; m < i; m++) Ne++;
    }
  }

  // compute connectivity matrix
  ivec vcon(Ne);
  ivec ccon(Ne);
  ivec vd(nvar);
  ivec cd(ncheck);
  int k = 0;
  int l = 0;
  for (int i = 0;i < C.length();i++) {
    for (int j = 0; j < C(i); j++) {
      for (int m = 0; m < i; m++) {
        vcon(k) = l;
        vd(l) = i;
        k++;
      }
      l++;
    }
  }
  k = 0;
  l = 0;
  for (int i = 0;i < R.length();i++) {
    for (int j = 0; j < R(i); j++) {
      for (int m = 0; m < i; m++) {
        ccon(k) = l;
        cd(l) = i;
        k++;
      }
      l++;
    }
  }
  it_assert(k == Ne, "C/R mismatch");

  // compute random permutations
  ivec ind = sort_index(randu(Ne));
  ivec cp = sort_index(randu(nvar));
  ivec rp = sort_index(randu(ncheck));

  // set the girth goal for various variable node degrees
  ivec Laim = zeros_i(Nmax);
  for (int i = 0; i < length(cycopt); i++) {
    Laim(i + 2) = cycopt(i);
  }
  for (int i = length(cycopt); i < Nmax - 2; i++) {
    Laim(i + 2) = cycopt(length(cycopt) - 1);
  }
  it_info_debug("Running with Laim=" << Laim.left(25));

  int failures = 0;
  const int Max_attempts = 100;
  const int apcl = 10;    // attempts before reducing girth target
  for (int k = 0; k < Ne; k++) {
    const int el = Ne - k - 2;
    if (k % 250 == 0) {
      it_info_debug("Processing edge: " << k << " out of " << Ne
                    << ". Variable node degree: " << vd(vcon(k))
                    << ". Girth target: " << Laim(vd(vcon(k)))
                    << ". Accumulated failures: " << failures);
    }
    const int c = cp(vcon(k));
    int L = Laim(vd(vcon(k)));
    int attempt = 0;
    while (true) {
      if (attempt > 0 && attempt % apcl == 0 && L >= 6) { L -= 2; };
      int r = rp(ccon(ind(k)));
      if (get(r, c)) { // double edge
        // set(r,c,0);
        if (el > 0) {
          //      int t=k+1+rand()%el;
          int t = k + 1 + randi(0, el - 1);
          int x = ind(t);
          ind(t) = ind(k);
          ind(k) = x;
          attempt++;
          if (attempt == Max_attempts) {
            failures++;
            break;
          }
        }
        else {  // almost at the last edge
          break;
        }
      }
      else {
        set(r, c, 1);
        if (L > 0) { // attempt to avoid cycles
          if (check_connectivity(r, c, r, c, 0, L) >= 0) {   // detected cycle
            set(r, c, 0);
            if (el > 0) {
              // make a swap in the index permutation
              //   int t=k+1+rand()%el;
              int t = k + 1 + randi(0, el - 1);
              int x = ind(t);
              ind(t) = ind(k);
              ind(k) = x;
              attempt++;
              if (attempt == Max_attempts) {  // give up
                failures++;
                set(r, c, 1);
                break;
              }
            }
            else {  // no edges left
              set(r, c, 1);
              break;
            }
          }
          else {
            break;
          }
        }
        else {
          break;
        }
      }
    }
  }
}

void LDPC_Parity_Unstructured::compute_CR(const vec& var_deg, const vec& chk_deg, const int Nvar,
    ivec &C, ivec &R)
{
  // compute the degree distributions from a node perspective
  vec Vi = linspace(1, length(var_deg), length(var_deg));
  vec Ci = linspace(1, length(chk_deg), length(chk_deg));
  // Compute number of cols with n 1's
  // C, R: Target number of columns/rows with certain number of ones
  C = to_ivec(round(Nvar * elem_div(var_deg, Vi)
                    / sum(elem_div(var_deg, Vi))));
  C = concat(0, C);
  int edges = sum(elem_mult(to_ivec(linspace(0, C.length() - 1,
                                    C.length())), C));
  R = to_ivec(round(edges * elem_div(chk_deg, Ci)));
  R = concat(0, R);
  vec Ri = linspace(0, length(R) - 1, length(R));
  vec Coli = linspace(0, length(C) - 1, length(C));

  // trim to deal with inconsistencies due to rounding errors
  if (sum(C) != Nvar) {
    ivec ind = find(C == max(C));
    C(ind(0)) = C(ind(0)) - (sum(C) - Nvar);
  }

  //the number of edges calculated from R must match the number of
  //edges calculated from C
  while (sum(elem_mult(to_vec(R), Ri)) !=
         sum(elem_mult(to_vec(C), Coli))) {
    //we're only changing R, this is probably(?) better for irac codes
    if (sum(elem_mult(to_vec(R), Ri)) > sum(elem_mult(to_vec(C), Coli))) {
      //remove an edge from R
      ivec ind = find(R == max(R));
      int old = R(ind(0));
      R.set(ind(0), old - 1);
      old = R(ind(0) - 1);
      R.set(ind(0) - 1, old + 1);
    }
    else {
      ivec ind = find(R == max(R));
      if (ind(0) == R.length() - 1) {
        R = concat(R, 0);
        Ri = linspace(0, length(R) - 1, length(R));
      }
      int old = R(ind(0));
      R.set(ind(0), old - 1);
      old = R(ind(0) + 1);
      R.set(ind(0) + 1, old + 1);
    }
  }

  C = concat(C, zeros_i(Nmax - length(C)));
  R = concat(R, zeros_i(Nmax - length(R)));

  it_info_debug("C=" << C << std::endl);
  it_info_debug("R=" << R << std::endl);

}

// ----------------------------------------------------------------------
// LDPC_Parity_Regular
// ----------------------------------------------------------------------

LDPC_Parity_Regular::LDPC_Parity_Regular(int Nvar, int k, int l,
    const std::string& method,
    const ivec& options)
{
  generate(Nvar, k, l, method, options);
}

void LDPC_Parity_Regular::generate(int Nvar, int k, int l,
                                   const std::string& method,
                                   const ivec& options)
{
  vec var_deg = zeros(k);
  vec chk_deg = zeros(l);
  var_deg(k - 1) = 1;
  chk_deg(l - 1) = 1;

  ivec C, R;
  compute_CR(var_deg, chk_deg, Nvar, C, R);
  it_info_debug("sum(C)=" << sum(C) << "  Nvar=" << Nvar);
  it_info_debug("sum(R)=" << sum(R) << "  approximate target=" << round_i(Nvar * k / static_cast<double>(l)));

  if (method == "rand") {
    generate_random_H(C, R, options);
  }
  else {
    it_error("not implemented");
  };
}


// ----------------------------------------------------------------------
// LDPC_Parity_Irregular
// ----------------------------------------------------------------------

LDPC_Parity_Irregular::LDPC_Parity_Irregular(int Nvar,
    const vec& var_deg,
    const vec& chk_deg,
    const std::string& method,
    const ivec& options)
{
  generate(Nvar, var_deg, chk_deg, method, options);
}

void LDPC_Parity_Irregular::generate(int Nvar, const vec& var_deg,
                                     const vec& chk_deg,
                                     const std::string& method,
                                     const ivec& options)
{
  ivec C, R;
  compute_CR(var_deg, chk_deg, Nvar, C, R);

  if (method == "rand") {
    generate_random_H(C, R, options);
  }
  else {
    it_error("not implemented");
  };
}

// ----------------------------------------------------------------------
// BLDPC_Parity
// ----------------------------------------------------------------------

BLDPC_Parity::BLDPC_Parity(const imat& base_matrix, int exp_factor)
{
  expand_base(base_matrix, exp_factor);
}

BLDPC_Parity::BLDPC_Parity(const std::string& filename, int exp_factor)
{
  load_base_matrix(filename);
  expand_base(H_b, exp_factor);
}

void BLDPC_Parity::expand_base(const imat& base_matrix, int exp_factor)
{
  Z = exp_factor;
  H_b = base_matrix;
  H_b_valid = true;

  initialize(H_b.rows() * Z, H_b.cols() * Z);
  for (int r = 0; r < H_b.rows(); r++) {
    for (int c = 0; c < H_b.cols(); c++) {
      int rz = r * Z;
      int cz = c * Z;
      switch (H_b(r, c)) {
      case -1:
        break;
      case 0:
        for (int i = 0; i < Z; ++i)
          set(rz + i, cz + i, 1);
        break;
      default:
        for (int i = 0; i < Z; ++i)
          set(rz + i, cz + (i + H_b(r, c)) % Z, 1);
        break;
      }
    }
  }
}


int BLDPC_Parity::get_exp_factor() const
{
  return Z;
}

imat BLDPC_Parity::get_base_matrix() const
{
  return H_b;
}

void BLDPC_Parity::set_exp_factor(int exp_factor)
{
  Z = exp_factor;
  if (H_b_valid) {
    expand_base(H_b, exp_factor);
  }
  else {
    calculate_base_matrix();
  }
}


void BLDPC_Parity::load_base_matrix(const std::string& filename)
{
  std::ifstream bm_file(filename.c_str());
  it_assert(bm_file.is_open(), "BLDPC_Parity::load_base_matrix(): Could not "
            "open file \"" << filename << "\" for reading");
  // delete old base matrix content
  H_b.set_size(0, 0);
  // parse text file content, row by row
  std::string line;
  int line_counter = 0;
  getline(bm_file, line);
  while (!bm_file.eof()) {
    line_counter++;
    std::stringstream ss(line);
    ivec row(0);
    while (ss.good()) {
      int val;
      ss >> val;
      row = concat(row, val);
    }
    if ((H_b.rows() == 0) || (row.size() == H_b.cols()))
      H_b.append_row(row);
    else
      it_warning("BLDPC_Parity::load_base_matrix(): Wrong size of "
                 "a parsed row number " << line_counter);
    getline(bm_file, line);
  }
  bm_file.close();

  // transpose parsed base matrix if necessary
  if (H_b.rows() > H_b.cols())
    H_b = H_b.transpose();

  H_b_valid = true;
  init_flag = false;
}

void BLDPC_Parity::save_base_matrix(const std::string& filename) const
{
  it_assert(H_b_valid, "BLDPC_Parity::save_base_matrix(): Base matrix is "
            "not valid");
  std::ofstream bm_file(filename.c_str());
  it_assert(bm_file.is_open(), "BLDPC_Parity::save_base_matrix(): Could not "
            "open file \"" << filename << "\" for writing");

  for (int r = 0; r < H_b.rows(); r++) {
    for (int c = 0; c < H_b.cols(); c++) {
      bm_file << std::setw(3) << H_b(r, c);
    }
    bm_file << "\n";
  }

  bm_file.close();
}


void BLDPC_Parity::calculate_base_matrix()
{
  std::string error_str = "BLDPC_Parity::calculate_base_matrix(): "
                          "Invalid BLDPC matrix. Cannot calculate base matrix from it.";
  int rows = H.rows() / Z;
  int cols = H.cols() / Z;
  it_assert((rows * Z == H.rows()) && (cols * Z == H.cols()), error_str);
  H_b.set_size(rows, cols);

  for (int r = 0; r < rows; ++r) {
    int rz = r * Z;
    for (int c = 0; c < cols; ++c) {
      int cz = c * Z;
      GF2mat_sparse H_Z = H.get_submatrix(rz, rz + Z - 1, cz, cz + Z - 1);

      if (H_Z.nnz() == 0) {
        H_b(r, c) = -1;
      }
      else if (H_Z.nnz() == Z) {
        // check for cyclic-shifted ZxZ matrix
        int shift = 0;
        while ((shift < Z) && (H_Z(0, shift) != 1))
          ++shift;
        it_assert(shift < Z, error_str);
        for (int i = 1; i < Z; ++i)
          it_assert(H_Z(0, shift) == H_Z(i, (i + shift) % Z), error_str);
        H_b(r, c) = shift;
      }
      else {
        it_error(error_str);
      }
    } // for (int c = 0; c < cols; ++c)
  } // for (int r = 0; r < rows; ++r)

  it_info("Base matrix calculated");
  H_b_valid = true;
}



// ----------------------------------------------------------------------
// LDPC_Generator_Systematic
// ----------------------------------------------------------------------

LDPC_Generator_Systematic::LDPC_Generator_Systematic(LDPC_Parity* const H,
    bool natural_ordering,
    const ivec& ind):
    LDPC_Generator("systematic"), G()
{
  ivec tmp;
  tmp = construct(H, natural_ordering, ind);
}


ivec LDPC_Generator_Systematic::construct(LDPC_Parity* const H,
    bool natural_ordering,
    const ivec& avoid_cols)
{
  int nvar = H->get_nvar();
  int ncheck = H->get_ncheck();

  // create dense representation of parity check matrix
  GF2mat Hd(H->get_H());

  // -- Determine initial column ordering --
  ivec col_order(nvar);
  if (natural_ordering) {
    for (int i = 0; i < nvar; i++) {
      col_order(i) = i;
    }
  }
  else {
    // take the columns in random order, but the ones to avoid at last
    vec col_importance = randu(nvar);
    for (int i = 0; i < length(avoid_cols); i++) {
      col_importance(avoid_cols(i)) = (-col_importance(avoid_cols(i)));
    }
    col_order = sort_index(-col_importance);
  }

  ivec actual_ordering(nvar);

  // Now partition P as P=[P1 P2]. Then find G so [P1 P2][I G]'=0.

  // -- Create P1 and P2 --
  GF2mat P1; //(ncheck,nvar-ncheck);      // non-invertible part
  GF2mat P2; //(ncheck,ncheck);           // invertible part

  it_info_debug("Computing a systematic generator matrix...");

  int j1 = 0, j2 = 0;
  ivec perm;
  GF2mat T, U;
  for (int k = 0; k < nvar; k++) {
    it_error_if(j1 >= nvar - ncheck, "LDPC_Generator_Systematic::construct(): "
                "Unable to obtain enough independent columns.");

    bvec c = Hd.get_col(col_order(k));
    if (j2 == 0) {     // first column in P2 is number col_order(k)
      P2 = GF2mat(c);
      P2.T_fact(T, U, perm);/* return value not used */
      actual_ordering(k) = nvar - ncheck;
      j2++;
    }
    else {
      if (j2 < ncheck) {
        if (P2.T_fact_update_addcol(T, U, perm, c)) {
          P2 = P2.concatenate_horizontal(c);
          actual_ordering(k) = nvar - ncheck + j2;
          j2++;
          continue;
        }
      }
      if (j1 == 0) {
        P1 = GF2mat(c);
        actual_ordering(k) = j1;
      }
      else {
        P1 = P1.concatenate_horizontal(c);
        actual_ordering(k) = j1;
      }
      j1++;
    }
  }

  it_info_debug("Rank of parity check matrix: " << j2);

  // -- Compute the systematic part of the generator matrix --
  G = (P2.inverse() * P1).transpose();

  // -- Permute the columns of the parity check matrix --
  GF2mat P = P1.concatenate_horizontal(P2);
  *H = LDPC_Parity(ncheck, nvar);
  // brute force copy from P to H
  for (int i = 0; i < ncheck; i++) {
    for (int j = 0; j < nvar; j++) {
      if (P.get(i, j)) {
        H->set(i, j, 1);
      }
    }
  }

  // -- Check that the result was correct --
  it_assert_debug((GF2mat(H->get_H())
                   * (gf2dense_eye(nvar - ncheck).concatenate_horizontal(G)).transpose()).is_zero(),
                  "LDPC_Generator_Systematic::construct(): Incorrect generator matrix G");

  G = G.transpose();  // store the generator matrix in transposed form
  it_info_debug("Systematic generator matrix computed.");

  mark_initialized();

  return actual_ordering;
}


void LDPC_Generator_Systematic::save(const std::string& filename) const
{
  it_file f(filename);
  int ver;
  f >> Name("Fileversion") >> ver;
  it_assert(ver == LDPC_binary_file_version,
            "LDPC_Generator_Systematic::save(): Unsupported file format");
  f << Name("G_type") << get_type();
  f << Name("G") << G;
  f.close();
}


void LDPC_Generator_Systematic::load(const std::string& filename)
{
  it_ifile f(filename);
  int ver;
  f >> Name("Fileversion") >> ver;
  it_assert(ver == LDPC_binary_file_version,
            "LDPC_Generator_Systematic::load(): Unsupported file format");
  std::string gen_type;
  f >> Name("G_type") >> gen_type;
  it_assert(gen_type == get_type(),
            "LDPC_Generator_Systematic::load(): Wrong generator type");
  f >> Name("G") >> G;
  f.close();

  mark_initialized();
}


void LDPC_Generator_Systematic::encode(const bvec &input, bvec &output)
{
  it_assert(is_initialized(), "LDPC_Generator_Systematic::encode(): Systematic "
            "generator not set up");
  it_assert(input.size() == G.cols(), "LDPC_Generator_Systematic::encode(): "
            "Improper input vector size (" << input.size() << " != "
            << G.cols() << ")");

  output = concat(input, G * input);
}


// ----------------------------------------------------------------------
// BLDPC_Generator
// ----------------------------------------------------------------------

BLDPC_Generator::BLDPC_Generator(const BLDPC_Parity* const H,
                                 const std::string type):
    LDPC_Generator(type), H_enc(), N(0), M(0), K(0), Z(0)
{
  construct(H);
}


void BLDPC_Generator::encode(const bvec &input, bvec &output)
{
  it_assert(is_initialized(), "BLDPC_Generator::encode(): Cannot encode with not "
            "initialized generator");
  it_assert(input.size() == K, "BLDPC_Generator::encode(): Input vector "
            "length is not equal to K");

  // copy systematic bits first
  output = input;
  output.set_size(N, true);

  // backward substitution to obtain the first Z parity check bits
  for (int k = 0; k < Z; k++) {
    for (int j = 0; j < K; j++) {
      output(K + k) += H_enc(M - 1 - k, j) * input(j);
    }
    for (int j = 0; j < k; j++) {
      output(K + k) += H_enc(M - 1 - k, K + j) * output(K + j);
    }
  }

  // forward substitution to obtain the next M-Z parity check bits
  for (int k = 0; k < M - Z; k++) {
    for (int j = 0; j < K; j++) {
      output(K + Z + k) += H_enc(k, j) * input(j);
    }
    for (int j = K; j < K + Z; j++) {
      output(K + Z + k) += H_enc(k, j) * output(j);
    }
    for (int j = K + Z; j < K + Z + k; j++) {
      output(K + Z + k) += H_enc(k, j) * output(j);
    }
  }
}


void BLDPC_Generator::construct(const BLDPC_Parity* const H)
{
  if (H != 0 && H->is_valid()) {
    H_enc = H->get_H(); // sparse to dense conversion
    Z = H->get_exp_factor();
    N = H_enc.cols();
    M = H_enc.rows();
    K = N - M;

    // ----------------------------------------------------------------------
    // STEP 1
    // ----------------------------------------------------------------------
    // loop over last M-Z columns of matrix H
    for (int i = 0; i < M - Z; i += Z) {
      // loop over last Z rows of matrix H
      for (int j = 0; j < Z; j++) {
        // Gaussian elimination by adding particular rows
        H_enc.add_rows(M - 1 - j, M - Z - 1 - j - i);
      }
    }

    // ----------------------------------------------------------------------
    // STEP 2
    // ----------------------------------------------------------------------
    // set first processed row index to M-Z
    int r1 = M - Z;
    // loop over columns with indexes K .. K+Z-1
    for (int c1 = K + Z - 1; c1 >= K; c1--) {
      int r2 = r1;
      // find the first '1' in column c1
      while (H_enc(r2, c1) == 0 && r2 < M - 1)
        r2++;
      // if necessary, swap rows r1 and r2
      if (r2 != r1)
        H_enc.swap_rows(r1, r2);

      // look for the other '1' in column c1 and get rid of them
      for (r2 = r1 + 1; r2 < M; r2++) {
        if (H_enc(r2, c1) == 1) {
          // Gaussian elimination by adding particular rows
          H_enc.add_rows(r2, r1);
        }
      }

      // increase first processed row index
      r1++;
    }

    mark_initialized();
  }
}


void BLDPC_Generator::save(const std::string& filename) const
{
  it_assert(is_initialized(),
            "BLDPC_Generator::save(): Can not save not initialized generator");
  // Every Z-th row of H_enc until M-Z
  GF2mat H_T(M / Z - 1, N);
  for (int i = 0; i < M / Z - 1; i++) {
    H_T.set_row(i, H_enc.get_row(i*Z));
  }
  // Last Z preprocessed rows of H_enc
  GF2mat H_Z = H_enc.get_submatrix(M - Z, 0, M - 1, N - 1);

  it_file f(filename);
  int ver;
  f >> Name("Fileversion") >> ver;
  it_assert(ver == LDPC_binary_file_version, "BLDPC_Generator::save(): "
            "Unsupported file format");
  f << Name("G_type") << get_type();
  f << Name("H_T") << H_T;
  f << Name("H_Z") << H_Z;
  f << Name("Z") << Z;
  f.close();
}

void BLDPC_Generator::load(const std::string& filename)
{
  GF2mat H_T, H_Z;

  it_ifile f(filename);
  int ver;
  f >> Name("Fileversion") >> ver;
  it_assert(ver == LDPC_binary_file_version, "BLDPC_Generator::load(): "
            "Unsupported file format");
  std::string gen_type;
  f >> Name("G_type") >> gen_type;
  it_assert(gen_type == get_type(),
            "BLDPC_Generator::load(): Wrong generator type");
  f >> Name("H_T") >> H_T;
  f >> Name("H_Z") >> H_Z;
  f >> Name("Z") >> Z;
  f.close();

  N = H_T.cols();
  M = (H_T.rows() + 1) * Z;
  K = N - M;

  H_enc = GF2mat(M - Z, N);
  for (int i = 0; i < H_T.rows(); i++) {
    for (int j = 0; j < Z; j++) {
      for (int k = 0; k < N; k++) {
        if (H_T(i, (k / Z)*Z + (k + Z - j) % Z)) {
          H_enc.set(i*Z + j, k, 1);
        }
      }
    }
  }
  H_enc = H_enc.concatenate_vertical(H_Z);

  mark_initialized();
}


// ----------------------------------------------------------------------
// LDPC_Code
// ----------------------------------------------------------------------

LDPC_Code::LDPC_Code(): H_defined(false), G_defined(false), dec_method(new std::string),
    max_iters(50), psc(true), pisc(false),
    llrcalc(LLR_calc_unit()) { set_decoding_method("BP");}

LDPC_Code::LDPC_Code(const LDPC_Parity* const H,
                     LDPC_Generator* const G_in,
                     bool perform_integrity_check):
    H_defined(false), G_defined(false), dec_method(new std::string), max_iters(50),
    psc(true), pisc(false), llrcalc(LLR_calc_unit())
{
    set_decoding_method("BP");
    set_code(H, G_in, perform_integrity_check);
}

LDPC_Code::LDPC_Code(const std::string& filename,
                     LDPC_Generator* const G_in):
    H_defined(false), G_defined(false), dec_method(new std::string), max_iters(50),
    psc(true), pisc(false), llrcalc(LLR_calc_unit())
{
  set_decoding_method("BP");
  load_code(filename, G_in);
}


void LDPC_Code::set_code(const LDPC_Parity* const H,
                         LDPC_Generator* const G_in,
                         bool perform_integrity_check)
{
  decoder_parameterization(H);
  setup_decoder();
  G = G_in;
  if (G != 0) {
    G_defined = true;
    if (perform_integrity_check) {
      integrity_check();
    } else {
      it_info_debug("LDPC_Code::set_code(): integrity check was not performed");
    }
  }
}

void LDPC_Code::load_code(const std::string& filename,
                          LDPC_Generator* const G_in)
{
  it_info_debug("LDPC_Code::load_code(): Loading LDPC codec from "
                << filename);

  it_ifile f(filename);
  int ver;
  f >> Name("Fileversion") >> ver;
  it_assert(ver == LDPC_binary_file_version, "LDPC_Code::load_code(): "
            "Unsupported file format");
  f >> Name("H_defined") >> H_defined;
  f >> Name("G_defined") >> G_defined;
  f >> Name("nvar") >> nvar;
  f >> Name("ncheck") >> ncheck;
  f >> Name("C") >> C;
  f >> Name("V") >> V;
  f >> Name("sumX1") >> sumX1;
  f >> Name("sumX2") >> sumX2;
  f >> Name("iind") >> iind;
  f >> Name("jind") >> jind;
  f.close();

  // load generator data
  if (G_defined) {
    it_assert(G_in != 0, "LDPC_Code::load_code(): Generator object is "
              "missing. Can not load the generator data from a file.");
    G = G_in;
    G->load(filename);
  }
  else {
    G = 0;
    it_info_debug("LDPC_Code::load_code(): Generator data not loaded. "
                  "Generator object will not be used.");
  }

  it_info_debug("LDPC_Code::load_code(): Successfully loaded LDPC codec "
                "from " << filename);

  setup_decoder();
}

void LDPC_Code::save_code(const std::string& filename) const
{
  it_assert(H_defined, "LDPC_Code::save_to_file(): There is no parity "
            "check matrix");
  it_info_debug("LDPC_Code::save_to_file(): Saving LDPC codec to "
                << filename);

  it_file f;
  f.open(filename, true);
  f << Name("Fileversion") << LDPC_binary_file_version;
  f << Name("H_defined") << H_defined;
  f << Name("G_defined") << G_defined;
  f << Name("nvar") << nvar;
  f << Name("ncheck") << ncheck;
  f << Name("C") << C;
  f << Name("V") << V;
  f << Name("sumX1") << sumX1;
  f << Name("sumX2") << sumX2;
  f << Name("iind") << iind;
  f << Name("jind") << jind;
  f.close();

  // save generator data;
  if (G_defined)
    G->save(filename);
  else
    it_info_debug("LDPC_Code::save_code(): Missing generator object - "
                  "generator data not saved");

  it_info_debug("LDPC_Code::save_code(): Successfully saved LDPC codec to "
                << filename);
}


void LDPC_Code::set_decoding_method(const std::string& method_in)
{
  it_assert((method_in == "bp") || (method_in == "BP"),
            "LDPC_Code::set_decoding_method(): Not implemented decoding method");
  *dec_method = method_in;
}

void LDPC_Code::set_exit_conditions(int max_iters_in,
                                    bool syndr_check_each_iter,
                                    bool syndr_check_at_start)
{
  it_assert(max_iters >= 0, "LDPC_Code::set_nrof_iterations(): Maximum "
            "number of iterations can not be negative");
  max_iters = max_iters_in;
  psc = syndr_check_each_iter;
  pisc = syndr_check_at_start;
}

void LDPC_Code::set_llrcalc(const LLR_calc_unit& llrcalc_in)
{
  llrcalc = llrcalc_in;
}


void LDPC_Code::encode(const bvec &input, bvec &output)
{
  it_assert(G_defined, "LDPC_Code::encode(): LDPC Generator is required "
            "for encoding");
  G->encode(input, output);
  it_assert_debug(syndrome_check(output), "LDPC_Code::encode(): Syndrome "
                  "check failed");
}

bvec LDPC_Code::encode(const bvec &input)
{
  bvec result;
  encode(input, result);
  return result;
}

void LDPC_Code::decode(const vec &llr_in, bvec &syst_bits)
{
  QLLRvec qllrin = llrcalc.to_qllr(llr_in);
  QLLRvec qllrout;
  bp_decode(qllrin, qllrout);
  syst_bits = (qllrout.left(nvar - ncheck) < 0);
}

bvec LDPC_Code::decode(const vec &llr_in)
{
  bvec syst_bits;
  decode(llr_in, syst_bits);
  return syst_bits;
}

void LDPC_Code::decode_soft_out(const vec &llr_in, vec &llr_out)
{
  QLLRvec qllrin = llrcalc.to_qllr(llr_in);
  QLLRvec qllrout;
  bp_decode(qllrin, qllrout);
  llr_out = llrcalc.to_double(qllrout);
}

vec LDPC_Code::decode_soft_out(const vec &llr_in)
{
  vec llr_out;
  decode_soft_out(llr_in, llr_out);
  return llr_out;
}

int LDPC_Code::bp_decode(const QLLRvec &LLRin, QLLRvec &LLRout)
{
  // Note the IT++ convention that a sure zero corresponds to
  // LLR=+infinity
  it_assert(H_defined, "LDPC_Code::bp_decode(): Parity check matrix not "
            "defined");
  it_assert((LLRin.size() == nvar) && (sumX1.size() == nvar)
            && (sumX2.size() == ncheck), "LDPC_Code::bp_decode(): Wrong "
            "input dimensions");

  if (pisc && syndrome_check(LLRin)) {
    LLRout = LLRin;
    return 0;
  }

  LLRout.set_size(LLRin.size());

  // allocate temporary variables used for the check node update
  ivec jj(max_cnd);
  QLLRvec m(max_cnd);
  QLLRvec ml(max_cnd);
  QLLRvec mr(max_cnd);
  
  // initial step
  for (int i = 0; i < nvar; i++) {
    int index = i;
    for (int j = 0; j < sumX1(i); j++) {
      mvc[index] = LLRin(i);
      index += nvar;
    }
  }

  bool is_valid_codeword = false;
  int iter = 0;
  do {
    iter++;
    if (nvar >= 100000) { it_info_no_endl_debug("."); }
    // --------- Step 1: check to variable nodes ----------
    for (int j = 0; j < ncheck; j++) {
      // The check node update calculations are hardcoded for degrees
      // up to 6.  For larger degrees, a general algorithm is used.
      switch (sumX2(j)) {
      case 0:
        it_error("LDPC_Code::bp_decode(): sumX2(j)=0");
      case 1:
        it_error("LDPC_Code::bp_decode(): sumX2(j)=1");
      case 2: {
        mcv[j+ncheck] = mvc[jind[j]];
        mcv[j] = mvc[jind[j+ncheck]];
        break;
      }
      case 3: {
        int j0 = j;
        QLLR m0 = mvc[jind[j0]];
        int j1 = j0 + ncheck;
        QLLR m1 = mvc[jind[j1]];
        int j2 = j1 + ncheck;
        QLLR m2 = mvc[jind[j2]];
        mcv[j0] = llrcalc.Boxplus(m1, m2);
        mcv[j1] = llrcalc.Boxplus(m0, m2);
        mcv[j2] = llrcalc.Boxplus(m0, m1);
        break;
      }
      case 4: {
        int j0 = j;
        QLLR m0 = mvc[jind[j0]];
        int j1 = j0 + ncheck;
        QLLR m1 = mvc[jind[j1]];
        int j2 = j1 + ncheck;
        QLLR m2 = mvc[jind[j2]];
        int j3 = j2 + ncheck;
        QLLR m3 = mvc[jind[j3]];
        QLLR m01 = llrcalc.Boxplus(m0, m1);
        QLLR m23 = llrcalc.Boxplus(m2, m3);
        mcv[j0] = llrcalc.Boxplus(m1, m23);
        mcv[j1] = llrcalc.Boxplus(m0, m23);
        mcv[j2] = llrcalc.Boxplus(m01, m3);
        mcv[j3] = llrcalc.Boxplus(m01, m2);
        break;
      }
      case 5: {
        int j0 = j;
        QLLR m0 = mvc[jind[j0]];
        int j1 = j0 + ncheck;
        QLLR m1 = mvc[jind[j1]];
        int j2 = j1 + ncheck;
        QLLR m2 = mvc[jind[j2]];
        int j3 = j2 + ncheck;
        QLLR m3 = mvc[jind[j3]];
        int j4 = j3 + ncheck;
        QLLR m4 = mvc[jind[j4]];
        QLLR m01 = llrcalc.Boxplus(m0, m1);
        QLLR m02 = llrcalc.Boxplus(m01, m2);
        QLLR m34 = llrcalc.Boxplus(m3, m4);
        QLLR m24 = llrcalc.Boxplus(m2, m34);
        mcv[j0] = llrcalc.Boxplus(m1, m24);
        mcv[j1] = llrcalc.Boxplus(m0, m24);
        mcv[j2] = llrcalc.Boxplus(m01, m34);
        mcv[j3] = llrcalc.Boxplus(m02, m4);
        mcv[j4] = llrcalc.Boxplus(m02, m3);
        break;
      }
      case 6: {
        int j0 = j;
        QLLR m0 = mvc[jind[j0]];
        int j1 = j0 + ncheck;
        QLLR m1 = mvc[jind[j1]];
        int j2 = j1 + ncheck;
        QLLR m2 = mvc[jind[j2]];
        int j3 = j2 + ncheck;
        QLLR m3 = mvc[jind[j3]];
        int j4 = j3 + ncheck;
        QLLR m4 = mvc[jind[j4]];
        int j5 = j4 + ncheck;
        QLLR m5 = mvc[jind[j5]];
        QLLR m01 = llrcalc.Boxplus(m0, m1);
        QLLR m23 = llrcalc.Boxplus(m2, m3);
        QLLR m45 = llrcalc.Boxplus(m4, m5);
        QLLR m03 = llrcalc.Boxplus(m01, m23);
        QLLR m25 = llrcalc.Boxplus(m23, m45);
        QLLR m0145 = llrcalc.Boxplus(m01, m45);
        mcv[j0] = llrcalc.Boxplus(m1, m25);
        mcv[j1] = llrcalc.Boxplus(m0, m25);
        mcv[j2] = llrcalc.Boxplus(m0145, m3);
        mcv[j3] = llrcalc.Boxplus(m0145, m2);
        mcv[j4] = llrcalc.Boxplus(m03, m5);
        mcv[j5] = llrcalc.Boxplus(m03, m4);
        break;
      }
      default: {
        int nodes = sumX2(j);
        if( nodes > max_cnd ) {
          std::ostringstream m_sout;
          m_sout << "check node degrees >" << max_cnd << " not supported in this version";
          it_error( m_sout.str() );
        }

        nodes--;
        jj[0] = j;
        m[0] = mvc[jind[jj[0]]];
        for(int i = 1; i <= nodes; i++ ) {
          jj[i] = jj[i-1] + ncheck;
          m[i] = mvc[jind[jj[i]]];
        }

	// compute partial sums from the left and from the right
        ml[0] = m[0];
        mr[0] = m[nodes];
        for(int i = 1; i < nodes; i++ ) {
          ml[i] = llrcalc.Boxplus( ml[i-1], m[i] );
          mr[i] = llrcalc.Boxplus( mr[i-1], m[nodes-i] );
        }

	// merge partial sums
        mcv[jj[0]] = mr[nodes-1];
        mcv[jj[nodes]] = ml[nodes-1];
        for(int i = 1; i < nodes; i++ )
          mcv[jj[i]] = llrcalc.Boxplus( ml[i-1], mr[nodes-1-i] );
      }
      }  // switch statement
    }
    
    // step 2: variable to check nodes
    for (int i = 0; i < nvar; i++) {
      switch (sumX1(i)) {
      case 0:
        it_error("LDPC_Code::bp_decode(): sumX1(i)=0");
      case 1: {
        /* This case is rare but apparently occurs for codes used in
	   the DVB-T2 standard.
	 */
	QLLR m0 = mcv[iind[i]];
        mvc[i] = LLRin(i);
        LLRout(i) = LLRin(i) + m0;
        break;
      }
      case 2: {
        QLLR m0 = mcv[iind[i]];
        int i1 = i + nvar;
        QLLR m1 = mcv[iind[i1]];
        mvc[i] = LLRin(i) + m1;
        mvc[i1] = LLRin(i) + m0;
        LLRout(i) = mvc[i1] + m1;
        break;
      }
      case 3: {
        int i0 = i;
        QLLR m0 = mcv[iind[i0]];
        int i1 = i0 + nvar;
        QLLR m1 = mcv[iind[i1]];
        int i2 = i1 + nvar;
        QLLR m2 = mcv[iind[i2]];
        LLRout(i) = LLRin(i) + m0 + m1 + m2;
        mvc[i0] = LLRout(i) - m0;
        mvc[i1] = LLRout(i) - m1;
        mvc[i2] = LLRout(i) - m2;
        break;
      }
      case 4: {
        int i0 = i;
        QLLR m0 = mcv[iind[i0]];
        int i1 = i0 + nvar;
        QLLR m1 = mcv[iind[i1]];
        int i2 = i1 + nvar;
        QLLR m2 = mcv[iind[i2]];
        int i3 = i2 + nvar;
        QLLR m3 = mcv[iind[i3]];
        LLRout(i) = LLRin(i) + m0 + m1 + m2 + m3;
        mvc[i0] = LLRout(i) - m0;
        mvc[i1] = LLRout(i) - m1;
        mvc[i2] = LLRout(i) - m2;
        mvc[i3] = LLRout(i) - m3;
        break;
      }
      default:   { // differential update
        QLLR mvc_temp = LLRin(i);
        int index_iind = i; // tracks i+jp*nvar
        for (int jp = 0; jp < sumX1(i); jp++) {
          mvc_temp +=  mcv[iind[index_iind]];
          index_iind += nvar;
        }
        LLRout(i) = mvc_temp;
        index_iind = i;  // tracks i+j*nvar
        for (int j = 0; j < sumX1[i]; j++) {
          mvc[index_iind] = mvc_temp - mcv[iind[index_iind]];
          index_iind += nvar;
        }
      }
      }
    }

    if (psc && syndrome_check(LLRout)) {
      is_valid_codeword = true;
      break;
    }
  }
  while (iter < max_iters);

  if (nvar >= 100000) { it_info_debug(""); }
  return (is_valid_codeword ? iter : -iter);
}


bool LDPC_Code::syndrome_check(const bvec &x) const
{
  QLLRvec llr = 1 - 2 * to_ivec(x);
  return syndrome_check(llr);
}

bool LDPC_Code::syndrome_check(const QLLRvec &LLR) const
{
  // Please note the IT++ convention that a sure zero corresponds to
  // LLR=+infinity
  int i, j, synd, vi;

  for (j = 0; j < ncheck; j++) {
    synd = 0;
    int vind = j; // tracks j+i*ncheck
    for (i = 0; i < sumX2(j); i++) {
      vi = V(vind);
      if (LLR(vi) < 0) {
        synd++;
      }
      vind += ncheck;
    }
    if ((synd&1) == 1) {
      return false;  // codeword is invalid
    }
  }
  return true;   // codeword is valid
};

QLLRvec LDPC_Code::soft_syndrome_check(const QLLRvec &LLR) const
{
  QLLRvec result(ncheck);
  int i,j,vi,vind;

  for (j=0; j<ncheck; j++) {
    result(j)=-QLLR_MAX;
    vind=j;
    for (i=0; i<sumX2(j); i++) {
      vi=V(vind);
      result(j) = llrcalc.Boxplus(LLR(vi),result(j));
      vind += ncheck;
    }
  }

  return result;
}

// ----------------------------------------------------------------------
// LDPC_Code private methods
// ----------------------------------------------------------------------

void LDPC_Code::decoder_parameterization(const LDPC_Parity* const Hmat)
{
  // copy basic parameters
  sumX1 = Hmat->sumX1;
  sumX2 = Hmat->sumX2;
  nvar = Hmat->nvar; //get_nvar();
  ncheck = Hmat->ncheck; //get_ncheck();
  int cmax = max(sumX1);
  int vmax = max(sumX2);

  // decoder parameterization
  V = zeros_i(ncheck * vmax);
  C = zeros_i(cmax * nvar);
  jind = zeros_i(ncheck * vmax);
  iind = zeros_i(nvar * cmax);

  it_info_debug("LDPC_Code::decoder_parameterization(): Computations "
                "- phase 1");
  for (int i = 0; i < nvar; i++) {
    ivec coli = Hmat->get_col(i).get_nz_indices();
    for (int j0 = 0; j0 < length(coli); j0++) {
      C(j0 + cmax*i) = coli(j0);
    }
  }

  it_info_debug("LDPC_Code::decoder_parameterization(): Computations "
                "- phase 2");
  it_info_debug("Computing decoder parameterization. Phase 2");
  for (int j = 0; j < ncheck; j++) {
    ivec rowj = Hmat->get_row(j).get_nz_indices();
    for (int i0 = 0; i0 < length(rowj); i0++) {
      V(j + ncheck*i0) = rowj(i0);
    }
  }

  it_info_debug("LDPC_Code::decoder_parameterization(): Computations "
                "- phase 3");
  it_info_debug("Computing decoder parameterization. Phase 3.");
  for (int j = 0; j < ncheck; j++) {
    for (int ip = 0; ip < sumX2(j); ip++) {
      int vip = V(j + ip * ncheck);
      int k = 0;
      while (1 == 1)  {
        if (C(k + vip*cmax) == j) {
          break;
        }
        k++;
      }
      jind(j + ip*ncheck) = vip + k * nvar;
    }
  }

  it_info_debug("LDPC_Code::decoder_parameterization(): Computations "
                "- phase 4");
  for (int i = 0; i < nvar; i++) {
    for (int jp = 0; jp < sumX1(i); jp++) {
      int cjp = C(jp + i * cmax);
      int k = 0;
      while (1 == 1) {
        if (V(cjp + k*ncheck) == i) {break; }
        k++;
      }
      iind(i + jp*nvar) = cjp + k * ncheck;
    }
  }

  H_defined = true;
}


void LDPC_Code::setup_decoder()
{
  if (H_defined) {
    mcv.set_size(max(sumX2) * ncheck);
    mvc.set_size(max(sumX1) * nvar);
  }
}


void LDPC_Code::integrity_check()
{
  if (G_defined) {
    it_info_debug("LDPC_Code::integrity_check(): Checking integrity of "
                  "the LDPC_Parity and LDPC_Generator data");
    bvec bv(nvar - ncheck), cw;
    bv.clear();
    bv(0) = 1;
    for (int i = 0; i < nvar - ncheck; i++) {
      G->encode(bv, cw);
      it_assert(syndrome_check(cw),
                "LDPC_Code::integrity_check(): Syndrome check failed");
      bv.shift_right(bv(nvar - ncheck - 1));
    }
  }
  else {
    it_info_debug("LDPC_Code::integrity_check(): No generator defined "
                  "- no check performed");
  }
}

// ----------------------------------------------------------------------
// Related functions
// ----------------------------------------------------------------------

std::ostream &operator<<(std::ostream &os, const LDPC_Code &C)
{
  ivec rdeg = zeros_i(max(C.sumX2) + 1);
  for (int i = 0; i < C.ncheck; i++)     {
    rdeg(C.sumX2(i))++;
  }

  ivec cdeg = zeros_i(max(C.sumX1) + 1);
  for (int j = 0; j < C.nvar; j++)     {
    cdeg(C.sumX1(j))++;
  }

  os << "--- LDPC codec ----------------------------------\n"
  << "Nvar : " << C.get_nvar() << "\n"
  << "Ncheck : " << C.get_ncheck() << "\n"
  << "Rate : " << C.get_rate() << "\n"
  << "Column degrees (node perspective): " << cdeg << "\n"
  << "Row degrees (node perspective): " << rdeg << "\n"
  << "-------------------------------------------------\n"
  << "Decoder parameters:\n"
  << " - method : " << C.get_decoding_method() << "\n"
  << " - max. iterations : " << C.max_iters << "\n"
  << " - syndrome check at each iteration : " << C.psc << "\n"
  << " - syndrome check at start : " << C.pisc << "\n"
  << "-------------------------------------------------\n"
  << C.llrcalc << "\n";
  return os;
}

} // namespace itpp
