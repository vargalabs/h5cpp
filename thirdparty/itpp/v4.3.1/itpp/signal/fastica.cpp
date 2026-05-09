/*!
 * \file
 * \brief Implementation of FastICA (Independent Component Analysis) for IT++
 * \author Francois Cayre and Teddy Furon
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
 *
 * This is IT++ implementation of the original Matlab package FastICA.
 *
 * This code is Copyright (C) 2004 by:
 *   Francois CAYRE and Teddy FURON
 *   TEMICS Project
 *   INRIA/Rennes (IRISA)
 *   Campus Universitaire de Beaulieu
 *   35042 RENNES cedex FRANCE
 *
 * Email : firstname.lastname@irisa.fr
 *
 * Matlab package is Copyright (C) 1998 by:
 *   Jarmo HURRI, Hugo GAVERT, Jaakko SARELA and Aapo HYVARINEN
 *   Laboratory of Information and Computer Science
 *   Helsinki University of Technology              *
 *
 * URL : http://www.cis.hut.fi/projects/ica/fastica/about.shtml
 *
 * If you use results given by this FastICA software in an article for
 * a scientific journal, conference proceedings or similar, please
 * include the following original reference in the bibliography :
 *
 *   A. Hyvarinen, Fast and Robust Fixed-Point Algorithms for
 *   Independent Component Analysis, IEEE Transactions on Neural
 *   Networks 10(3):626-634, 1999
 *
 * Differences with the original Matlab implementation:
 * - no GUI
 * - return something even in the case of a convergence problem
 * - optimization of SVD decomposition (performed 2 times in Matlab,
 *   only 1 time in IT++)
 * - default approach is SYMM with non-linearity POW3
 */

#include <itpp/signal/fastica.h>
#include <itpp/signal/sigfun.h>
#include <itpp/signal/resampling.h>
#include <itpp/base/algebra/eigen.h>
#include <itpp/base/algebra/svd.h>
#include <itpp/base/math/trig_hyp.h>
#include <itpp/base/matfunc.h>
#include <itpp/base/random.h>
#include <itpp/base/sort.h>
#include <itpp/base/specmat.h>
#include <itpp/base/svec.h>
#include <itpp/base/math/min_max.h>
#include <itpp/stat/misc_stat.h>


using namespace itpp;


/*!
  \brief Local functions for FastICA
  @{
*/
static void selcol(const mat oldMatrix, const vec maskVector, mat & newMatrix);
static int pcamat(const mat vectors, const int numOfIC, int firstEig, int lastEig, mat & Es, vec & Ds);
static void remmean(mat inVectors, mat & outVectors, vec & meanValue);
static void whitenv(const mat vectors, const mat E, const mat D, mat & newVectors, mat & whiteningMatrix, mat & dewhiteningMatrix);
static mat orth(const mat A);
static mat mpower(const mat A, const double y);
static ivec getSamples(const int max, const double percentage);
static vec sumcol(const mat A);
static bool fpica(const mat X, const mat whiteningMatrix, const mat dewhiteningMatrix, const int approach, const int numOfIC, const int g, const int finetune, const double a1, const double a2, double myy, const int stabilization, const double epsilon, const int maxNumIterations, const int maxFinetune, const int initState, mat guess, double sampleSize, mat & A, mat & W);
/*! @} */

namespace itpp
{

// Constructor, init default values
Fast_ICA::Fast_ICA(mat ma_mixedSig)
{

  // Init default values
  approach = FICA_APPROACH_SYMM;
  g = FICA_NONLIN_POW3;
  finetune = true;
  a1 = 1.0;
  a2 = 1.0;
  mu = 1.0;
  epsilon = 0.0001;
  sampleSize = 1.0;
  stabilization = false;
  maxNumIterations = 100000;
  maxFineTune = 100;
  firstEig = 1;

  mixedSig = ma_mixedSig;

  lastEig = mixedSig.rows();
  numOfIC = mixedSig.rows();
  PCAonly = false;
  initState = FICA_INIT_RAND;

}

// Call main function
bool Fast_ICA::separate(void)
{

  int Dim = numOfIC;

  mat mixedSigC;
  vec mixedMean;

  mat guess;
  if (initState == FICA_INIT_RAND)
    guess = zeros(Dim, Dim);
  else
    guess = mat(initGuess);

  VecPr = zeros(mixedSig.rows(), numOfIC);

  icasig = zeros(numOfIC, mixedSig.cols());

  remmean(mixedSig, mixedSigC, mixedMean);

  if (pcamat(mixedSigC, numOfIC, firstEig, lastEig, E, D) < 1) {
    // no principal components could be found (e.g. all-zero data): return the unchanged input
    icasig = mixedSig;
    return false;
  }

  whitenv(mixedSigC, E, diag(D), whitesig, whiteningMatrix, dewhiteningMatrix);

  Dim = whitesig.rows();
  if (numOfIC > Dim) numOfIC = Dim;

  ivec NcFirst = to_ivec(zeros(numOfIC));
  vec NcVp = D;
  for (int i = 0; i < NcFirst.size(); i++) {

    NcFirst(i) = max_index(NcVp);
    NcVp(NcFirst(i)) = 0.0;
    VecPr.set_col(i, dewhiteningMatrix.get_col(i));

  }

  bool result = true;
  if (PCAonly == false) {

    result = fpica(whitesig, whiteningMatrix, dewhiteningMatrix, approach, numOfIC, g, finetune, a1, a2, mu, stabilization, epsilon, maxNumIterations, maxFineTune, initState, guess, sampleSize, A, W);

    icasig = W * mixedSig;

  }

  else { // PCA only : returns E as IcaSig
    icasig = VecPr;
  }
  return result;
}

void Fast_ICA::set_approach(int in_approach) { approach = in_approach; if (approach == FICA_APPROACH_DEFL) finetune = true; }

void Fast_ICA::set_nrof_independent_components(int in_nrIC) { numOfIC = in_nrIC; }

void Fast_ICA::set_non_linearity(int in_g) { g = in_g; }

void Fast_ICA::set_fine_tune(bool in_finetune) { finetune = in_finetune; }

void Fast_ICA::set_a1(double fl_a1) { a1 = fl_a1; }

void Fast_ICA::set_a2(double fl_a2) { a2 = fl_a2; }

void Fast_ICA::set_mu(double fl_mu) { mu = fl_mu; }

void Fast_ICA::set_epsilon(double fl_epsilon) { epsilon = fl_epsilon; }

void Fast_ICA::set_sample_size(double fl_sampleSize) { sampleSize = fl_sampleSize; }

void Fast_ICA::set_stabilization(bool in_stabilization) { stabilization = in_stabilization; }

void Fast_ICA::set_max_num_iterations(int in_maxNumIterations) { maxNumIterations = in_maxNumIterations; }

void Fast_ICA::set_max_fine_tune(int in_maxFineTune) { maxFineTune = in_maxFineTune; }

void Fast_ICA::set_first_eig(int in_firstEig) { firstEig = in_firstEig; }

void Fast_ICA::set_last_eig(int in_lastEig) { lastEig = in_lastEig; }

void Fast_ICA::set_pca_only(bool in_PCAonly) { PCAonly = in_PCAonly; }

void Fast_ICA::set_init_guess(mat ma_initGuess)
{
  initGuess = ma_initGuess;
  initState = FICA_INIT_GUESS;
}

mat Fast_ICA::get_mixing_matrix() { if (PCAonly) { it_warning("No ICA performed."); return (zeros(1, 1));} else return A; }

mat Fast_ICA::get_separating_matrix() { if (PCAonly) { it_warning("No ICA performed."); return(zeros(1, 1)); } else return W; }

mat Fast_ICA::get_independent_components() { if (PCAonly) { it_warning("No ICA performed."); return(zeros(1, 1)); } else return icasig; }

int Fast_ICA::get_nrof_independent_components() { return numOfIC; }

mat Fast_ICA::get_principal_eigenvectors() { return VecPr; }

mat Fast_ICA::get_whitening_matrix() { return whiteningMatrix; }

mat Fast_ICA::get_dewhitening_matrix() { return dewhiteningMatrix; }

mat Fast_ICA::get_white_sig() { return whitesig; }

} // namespace itpp


static void selcol(const mat oldMatrix, const vec maskVector, mat & newMatrix)
{

  int numTaken = 0;

  for (int i = 0; i < size(maskVector); i++) if (maskVector(i) == 1) numTaken++;

  newMatrix = zeros(oldMatrix.rows(), numTaken);

  numTaken = 0;

  for (int i = 0; i < size(maskVector); i++) {

    if (maskVector(i) == 1) {

      newMatrix.set_col(numTaken, oldMatrix.get_col(i));
      numTaken++;

    }
  }

  return;

}

static int pcamat(const mat vectors, const int numOfIC, int firstEig, int lastEig, mat & Es, vec & Ds)
{

  mat Et;
  vec Dt;
  cmat Ec;
  cvec Dc;
  double lowerLimitValue = 0.0,
                           higherLimitValue = 0.0;

  int oldDimension = vectors.rows();

  mat covarianceMatrix = cov(transpose(vectors), 0);

  eig_sym(covarianceMatrix, Dt, Et);

  int maxLastEig = 0;

  // Compute rank
  for (int i = 0; i < Dt.length(); i++) if (Dt(i) > FICA_TOL) maxLastEig++;
  if (maxLastEig < 1) return 0;

  // Force numOfIC components
  if (maxLastEig > numOfIC) maxLastEig = numOfIC;

  vec eigenvalues = zeros(size(Dt));
  vec eigenvalues2 = zeros(size(Dt));

  eigenvalues2 = Dt;

  sort(eigenvalues2);

  vec lowerColumns = zeros(size(Dt));

  for (int i = 0; i < size(Dt); i++) eigenvalues(i) = eigenvalues2(size(Dt) - i - 1);

  if (lastEig > maxLastEig) lastEig = maxLastEig;

  if (lastEig < oldDimension) lowerLimitValue = (eigenvalues(lastEig - 1) + eigenvalues(lastEig)) / 2;
  else lowerLimitValue = eigenvalues(oldDimension - 1) - 1;

  for (int i = 0; i < size(Dt); i++) if (Dt(i) > lowerLimitValue) lowerColumns(i) = 1;

  if (firstEig > 1) higherLimitValue = (eigenvalues(firstEig - 2) + eigenvalues(firstEig - 1)) / 2;
  else higherLimitValue = eigenvalues(0) + 1;

  vec higherColumns = zeros(size(Dt));
  for (int i = 0; i < size(Dt); i++) if (Dt(i) < higherLimitValue) higherColumns(i) = 1;

  vec selectedColumns = zeros(size(Dt));
  for (int i = 0; i < size(Dt); i++) selectedColumns(i) = (lowerColumns(i) == 1 && higherColumns(i) == 1) ? 1 : 0;

  selcol(Et, selectedColumns, Es);

  int numTaken = 0;

  for (int i = 0; i < size(selectedColumns); i++) if (selectedColumns(i) == 1) numTaken++;

  Ds = zeros(numTaken);

  numTaken = 0;

  for (int i = 0; i < size(Dt); i++)
    if (selectedColumns(i) == 1) {
      Ds(numTaken) = Dt(i);
      numTaken++;
    }

  return lastEig;

}


static void remmean(mat inVectors, mat & outVectors, vec & meanValue)
{

  outVectors = zeros(inVectors.rows(), inVectors.cols());
  meanValue = zeros(inVectors.rows());

  for (int i = 0; i < inVectors.rows(); i++) {

    meanValue(i) = mean(inVectors.get_row(i));

    for (int j = 0; j < inVectors.cols(); j++) outVectors(i, j) = inVectors(i, j) - meanValue(i);

  }

}

static void whitenv(const mat vectors, const mat E, const mat D, mat & newVectors, mat & whiteningMatrix, mat & dewhiteningMatrix)
{

  whiteningMatrix = zeros(E.cols(), E.rows());
  dewhiteningMatrix = zeros(E.rows(), E.cols());

  for (int i = 0; i < D.cols(); i++) {
    whiteningMatrix.set_row(i, std::pow(std::sqrt(D(i, i)), -1)*E.get_col(i));
    dewhiteningMatrix.set_col(i, std::sqrt(D(i, i))*E.get_col(i));
  }

  newVectors = whiteningMatrix * vectors;

  return;

}

static mat orth(const mat A)
{

  mat Q;
  mat U, V;
  vec S;
  double eps = 2.2e-16;
  double tol = 0.0;
  int mmax = 0;
  int r = 0;

  svd(A, U, S, V);
  if (A.rows() > A.cols()) {

    U = U(0, U.rows() - 1, 0, A.cols() - 1);
    S = S(0, A.cols() - 1);
  }

  mmax = (A.rows() > A.cols()) ? A.rows() : A.cols();

  tol = mmax * eps * max(S);

  for (int i = 0; i < size(S); i++) if (S(i) > tol) r++;

  Q = U(0, U.rows() - 1, 0, r - 1);

  return (Q);
}

static mat mpower(const mat A, const double y)
{

  mat T = zeros(A.rows(), A.cols());
  mat dd = zeros(A.rows(), A.cols());
  vec d = zeros(A.rows());
  vec dOut = zeros(A.rows());

  eig_sym(A, d, T);

  dOut = pow(d, y);

  diag(dOut, dd);

  for (int i = 0; i < T.cols(); i++) T.set_col(i, T.get_col(i) / norm(T.get_col(i)));

  return (T*dd*transpose(T));

}

static ivec getSamples(const int max, const double percentage)
{

  vec rd = randu(max);
  sparse_vec sV;
  ivec out;
  int sZ = 0;

  for (int i = 0; i < max; i++) if (rd(i) < percentage) { sV.add_elem(sZ, i); sZ++; }

  out = to_ivec(full(sV));

  return (out);

}

static vec sumcol(const mat A)
{

  vec out = zeros(A.cols());

  for (int i = 0; i < A.cols(); i++) { out(i) = sum(A.get_col(i)); }

  return (out);

}

static bool fpica(const mat X, const mat whiteningMatrix, const mat dewhiteningMatrix, const int approach, const int numOfIC, const int g, const int finetune, const double a1, const double a2, double myy, const int stabilization, const double epsilon, const int maxNumIterations, const int maxFinetune, const int initState, mat guess, double sampleSize, mat & A, mat & W)
{

  int vectorSize = X.rows();
  int numSamples = X.cols();
  int gOrig = g;
  int gFine = finetune + 1;
  double myyOrig = myy;
  double myyK = 0.01;
  int failureLimit = 5;
  int usedNlinearity = 0;
  double stroke = 0.0;
  int notFine = 1;
  int loong = 0;
  int initialStateMode = initState;
  double minAbsCos = 0.0, minAbsCos2 = 0.0;

  if (sampleSize * numSamples < 1000) sampleSize = (1000 / (double)numSamples < 1.0) ? 1000 / (double)numSamples : 1.0;

  if (sampleSize != 1.0) gOrig += 2;
  if (myy != 1.0) gOrig += 1;

  int fineTuningEnabled = 1;

  if (!finetune) {
    if (myy != 1.0) gFine = gOrig;
    else gFine = gOrig + 1;
    fineTuningEnabled = 0;
  }

  int stabilizationEnabled = stabilization;

  if (!stabilization && myy != 1.0) stabilizationEnabled = true;

  usedNlinearity = gOrig;

  if (initState == FICA_INIT_GUESS && guess.rows() != whiteningMatrix.cols()) {
    initialStateMode = 0;

  }
  else if (guess.cols() < numOfIC) {

    mat guess2 = randu(guess.rows(), numOfIC - guess.cols()) - 0.5;
    guess = concat_horizontal(guess, guess2);
  }
  else if (guess.cols() > numOfIC) guess = guess(0, guess.rows() - 1, 0, numOfIC - 1);

  if (approach == FICA_APPROACH_SYMM) {

    usedNlinearity = gOrig;
    stroke = 0;
    notFine = 1;
    loong = 0;

    A = zeros(vectorSize, numOfIC);
    mat B = zeros(vectorSize, numOfIC);

    if (initialStateMode == 0) B = orth(randu(vectorSize, numOfIC) - 0.5);
    else B = whiteningMatrix * guess;

    mat BOld = zeros(B.rows(), B.cols());
    mat BOld2 = zeros(B.rows(), B.cols());

    for (int round = 0; round < maxNumIterations; round++) {

      if (round == maxNumIterations - 1) {

        // If there is a convergence problem,
        // we still want ot return something.
        // This is difference with original
        // Matlab implementation.
        A = dewhiteningMatrix * B;
        W = transpose(B) * whiteningMatrix;

        return false;
      }

      B = B * mpower(transpose(B) * B , -0.5);

      minAbsCos = min(abs(diag(transpose(B) * BOld)));
      minAbsCos2 = min(abs(diag(transpose(B) * BOld2)));

      if (1 - minAbsCos < epsilon) {

        if (fineTuningEnabled && notFine) {

          notFine = 0;
          usedNlinearity = gFine;
          myy = myyK * myyOrig;
          BOld = zeros(B.rows(), B.cols());
          BOld2 = zeros(B.rows(), B.cols());

        }

        else {

          A = dewhiteningMatrix * B;
          break;

        }

      } // IF epsilon

      else if (stabilizationEnabled) {
        if (!stroke && (1 - minAbsCos2 < epsilon)) {

          stroke = myy;
          myy /= 2;
          if (mod(usedNlinearity, 2) == 0) usedNlinearity += 1 ;

        }
        else if (stroke) {

          myy = stroke;
          stroke = 0;
          if (myy == 1 && mod(usedNlinearity, 2) != 0) usedNlinearity -= 1;

        }
        else if (!loong && (round > maxNumIterations / 2)) {

          loong = 1;
          myy /= 2;
          if (mod(usedNlinearity, 2) == 0) usedNlinearity += 1;

        }

      } // stabilizationEnabled

      BOld2 = BOld;
      BOld = B;

      switch (usedNlinearity) {

        // pow3
      case FICA_NONLIN_POW3 : {
        B = (X * pow(transpose(X) * B, 3)) / numSamples - 3 * B;
        break;
      }
      case(FICA_NONLIN_POW3+1) : {
        mat Y = transpose(X) * B;
        mat Gpow3 = pow(Y, 3);
        vec Beta = sumcol(pow(Y, 4));
        mat D = diag(pow(Beta - 3 * numSamples , -1));
        B = B + myy * B * (transpose(Y) * Gpow3 - diag(Beta)) * D;
        break;
      }
      case(FICA_NONLIN_POW3+2) : {
        mat Xsub = X.get_cols(getSamples(numSamples, sampleSize));
        B = (Xsub * pow(transpose(Xsub) * B, 3)) / Xsub.cols() - 3 * B;
        break;
      }
      case(FICA_NONLIN_POW3+3) : {
        mat Ysub = transpose(X.get_cols(getSamples(numSamples, sampleSize))) * B;
        mat Gpow3 = pow(Ysub, 3);
        vec Beta = sumcol(pow(Ysub, 4));
        mat D = diag(pow(Beta - 3 * Ysub.rows() , -1));
        B = B + myy * B * (transpose(Ysub) * Gpow3 - diag(Beta)) * D;
        break;
      }

      // TANH
      case FICA_NONLIN_TANH : {
        mat hypTan = tanh(a1 * transpose(X) * B);
        B = (X * hypTan) / numSamples - elem_mult(reshape(repeat(sumcol(1 - pow(hypTan, 2)), B.rows()), B.rows(), B.cols()), B) / numSamples * a1;
        break;
      }
      case(FICA_NONLIN_TANH+1) : {
        mat Y = transpose(X) * B;
        mat hypTan = tanh(a1 * Y);
        vec Beta = sumcol(elem_mult(Y, hypTan));
        vec Beta2 = sumcol(1 - pow(hypTan, 2));
        mat D = diag(pow(Beta - a1 * Beta2 , -1));
        B = B + myy * B * (transpose(Y) * hypTan - diag(Beta)) * D;
        break;
      }
      case(FICA_NONLIN_TANH+2) : {
        mat Xsub = X.get_cols(getSamples(numSamples, sampleSize));
        mat hypTan = tanh(a1 * transpose(Xsub) * B);
        B = (Xsub * hypTan) / Xsub.cols() -  elem_mult(reshape(repeat(sumcol(1 - pow(hypTan, 2)), B.rows()), B.rows(), B.cols()), B) / Xsub.cols() * a1;
        break;
      }
      case(FICA_NONLIN_TANH+3) : {
        mat Ysub = transpose(X.get_cols(getSamples(numSamples, sampleSize))) * B;
        mat hypTan = tanh(a1 * Ysub);
        vec Beta = sumcol(elem_mult(Ysub, hypTan));
        vec Beta2 = sumcol(1 - pow(hypTan, 2));
        mat D = diag(pow(Beta - a1 * Beta2 , -1));
        B = B + myy * B * (transpose(Ysub) * hypTan - diag(Beta)) * D;
        break;
      }

      // GAUSS
      case FICA_NONLIN_GAUSS : {
        mat U = transpose(X) * B;
        mat Usquared = pow(U, 2);
        mat ex = exp(-a2 * Usquared / 2);
        mat gauss = elem_mult(U, ex);
        mat dGauss = elem_mult(1 - a2 * Usquared, ex);
        B = (X * gauss) / numSamples - elem_mult(reshape(repeat(sumcol(dGauss), B.rows()), B.rows(), B.cols()), B) / numSamples;
        break;
      }
      case(FICA_NONLIN_GAUSS+1) : {
        mat Y = transpose(X) * B;
        mat ex = exp(-a2 * pow(Y, 2) / 2);
        mat gauss = elem_mult(Y, ex);
        vec Beta = sumcol(elem_mult(Y, gauss));
        vec Beta2 = sumcol(elem_mult(1 - a2 * pow(Y, 2), ex));
        mat D = diag(pow(Beta - Beta2 , -1));
        B = B + myy * B * (transpose(Y) * gauss - diag(Beta)) * D;
        break;
      }
      case(FICA_NONLIN_GAUSS+2) : {
        mat Xsub = X.get_cols(getSamples(numSamples, sampleSize));
        mat U = transpose(Xsub) * B;
        mat Usquared = pow(U, 2);
        mat ex = exp(-a2 * Usquared / 2);
        mat gauss = elem_mult(U, ex);
        mat dGauss = elem_mult(1 - a2 * Usquared, ex);
        B = (Xsub * gauss) / Xsub.cols() - elem_mult(reshape(repeat(sumcol(dGauss), B.rows()), B.rows(), B.cols()), B) / Xsub.cols();
        break;
      }
      case(FICA_NONLIN_GAUSS+3) : {
        mat Ysub = transpose(X.get_cols(getSamples(numSamples, sampleSize))) * B;
        mat ex = exp(-a2 * pow(Ysub, 2) / 2);
        mat gauss = elem_mult(Ysub, ex);
        vec Beta = sumcol(elem_mult(Ysub, gauss));
        vec Beta2 = sumcol(elem_mult(1 - a2 * pow(Ysub, 2), ex));
        mat D = diag(pow(Beta - Beta2 , -1));
        B = B + myy * B * (transpose(Ysub) * gauss - diag(Beta)) * D;
        break;
      }

      // SKEW
      case FICA_NONLIN_SKEW : {
        B = (X * (pow(transpose(X) * B, 2))) / numSamples;
        break;
      }
      case(FICA_NONLIN_SKEW+1) : {
        mat Y = transpose(X) * B;
        mat Gskew = pow(Y, 2);
        vec Beta = sumcol(elem_mult(Y, Gskew));
        mat D = diag(pow(Beta , -1));
        B = B + myy * B * (transpose(Y) * Gskew - diag(Beta)) * D;
        break;
      }
      case(FICA_NONLIN_SKEW+2) : {
        mat Xsub = X.get_cols(getSamples(numSamples, sampleSize));
        B = (Xsub * (pow(transpose(Xsub) * B, 2))) / Xsub.cols();
        break;
      }
      case(FICA_NONLIN_SKEW+3) : {
        mat Ysub = transpose(X.get_cols(getSamples(numSamples, sampleSize))) * B;
        mat Gskew = pow(Ysub, 2);
        vec Beta = sumcol(elem_mult(Ysub, Gskew));
        mat D = diag(pow(Beta , -1));
        B = B + myy * B * (transpose(Ysub) * Gskew - diag(Beta)) * D;
        break;
      }


      } // SWITCH usedNlinearity

    } // FOR maxIterations

    W = transpose(B) * whiteningMatrix;


  } // IF FICA_APPROACH_SYMM APPROACH

  // DEFLATION
  else {

    // FC 01/12/05
    A = zeros(whiteningMatrix.cols(), numOfIC);
    //    A = zeros( vectorSize, numOfIC );
    mat B = zeros(vectorSize, numOfIC);
    W = transpose(B) * whiteningMatrix;
    int round = 1;
    int numFailures = 0;

    while (round <= numOfIC) {

      myy = myyOrig;

      usedNlinearity = gOrig;
      stroke = 0;

      notFine = 1;
      loong = 0;
      int endFinetuning = 0;

      vec w = zeros(vectorSize);

      if (initialStateMode == 0)

        w = randu(vectorSize) - 0.5;

      else w = whiteningMatrix * guess.get_col(round);

      w = w - B * transpose(B) * w;

      w /= norm(w);

      vec wOld = zeros(vectorSize);
      vec wOld2 = zeros(vectorSize);

      int i = 1;
      int gabba = 1;

      while (i <= maxNumIterations + gabba) {

        w = w - B * transpose(B) * w;

        w /= norm(w);

        if (notFine) {

          if (i == maxNumIterations + 1) {

            round--;

            numFailures++;

            if (numFailures > failureLimit) {

              if (round == 0) {

                A = dewhiteningMatrix * B;
                W = transpose(B) * whiteningMatrix;

              } // IF round

              return false;

            } // IF numFailures > failureLimit

            break;

          } // IF i == maxNumIterations+1

        } // IF NOTFINE

        else if (i >= endFinetuning) wOld = w;

        if (norm(w - wOld) < epsilon || norm(w + wOld) < epsilon) {

          if (fineTuningEnabled && notFine) {

            notFine = 0;
            gabba = maxFinetune;
            wOld = zeros(vectorSize);
            wOld2 = zeros(vectorSize);
            usedNlinearity = gFine;
            myy = myyK * myyOrig;
            endFinetuning = maxFinetune + i;

          } // IF finetuning

          else {

            numFailures = 0;

            B.set_col(round - 1, w);

            A.set_col(round - 1, dewhiteningMatrix*w);

            W.set_row(round - 1, transpose(whiteningMatrix)*w);

            break;

          } // ELSE finetuning

        } // IF epsilon

        else if (stabilizationEnabled) {

          if (stroke == 0.0 && (norm(w - wOld2) < epsilon || norm(w + wOld2) < epsilon)) {

            stroke = myy;
            myy /= 2.0 ;

            if (mod(usedNlinearity, 2) == 0) {

              usedNlinearity++;

            } // IF MOD

          }// IF !stroke

          else if (stroke != 0.0) {

            myy = stroke;
            stroke = 0.0;

            if (myy == 1 && mod(usedNlinearity, 2) != 0) {
              usedNlinearity--;
            }

          } // IF Stroke

          else if (notFine && !loong && i > maxNumIterations / 2) {

            loong = 1;
            myy /= 2.0;

            if (mod(usedNlinearity, 2) == 0) {

              usedNlinearity++;

            } // IF MOD

          } // IF notFine

        } // IF stabilization


        wOld2 = wOld;
        wOld = w;

        switch (usedNlinearity) {

          // pow3
        case FICA_NONLIN_POW3 : {
          w = (X * pow(transpose(X) * w, 3)) / numSamples - 3 * w;
          break;
        }
        case(FICA_NONLIN_POW3+1) : {
          vec Y = transpose(X) * w;
          vec Gpow3 = X * pow(Y, 3) / numSamples;
          double Beta = dot(w, Gpow3);
          w = w - myy * (Gpow3 - Beta * w) / (3 - Beta);
          break;
        }
        case(FICA_NONLIN_POW3+2) : {
          mat Xsub = X.get_cols(getSamples(numSamples, sampleSize));
          w = (Xsub * pow(transpose(Xsub) * w, 3)) / Xsub.cols() - 3 * w;
          break;
        }
        case(FICA_NONLIN_POW3+3) : {
          mat Xsub = X.get_cols(getSamples(numSamples, sampleSize));
          vec Gpow3 = Xsub * pow(transpose(Xsub) * w, 3) / (Xsub.cols());
          double Beta = dot(w, Gpow3);
          w = w - myy * (Gpow3 - Beta * w) / (3 - Beta);
          break;
        }

        // TANH
        case FICA_NONLIN_TANH : {
          vec hypTan = tanh(a1 * transpose(X) * w);
          w = (X * hypTan - a1 * sum(1 - pow(hypTan, 2)) * w) / numSamples;
          break;
        }
        case(FICA_NONLIN_TANH+1) : {
          vec Y = transpose(X) * w;
          vec hypTan = tanh(a1 * Y);
          double Beta = dot(w, X * hypTan);
          w = w - myy * ((X * hypTan - Beta * w) / (a1 * sum(1 - pow(hypTan, 2)) - Beta));
          break;
        }
        case(FICA_NONLIN_TANH+2) : {
          mat Xsub = X.get_cols(getSamples(numSamples, sampleSize));
          vec hypTan = tanh(a1 * transpose(Xsub) * w);
          w = (Xsub * hypTan - a1 * sum(1 - pow(hypTan, 2)) * w) / Xsub.cols();
          break;
        }
        case(FICA_NONLIN_TANH+3) : {
          mat Xsub = X.get_cols(getSamples(numSamples, sampleSize));
          vec hypTan = tanh(a1 * transpose(Xsub) * w);
          double Beta = dot(w, Xsub * hypTan);
          w = w - myy * ((Xsub * hypTan - Beta * w) / (a1 * sum(1 - pow(hypTan, 2)) - Beta));
          break;
        }

        // GAUSS
        case FICA_NONLIN_GAUSS : {
          vec u = transpose(X) * w;
          vec Usquared = pow(u, 2);
          vec ex = exp(-a2 * Usquared / 2);
          vec gauss = elem_mult(u, ex);
          vec dGauss = elem_mult(1 - a2 * Usquared, ex);
          w = (X * gauss - sum(dGauss) * w) / numSamples;
          break;
        }
        case(FICA_NONLIN_GAUSS+1) : {
          vec u = transpose(X) * w;
          vec Usquared = pow(u, 2);

          vec ex = exp(-a2 * Usquared / 2);
          vec gauss = elem_mult(u, ex);
          vec dGauss = elem_mult(1 - a2 * Usquared, ex);
          double Beta = dot(w, X * gauss);
          w = w - myy * ((X * gauss - Beta * w) / (sum(dGauss) - Beta));
          break;
        }
        case(FICA_NONLIN_GAUSS+2) : {
          mat Xsub = X.get_cols(getSamples(numSamples, sampleSize));
          vec u = transpose(Xsub) * w;
          vec Usquared = pow(u, 2);
          vec ex = exp(-a2 * Usquared / 2);
          vec gauss = elem_mult(u, ex);
          vec dGauss = elem_mult(1 - a2 * Usquared, ex);
          w = (Xsub * gauss - sum(dGauss) * w) / Xsub.cols();
          break;
        }
        case(FICA_NONLIN_GAUSS+3) : {
          mat Xsub = X.get_cols(getSamples(numSamples, sampleSize));
          vec u = transpose(Xsub) * w;
          vec Usquared = pow(u, 2);
          vec ex = exp(-a2 * Usquared / 2);
          vec gauss = elem_mult(u, ex);
          vec dGauss = elem_mult(1 - a2 * Usquared, ex);
          double Beta = dot(w, Xsub * gauss);
          w = w - myy * ((Xsub * gauss - Beta * w) / (sum(dGauss) - Beta));
          break;
        }

        // SKEW
        case FICA_NONLIN_SKEW : {
          w = (X * (pow(transpose(X) * w, 2))) / numSamples;
          break;
        }
        case(FICA_NONLIN_SKEW+1) : {
          vec Y = transpose(X) * w;
          vec Gskew = X * pow(Y, 2) / numSamples;
          double Beta = dot(w, Gskew);
          w = w - myy * (Gskew - Beta * w / (-Beta));
          break;
        }
        case(FICA_NONLIN_SKEW+2) : {
          mat Xsub = X.get_cols(getSamples(numSamples, sampleSize));
          w = (Xsub * (pow(transpose(Xsub) * w, 2))) / Xsub.cols();
          break;
        }
        case(FICA_NONLIN_SKEW+3) : {
          mat Xsub = X.get_cols(getSamples(numSamples, sampleSize));
          vec Gskew = Xsub * pow(transpose(Xsub) * w, 2) / Xsub.cols();
          double Beta = dot(w, Gskew);
          w = w - myy * (Gskew - Beta * w) / (-Beta);
          break;
        }

        } // SWITCH nonLinearity

        w /= norm(w);
        i++;

      } // WHILE i<= maxNumIterations+gabba

      round++;

    } // While round <= numOfIC

  } // ELSE Deflation
  return true;
} // FPICA
