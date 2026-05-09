/*!
 * \file
 * \brief Sorting functions
 * \author Tony Ottosson, Mark Dobossy and Adam Piatyszek
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

#ifndef SORT_H
#define SORT_H

#include <itpp/base/vec.h>
#include <itpp/base/converters.h>
#include <itpp/base/math/log_exp.h>


namespace itpp
{

/*!
 * \brief Sorting algorithms that can be used in a \a Sort class
 *
 * - Introsort (the default and the fastest method in most cases)
 * - Quicksort
 * - Heapsort
 * - Insertion Sort (suitable for very short vectors)
 */
enum SORTING_METHOD { INTROSORT = 0, QUICKSORT = 1, HEAPSORT = 2,
                      INSERTSORT = 3
                    };

/*!
 * \brief Class for sorting of vectors
 *
 * A class which takes a vector, and sorts its values descending. There
 * are two types of sort: a normal sort (accessed via the Sort::sort()
 * function) which sorts the vector passed in the argument, and an index
 * sort (accessed via the Sort::sort_index() function) which leaves the
 * passed vector intact, but returns an index vector describing the sorted
 * order.
 *
 * The Sort class has four sorting methods implemented:
 * - Introsort [1,2]: It is a sorting algorithm developed by David Musser.
 *   I starts as a quicksort, but switches to a heapsort in cases where
 *   the recursion becomes too deep. Additionally, when sub vectors become
 *   smaller than 16 elements, it switches to an insertion sort. Introsort
 *   has a worst-case of \f$\Theta(n\log n)\f$ comparisons, bu has the
 *   efficiency of the quick sort algorithm in cases where the data is
 *   well conditioned for quicksort.
 * - Quicksort [3]: It is a comparison sorting algorithm that has an
 *   average complexity of \f$\Theta(n\log n)\f$ comparisons. For most
 *   data sets, the quicksort will be significantly more efficient than
 *   this average. However for data sets not well suited to it, quicksort
 *   may require as many as \f$\Theta(n^2)\f$ comparisons. Example of such
 *   ill-suited data sets are those which are nearly in order, and data
 *   sets with multiple elements of the same value.
 * - Heapsort [4]: It is a comparison sorting algorithm. While it is
 *   usually seen to be slower than quicksort routines, its worst-case
 *   requires only \f$\Theta(n\log n)\f$ comparisons. This makes it an
 *   ideal quicksort replacement for data sets that that are
 *   ill-conditioned for quicksorting.
 * - Insertion sort [5]: An insertion sort is a simple comparison sort,
 *   which is widely considered to be one of the most efficient algorithms
 *   for very small data sets (10-20 elements).
 *   http://en.wikipedia.org/wiki/Insertion_sort
 *
 * References:
 * - [1] http://en.wikipedia.org/wiki/Introsort
 * - [2] http://www.cs.rpi.edu/~musser/gp/introsort.ps
 * - [3] http://en.wikipedia.org/wiki/Quicksort
 * - [4] http://en.wikipedia.org/wiki/Heapsort
 * - [5] http://en.wikipedia.org/wiki/Insertion_sort
 *
 * \author Tony Ottosson (Quicksort), Mark Dobossy (Introsort, Heapsort
 *         and Insertion Sort) and Adam Piatyszek (Sort class design, code
 *         clean-up)
 */
template<class T>
class Sort
{
public:
  //! Constructor that sets Intro Sort method by default
  Sort(SORTING_METHOD method = INTROSORT): sort_method(method) {}

  //! Set sorting method
  void set_method(SORTING_METHOD method) { sort_method = method; }

  //! Get current sorting method
  SORTING_METHOD get_method() const { return sort_method; }

  /*!
   * \brief Sorting function of a subset of a vector \a data
   *
   * \param low Start index of a subvector to be sorted
   * \param high End index of a subvector to be sorted
   * \param data Data vector, in which a part of it is to be sorted
   */
  void sort(int low, int high, Vec<T> &data);

  /*!
   * \brief Sorting function that returns a sorted index vector
   *
   * \param low Start index of a subvector to be sorted
   * \param high End index of a subvector to be sorted
   * \param data Data vector, in which a part of it is to be sorted
   */
  ivec sort_index(int low, int high, const Vec<T> &data);

  /*!
   * \brief Introsort function of a subset of a vector \c data
   *
   * \param low Start index of a subvector to be sorted
   * \param high End index of a subvector to be sorted
   * \param max_depth Maximum recursion depth before switching to heap sort
   *        recommended value: log2 of the length of the data vector
   * \param data Data vector, in which a part of it is to be sorted
   *
   * \note An introsort is not a stable sort (i.e. it may not maintain
   *       the relative order of elements with equal value.)
   * \note This function uses recurrence.
   */
  void intro_sort(int low, int high, int max_depth, Vec<T> &data);

  /*!
   * \brief Introsort function, which returns a sorted index vector
   *
   * \param low Start index of a subvector to be sorted
   * \param high End index of a subvector to be sorted
   * \param max_depth Maximum recursion depth before switching to heap sort
   *        recommended value: log2 of the length of the data vector
   * \param data Data vector, in which a part of it is to be sorted
   *
   * \note An Introsort is not a stable sort (i.e. it may not maintain
   *       the relative order of elements with equal value.)
   * \note This function uses recurrence.
   */
  ivec intro_sort_index(int low, int high, int max_depth,
                        const Vec<T> &data);

private:
  SORTING_METHOD sort_method;

  void IntroSort(int low, int high, int max_depth, T data[]);
  void IntroSort_Index(int low, int high, int max_depth, int indexlist[],
                       const T data[]);

  void QuickSort(int low, int high, T data[]);
  void QuickSort_Index(int low, int high, int indexlist[], const T data[]);

  void HeapSort(int low, int high, T data[]);
  void HeapSort_Index(int low, int high, int indexlist[], const T data[]);

  void InsertSort(int low, int high, T data[]);
  void InsertSort_Index(int low, int high, int indexlist[], const T data[]);
};


/*!
 * \relates Vec
 * \brief Sort the \a data vector in increasing order
 *
 * \param data Vector to be sorted
 * \param method Sorting method: INTROSORT (default), QUICKSORT, HEAPSORT
 * or INSERTSORT
 */
template<class T>
void sort(Vec<T> &data, SORTING_METHOD method = INTROSORT)
{
  Sort<T> s(method);
  s.sort(0, data.size() - 1, data);
}

/*!
 * \relates Vec
 * \brief Return an index vector corresponding to a sorted vector \a data
 * in increasing order
 *
 * \param data Vector for which to return a sorted index vector
 * \param method Sorting method: INTROSORT (default), QUICKSORT, HEAPSORT
 * or INSERTSORT
 */
template<class T>
ivec sort_index(const Vec<T> &data, SORTING_METHOD method = INTROSORT)
{
  Sort<T> s(method);
  return s.sort_index(0, data.size() - 1, data);
}


// ----------------------------------------------------------------------
// Public functions for various sorting methods
// ----------------------------------------------------------------------

template<class T>
void Sort<T>::sort(int low, int high, Vec<T> &data)
{
  int N = data.size();
  // Nothing to sort if data vector has only one or zero elements
  if (N < 2)
    return;

  it_assert((low >= 0) && (high > low) && (high < N), "Sort::sort(): "
            "low or high out of bounds");

  switch (sort_method) {
  case INTROSORT:
    IntroSort(low, high, levels2bits(N), data._data());
    break;
  case QUICKSORT:
    QuickSort(low, high, data._data());
    break;
  case HEAPSORT:
    HeapSort(low, high, data._data());
    break;
  case INSERTSORT:
    InsertSort(low, high, data._data());
    break;
  default:
    it_error("Sort<T>::sort(): Unknown sorting method");
  }
}


template<class T>
ivec Sort<T>::sort_index(int low, int high, const Vec<T> &data)
{
  int N = data.size();
  // Nothing to sort if data vector has only one or zero elements
  if (N == 1)
    return ivec("0");
  else if (N == 0)
    return ivec();

  it_assert((low >= 0) && (high > low) && (high < N), "Sort::sort(): "
            "low or high out of bounds");

  ivec indexlist(N);
  for (int i = 0; i < N; ++i) {
    indexlist(i) = i;
  }

  switch (sort_method) {
  case INTROSORT:
    IntroSort_Index(low, high, levels2bits(N), indexlist._data(),
                    data._data());
    break;
  case QUICKSORT:
    QuickSort_Index(low, high, indexlist._data(), data._data());
    break;
  case HEAPSORT:
    HeapSort_Index(low, high, indexlist._data(), data._data());
    break;
  case INSERTSORT:
    InsertSort_Index(low, high, indexlist._data(), data._data());
    break;
  default:
    it_error("Sort<T>::sort_index(): Unknown sorting method");
  }

  return indexlist;
}


// INTRO SORT
template<class T>
void Sort<T>::intro_sort(int low, int high, int max_depth, Vec<T> &data)
{
  it_assert((low >= 0) && (high > low) && (high < data.size()),
            "Sort::sort(): low or high out of bounds");
  IntroSort(low, high, max_depth, data._data());
}

// INTRO SORT INDEX
template<class T>
ivec Sort<T>::intro_sort_index(int low, int high, int max_depth,
                               const Vec<T> &data)
{
  int N = data.size();
  it_assert((low >= 0) && (high > low) && (high < N),
            "Sort::sort(): low or high out of bounds");

  ivec indexlist(N);
  for (int i = 0; i < N; ++i) {
    indexlist(i) = i;
  }

  IntroSort_Index(low, high, max_depth, indexlist._data(), data._data());

  return indexlist;
}


// ----------------------------------------------------------------------
// Private functions for sorting methods
// ----------------------------------------------------------------------

template<class T>
void Sort<T>::IntroSort(int low, int high, int max_depth, T data[])
{
  if (high - low > 16) {
    max_depth--;
    if (max_depth == 0) {
      HeapSort(low, high, data);
      return;
    }

    if (high > low) {
      T a = data[low];
      int plow = low;
      int phigh = high;
      T test = data[phigh];
      while (plow < phigh) {
        if (test < a) {
          data[plow] = test;
          plow++;
          test = data[plow];
        }
        else {
          data[phigh] = test;
          phigh--;
          test = data[phigh];
        }
      }
      data[plow] = a;
      IntroSort(low, plow - 1, max_depth, data);
      IntroSort(plow + 1, high, max_depth, data);
      return;
    }
  }
  else {
    InsertSort(low, high, data);
    return;
  }
}

template<class T>
void Sort<T>::IntroSort_Index(int low, int high, int max_depth,
                              int indexlist[], const T data[])
{
  if (high - low > 16) {
    max_depth--;
    if (max_depth == 0) {
      HeapSort_Index(low, high, indexlist, data);
      return;
    }

    if (high > low) {
      int aindex = indexlist[low];
      T a = data[aindex];
      int plow = low;
      int phigh = high;
      int testindex = indexlist[phigh];
      T test = data[testindex];
      while (plow < phigh) {
        if (test < a) {
          indexlist[plow] = testindex;
          plow++;
          testindex = indexlist[plow];
          test = data[testindex];
        }
        else {
          indexlist[phigh] = testindex;
          phigh--;
          testindex = indexlist[phigh];
          test = data[testindex];
        }
      }
      indexlist[plow] = aindex;
      IntroSort_Index(low, plow - 1, max_depth, indexlist, data);
      IntroSort_Index(plow + 1, high, max_depth, indexlist, data);
    }
  }
  else {
    InsertSort_Index(low, high, indexlist, data);
    return;
  }
}

template <class T>
void Sort<T>::QuickSort(int low, int high, T data[])
{
  if (high > low) {
    T a = data[low];
    int plow = low;
    int phigh = high;
    T test = data[phigh];
    while (plow < phigh) {
      if (test < a) {
        data[plow] = test;
        plow++;
        test = data[plow];
      }
      else {
        data[phigh] = test;
        phigh--;
        test = data[phigh];
      }
    }
    data[plow] = a;
    QuickSort(low, plow - 1, data);
    QuickSort(plow + 1, high, data);
  }
}

template<class T>
void Sort<T>::QuickSort_Index(int low, int high, int indexlist[],
                              const T data[])
{
  if (high > low) {
    int aindex = indexlist[low];
    T a = data[aindex];
    int plow = low;
    int phigh = high;
    int testindex = indexlist[phigh];
    T test = data[testindex];
    while (plow < phigh) {
      if (test < a) {
        indexlist[plow] = testindex;
        plow++;
        testindex = indexlist[plow];
        test = data[testindex];
      }
      else {
        indexlist[phigh] = testindex;
        phigh--;
        testindex = indexlist[phigh];
        test = data[testindex];
      }
    }
    indexlist[plow] = aindex;
    QuickSort_Index(low, plow - 1, indexlist, data);
    QuickSort_Index(plow + 1, high, indexlist, data);
  }
}

template<class T>
void Sort<T>::HeapSort(int low, int high, T data[])
{
  int size = (high + 1) - low;
  int i = size / 2;
  T temp;
  while (1) {
    if (i > 0)
      temp = data[--i + low];
    else {
      if (size-- == 0)
        break;
      temp = data[size + low];
      data[size + low] = data[low];
    }

    int parent = i;
    int child = i * 2 + 1;

    while (child < size) {
      if (child + 1 < size  &&  data[child + 1 + low] > data[child + low])
        child++;
      if (data[child + low] > temp) {
        data[parent + low] = data[child + low];
        parent = child;
        child = parent * 2 + 1;
      }
      else
        break;
    }
    data[parent + low] = temp;
  }
}

template<class T>
void Sort<T>::HeapSort_Index(int low, int high, int indexlist[],
                             const T data[])
{
  int size = (high + 1) - low;
  int i = size / 2;

  while (1) {
    T tempValue;
    int tempIndex;
    if (i > 0) {
      i--;
      tempValue = data[indexlist[i + low]];
      tempIndex = indexlist[i + low];
    }
    else {
      if (size-- == 0)
        break;
      tempValue = data[indexlist[size + low]];
      tempIndex = indexlist[size + low];
      indexlist[size+low] = indexlist[low];
    }

    int parent = i;
    int child = i * 2 + 1;

    while (child < size) {
      if ((child + 1 < size)
          && data[indexlist[child + 1 + low]] > data[indexlist[child + low]])
        child++;
      if (data[indexlist[child + low]] > tempValue) {
        indexlist[parent + low] = indexlist[child + low];
        parent = child;
        child = parent * 2 + 1;
      }
      else
        break;
    }
    indexlist[parent + low] = tempIndex;
  }
}

template<class T>
void Sort<T>::InsertSort(int low, int high, T data[])
{
  for (int i = low + 1; i <= high; i++) {
    T value = data[i];
    int j;
    for (j = i - 1; j >= low && data[j] > value; j--) {
      data[j + 1] = data[j];
    }
    data[j + 1] = value;
  }
}

template<class T>
void Sort<T>::InsertSort_Index(int low, int high, int indexlist[],
                               const T data[])
{
  for (int i = low + 1; i <= high; i++) {
    T value = data[indexlist[i]];
    int tempIndex = indexlist[i];
    int j;
    for (j = i - 1; j >= low && data[indexlist[j]] > value; j--) {
      indexlist[j + 1] = indexlist[j];
    }
    indexlist[j + 1] = tempIndex;
  }
}


} // namespace itpp

#endif // #ifndef SORT_H
