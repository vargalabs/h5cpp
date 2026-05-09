#include <iostream>

#include <itpp/base/circular_buffer.h>
#include <itpp/base/random.h>

using namespace std;
using namespace itpp;

int main()
{
//Set a fixed seed to the random number generator
  RNG_reset(12345);

  int  nrof_elements;
  ivec index_vec;
  vec  a = randn(3);
  vec  b = randn(8);
  vec  out_vec;
  double   out;
  Array <double> out_array;

  Circular_Buffer<double> cb1(10);

//Put the elements of a to the buffer
  cb1.put(a);

//Vector output: Peek at the two oldest elements of the buffer (the two first extracted)
  cb1.peek(out_vec, 2);
  cout << "peek(out_vec,2) = " << out_vec << ": display 2 first elements of the buffer, without affecting the content"  << endl;

//Vector output: Peek at all elements of the buffer in reverse order
  cb1.peek_reverse(out_vec);
  cout << "peek_reverse(out_vec,-1) = " << out_vec << ": display buffer, without affecting the content"  << endl;

//Array output: Peek at all elements of the buffer in reverse order
  cb1.peek_reverse(out_array);
  cout << "peek_reverse(out_array,-1) = " << out_array << ": display buffer, without affecting the content"  << endl;

//Put the elements of \c b to the buffer
  cb1.put(b);

//Extract the oldest element of the buffer
  cb1.get(out_vec, 1);
  cout << "get(out_vec,1) = " << out_vec << endl ;

//Extract all element of the buffer
  cb1.get(out_vec);
  cout << "get(out_vec) = " << out_vec << endl ;

  cb1.get(out_vec);
  cb1.get(out_array, 0);
  cout << "get(out_vec) = " << out_vec;
  cout << "; get(out_array,0) = " << out_array << ": empty buffer, no content" << endl;

  for (int i = 0;i < a.length();i++) { cb1.put(a(i)); }
  for (int i = 0;i < b.length();i++) { cb1.put(b(i)); }

  cout << "buffer size = " << cb1.size() << endl;
  cout << "nrof_elements = " << cb1.nrof_elements() << endl;

  nrof_elements = cb1.nrof_elements();
  for (int i = 0;i < nrof_elements;i++) {
    cout << "i = " << i;
    cout << ": peek() = " << cb1.peek();
    cout << "; peek_reverse() = " << cb1.peek_reverse();
    cout << "; get() = " << cb1.get();
    cout << endl;
  }

//Test of error handling
//cout << "get() = " << cb1.get() << endl;

  cb1.put(b);
  cout << "buffer size = " << cb1.size() << endl;
  cout << "nrof_elements = " << cb1.nrof_elements() << endl;

  cb1.peek(out_vec, -1);
  cout << "peek(out_vec) = " << out_vec << endl;

  index_vec = "5 3 7 1";
  for (int i = 0;i < index_vec.size();i++) {
    cb1.peek(index_vec(i), out);
    cout << "peek at index " << index_vec(i);
    cout << ": peek(" << index_vec(i) << ") = " << out << endl;
  }

  cb1.peek(index_vec, out_vec);
  cout << "peek at indices " << index_vec << ":" << endl;
  cout << "peek(index_vec,out_vec) = " << out_vec << endl;

  cb1.set_size(15, true);
  cout << "buffer size = " << cb1.size() << endl;
  cout << "nrof_elements = " << cb1.nrof_elements() << endl;

  cb1.peek(out_vec, -1);
  cout << "peek(out_vec) = " << out_vec << endl;

  cb1.set_size(5, true);
  cout << "buffer size = " << cb1.size() << endl;
  cout << "nrof_elements = " << cb1.nrof_elements() << endl;

  cb1.peek(out_vec, -1);
  cout << "peek(out_vec) = " << out_vec << endl;

  cout << "Copy constructor: " << endl;
  Circular_Buffer<double> cb2(cb1);
  cb2.peek(out_vec, -1);
  cout << "Circular_Buffer<double> cb2(cb1): cb2.peek(out_vec) = " << out_vec << endl;

  cout << "Copy operator: " << endl;
  Circular_Buffer<double> cb3 = cb1;
  cb3.peek(out_array, -1);
  cout << "Circular_Buffer<double> cb3=cb1: cb3.peek(out_array) = " << out_array << endl;

  return 0;
}
