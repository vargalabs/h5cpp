#include <itpp/itcomm.h>

using namespace itpp;

int main()
{
  // Declare the it_file class
  it_file ff;

  // Open a file with the name "it_file_test.it"
  ff.open("it_file_test.it");

  // Create some data to put into the file
  vec a = linspace(1, 20, 20);

  // Put the variable a into the file. The Name("a") tells the file class
  // that the next variable shall be named "a".
  ff << Name("a") << a;

  // Force the file to be written to disc. This is useful when performing
  // iterations and ensures that the information is not stored in any cache
  // memory. In this simple example it is not necessary to flush the file.
  ff.flush();

  // Close the file
  ff.close();

  // Exit program
  return 0;
}
