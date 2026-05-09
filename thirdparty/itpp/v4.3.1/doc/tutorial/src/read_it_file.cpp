#include <itpp/itcomm.h>

using namespace itpp;

int main()
{
  // Declare the it_file class
  it_file ff;

  // Open the file "it_file_test.it" for reading
  ff.open("it_file_test.it");

  // Read the variable a from the file. Put result in vector a.
  vec a;
  ff >> Name("a") >> a;

  // Print the result
  std::cout << "a = " << a << std::endl;

  // Exit the program
  return 0;
}
