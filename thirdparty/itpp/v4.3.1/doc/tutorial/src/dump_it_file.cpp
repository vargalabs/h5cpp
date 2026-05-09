#include <itpp/itbase.h>
#include <iostream>
#include <cstdlib>

using namespace itpp;
using namespace std;

int main(int argc, char *argv[])
{
	it_ifile ff;

	if (1 >= argc)
	{
		cout << "Usage: " << argv[0] << " file.it" << endl;
		return EXIT_FAILURE;
	}

	for (int i = 1; i < argc; ++i)
	{
		ff.open(argv[i]);
		cout << "Opening " << argv[i] << endl;
		int n = 0;
		string name;
		string type;
		string desc;
		uint64_t bytes;
		while(true == ff.seek(n++))
		{
			ff.info(name, type, desc, bytes);
			cout << name << ": ";
			if ("dmat" == type)
			{
				mat tmp;
				ff >> tmp;
				cout << tmp;
			} else if ("imat" == type)
			{
				imat tmp;
				ff >> tmp;
				cout << tmp;
			} else if ("dvec" == type)
			{
				vec tmp;
				ff >> tmp;
				cout << tmp;
			} else if ("ivec" == type)
			{
				ivec tmp;
				ff >> tmp;
				cout << tmp;
			} else if ("float64" == type)
			{
				double tmp;
				ff >> tmp;
				cout << tmp;
			} else if ("int32" == type)
			{
				int tmp;
				ff >> tmp;
				cout << tmp;
			} else
			{
				cout << type;
			}
			cout << endl;
		}
		ff.close();
	}
	return EXIT_SUCCESS;
}
