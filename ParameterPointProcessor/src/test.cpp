#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <cmath>

using namespace std;

int main(int argc, char* argv[]) {


  if ( argc != 2)
  {
	 printf("test usage:\n");
	 printf("test <a>");
  }

//   double a = (double)atof(argv[1]);
//	double out = 2.0 * a;
//	std::cout << out << std::endl;
//	std::cout << 4.0 << " " << 5.0  << std::endl;
	
  string fname = argv[1];
  fstream file;
  file.open(fname.c_str(), ios::in);
  double a;
  file >> a;
  file.close();

	double out = 2.0 * a;
	std::cout << out << std::endl;

	return 0;
  
}
