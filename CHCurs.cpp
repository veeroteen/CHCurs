#include <iostream>
#include "Primitives.h"
#include "FEM.h"

int main()
{

   std::string parasha("cfg.txt");
   NNFEM<double,Tetrahedron<double>> fm(parasha);
   fm.Solve(3);
   fm.printU();

}

