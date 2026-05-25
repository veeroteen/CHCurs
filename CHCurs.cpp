#include <iostream>
#include "Primitives.h"
#include "FEM.h"

int main()
{

   std::string cfg("cfg.txt");
   NNFEM<double,Tetrahedron<double>> fm(cfg);
   fm.Solve(1);
   fm.printU();

}

