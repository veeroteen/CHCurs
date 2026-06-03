#include <iostream>
#include "Primitives.h"
#include "FEM.h"
#include "NNFEM.h"

int main()
{

   std::string cfg("cfg.txt");
   LFEM<double,Tetrahedron<double>> fm(cfg);
   fm.Solve(1);

}

