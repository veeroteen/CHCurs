
#include <iostream>
#include "Translator.h"
#include "Functions.h"
int main()
{
   for (size_t i = 1; i <= 1; i++)
   {
      std::string path = "tetr" + std::to_string(i) + ".msh";

      std::string outDir = "../test" + std::to_string(i);
      Translator trans(path, outDir);
      trans.setDirih(outDir, u);
      trans.setNeuman(outDir, gu);
      //trans.setRobin(outDir, gu, ustr);
      trans.setRobinZero();
      trans.setNodes(outDir, dgu);
      trans.setElements(outDir);
   }
}

