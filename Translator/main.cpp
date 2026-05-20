
#include <iostream>
#include "Translator.h"
#include "Functions.h"
int main()
{

      std::string path = "tetr" + std::to_string(2) + ".msh";

      std::string outDir = "../test" + std::to_string(2);
      Translator trans(path, outDir);
      trans.setDirih(outDir, u);
      trans.setNeumanZero();
      //trans.setNeuman(outDir, gu);
      trans.setRobinZero();
      //trans.setRobin(outDir, gu, ustr);
      trans.setNodes(outDir, fu);
      trans.setElements(outDir);
   
}

