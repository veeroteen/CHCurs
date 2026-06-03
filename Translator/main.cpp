
#include <iostream>
#include "Translator.h"
#include "Functions.h"
int main()
{

      std::string path = "tetr" + std::to_string(1) + ".msh";

      std::string outDir = "../test" + std::to_string(3);
      Translator trans(path, outDir);
      trans.setDirih(outDir, ustr);
      trans.setNeumanZero();
      //trans.setNeuman(outDir, gu);
      trans.setRobinZero();
      //trans.setRobin(outDir, gu, ustr);
      trans.setNodes(outDir, fu,ustr);
      trans.setElements(outDir);
   
}

