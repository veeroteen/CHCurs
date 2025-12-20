
#include <iostream>
#include "Translator.h"
#include "Functions.h"
int main()
{
   std::string path = "tetr.msh";
   Translator trans(path);
   std::string outDir = "test4";
   trans.setDirih(outDir,u2);
   trans.setNewman(outDir,gu2);
   trans.setNodes(outDir,dgu2);
   trans.setElements(outDir);
}

