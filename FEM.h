//FEM.h
#pragma once
#include"MISC.h"
#include "Primitives.h"
#include <fstream>
#include <string>
#include "CS.h"
#include <set>
#include "CGM.h"
template<Field T, ElemType El>
class FEM
{
protected:
   std::vector<El> elements;
   std::vector<Node<T, El::GetDim()>> nodes;
   std::vector<T> *f;
   CGM<T> *SLAU;
public:
   FEM(std::string &confPath)
   {
      std::ifstream conf(confPath);
      std::string DIR;
      conf >> DIR;
      conf.close();
      size_t size = 0;
      size_t elemCount = 0;
      size_t Dim = El::GetDim();
      conf.open(DIR + "/config.txt");
      size_t dirih = 0, neumann = 0, robin = 0;
      T gamma = 0;
      conf >> gamma;
      conf >> size;
      conf >> elemCount;
      conf >> dirih;
      conf >> neumann;
      conf >> robin;
      
      std::string
         dirihPath = DIR + "/dirih.txt",
         newumannPath = DIR + "/neumann.txt",
         robinPath = DIR + "/robin.txt",
         nodesPath = DIR + "/nodes.txt",
         elemntsPath = DIR + "/elems.txt",
         fPath = DIR + "/f.txt";
      f = new std::vector<T>(size,0);
      std::vector<T> fl(size, 0);
      conf.close();
      std::ifstream file(fPath);
      {
         T val = T();
         for(size_t i = 0; i < size; i++)
         {
            file >> val;
            fl[i] = val;
         }
      }
      file.close();

      nodes.resize(size);
      file.open(nodesPath);
      {
         T val = T();
         std::vector<T> buff(Dim);
         for (size_t i = 0; i < size; i++)
         {
            for (size_t j = 0; j < Dim; j++)
            {
               file >> val;
               buff[j] = val;
            }
            nodes[i] = Node<T, El::GetDim()>(buff.begin());
         }
      }
      file.close();
      

      std::vector<std::set<size_t>> dict(size);
      
      size_t llcount = 0;
      file.open(elemntsPath);
      {
         elements.resize(elemCount);
         T val = T();
         size_t it = 0;
         std::vector<size_t> buff(El::GetNodesCount());

         for (size_t i = 0; i < elemCount; i++)
         {
            for (size_t j = 0; j < El::GetNodesCount(); j++)
            {
               file >> it;
               buff[j] = it;
            }

            for(size_t h = 0; h < buff.size();h++)
            {
               for (size_t a = 0; a < h; a++)
               {
                  if(buff[a] < buff[h])
                  {
                     auto flag = dict[buff[h]].insert(buff[a]);
                     if (flag.second)
                     {
                        llcount++;
                     }
                  }
               }

               for(size_t a = h+1; a < buff.size();a++)
               {
                  if (buff[a] < buff[h])
                  {
                     auto flag = dict[buff[h]].insert(buff[a]);
                     if (flag.second)
                     {
                        llcount++;
                     }
                  }
               }
            
            }

            file >> val;
            elements[i] = El(buff.begin(), val, nodes);
            elements[i].setBasis(nodes);
         }
      }
      file.close();


      std::vector<size_t> *il = new std::vector<size_t>(size + 1,0);
      std::vector<size_t> *jl = new std::vector<size_t>(llcount, 0);

      auto &ilr = *il;
      auto &jlr = *jl;
      auto &fr = *f;
      size_t jli = 0;
      
      for(size_t a = 0; a < dict.size();a++)
      {
         ilr[a+1] = dict[a].size() + ilr[a];
         for(auto &n : dict[a])
         {
            jlr[jli] = n;
            jli++;
         }
      }
      std::vector<T> *di = new std::vector<T>(size,0);
      std::vector<T> *ll = new std::vector<T>(llcount, 0);
      auto &llr = *ll;
      auto &dir = *di;

      for (size_t e = 0; e < elements.size(); e++) 
      {
         std::vector<T> local(El::GetNodesCount() * (El::GetNodesCount() + 1) / 2 ,0);
         for (size_t i = 0; i < El::GetNodesCount(); i++)
         {
            for (size_t j = 0; j <= i; j++)
            {
               local[((i * (i + 1)) / 2) + j] += elements[e].scalGrad(i, j) * elements[e].V * elements[e].h;
               local[((i * (i + 1)) / 2) + j] += gamma * elements[e].V / (i == j ? 10 : 20);
            }
            fr[elements[e].vertices[i]] += fl[elements[e].vertices[i]] * elements[e].V / 4;
         }

         for(size_t i = 0; i < elements[e].vertices.size(); i++)
         {
            dir[elements[e].vertices[i]] += local[((i*(i+1))/2) + i];
            for(size_t j = 0; j <= i;j++)
            {
               if (elements[e].vertices[j] < elements[e].vertices[i])
               {
                  for (size_t it = ilr[elements[e].vertices[i]]; it < ilr[elements[e].vertices[i] + 1]; it++)
                  {
                     if(jlr[it] == elements[e].vertices[j])
                     {
                        llr[it] += local[((i * (i + 1)) / 2) + j];
                     }
                  }
               }
            }

         }
      }

      Neumann(newumannPath,neumann);
      
      file.open(robinPath);
      if (file.is_open())
      {
         T beta,h;
         for (size_t i = 0; i < robin; i++)
         {
            std::array<size_t, 3> triangle;
            size_t tetr = 0;
            file >> triangle[0] >> triangle[1] >> triangle[2] >> tetr >> beta;

            h = elements[tetr].h;
            T S = triangleArea(nodes[triangle[0]], nodes[triangle[1]], nodes[triangle[2]]);
            std::string fun;
            std::getline(file, fun);

            for (size_t j = 0; j < triangle.size(); j++)
            {
               std::string bas;
               polyToStr(bas, (elements[tetr].getBasis(triangle[j]))->basis);
               
               std::vector<std::string> poly = { bas,fun };
               fr[triangle[j]] += integrate(elements[tetr], triangle, nodes, S, poly);

               dir[triangle[j]] += beta * S / 6;

               for (size_t k = 0; k < j; k++)
               {
                  for (size_t it = ilr[triangle[j]]; it < ilr[triangle[j] + 1]; it++)
                  {
                     if (jlr[it] == triangle[k])
                     {
                        llr[it] += beta * S / 12;
                        break;
                     }
                  }
               }
            }
         }
         file.close();
      }
      file.open(dirihPath);
      if(file.is_open())
      {
         size_t node;
         T value;
         for(size_t i = 0; i < dirih; i++)
         {
            file >> node >> value;
            fr[node] = value;
            dir[node] = 1;
            
            for (size_t it = ilr[node]; it < ilr[node + 1]; it++)
            {
               fr[jlr[it]] -= value * llr[it];
               llr[it] = 0;
            }
            
            for (size_t p = node + 1; p < size; p++)
            {
               for (size_t it = ilr[p]; it < ilr[p + 1] && jlr[it] <= node; it++)
               {
                  if (jlr[it] == node)
                  {

                     fr[p] = fr[p] - value * llr[it];
                     llr[it] = 0;

                  }
               }
                  
            }

         }
         file.close();
      }

      SLAU = new CGM<T>(il, jl, nullptr, nullptr, ll, nullptr, di, f);
   }
   
   void Neumann(std::string &path,size_t count)
   {
      std::ifstream file(path);
      if (file.is_open())
      {
         T h, teta;
         auto &fr = *f;
         for (size_t i = 0; i < count; i++)
         {
            std::array<size_t, 3> triangle;
            file >> triangle[0] >> triangle[1] >> triangle[2];
            size_t tetr = 0;
            file >> tetr;
            std::string fun;
            std::getline(file, fun);
            h = elements[tetr].h;
            T S = triangleArea(nodes[triangle[0]], nodes[triangle[1]], nodes[triangle[2]]);
            for (size_t j = 0; j < triangle.size(); j++)
            {
               std::string bas;
               polyToStr(bas, (elements[tetr].getBasis(triangle[j]))->basis);
               std::vector<std::string> poly{ bas , fun };
               fr[triangle[j]] += integrate(elements[tetr], triangle, nodes, S, poly);
            }

         }
      
         file.close();
      }
   }

   void Solve(unsigned i)
   {
      SLAU->Solve(i);
      SLAU->outX();
   }

};