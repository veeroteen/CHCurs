//FEM.h
#pragma once
#include"MISC.h"
#include "Primitives.h"
#include <fstream>
#include <string>
#include "CS.h"
#include <set>
#include "CGM.h"
#include "Functions.h"
template<Field T, ElemType El>
class FEM
{
private:
   std::vector<El> elements;
   std::vector<Node<T, El::GetDim()>> nodes;
   std::vector<T> u;
   std::vector<T> *f;
   CGM<T> *SLAU;
   std::string
      dirihPath,
      newumannPath,
      robinPath,
      nodesPath,
      elemntsPath,
      fPath;
   size_t dirih = 0, neumann = 0, robin = 0;
   T gamma = 0;
   size_t size = 0;
   size_t elemCount = 0;

   void Neumann(size_t count)
   {
      std::ifstream file(newumannPath);
      if (file.is_open())
      {
         auto &fr = *f;
         for (size_t i = 0; i < count; i++)
         {
            std::array<size_t, 3> triangle;
            file >> triangle[0] >> triangle[1] >> triangle[2];
            size_t tetr = 0;
            file >> tetr;
            std::string fun;
            std::getline(file, fun);
            Triangle<T> triag = elements[tetr].getSubElement(triangle, nodes);
            T S = triag.S;
            for (size_t j = 0; j < triangle.size(); j++)
            {
               std::string bas;
               polyToStr(bas, (elements[tetr].getBasis(triangle[j]))->basis);
               std::vector<StringFun<T>> poly{ bas , fun };
               fr[triangle[j]] += integrate(elements[tetr], triangle, nodes, S, poly);
            }

         }

         file.close();
      }
   }
   void Robin( size_t count, std::vector<T> &dir, std::vector<T> &llr,std::vector<size_t> ilr, std::vector<size_t> jlr)
   {
      std::ifstream file(robinPath);
      if (file.is_open())
      {
         auto &fr = *f;
         std::string beta;
         for (size_t i = 0; i < count; i++)
         {
            std::array<size_t, 3> triangle;
            size_t tetr = 0;
            file >> triangle[0] >> triangle[1] >> triangle[2] >> tetr >> beta;

            Triangle<T> triag = elements[tetr].getSubElement(triangle, nodes);
            T S = triag.S;
            std::string fun;
            std::getline(file, fun);

            for (size_t j = 0; j < triangle.size(); j++)
            {
               std::string bas;
               polyToStr(bas, (elements[tetr].getBasis(triangle[j]))->basis);

               std::vector<StringFun<T>> poly = { bas,fun };
               fr[triangle[j]] += integrate(elements[tetr], triangle, nodes, S, poly);
               poly.clear();
               poly.emplace_back(bas);
               poly.emplace_back(bas);
               poly.emplace_back(beta);

               dir[triangle[j]] += integrate(elements[tetr], triangle, nodes, S, poly);
               for (size_t k = 0; k < j; k++)
               {
                  for (size_t it = ilr[triangle[j]]; it < ilr[triangle[j] + 1]; it++)
                  {
                     if (jlr[it] == triangle[k])
                     {
                        polyToStr(bas, (elements[tetr].getBasis(triangle[k]))->basis);
                        poly[1] = StringFun<T>(bas);
                        llr[it] += integrate(elements[tetr], triangle, nodes, S, poly);
                        break;
                     }
                  }
               }
            }
         }
         file.close();
      }
   }
   void Dirih( size_t count, std::vector<T> &dir, std::vector<T> &llr, std::vector<size_t> ilr, std::vector<size_t> jlr)
   {
      std::ifstream file(dirihPath);
      if (file.is_open())
      {
         auto &fr = *f;
         size_t node;
         T value;
         for (size_t i = 0; i < count; i++)
         {
            file >> node >> value;
            fr[node] = value;
            dir[node] = 1;

            for (size_t it = ilr[node]; it < ilr[node + 1]; it++)
            {
               fr[jlr[it]] -= value * llr[it];
               llr[it] = 0;
            }

            for (size_t p = node + 1; p < dir.size(); p++)
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
   }
   void setF( StringFun<T> &fl)
   {
      std::ifstream file(fPath);
      std::string fun = "";
      std::getline(file, fun);
      fl = StringFun<T>(fun);
      file.close();
   }
   void setNodes()
   {
      std::ifstream file(nodesPath);
      {
         std::vector<T> buff(El::GetDim());
         for (size_t i = 0; i < nodes.size(); i++)
         {
            T val = T(0);
            for (size_t j = 0; j < El::GetDim(); j++)
            {
               file >> val;
               buff[j] = val;
            }
            file >> val;
            u[i] = val;
            nodes[i] = Node<T, El::GetDim()>(buff.begin(), u.begin()+i);
         }
      }
      file.close();
   }
   void setElements(size_t count, std::vector<std::set<size_t>> &dict, size_t &llcount)
   {
      
      std::ifstream file(elemntsPath);
      {
         elements.resize(count);
         size_t it = 0;
         std::vector<size_t> buff(El::GetNodesCount());

         for (size_t i = 0; i < count; i++)
         {
            for (size_t j = 0; j < El::GetNodesCount(); j++)
            {
               file >> it;
               buff[j] = it;
            }

            for (size_t h = 0; h < buff.size(); h++)
            {
               for (size_t a = 0; a < h; a++)
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

               for (size_t a = h + 1; a < buff.size(); a++)
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
            std::string lambda;
            file >> lambda;

            elements[i] = El(buff.begin(), lambda, nodes);
            elements[i].setBasis(nodes);
         }
      }
      file.close();
   
   }
public:
   FEM(std::string &confPath)
   {
      std::ifstream conf(confPath);
      std::string DIR;
      conf >> DIR;
      conf.close();

      conf.open(DIR + "/config.txt");

      conf >> gamma;
      conf >> size;
      conf >> elemCount;
      conf >> dirih;
      conf >> neumann;
      conf >> robin;

      dirihPath = DIR + "/dirih.txt";
      newumannPath = DIR + "/neumann.txt";
      robinPath = DIR + "/robin.txt";
      nodesPath = DIR + "/nodes.txt";
      elemntsPath = DIR + "/elems.txt";
      fPath = DIR + "/f.txt";
      conf.close();



     
   }
   


   void Solve(unsigned i)
   {
      StringFun<T> fl;
      setF(fl);

      f = new std::vector<T>(size);
      nodes.resize(size);
      u.resize(size);
      setNodes();

      std::vector<std::set<size_t>> dict(size);
      size_t llcount = 0;
      setElements(elemCount, dict, llcount);


      std::vector<size_t> *il = new std::vector<size_t>(size + 1, 0);
      std::vector<size_t> *jl = new std::vector<size_t>(llcount, 0);
      auto &ilr = *il;
      auto &jlr = *jl;
      auto &fr = *f;
      size_t jli = 0;
      for (size_t a = 0; a < dict.size(); a++)
      {
         ilr[a + 1] = dict[a].size() + ilr[a];
         for (auto &n : dict[a])
         {
            jlr[jli] = n;
            jli++;
         }
      }
      std::vector<T> *di = new std::vector<T>(size, 0);
      std::vector<T> *ll = new std::vector<T>(llcount, 0);
      auto &llr = *ll;
      auto &dir = *di;

      T delta = T(0);
      size_t count = 0;
      do 
      {
         for (size_t e = 0; e < elements.size(); e++)
         {
            std::vector<T> local(El::GetNodesCount() * (El::GetNodesCount() + 1) / 2, 0);
            for (size_t i = 0; i < El::GetNodesCount(); i++)
            {
               std::vector<StringFun<T>> polys;
               std::string base = "";
               polyToStr(base, elements[e].basis[i].basis);
               polys.emplace_back(base);

               for (size_t j = 0; j <= i; j++)
               {
                  polyToStr(base, elements[e].basis[j].basis);
                  polys.emplace_back(base);
                  std::vector<StringFun<T>> grads;
                  std::string base = std::to_string(elements[e].scalGrad(i, j));
                  grads.emplace_back(base);
                  grads.emplace_back(elements[e].h);
                  local[((i * (i + 1)) / 2) + j] += integrate(elements[e], elements[e].vertices, nodes, elements[e].V, grads);

                  local[((i * (i + 1)) / 2) + j] += gamma * integrate(elements[e], elements[e].vertices, nodes, elements[e].V, polys);
                  polys.pop_back();
               }
               polys.push_back(fl);
               fr[elements[e].vertices[i]] += integrate(elements[e], elements[e].vertices, nodes, elements[e].V, polys);
            }

            for (size_t i = 0; i < elements[e].vertices.size(); i++)
            {
               dir[elements[e].vertices[i]] += local[((i * (i + 1)) / 2) + i];
               for (size_t j = 0; j <= i; j++)
               {
                  if (elements[e].vertices[j] < elements[e].vertices[i])
                  {
                     for (size_t it = ilr[elements[e].vertices[i]]; it < ilr[elements[e].vertices[i] + 1]; it++)
                     {
                        if (jlr[it] == elements[e].vertices[j])
                        {
                           llr[it] += local[((i * (i + 1)) / 2) + j];
                        }
                     }
                  }
               }

            }
         }

         Neumann(neumann);
         Robin(robin, dir, llr, ilr, jlr);
         Dirih(dirih, dir, llr, ilr, jlr);


         SLAU = new CGM<T>(il, jl, nullptr, nullptr, ll, nullptr, di, f);

         SLAU->Solve(i);
         auto next = SLAU->getX();
         delta = mod(u, next);

         for (size_t i = 0; i < nodes.size(); i++)
         {
            u[i] = next[i];
         }
         for(auto &a : llr)
         {
            a = T(0);
         }
         for(auto &a: dir)
         {
            a = T(0);
         }
         for(auto &a : fr)
         {
            a = T(0);
         }
         count++;
      } while (delta > 1E-14);

      std::cout << count << std::endl;
   }
   void printU()
   {
      for(size_t i = 0; i < u.size(); i++)
      {
         std::cout << std::setprecision(16) << u[i] << std::endl;
      }
   
   }


};