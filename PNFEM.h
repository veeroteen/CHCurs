#pragma once
#include "FEM.h"


template<Field T, ElemType El>
class PNFEM : public FEM <T, El >
{
   using FEM<T, El>::elements;
   using FEM<T, El>::nodes;
   using FEM<T, El>::u;
   using FEM<T, El>::f;
   using FEM<T, El>::size;
   using FEM<T, El>::gamma;

   using FEM<T, El>::SLAU;
   using FEM<T, El>::Neumann;
   using FEM<T, El>::Robin;
   using FEM<T, El>::Dirih;
   using FEM<T, El>::setNodes;
   using FEM<T, El>::setElements;
   using FEM<T, El>::setF;

public:
   using FEM<T, El>::printU;
   PNFEM(std::string &confPath) : FEM<T, El>(confPath) {};
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
      setElements(this->elemCount, dict, llcount);


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
      while (true)
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
               polys.push_back(gamma);
               for (size_t j = 0; j <= i; j++)
               {
                  polyToStr(base, elements[e].basis[j].basis);
                  polys.emplace_back(base);
                  std::vector<StringFun<T>> grads;
                  std::string base = std::format("{:.15f}", elements[e].scalGrad(i, j));
                  grads.emplace_back(base);
                  grads.emplace_back(elements[e].h);
                  local[((i * (i + 1)) / 2) + j] += integrate(elements[e], elements[e].vertices, nodes, elements[e].V, grads);
                  local[((i * (i + 1)) / 2) + j] += integrate(elements[e], elements[e].vertices, nodes, elements[e].V, polys);
                  polys.pop_back();
               }
               polys.pop_back();
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

         Neumann(this->neumann);
         Robin(this->robin, dir, llr, ilr, jlr);
         Dirih(this->dirih, dir, llr, ilr, jlr);

         CMatrix<T> matrix(il, jl, nullptr, nullptr, ll, nullptr, di, true);
         std::vector<T> Re(u.size(), 0);
         matrix.multiplyA(u, Re);
         diff(Re, fr, Re);
         if (sqrt(scalar(Re, Re)) < 1E-14)
         {
            break;
         }

         SLAU = new CGM<T>(il, jl, nullptr, nullptr, ll, nullptr, di, f);

         SLAU->Solve(i);
         auto next = SLAU->getX();
         delta = mod(u, next);
         for (size_t i = 0; i < nodes.size(); i++)
         {
            u[i] = next[i];
         }
         for (auto &a : llr)
         {
            a = T(0);
         }
         for (auto &a : dir)
         {
            a = T(0);
         }
         for (auto &a : fr)
         {
            a = T(0);
         }
         count++;
      }

      std::cout << count << std::endl;
   }
};


