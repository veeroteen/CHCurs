#pragma once
#include "CS.h"

template <Field T>
class LOS :public ThreeStageBase<T>
{
	using ThreeStageBase<T>::matrix;
	using ThreeStageBase<T>::f;
	using ThreeStageBase<T>::x;
	using ThreeStageBase<T>::symmetry;
	using ThreeStageBase<T>::A;
	using ThreeStageBase<T>::incompLU;
	using ThreeStageBase<T>::incompChol;
	using ThreeStageBase<T>::eps;
	using ThreeStageBase<T>::iterC;

private:
	void iteration(std::vector<T> &r, std::vector<T> &z, std::vector<T> &p)
	{
		T b = scalar(p, p);
		T a = scalar(p, r) / b;

		for (size_t i = 0; i < x->size(); i++)
		{
			(*x)[i] = (*x)[i] + a * z[i];
			r[i] = r[i] - a * p[i];

		}
		std::vector<T> temp(r);
		A->multiplyA(temp);

		b = -scalar(p, temp) / b;


		for (size_t i = 0; i < z.size(); i++)
		{
			z[i] = r[i] + b * z[i];
			p[i] = temp[i] + b * p[i];

		}
	}
	void iterationDiag(std::vector<T> &r, std::vector<T> &z, std::vector<T> &p, std::vector<T> &diag)
	{
		T b = scalar(p, p);
		T a = scalar(p, r) / b;


		for (size_t i = 0; i < x->size(); i++)
		{
			(*x)[i] = (*x)[i] + a * z[i];
			r[i] = r[i] - a * p[i];

		}
		std::vector<T> temp(r);

		matrix.diagSolve(diag, r, temp);
		A->multiplyA(temp);
		matrix.diagSolve(diag, temp, temp);


		b = -scalar(p, temp) / b;

		for (size_t i = 0; i < z.size(); i++)
		{
			z[i] = r[i] / (diag)[i] + b * z[i];
			p[i] = temp[i] + b * p[i];

		}
	}
	void iterationLU(std::vector<T> &r, std::vector<T> &z, std::vector<T> &p)
	{
		T b = scalar(p, p);
		T a = scalar(p, r) / b;


		for (size_t i = 0; i < x->size(); i++)
		{
			(*x)[i] = (*x)[i] + (a * z[i]);
			r[i] = r[i] - (a * p[i]);

		}
		std::vector<T> temp(r);

		matrix.revrsSolutionCSC(*matrix.iu, *matrix.ju, temp, temp, *(matrix.lu), false);
		A->multiplyA(temp);
		matrix.forwSolutionCSR(*matrix.il, *matrix.jl, temp, temp, *(matrix.ll));


		b = -scalar(p, temp) / b;



		std::vector<T> Ur(r);
		matrix.revrsSolutionCSC(*matrix.iu, *matrix.ju, r, Ur, *(matrix.lu), false);
		for (size_t i = 0; i < z.size(); i++)
		{
			z[i] = Ur[i] + b * z[i];
			p[i] = temp[i] + b * p[i];

		}
	}

public:
	LOS(std::vector<size_t> *il, std::vector<size_t> *jl, std::vector<size_t> *iu, std::vector<size_t> *ju, std::vector<T> *ll, std::vector<T> *lu, std::vector<T> *di, std::vector<T> *f) : ThreeStageBase<T>(il, jl, iu, ju, ll, lu, di, f)
	{
	}
	LOS(std::istream &input) : ThreeStageBase<T>(input)
	{
	}
	LOS(std::string &path) : ThreeStageBase<T>(path,false)
	{
	}

	void Solve(unsigned par)
	{
		size_t iter = 0;
		switch (par)
		{
			case 1:
			{

				std::vector<T> r(*f);
				std::vector<T> z(r);
				std::vector<T> p(r);
				A->multiplyA(p);

				for (size_t i = 0; i < iterC && sqrt(scalar(r, r) / scalar(*f, *f)) > eps; i++)
				{
					iteration(r, z, p);

					iter++;
				}
				break;
			}
			case 2:
			{
				std::vector<T> diag(*(A->di));
				for (auto &a : diag)
				{
					a = sqrt(a);
				}


				std::vector<T> r(*f);
				matrix.diagSolve(diag, r, r);
				std::vector<T> z(r);
				matrix.diagSolve(diag, z, z);

				std::vector<T> p(z);
				A->multiplyA(p);
				matrix.diagSolve(diag, p, p);

				for (size_t i = 0; i < iterC && sqrt(scalar(r, r) / scalar(*f, *f)) > eps; i++)
				{
					iterationDiag(r, z, p, diag);

					iter++;
				}
				break;
			}
			case 3:
			{
				incompLU();
				std::vector<T> r(*f);
				matrix.forwSolutionCSR(*matrix.il, *matrix.jl, r, r, *matrix.ll);
				std::vector<T> z(r.size(), 0);
				matrix.revrsSolutionCSC(*matrix.iu, *matrix.ju, r, z, *matrix.lu, false);
				std::vector<T> p(z);
				A->multiplyA(p);
				matrix.forwSolutionCSR(*matrix.il, *matrix.jl, p, p, *matrix.ll);

				for (size_t i = 0; i < iterC && sqrt(scalar(r, r) / scalar(*f, *f)) > eps; i++)
				{
					iterationLU(r, z, p);
					iter++;

				}

				break;
			}

		}
	}
	void printTL()
	{
		for (size_t i = 0; i < A->di->size(); i++)
		{
			size_t count = i;
			size_t j = 0;
			for (size_t il = (*matrix.il)[i]; il < (*matrix.il)[i + 1]; il++)
			{
				for (; j < (*matrix.jl)[il]; j++)
				{
					std::cout << 0 << " ";
					count--;
				}
				std::cout << std::setprecision(8) << (*matrix.ll)[il] << " ";
				j++;
				count--;

			}
			for (size_t t = 0; t < count; t++)
			{
				std::cout << 0 << " ";
			}
			std::cout << std::setprecision(8) << 1 << " ";
			std::cout << std::endl;
		}

	}
	void printV(std::vector<T> &p)
	{
		for (auto &a : p)
		{
			std::cout << std::setprecision(16) << a << std::endl;
		}
		std::cout << std::endl;
	}

};