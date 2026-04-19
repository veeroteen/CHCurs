//CSR.h
#pragma once
#include <vector>
#include <fstream>
#include <numeric>
#include <algorithm>
#include <iomanip>
#include "MISC.h"

template <Field T>
T scalar(const std::vector<T> &a, const std::vector<T> b)
{
	T res = T();

	for (size_t i = 0; i < a.size(); i++)
	{
		res += a[i] * b[i];
	}
	return res;

}

template <Field T>
void diff(std::vector<T> &a, std::vector<T> &b,std::vector<T> &res)
{
	for(size_t i = 0; i < a.size();i++)
	{
		res[i] = a[i] - b[i];
	}
}


template <Field T>
void Vdiff(std::vector<T> &a, std::vector<T> &b, std::vector<T> &res)
{
	res.resize(a.size());
	for (size_t i = 0; i < a.size(); i++)
	{
		res[i] = a[i] - b[i];
	}

}


template <Field T>
struct CMatrix
{
	std::vector<T> *ll, *lu, *di;
	std::vector<size_t> *il, *jl, *iu, *ju;
	CMatrix(std::vector<size_t> *il, std::vector<size_t> *jl, std::vector<size_t> *iu, std::vector<size_t> *ju, std::vector<T> *ll, std::vector<T> *lu, std::vector<T> *di) :
		il(il), jl(jl), iu(iu), ju(ju), ll(ll), lu(lu), di(di)
	{
	}
	void load(std::vector<size_t> *_il, std::vector<size_t> *_jl, std::vector<size_t> *_iu, std::vector<size_t> *_ju, std::vector<T> *_ll, std::vector<T> *_lu, std::vector<T> *_di)
	{
		il = _il;
		jl = _jl;
		iu = _iu;
		ju = _ju;
		ll = _ll;
		lu = _lu;
		di = _di;
	}
	CMatrix()
	{
		il = new std::vector<size_t>();
		jl = new std::vector<size_t>();
		iu = new std::vector<size_t>();
		ju = new std::vector<size_t>();
		ll = new std::vector<T>();
		lu = new std::vector<T>();
		di = new std::vector<T>();
	}
	~CMatrix()
	{
		if (il != nullptr)
		{
			delete il;
		}
		if (jl != nullptr)
		{
			delete jl;
		}
		if (ll != nullptr)
		{
			delete ll;
		}
		if (lu != nullptr)
		{
			delete lu;
		}
		if (di != nullptr)
		{
			delete di;
		}
	}
};


template <Field T>
class ThreeStageBase
{
private:

protected:
	CMatrix<T> matrix;
	std::vector <T> *f;
	std::vector <T> *x;
	CMatrix<T> *A;
	bool symmetry = true;
	T eps;
	size_t iterC;

	ThreeStageBase(std::vector<size_t> *il, std::vector<size_t> *jl, std::vector<size_t> *iu, std::vector<size_t> *ju, std::vector<T> *ll, std::vector<T> *lu, std::vector<T> *di, std::vector<T> *f)
	{
		if(ju != nullptr)
		{
			symmetry = false;
		}
		matrix.load(il, jl, iu, ju, ll, lu, di);
		this->f = f;
		x = new std::vector<T>(f->size(), 0);
		A = new CMatrix<T>(matrix.il, matrix.jl, nullptr, nullptr, matrix.ll, nullptr, matrix.di);
		eps = 1E-30;
		iterC = 10000;
	}
	ThreeStageBase(std::string &path)
	{
		std::ifstream file(path + "/kuslau.txt");
		size_t size;
		file >> size;
		file >> iterC;
		file >> eps;
		file.close();

		file.open(path + "/di.txt");
		matrix.di->resize(size);
		f = new std::vector<T>(size);
		for (size_t i = 0; i < size; i++)
		{
			file >> (*matrix.di)[i];
		}
		file.close();

		file.open(path + "/ggl.txt");
		if (file.is_open())
		{
			symmetry = false;
			T buff = T();
			while (file >> buff)
			{
				file >> buff;
				matrix.ll->push_back(buff);
			}
			file.close();

			file.open(path + "/ggu.txt");
			while (file >> buff)
			{
				file >> buff;
				matrix.lu->push_back(buff);
			}
		}
		else
		{
			T buff = T();
			std::string str;
			file.open(path + "/gg.txt");
			char a;
			while (file >> buff)
			{
				matrix.ll->push_back(buff);
			}
			file.close();
		}

		size_t buff;
		file.open(path + "/ig.txt");
		while (file >> buff)
		{
			matrix.il->push_back(buff - 1);
		}
		file.close();

		file.open(path + "/jg.txt");
		while (file >> buff)
		{
			matrix.jl->push_back(buff - 1);
		}
		file.close();

		f->resize(size);
		file.open(path + "/pr.txt");
		for (size_t i = 0; i < size; i++)
		{
			file >> (*f)[i];
		}
		file.close();
		x = new std::vector<T>(size, 0);
		if (symmetry)
		{
			
			A = new CMatrix<T>(matrix.il, matrix.jl,nullptr , nullptr, matrix.ll,nullptr , matrix.di);
		}
		else
		{
			A = new CMatrix<T>(matrix.il, matrix.jl, matrix.iu, matrix.ju, matrix.ll, matrix.lu, matrix.di);
		}
	}

	/// <summary>
	/// input stream with raw SLE. First line contains size
	/// </summary>
	/// <param name="input"></param>
	ThreeStageBase(std::istream &input)
	{
		T a = T();
		size_t i = 0, j = 0, size = 0;
		input >> size;
		f = new std::vector<T>;
		x = new std::vector<T>;
		f->reserve(size);
		x->resize(size, 0);
		matrix.il->reserve(size + 1);
		matrix.di->reserve(size);
		matrix.il->push_back(0);
		matrix.iu->resize(size + 1, 0);
		std::vector<T> *lu = new std::vector<T>();
		std::vector<size_t> temp;

		while (!input.eof())
		{
			input >> a;
			if (input.peek() != '\n' && input.peek() != -1)
			{
				if (i == j)//diag
				{
					matrix.di->push_back(a);
				}
				else if (a != 0)
				{
					if (j > i)//lu
					{
						lu->push_back(a);
						matrix.ju->push_back(i);
						temp.push_back(j);
						(*matrix.iu)[j+1]++;
					}
					else//ll
					{
						matrix.ll->push_back(a);
						matrix.jl->push_back(j);
					}
				}
				j++;
			}
			else
			{
				f->push_back(a);
				i++;
				j = 0;
				matrix.il->push_back(matrix.ll->size());
			}
		}

		for(size_t i = 1; i < matrix.iu->size();i++)
		{
			(*matrix.iu)[i] += (*matrix.iu)[i - 1];
		}


		std::vector<size_t> idx(temp.size());
		std::iota(idx.begin(), idx.end(), 0);

		std::stable_sort(idx.begin(), idx.end(),
			[&](size_t a, size_t b) {
				return (temp)[a] < (temp)[b];
			});

		for (size_t i = 0; i < idx.size(); ++i)
		{
			if (idx[i] != i)
			{
				for (size_t j = i + 1; j < idx.size(); j++)
				{
					if (idx[j] == i)
					{
						std::swap((*lu)[idx[i]], (*lu)[idx[j]]);
						std::swap((*matrix.ju)[idx[i]], (*matrix.ju)[idx[j]]);
						std::swap(idx[i], idx[j]);
						continue;
					}

				}
			}

		}

		for (size_t i = 0; i < lu->size(); i++)
		{
			if (lu->size() != matrix.ll->size() || (*lu)[i] != (*matrix.ll)[i])
			{
				symmetry = false;
				break;
			}
		}
		if (symmetry)
		{
			delete lu;
			delete matrix.lu;
			matrix.lu = nullptr;
			delete matrix.ju;
			matrix.ju = nullptr;
			delete matrix.iu;
			matrix.iu = nullptr;
			A = new CMatrix<T>(matrix.il, matrix.jl, nullptr, nullptr, matrix.ll, nullptr, matrix.di);
		}
		else
		{
			delete matrix.lu;
			matrix.lu = lu;
			A = new CMatrix<T>(matrix.il, matrix.jl, matrix.iu, matrix.ju, matrix.ll, matrix.lu, matrix.di);
		}
		eps = 1E-15;
		iterC = 10000;
	}

	void forwSolutionCSR(std::vector<size_t> &ia, std::vector<size_t> &ja, std::vector<T> &f, std::vector<T> &x, std::vector<T> &al, bool loverD = true)
	{
		if (symmetry)
		{
			for (size_t i = 0; i < f.size(); i++)
			{
				T buff = T();
				for (size_t k = (ia)[i]; k < (ia)[i + 1]; k++)
				{
					buff += al[k] * x[(ja)[k]];
					
				}
				x[i] = (f[i] - buff) / (*matrix.di)[i];

			}
		}
		else
		{
			for (size_t i = 0; i < f.size(); i++)
			{
				T buff = T();
				for (size_t k = ia[i]; k < ia[i + 1]; k++)
				{
					buff += al[k] * x[ja[k]];

				}
				x[i] = (f[i] - buff) / (loverD ? 1 : (*matrix.di)[i]);

			}
		}

	}
	void revrsSolutionCSC(std::vector<size_t> &ia, std::vector<size_t> &ja, std::vector<T> &f, std::vector<T> &x, std::vector<T> &lu, bool loverD = true)
	{
		if (symmetry)
		{
			for (size_t i = 0; i < x.size(); i++)
			{
				x[i] = f[i] / (*matrix.di)[i];

			}
			for (size_t j = f.size() - 1; j < f.size(); j--)
			{

				for (size_t k = (ia)[j]; k < (ia)[j + 1]; k++)
				{
					x[ja[k]] -= lu[k] * x[j] / (*matrix.di)[ja[k]];

				}

			}
		}
		else
		{
			for (size_t i = 0; i < x.size(); i++)
			{
				x[i] = f[i] / (loverD ? 1 : (*matrix.di)[i]);

			}
			for (size_t j = f.size() - 1; j < f.size(); j--)
			{
				for (size_t k = (ia)[j]; k < (ia)[j + 1]; k++)
				{
					x[ja[k]] -= lu[k] * x[j] / (loverD ? 1 : (*matrix.di)[ja[k]]);

				}

			}

		}
	}
	


	void diagMult(std::vector<T> &diag, std::vector<T> &f, std::vector<T> &x)
	{
		for (size_t i = 0; i < diag.size(); i++)
		{
			x[i] = f[i] * diag[i];
		}
	}
	void diagSolve(std::vector<T> &diag, std::vector<T> &f, std::vector<T> &x)
	{
		for (size_t i = 0; i < diag.size(); i++)
		{
			x[i] = f[i] / diag[i];
		}
	}

	void incompLU()
	{

		matrix.di = new std::vector<T>(*A->di);
		matrix.ll = new std::vector<T>(*A->ll);
		matrix.lu = new std::vector<T>(*A->lu);

		for (size_t i = 0; i < f->size(); i++)
		{
			//L
			for (size_t k = (*matrix.il)[i]; k < (*matrix.il)[i + 1]; k++)
			{
				size_t jl = (*matrix.jl)[k];
				T buff = T();
				for (size_t il = (*matrix.il)[i]; il < k; il++)
				{
					size_t tjl = (*matrix.jl)[il];

					for (size_t kl = (*matrix.iu)[jl]; kl < (*matrix.iu)[jl + 1]; kl++)
					{
						if ((*matrix.ju)[kl] > jl)
						{
							break;
						}
						if ((*matrix.ju)[kl] == tjl)
						{
							buff += (*matrix.lu)[kl] * (*matrix.ll)[il];
							break;
						}
					}
				}

				(*matrix.ll)[k] = ((*matrix.ll)[k] - buff) / (*matrix.di)[jl];
			}

			//U
			for (size_t k = (*matrix.iu)[i]; k < (*matrix.iu)[i + 1]; k++)
			{
				T buff = T();
				size_t ju = (*matrix.ju)[k];
				for (size_t iu = (*matrix.iu)[i]; iu < k; iu++)
				{
					size_t tju = (*matrix.ju)[iu];
					for (size_t il = (*matrix.il)[ju]; il < (*matrix.il)[ju + 1]; il++)
					{
						if ((*matrix.jl)[il] == tju)
						{
							buff += (*matrix.lu)[iu] * (*matrix.ll)[il];
							break;
						}
					}
				}
				(*matrix.lu)[k] = (*matrix.lu)[k] - buff;
			}

			//diag
			T diag = T();
			for (size_t il = (*matrix.il)[i]; il < (*matrix.il)[i + 1]; il++)
			{
				size_t jl = (*matrix.jl)[il];
				for (size_t iu = (*matrix.iu)[i]; iu < (*matrix.iu)[i + 1]; iu++)
				{
					if ((*matrix.ju)[iu] == jl)
					{
						diag += (*matrix.lu)[iu] * (*matrix.ll)[il];
						break;
					}
				}
			}
			(*matrix.di)[i] = (*matrix.di)[i] - diag;

		}
	}
	void incompChol()
	{
		matrix.di = new std::vector<T>(*(A->di));
		matrix.ll = new std::vector<T>(*(A->ll));


		auto it = matrix.il->begin();

		for (size_t i = 0; i < f->size(); i++)
		{
			T diag = 0;

			for (size_t k = *it; k < *(it + 1); k++)
			{
				size_t j = (*matrix.jl)[k];
				T buff = T();
				size_t p1 = *it, p2 = (*matrix.il)[j];
				while (p1 < *(it + 1) && p2 < (*matrix.il)[j + 1])
				{
					size_t c1 = (*matrix.jl)[p1];
					size_t c2 = (*matrix.jl)[p2];
					if (c1 == c2)
					{
						buff += (*matrix.ll)[p1] * (*matrix.ll)[p2];
						p1++;
						p2++;
					}
					else if (c1 < c2)
					{
						p1++;
					}
					else
					{
						p2++;
					}
				}

				(*matrix.ll)[k] = ((*matrix.ll)[k] - buff) / (*matrix.di)[j];
				diag += (*matrix.ll)[k] * (*matrix.ll)[k];
			}
			it++;
			(*matrix.di)[i] = sqrt((*matrix.di)[i] - diag);

		}
	}
	void multiplyA(std::vector<T> &x)
	{
		if (symmetry)
		{
			std::vector<T> tx(x.size(), 0);
			for (size_t i = 0; i < x.size(); i++)
			{
				tx[i] += (*A->di)[i] * x[i];
				for (size_t k = (*A->il)[i]; k < (*A->il)[i + 1]; ++k)
				{
					tx[i] += (*A->ll)[k] * x[(*A->jl)[k]];
				}
			}
			for (int j = 0; j < x.size(); ++j)
			{
				for (int k = (*A->il)[j]; k < (*A->il)[j + 1]; ++k)
				{
					tx[(*A->jl)[k]] += (*A->ll)[k] * x[j];
				}
			}
			x = tx;
		}
		else
		{
			std::vector<T> tx(x);
			x.clear();
			x.resize(tx.size());
			for (size_t i = 0; i < x.size(); i++)
			{
				x[i] += (*A->di)[i] * tx[i];
				for (size_t k = (*A->il)[i]; k < (*A->il)[i + 1]; k++)
				{
					x[i] += (*A->ll)[k] * tx[(*A->jl)[k]];
				}

				for (size_t k = (*A->iu)[i]; k < (*A->iu)[i + 1]; ++k)
				{
					x[(*A->ju)[k]] += (*A->lu)[k] * tx[i];
				}

			}
		}
	}
	void multiplyA(const std::vector<T> &x, std::vector<T> &res)
	{
		res.clear();
		res.resize(x.size());
		if (symmetry)
		{
			for (size_t i = 0; i < x.size(); i++)
			{
				res[i] += (*A->di)[i] * x[i];
				for (size_t k = (*A->il)[i]; k < (*A->il)[i + 1]; ++k)
				{
					res[i] += (*A->ll)[k] * x[(*A->jl)[k]];
				}
			}
			for (int j = 0; j < x.size(); ++j)
			{
				for (int k = (*A->il)[j]; k < (*A->il)[j + 1]; ++k)
				{
					res[(*A->jl)[k]] += (*A->ll)[k] * x[j];
				}
			}
		}
		else
		{

			for (size_t i = 0; i < x.size(); i++)
			{
				res[i] += (*A->di)[i] * x[i];
				for (size_t k = (*A->il)[i]; k < (*A->il)[i + 1]; k++)
				{
					res[i] += (*A->ll)[k] * x[(*A->jl)[k]];
				}

				for (size_t k = (*A->iu)[i]; k < (*A->iu)[i + 1]; ++k)
				{
					res[(*A->ju)[k]] += (*A->lu)[k] * x[i];
				}

			}
		}
	}
	void multiplyAT(std::vector<T> &x)
	{

		if (symmetry)
		{
			multiplyA(x);
		}
		else
		{
			std::vector<T> tx(x.size(), 0);
			for (size_t i = 0; i < x.size(); i++)
			{
				tx[i] += (*A->di)[i] * x[i];

				for (size_t k = (*A->il)[i]; k < (*A->il)[i + 1]; k++)
				{
					tx[(*A->jl)[k]] += (*A->ll)[k] * x[i];
				}

				for (size_t k = (*matrix.iu)[i]; k < (*matrix.iu)[i + 1]; ++k)
				{
					tx[i] += (*A->lu)[k] * x[(*A->ju)[k]];
				}

			}
			x = tx;
		}
	}

public:
	void setF()
	{
		multiplyA(*f);
	}
	void OutCSR(std::ostream &stream)
	{
		stream << "il: ";
		for (auto a = matrix.il->begin(); a != matrix.il->end(); a++)
		{
			stream << *a << " ";

		}
		stream << std::endl;

		stream << "jl: ";
		for (auto a = matrix.jl->begin(); a != matrix.jl->end(); a++)
		{
			stream << *a << " ";

		}
		stream << std::endl;

		stream << "iu: ";
		if (!symmetry)
		{
			for (auto a = matrix.iu->begin(); a != matrix.iu->end(); a++)
			{
				stream << *a << " ";

			}
		}
		else
		{
			std::cout << "nullptr";
		}
		stream << std::endl;

		if (!symmetry)
		{
			stream << "ju: ";
			for (auto a = matrix.ju->begin(); a != matrix.ju->end(); a++)
			{
				stream << *a << " ";

			}
		}
		else
		{
			std::cout << "nullptr";
		}
		stream << std::endl;

		stream << "ll: ";
		for (auto a = matrix.ll->begin(); a != matrix.ll->end(); a++)
		{
			stream << *a << " ";

		}
		stream << std::endl;

		if (!symmetry)
		{
			stream << "lu: ";
			for (auto a = matrix.lu->begin(); a != matrix.lu->end(); a++)
			{
				stream << *a << " ";

			}
		}
		else
		{
			std::cout << "nullptr";
		}
		stream << std::endl;

		stream << "di: ";
		for (auto a = matrix.di->begin(); a != matrix.di->end(); a++)
		{
			stream << *a << " ";

		}
		stream << std::endl;

		stream << "f: ";
		for (auto a = f->begin(); a != f->end(); a++)
		{
			stream << *a << " ";

		}
		stream << std::endl;
		stream << "x: ";
		for (auto a = x->begin(); a != x->end(); a++)
		{
			stream << *a << " ";

		}
		stream << std::endl;


	}
	void outX()
	{
		for (auto &a : *x)
		{
			std::cout << std::setprecision(16) << a << std::endl;
		}
	}
};



