#pragma once
#include "MISC.h"
#include <array>
#include <vector>
#include "Gaus.h"
template <Field T,size_t n>
struct Node
{
	std::array<T, n> node;
	Node() = default;
	template<std::input_iterator It>
	Node(It begin)
	{
		std::copy_n(begin, n, node.begin());
	}

	T&  operator [](size_t i) 
	{
		return node[i];
	}

};

template <Field T>
struct Elem
{
protected:
	Elem() = default;
	Elem(T h) : h(h) {}
public:
	T h;
	T V;

};

template <Field T>
struct Basis
{
	std::vector<T> basis;
	Basis() = default;
	Basis(const std::vector<T> &vec)
	{
		basis = std::vector<T>(vec);
	}
	T& operator[](size_t i)
	{
		return basis[i];
	}

};


template <Field T>
class Tetrahedron : public Elem<T>
{
public:
	using Elem<T>::V;
private:
	T Volume(std::vector<Node<T, 3>> &nodes)
	{
		Node<T, 3> head = nodes[vertices[0]];
		std::vector<Node<T, 3>> directions(3);
		for (uint8_t i = 0; i < 3; i++)
		{
			directions[i] = nodes[vertices[i+1]];
			for (size_t j = 0; j < 3; j++)
			{
				directions[i][j] = directions[i][j] - head[j];
			}
		}
		directions[0][0] = directions[0][0] * ((directions[1][1] * directions[2][2]) - (directions[2][1] * directions[1][2]));
		directions[0][1] = directions[0][1] * ((directions[1][0] * directions[2][2]) - (directions[2][0] * directions[1][2]));
		directions[0][2] = directions[0][2] * ((directions[1][0] * directions[2][1]) - (directions[2][0] * directions[1][1]));
		T res = abs((directions[0][0] - directions[0][1] + directions[0][2])) / 6;
		if(res == 0)
		{
			throw std::invalid_argument("V == 0");
		}
		return res;
	}

public:
	
	std::array<size_t,4> vertices;
	std::array<Basis<T>,4> basis;
	Tetrahedron() = default;

	template<std::input_iterator It>
	Tetrahedron(It it, T h, std::vector<Node<T, 3>> &nodes) : Elem<T>(h)
	{
		std::copy_n(it, vertices.size(), vertices.begin());
		V = Volume(nodes);
	}
	static constexpr size_t GetNodesCount()
	{
		return 4;
	}
	static constexpr size_t GetDim()
	{
		return 3;
	}

	auto getBasis(size_t node)
	{
		for(size_t i = 0; i < vertices.size();i++)
		{
			if(node == vertices[i])
			{
				return basis.begin() + i;
			}
		}
		return basis.end();
	}
	void setBasis(std::vector<Node<T, 3>> &nodes)
	{
		std::vector<std::vector<T>*> *M = new std::vector<std::vector<T>*>;
		M->reserve(4);
		for(size_t i = 0; i < vertices.size();i++)
		{
			std::vector<T> *buff = new std::vector<T>;
			buff->reserve(4);
			buff->push_back(1);
			for (size_t j = 0; j < 3; j++)
			{
				buff->push_back(nodes[vertices[i]][j]);
			}
			M->push_back(buff);
		}
		GaussAlgoth<T> SLAU(M);
		std::vector<std::vector<T>> B = { { 1,0,0,0 },{ 0,1,0,0 },{ 0,0,1,0 }, { 0,0,0,1 } };
		SLAU.transform(B);
		
		for(size_t i = 0; i < 4; i++)
		{
			SLAU.calcX(B[i]);
			basis[i] = Basis(B[i]);
		}
		if (M != nullptr)
		{
			for (size_t i = 0; i < M->size(); i++)
			{
				if((*M)[i] != nullptr)
				{
					delete (*M)[i];
				}

			}
			delete M;
		}
	
	}
	T scalGrad(const size_t a,const size_t b)
	{
		T res = 0;
		for(size_t i = 1; i < basis.size(); i++)
		{
			res += basis[a][i] * basis[b][i];
		}
		return res;
	}

};

template<typename T>
concept ElemType = requires(T e)
{
	{ T::GetNodesCount() } -> std::convertible_to<size_t>;
	{ T::GetDim() } -> std::convertible_to<size_t>;
};

template <Field T>
struct Neumann
{
	std::vector<size_t> nodes;
	T teta;
	size_t elemN;
};

template <Field T>
T triangleArea(Node<T,3> &A, Node<T, 3> &B, Node<T, 3> &C)
{
	std::array<T, 3> a, b,c;
	for(size_t i = 0; i < a.size();i++)
	{
		a[i] = B[i] - A[i];
		b[i] = C[i] - A[i];
	}

	c[0] = a[1] * b[2] - a[2] * b[1];
	c[1] = a[2] * b[0] - a[0] * b[2];
	c[2] = a[0] * b[1] - a[1] * b[0];

	return std::sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]) / 2;
}

template <Field T, ElemType El>
T integrate(const El &element,const std::array<size_t,3> &triangle, std::vector<Node<T, El::GetDim()>> &nodes,const T &S, std::vector<std::vector<double>> &polys)
{
	constexpr std::array<T, 4> xi = { 1.0 / 3.0, 0.2, 0.6, 0.2 };
	constexpr std::array<T, 4> eta = { 1.0 / 3.0, 0.2, 0.2, 0.6 };
	constexpr std::array<T, 4> w = { -27.0 / 96.0, 25.0 / 96.0, 25.0 / 96.0, 25.0 / 96.0 };

	T result = T();
	for(size_t i = 0; i < 4; i++)
	{
		T l2 = xi[i];
		T l3 = eta[i];
		T l1 = 1 - l2 - l3;
		T x = l1 * nodes[triangle[0]][0] + l2 * nodes[triangle[1]][0] + l3 * nodes[triangle[2]][0];
		T y = l1 * nodes[triangle[0]][1] + l2 * nodes[triangle[1]][1] + l3 * nodes[triangle[2]][1];
		T z = l1 * nodes[triangle[0]][2] + l2 * nodes[triangle[1]][2] + l3 * nodes[triangle[2]][2];

		T g = T(1);
		
		for(auto &a : polys)
		{
			g *= calcPoly(a,x, y, z);
		}
		

		result += w[i] * g;
	}



	return result * S * 2;
}

