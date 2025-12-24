#pragma once
#include <string>
#include <fstream>
#include <vector>
#include <array>
#include <algorithm>
#include <iomanip>
#include "Functions.h"
struct Node
{
	std::array<double,3> cords;

	double& operator[] (size_t i)
	{
		return cords[i];
	}
	Node() = default;
	Node(double x, double y, double z)
	{
		cords[0] = x;
		cords[1] = y;
		cords[2] = z;
	}

	bool operator== (Node &a)
	{
		for(size_t i = 0; i < 3; i++)
		{
			if (cords[i] != a[i])
			{
				return false;
			}
		}

		return true;
	}

};
struct Elems
{
	std::array<size_t, 4> cords;

	size_t &operator[] (size_t i)
	{
		return cords[i];
	}

};

class Translator
{
	std::vector<Node> nodes;
	std::vector<Elems> elems;
	double lambda = 1;
	double beta = 1;
	std::ofstream config;
	bool onEdge(std::array<size_t,3> &heads,std::array<double,3> &normal)
	{
		for (size_t i = 0; i < 3; i++)
		{
			if (nodes[heads[0]][i] == 0.5 && nodes[heads[1]][i] == 0.5 && nodes[heads[2]][i] == 0.5)
			{
				normal[i] = 1;
				return true;
			}
			else if(nodes[heads[0]][i] == 0 && nodes[heads[1]][i] == 0 && nodes[heads[2]][i] == 0)
			{
				normal[i] = -1;
				return true;
			}
		}
		return false;
	}


public:
	Translator(std::string &in,std::string &outDir)
	{
		config.open(outDir + "/config.txt");
		std::ifstream input(in);
		std::string buff;
		for (size_t i = 0; i < 8; i++)
		{
			std::getline(input, buff);
		}
		size_t count = 0;
		input >> count;
		nodes.resize(count);

		for(size_t i = 0; i < count; i++)
		{
			double buff = 0;
			input >> buff;
			input >> nodes[i][0];
			input >> nodes[i][1];
			input >> nodes[i][2];
		}

		for (size_t i = 0; i < 3; i++)
		{
			std::getline(input, buff);
		}

		input >> count;
		elems.resize(count);
		for (size_t i = 0; i < count; i++)
		{
			size_t node = 0;
			input >> node;
			input >> node;
			input >> node;
			input >> node;
			input >> node;

			input >> elems[i][0];
			elems[i][0]--;
			input >> elems[i][1];
			elems[i][1]--;
			input >> elems[i][2];
			elems[i][2]--;
			input >> elems[i][3];
			elems[i][3]--;

			std::sort(elems[i].cords.begin(), elems[i].cords.end());
		}
		input.close();
		double gamma = 0;
		config << gamma << std::endl;
		config << nodes.size() << std::endl;
		config << elems.size() << std::endl;
	}

	template <typename F>
	void setDirih(std::string &outDir,F f)
	{
		size_t count = 0;
		std::ofstream out(outDir + "/dirih.txt");
		for(size_t i = 0; i < nodes.size();i++)
		{
			if(nodes[i][2] == 0)
			{
				out << i << " " << std::setprecision(16) << f(nodes[i].cords) << std::endl;
				count++;
			}
		}
		out.close();
		config << count << std::endl;
	}
	template <typename F>
	void setNeuman(std::string &outDir,F guf)
	{
		size_t count = 0;
		std::ofstream file(outDir + "/neumann.txt");
		for(size_t i = 0; i < elems.size();i++)
		{
			std::array<size_t, 3> chain;
			std::array<double, 3> normal = {0,0,0};
			auto &element = elems[i];
			chain = { element[0],element[1],element[2] };
			if(onEdge(chain,normal))
			{
				file << chain[0] << " " << chain[1] << " " << chain[2] << " " << i << " ";
				std::string fun;
				guf(normal,fun);
				fun = std::to_string(lambda) + "*(" + fun + ")";
				file << fun;
				file << std::endl;
				count++;
			}
			normal = { 0,0,0 };
			chain = { element[0],element[1],element[3] };
			if (onEdge(chain, normal))
			{
				file << chain[0] << " " << chain[1] << " " << chain[2] << " " << i << " ";
				std::string fun;
				guf(normal, fun);
				fun = std::to_string(lambda) + "*(" + fun + ")";
				file << fun;
				file << std::endl;
				count++;
			}
			normal = { 0,0,0 };
			chain = { element[0],element[2],element[3] };
			if (onEdge(chain, normal))
			{
				file << chain[0] << " " << chain[1] << " " << chain[2] << " " << i << " ";
				std::string fun;
				guf(normal, fun);
				fun = std::to_string(lambda) + "*(" + fun + ")";
				file << fun;
				file << std::endl;
				count++;
			}
			normal = { 0,0,0 };
			chain = { element[1],element[2],element[3] };
			if (onEdge(chain, normal))
			{
				file << chain[0] << " " << chain[1] << " " << chain[2] << " " << i << " ";
				std::string fun;
				guf(normal, fun);
				fun = std::to_string(lambda) + "*(" + fun + ")";
				file << fun;
				file << std::endl;
				count++;
			}
		}
		file.close();
		config << count << std::endl;
	}
	template <typename F, typename Fs>
	void setRobin(std::string &outDir, F guf,Fs uStr)
	{ 
		size_t count = 0;
		std::ofstream file(outDir + "/robin.txt");
		for (size_t i = 0; i < elems.size(); i++)
		{
			std::array<size_t, 3> chain;
			std::array<double, 3> normal = { 0,0,0 };
			auto &element = elems[i];
			chain = { element[0],element[1],element[2] };
			if (onEdge(chain, normal))
			{
				file << chain[0] << " " << chain[1] << " " << chain[2] << " " << i << " " << beta << " ";
				std::string fun;
				guf(normal, fun);
				fun = std::to_string(lambda) + "*(" + fun + ")+";
				file << fun;
				uStr(fun);
				fun = std::to_string(beta) + "*(" + fun + ")";
				file << fun;
				file << std::endl;
				count++;
			}
			normal = { 0,0,0 };
			chain = { element[0],element[1],element[3] };
			if (onEdge(chain, normal))
			{
				file << chain[0] << " " << chain[1] << " " << chain[2] << " " << i << " ";
				std::string fun;
				guf(normal, fun);
				fun = std::to_string(lambda) + "*(" + fun + ")+";
				file << fun;
				uStr(fun);
				fun = std::to_string(beta) + "*(" + fun + ")";
				file << fun;
				file << std::endl;
				count++;
			}
			normal = { 0,0,0 };
			chain = { element[0],element[2],element[3] };
			if (onEdge(chain, normal))
			{
				file << chain[0] << " " << chain[1] << " " << chain[2] << " " << i << " ";
				std::string fun;
				guf(normal, fun);
				fun = std::to_string(lambda) + "*(" + fun + ")+";
				file << fun;
				uStr(fun);
				fun = std::to_string(beta) + "*(" + fun + ")";
				file << fun;
				file << std::endl;
				count++;
			}
			normal = { 0,0,0 };
			chain = { element[1],element[2],element[3] };
			if (onEdge(chain, normal))
			{
				file << chain[0] << " " << chain[1] << " " << chain[2] << " " << i << " ";
				std::string fun;
				guf(normal, fun);
				fun = std::to_string(lambda) + "*(" + fun + ")+";
				file << fun;
				uStr(fun);
				fun = std::to_string(beta) + "*(" + fun + ")";
				file << fun;
				file << std::endl;
				count++;
			}

		}
		file.close();
		config << count << std::endl;
	}
	void setRobinZero()
	{
		config << 0 << std::endl;
	}


	template <typename F>
	void setNodes(std::string &outDir,F dguf)
	{
		std::ofstream f(outDir + "/f.txt");
		std::ofstream file(outDir + "/nodes.txt");
		for(auto &a : nodes)
		{
			file << a[0] << " " << a[1] << " " << a[2] << std::endl;
			f <<std::setprecision(16) << dguf(a.cords) << std::endl;
		}
		file.close();
		f.close();
	}
	void setElements(std::string &outDir)
	{
		std::ofstream file(outDir + "/elems.txt");
		for(auto &a : elems)
		{
			file << a[0] << " " << a[1] << " " << a[2] << " " << a[3] << " " << 2 << std::endl;
		
		}
		file.close();
	}


	~Translator()
	{
		config.close();
	}
};