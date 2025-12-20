#pragma once
#include <string>
#include <fstream>
#include <vector>
#include <array>
#include <algorithm>
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

	bool onEdge(std::array<size_t,3> &heads,std::array<double,3> &normal)
	{
		for (size_t i = 0; i < 3; i++)
		{
			if (nodes[heads[0]][i] == 2 && nodes[heads[1]][i] == 2 && nodes[heads[2]][i] == 2)
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
	Translator(std::string &in)
	{
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
	}

	template <typename F>
	void setDirih(std::string &outDir,F f)
	{
		std::vector<Node> heads =
		{
			Node(0,0,0),
			Node(0,0,2),
			Node(0,2,0),
			Node(0,2,2),
			Node(2,0,0),
			Node(2,0,2),
			Node(2,2,0),
			Node(2,2,2)
		};

		std::ofstream out(outDir + "/dirih.txt");
		for(size_t i = 0; i < nodes.size();i++)
		{
			for(size_t j = 0; j < heads.size();j++)
			{
				if(heads[j] == nodes[i])
				{
					out << i << " " << f(nodes[i].cords) << std::endl;
					break;
				}
			}
		}
		out.close();

	}
	template <typename F>
	void setNewman(std::string &outDir,F gu)
	{
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
				std::vector<double> poly;
				size_t count = gu(normal,poly);
				file << count;
				for(size_t k = 0; k < count; k++)
				{
					file << " " << poly[k];
				}
				file << std::endl;
			}
			normal = { 0,0,0 };
			chain = { element[0],element[1],element[3] };
			if (onEdge(chain, normal))
			{
				file << chain[0] << " " << chain[1] << " " << chain[2] << " " << i << " ";
				std::vector<double> poly;
				size_t count = gu(normal, poly);
				file << count;
				for (size_t k = 0; k < count; k++)
				{
					file << " " << poly[k];
				}
				file << std::endl;
			}
			normal = { 0,0,0 };
			chain = { element[0],element[2],element[3] };
			if (onEdge(chain, normal))
			{
				file << chain[0] << " " << chain[1] << " " << chain[2] << " " << i << " ";
				std::vector<double> poly;
				size_t count = gu(normal, poly);
				file << count;
				for (size_t k = 0; k < count; k++)
				{
					file << " " << poly[k];
				}
				file << std::endl;
			}
			normal = { 0,0,0 };
			chain = { element[1],element[2],element[3] };
			if (onEdge(chain, normal))
			{
				file << chain[0] << " " << chain[1] << " " << chain[2] << " " << i << " ";
				std::vector<double> poly;
				size_t count = gu(normal, poly);
				file << count;
				for (size_t k = 0; k < count; k++)
				{
					file << " " << poly[k];
				}
				file << std::endl;
			}



		}
		file.close();
	
	}
	template <typename F>
	void setNodes(std::string &outDir,F dgu)
	{
		std::ofstream f(outDir + "/f.txt");
		std::ofstream file(outDir + "/nodes.txt");
		for(auto &a : nodes)
		{
			file << a[0] << " " << a[1] << " " << a[2] << std::endl;
			f << dgu(a.cords) << std::endl;
		}
		file.close();
		f.close();
	}
	void setElements(std::string &outDir)
	{
		std::ofstream file(outDir + "/elems.txt");
		for(auto &a : elems)
		{
			file << a[0] << " " << a[1] << " " << a[2] << " " << a[3] << " " << 1 << std::endl;
		
		}
		file.close();
	}


};