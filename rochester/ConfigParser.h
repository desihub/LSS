#ifndef CONFIGPARSER
#define CONFIGPARSER
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <typeinfo>

using namespace std;

class ConfigParser
{
	private:

		map< string, string > info;

	public:
		ConfigParser(string filename);

		template<typename T>  T Get(string name)
		{
			T i;
			istringstream(info[name]) >> i;
			return(i);
		}

		template<typename T> vector<T> GetVector(string name)
		{
			vector<T> res;
			string buff;
			for(char l : info[name])
			{
				if(l == ',')
				{
					res.push_back(0);
					istringstream(buff) >> res.back();
					buff = "";
					continue;
				}
				buff.push_back(l);
			}
			res.push_back(0);
			istringstream(buff) >> res.back();
			return(res);
		}

};

#endif
