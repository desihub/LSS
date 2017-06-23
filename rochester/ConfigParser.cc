#include <fstream>
#include <iostream>
#include <sstream>
#include "ConfigParser.h"

string cleanspaces(const string& str)
{
	string nstr;
	for(char l : str)
	{
		if(!isspace(l))
		{
			nstr.push_back(l);
		}
	}
	return(nstr);
}

ConfigParser::ConfigParser(string filename)
{
	fstream infile(filename, ios_base::in);
	string line;
	while(getline(infile, line))
	{
		size_t commentpos = line.find("#");
		if(commentpos != string::npos)
		{
			line = line.substr(0, commentpos);
		}

		size_t eqpos = line.find("=");
		if(eqpos == string::npos) {continue;}

		string parameter = line.substr(0, eqpos);
		string val = line.substr(eqpos+1);

		parameter = cleanspaces(parameter);
		val = cleanspaces(val);
		//cout << parameter << ":" << val << endl;

		info[parameter] = val;
	}
}

