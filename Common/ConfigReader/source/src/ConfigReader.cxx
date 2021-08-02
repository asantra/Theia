#include "ConfigReader.h"

using namespace std;

ConfigReader::ConfigReader() {}
ConfigReader::~ConfigReader() {}

void ConfigReader::readConfigFile(string configFileName)
{
    // std::ifstream is RAII, i.e. no need to call close
    std::ifstream cFile (configFileName);
    if (cFile.is_open()){
        std::string line;
        while(getline(cFile, line)){
            line.erase(std::remove_if(line.begin(), line.end(), [](unsigned char x){return std::isspace(x);}), line.end());
            while (line.length()==0 )                 
                getline(cFile,line);
            
            if( line.empty() || line[0] == '#' )
            {
                continue;
            }
            auto delimiterPos = line.find("=");
            auto name = line.substr(0, delimiterPos);
            auto value = line.substr(delimiterPos + 1);
            std::cout << name << " " << value << '\n';
        }
    }
    else {
        std::cerr << "Couldn't open config file for reading.\n";
    }
}
