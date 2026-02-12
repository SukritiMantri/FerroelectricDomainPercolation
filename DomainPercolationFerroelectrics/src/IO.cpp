//
//  IO.cpp
//  DomainPercolationFerroelectrics
//
//  Created by Sukriti Mantri on 2/10/26.
//

#include "IO.h"
#include <fstream>

using namespace std;

void writeBinaryFile(
    const std::string& filename,
    double* data,
    int count
)
{
    ofstream file(filename, ios::out | ios::binary);
    if(file.is_open())
    {
        file.write((char*)data, count * sizeof(double));
        file.close();
    }
}
