//
//  IO.h.h
//  DomainPercolationFerroelectrics
//
//  Created by Sukriti Mantri on 2/10/26.
//

#pragma once

#include <string>

// Write binary array to file
void writeBinaryFile(
    const std::string& filename,
    double* data,
    int count
);
