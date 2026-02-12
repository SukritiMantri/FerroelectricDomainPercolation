//
//  EulerAngles.h
//  DomainPercolationFerroelectrics
//
//  Created by Sukriti Mantri on 2/10/26.
//

#pragma once

// Random texture
void assignRandomEulerAngles(
    int grainNo,
    double* EA1,
    double* EA2,
    double* EA3
);

// MD texture
void assignMDEulerAngles(
    int grainNo,
    double r,
    double f,
    int htex,
    int ktex,
    int ltex,
    int maxno,
    double* EA1,
    double* EA2,
    double* EA3
);

// Σ3 grain-boundary weighting
void applySigma3Weighting(
    int gbweightno,
    int countGB,
    double a1,
    double a2,
    double a3,
    double theta2,
    double* GBID1,
    double* GBID2,
    double* EA1,
    double* EA2,
    double* EA3
);
