//
//  Percolation.h
//  DomainPercolationFerroelectrics
//
//  Created by Sukriti Mantri on 2/10/26.
//

#pragma once

void runPercolationForMicrostructure(
    int grainNo,
    int Iterationseclevel,
    int countGB,
    double chTol,
    double tolerance,
    double* GBx,
    double* GBy,
    double* GBz,
    double* GBID1,
    double* GBID2,
    double* EA1,
    double* EA2,
    double* EA3,
    double& meanvaluenet,
    double& meanvaluetotalnet,
    double& standerddeviationmax,
    double& standerddeviationtot
);
