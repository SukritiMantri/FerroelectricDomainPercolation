//
//  GrainBoundaries.h
//  DomainPercolationFerroelectrics
//
//  Created by Sukriti Mantri on 2/10/26.
//

#pragma once

// Detect grain boundary points and assign IDs
int findGrainBoundaries(
    int steps,
    double* belong,
    double* xCoor,
    double* yCoor,
    double* zCoor,
    double* gbX,
    double* gbY,
    double* gbZ,
    double* gb_GID1,
    double* gb_GID2,
    double* gb_ID
);

// Fit grain-boundary planes (A x + B y + z + C = 0)
int fitGrainBoundaryPlanes(
    int gbPOINTno,
    double* gbX,
    double* gbY,
    double* gbZ,
    double* gb_GID1,
    double* gb_GID2,
    double* gb_ID,
    double* A,
    double* B,
    double* C,
    double* meangbX,
    double* meangbY,
    double* meangbZ
);

// Compute GB normals from plane coefficients
int computeGrainBoundaryNormals(
    int gbPOINTno,
    double* A,
    double* B,
    double* C,
    double* gb_GID1,
    double* gb_GID2,
    double* GBx,
    double* GBy,
    double* GBz,
    double* GBID1,
    double* GBID2
);
