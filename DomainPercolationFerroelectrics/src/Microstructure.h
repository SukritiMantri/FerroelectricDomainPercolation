//
//  Microstructure.h
//  DomainPercolationFerroelectrics
//
//  Created by Sukriti Mantri on 2/10/26.
//

#pragma once

// Coordinate grid
void generateCoordinates(
    int steps,
    double* xCoor,
    double* yCoor,
    double* zCoor
);

// Initial grain center selection
void generateGrainCenters(
    int grainNo,
    int steps,
    double* grainCX,
    double* grainCY,
    double* grainCZ
);

// Assign each point to nearest grain
void assignPointsToGrains(
    int steps,
    int grainNo,
    double* xCoor,
    double* yCoor,
    double* zCoor,
    double* grainCX,
    double* grainCY,
    double* grainCZ,
    double* belong
);

// Recompute grain centers from assigned points
void computeGrainCentersFromPoints(
    int steps,
    int grainNo,
    double* xCoor,
    double* yCoor,
    double* zCoor,
    double* belong,
    double* grainCenterCx,
    double* grainCenterCy,
    double* grainCenterCz,
    double* pointbelongtograin
);

// Grain size and radius
void computeGrainSizeAndRadius(
    int steps,
    int grainNo,
    double* belong,
    double* grainSize,
    double* grainRad
);
