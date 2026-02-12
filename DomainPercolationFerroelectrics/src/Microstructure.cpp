//
//  Untitled.cpp
//  DomainPercolationFerroelectrics
//
//  Created by Sukriti Mantri on 2/10/26.
//
#include "Microstructure.h"
#include <cmath>
#include <cstdlib>
#include <iostream>

using namespace std;

void generateCoordinates(
    int steps,
    double* xCoor,
    double* yCoor,
    double* zCoor
)
{
    int coor = 0;

    for(int i = 0; i < steps; i++)
    {
        for(int j = 0; j < steps; j++)
        {
            for(int k = 0; k < steps; k++)
            {
                xCoor[coor] = k;
                yCoor[coor] = j;
                zCoor[coor] = i;
                coor++;
            }
        }
    }
}

void generateGrainCenters(
    int grainNo,
    int steps,
    double* grainCX,
    double* grainCY,
    double* grainCZ
)
{
    double vol = cbrt((steps * steps * steps) / grainNo);

    int minimum = int(steps * 0.1);
    int total = steps - 2 * minimum;

    int i = 0;
    do
    {
        grainCX[i] = rand() % total + minimum;
        grainCY[i] = rand() % total + minimum;
        grainCZ[i] = rand() % total + minimum;

        int judge = 0;

        if(i > 0)
        {
            for(int j = 0; j < i; j++)
            {
                double disG =
                    sqrt(pow(grainCX[j] - grainCX[i], 2) +
                         pow(grainCY[j] - grainCY[i], 2) +
                         pow(grainCZ[j] - grainCZ[i], 2));

                if(disG < (vol / 5))
                {
                    judge = 0;
                    break;
                }
                else
                    judge = 1;
            }
        }
        else
            judge = 1;

        if(judge == 1)
            i++;

    } while(i < grainNo);
}


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
)
{
    double* distance = new double[grainNo];

    for(int coor = 0; coor < steps * steps * steps; coor++)
    {
        for(int grID = 0; grID < grainNo; grID++)
        {
            distance[grID] =
                sqrt(pow(xCoor[coor] - grainCX[grID], 2) +
                     pow(yCoor[coor] - grainCY[grID], 2) +
                     pow(zCoor[coor] - grainCZ[grID], 2));
        }

        double smallest = distance[0];

        for(int grID = 0; grID < grainNo; grID++)
        {
            if(distance[grID] <= smallest)
            {
                smallest = distance[grID];
                belong[coor] = grID;
            }
        }
    }

    delete[] distance;
}


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
)
{
    for(int i = 0; i < grainNo; i++)
    {
        grainCenterCx[i] = 0;
        grainCenterCy[i] = 0;
        grainCenterCz[i] = 0;
        pointbelongtograin[i] = 0;
    }

    for(int coor = 0; coor < steps * steps * steps; coor++)
    {
        int id = belong[coor];
        grainCenterCx[id] += xCoor[coor];
        grainCenterCy[id] += yCoor[coor];
        grainCenterCz[id] += zCoor[coor];
        pointbelongtograin[id]++;
    }

    for(int grID = 0; grID < grainNo; grID++)
    {
        grainCenterCx[grID] =
            round(grainCenterCx[grID] / pointbelongtograin[grID]);
        grainCenterCy[grID] =
            round(grainCenterCy[grID] / pointbelongtograin[grID]);
        grainCenterCz[grID] =
            round(grainCenterCz[grID] / pointbelongtograin[grID]);
    }
}


void computeGrainSizeAndRadius(
    int steps,
    int grainNo,
    double* belong,
    double* grainSize,
    double* grainRad
)
{
    for(int grID = 0; grID < grainNo; grID++)
    {
        grainSize[grID] = 0;
        for(int i = 0; i < steps * steps * steps; i++)
        {
            if(belong[i] == grID)
                grainSize[grID]++;
        }
    }

    for(int i = 0; i < grainNo; i++)
    {
        grainRad[i] = cbrt(grainSize[i]);
    }
}



