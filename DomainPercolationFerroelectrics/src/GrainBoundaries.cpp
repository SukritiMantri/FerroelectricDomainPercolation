//
//  GrainBoundaries.cpp
//  DomainPercolationFerroelectrics
//
//  Created by Sukriti Mantri on 2/10/26.
//

#include "GrainBoundaries.h"
#include <cmath>
#include <cstdio>

using namespace std;
// finding neighbouring grains
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
)
{
    int coor = -1;

    int coor1, coor2, coor3, coor4;

    for(int i = 0; i < steps - 1; i++)
    {
        for(int j = 0; j < steps - 1; j++)
        {
            for(int k = 0; k < steps - 1; k++)
            {
                coor1 = i * steps * steps + j * steps + k;
                coor2 = (i + 1) * steps * steps + j * steps + k;
                coor3 = i * steps * steps + (j + 1) * steps + k;
                coor4 = i * steps * steps + j * steps + (k + 1);

                if(belong[coor1] != belong[coor2])
                {
                    coor++;
                    gbX[coor] = (xCoor[coor1] + xCoor[coor2]) / 2;
                    gbY[coor] = (yCoor[coor1] + yCoor[coor2]) / 2;
                    gbZ[coor] = (zCoor[coor1] + zCoor[coor2]) / 2;
                    gb_GID1[coor] = belong[coor1];
                    gb_GID2[coor] = belong[coor2];

                    double intpart, fracpart;
                    double p1 = fmax(gb_GID1[coor], gb_GID2[coor]) / 10;
                    fracpart = modf(p1, &intpart);
                    double g1i = intpart, g1f = fracpart;

                    double p2 = fmin(gb_GID1[coor], gb_GID2[coor]) / 10;
                    fracpart = modf(p2, &intpart);
                    double g2i = intpart, g2f = fracpart;

                    gb_ID[coor] =
                        (g2i * 1000) + (g2f * 1000) +
                        (g1i * 10) + (g1f * 10);
                }

                if(belong[coor1] != belong[coor3])
                {
                    coor++;
                    gbX[coor] = (xCoor[coor1] + xCoor[coor3]) / 2;
                    gbY[coor] = (yCoor[coor1] + yCoor[coor3]) / 2;
                    gbZ[coor] = (zCoor[coor1] + zCoor[coor3]) / 2;
                    gb_GID1[coor] = belong[coor1];
                    gb_GID2[coor] = belong[coor3];

                    double intpart, fracpart;
                    double p1 = fmax(gb_GID1[coor], gb_GID2[coor]) / 10;
                    fracpart = modf(p1, &intpart);
                    double g1i = intpart, g1f = fracpart;

                    double p2 = fmin(gb_GID1[coor], gb_GID2[coor]) / 10;
                    fracpart = modf(p2, &intpart);
                    double g2i = intpart, g2f = fracpart;

                    gb_ID[coor] =
                        (g2i * 1000) + (g2f * 1000) +
                        (g1i * 10) + (g1f * 10);
                }

                if(belong[coor1] != belong[coor4])
                {
                    coor++;
                    gbX[coor] = (xCoor[coor1] + xCoor[coor4]) / 2;
                    gbY[coor] = (yCoor[coor1] + yCoor[coor4]) / 2;
                    gbZ[coor] = (zCoor[coor1] + zCoor[coor4]) / 2;
                    gb_GID1[coor] = belong[coor1];
                    gb_GID2[coor] = belong[coor4];

                    double intpart, fracpart;
                    double p1 = fmax(gb_GID1[coor], gb_GID2[coor]) / 10;
                    fracpart = modf(p1, &intpart);
                    double g1i = intpart, g1f = fracpart;

                    double p2 = fmin(gb_GID1[coor], gb_GID2[coor]) / 10;
                    fracpart = modf(p2, &intpart);
                    double g2i = intpart, g2f = fracpart;

                    gb_ID[coor] =
                        (g2i * 1000) + (g2f * 1000) +
                        (g1i * 10) + (g1f * 10);
                }
            }
        }
    }

    return coor + 1;
}

/////Finding the coefficients to calculate the grain boundaries plane normal

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
)
{
    int count = -1;

    for(int i = 0; i < gbPOINTno; i++)
    {
        if(!isnan(gb_ID[i]))
        {
            count++;
            double gbid = gb_ID[i];

            double meanX = gbX[i];
            double meanY = gbY[i];
            double meanZ = gbZ[i];
            double trial = 1;

            for(int j = i + 1; j < gbPOINTno; j++)
            {
                if(gb_ID[j] == gbid)
                {
                    trial++;
                    meanX += gbX[j];
                    meanY += gbY[j];
                    meanZ += gbZ[j];
                }
            }

            meanX /= trial;
            meanY /= trial;
            meanZ /= trial;

            meangbX[i] = meanX;
            meangbY[i] = meanY;
            meangbZ[i] = meanZ;

            double squareX = pow(gbX[i] - meanX, 2);
            double squareY = pow(gbY[i] - meanY, 2);
            double squareZ = pow(gbZ[i] - meanZ, 2);
            double multXY = (gbX[i] - meanX) * (gbY[i] - meanY);
            double multXZ = (gbX[i] - meanX) * (gbZ[i] - meanZ);
            double multYZ = (gbY[i] - meanY) * (gbZ[i] - meanZ);

            for(int j = i + 1; j < gbPOINTno; j++)
            {
                if(gb_ID[j] == gbid)
                {
                    squareX += pow(gbX[j] - meanX, 2);
                    squareY += pow(gbY[j] - meanY, 2);
                    squareZ += pow(gbZ[j] - meanZ, 2);
                    multXY += (gbX[j] - meanX) * (gbY[j] - meanY);
                    multXZ += (gbX[j] - meanX) * (gbZ[j] - meanZ);
                    multYZ += (gbY[j] - meanY) * (gbZ[j] - meanZ);

                    gb_ID[j] = NAN;
                    gb_GID1[j] = NAN;
                    gb_GID2[j] = NAN;
                }
            }

            B[i] = (multXY * multXZ - squareX * multYZ) /
                   (multXY * multXY - squareX * squareY);
            A[i] = (squareY * multXZ - multXY * multYZ) /
                   (squareX * squareY - multXY * multXY);

            C[i] = meanZ - A[i] * meanX - B[i] * meanY;
        }
        else
        {
            A[i] = B[i] = C[i] = NAN;
            meangbX[i] = meangbY[i] = meangbZ[i] = NAN;
        }
    }

    return count + 1;
}



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
)
{
    int count = 0;

    for(int i = 0; i < gbPOINTno; i++)
    {
        if(!isnan(A[i]))
        {
            GBx[count] = A[i];
            GBy[count] = B[i];
            GBz[count] = -1;
            GBID1[count] = gb_GID1[i];
            GBID2[count] = gb_GID2[i];
            count++;
        }
    }

    return count;
}

