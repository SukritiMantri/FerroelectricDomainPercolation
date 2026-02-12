//
//  EulerAngles.cpp
//  DomainPercolationFerroelectrics
//
//  Created by Sukriti Mantri on 2/10/26.
//

#include "EulerAngles.h"
#include <cmath>
#include <cstdlib>
#include <iostream>

using namespace std;

void assignRandomEulerAngles(
    int grainNo,
    double* EA1,
    double* EA2,
    double* EA3
)
{
    for(int i = 0; i < grainNo; i++)
    {
        EA1[i] = rand() % 101;
        EA2[i] = rand() % 101;
        EA3[i] = rand() % 101;
    }

    for(int i = 0; i < grainNo; i++)
    {
        EA1[i] = (M_PI / 2) * (EA1[i] / 100);
        EA2[i] = acos(EA2[i] / 100);
        EA3[i] = (M_PI / 2) * (EA3[i] / 100);
    }
}

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
)
{
    int i = 0;
    int count = 0;

    do
    {
        EA1[i] = rand() % 101;
        EA2[i] = rand() % 101;
        EA3[i] = rand() % 101;

        EA1[i] = (M_PI / 2) * (EA1[i] / 100);
        EA2[i] = acos(EA2[i] / 100);
        EA3[i] = (M_PI / 2) * (EA3[i] / 100);

        double phi = EA1[i];
        double psi = EA3[i];
        double theta = EA2[i];

        double Matrixcheck[3][3];

        Matrixcheck[0][0] = cos(psi)*cos(phi)-cos(theta)*sin(phi)*sin(psi);
        Matrixcheck[0][1] = cos(psi)*sin(phi)+cos(theta)*cos(phi)*sin(psi);
        Matrixcheck[0][2] = sin(psi)*sin(theta);
        Matrixcheck[1][0] = -sin(psi)*cos(phi)-cos(theta)*sin(phi)*cos(psi);
        Matrixcheck[1][1] = -sin(psi)*sin(phi)+cos(theta)*cos(phi)*cos(psi);
        Matrixcheck[1][2] = cos(psi)*sin(theta);
        Matrixcheck[2][0] = sin(theta)*sin(phi);
        Matrixcheck[2][1] = -sin(theta)*cos(phi);
        Matrixcheck[2][2] = cos(theta);

        double plane1[6][3] =
        {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
        double plane2[6][3] = {0};

        for(int i2 = 0; i2 < 6; i2++)
            for(int k = 0; k < 3; k++)
                for(int j = 0; j < 3; j++)
                    plane2[i2][k] += plane1[i2][j] * Matrixcheck[k][j];

        double angleMin = 90;

        for(int l = 0; l < 6; l++)
        {
            double dotproduct =
                (plane2[l][0]*htex +
                 plane2[l][1]*ktex +
                 plane2[l][2]*ltex) /
                (sqrt(htex*htex+ktex*ktex+ltex*ltex) *
                 sqrt(plane2[l][0]*plane2[l][0] +
                      plane2[l][1]*plane2[l][1] +
                      plane2[l][2]*plane2[l][2]));

            if(dotproduct < 0) dotproduct = -dotproduct;
            if(dotproduct > 1) dotproduct = 1;

            double angle = 180 * acos(dotproduct) / M_PI;
            if(angle < angleMin) angleMin = angle;
        }

        angleMin *= M_PI / 180;

        double Prob =
            f * (1 / sqrt(pow(
                (pow(r,2)*pow(cos(angleMin),2) +
                 (pow(sin(angleMin),2)/r)),3))) + (1 - f);

        int randnumber = rand() % maxno;

        if(double(randnumber) < Prob)
            i++;

        count++;

    } while(i < grainNo);
}

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
)
{
    double* gblistweight = new double[countGB];

    for(int i = 0; i < countGB; i++)
        gblistweight[i] = 1;

    for(int weighind = 0; weighind < gbweightno; weighind++)
    {
        int randsel = rand() % countGB;

        if(gblistweight[randsel] == 1)
        {
            int grainint1 = GBID1[randsel];
            int grainint2 = GBID2[randsel];

            double phi = EA1[grainint1];
            double psi = EA3[grainint1];
            double theta = EA2[grainint1];

            double Matrix1[3][3];
            double Matrix2[3][3];
            double Matrixmult[3][3];

            Matrix1[0][0] = cos(psi)*cos(phi)-cos(theta)*sin(phi)*sin(psi);
            Matrix1[0][1] = cos(psi)*sin(phi)+cos(theta)*cos(phi)*sin(psi);
            Matrix1[0][2] = sin(psi)*sin(theta);
            Matrix1[1][0] = -sin(psi)*cos(phi)-cos(theta)*sin(phi)*cos(psi);
            Matrix1[1][1] = -sin(psi)*sin(phi)+cos(theta)*cos(phi)*cos(psi);
            Matrix1[1][2] = cos(psi)*sin(theta);
            Matrix1[2][0] = sin(theta)*sin(phi);
            Matrix1[2][1] = -sin(theta)*cos(phi);
            Matrix1[2][2] = cos(theta);

            Matrix2[0][0] = cos(theta2)+(a1*a1*(1-cos(theta2)));
            Matrix2[0][1] = a1*a2*(1-cos(theta2))+a3*sin(theta2);
            Matrix2[0][2] = a1*a3*(1-cos(theta2))-a2*sin(theta2);
            Matrix2[1][0] = a2*a1*(1-cos(theta2))-a3*sin(theta2);
            Matrix2[1][1] = cos(theta2)+(a2*a2*(1-cos(theta2)));
            Matrix2[1][2] = a2*a3*(1-cos(theta2))+a1*sin(theta2);
            Matrix2[2][0] = a3*a1*(1-cos(theta2))+a2*sin(theta2);
            Matrix2[2][1] = a3*a2*(1-cos(theta2))-a1*sin(theta2);
            Matrix2[2][2] = cos(theta2)+a3*a3*(1-cos(theta2));

            for(int i = 0; i < 3; i++)
                for(int j = 0; j < 3; j++)
                {
                    Matrixmult[i][j] = 0;
                    for(int k = 0; k < 3; k++)
                        Matrixmult[i][j] += Matrix2[i][k] * Matrix1[k][j];
                }

            double euler[3];

            euler[1] = acos(Matrixmult[2][2]);

            if(sin(euler[1]) == 0)
            {
                euler[0] = acos(Matrixmult[0][0]);
                euler[2] = 0;
            }
            else
            {
                euler[0] = asin(Matrixmult[2][0]/sin(euler[1]));
                euler[2] = asin(Matrixmult[0][2]/sin(euler[1]));
            }

            EA1[grainint2] = euler[0]*180/M_PI;
            EA2[grainint2] = euler[1]*180/M_PI;
            EA3[grainint2] = euler[2]*180/M_PI;

            gblistweight[randsel] = 0;
        }
        else
            weighind--;
    }

    delete[] gblistweight;
}
