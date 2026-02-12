//
//  main.cpp
//  DomainPercolationFerroelectrics
//
//  Created by Sukriti Mantri on 2/10/26.
//  Original code created in April 2021
//  This is a refactored code to improve readability, reproducibility, and extensibility

#include <iostream>

#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <fstream>
#include <ctime>
#include <stdlib.h>     /* srand, rand */
#include <vector>
#include <algorithm>
#include <string>

#include "Microstructure.h"
#include "GrainBoundaries.h"
#include "EulerAngles.h"
#include "Percolation.h"
#include "IO.h"





#define pi 3.141592653589793
using namespace std;

int sign(double i)
{
    if (i>0)
        return 1;
    else if(i<0)
        return -1;
    else
        return 0;
}


//CODE BEGINS//
int main()
{
    
    
    cout << "are you doing random: 1 for yes and 0 for no\n";
    int randomweigh;
    cin >> randomweigh;
    string random = std::to_string(randomweigh);
    
    cout << "are you doing Sigma 3 weighting: 1 for yes and 0 for no\n";
    int s3;
    cin >> s3;
    int gbweightno = 25;
    double a1 = 1/sqrt(3);
    double a2 = 1/sqrt(3);
    double a3 = 1/sqrt(3);
    double theta2 = pi/3;
    string sigma = std::to_string(s3);
    sigma += "_perc";
    sigma += std::to_string(gbweightno*100/50);

   
    
    cout << "are you doing MD weighting: 1 for yes and 0 for no\n";
    int MD;
    cin >> MD;
    double r=0.1;  //this is the r index for MD function
    double f=1;
    int htex=1;
    int ktex=0;
    int ltex=0;
    int maxno = (f*(1/(pow(r,3)))+1-f)+1;
    string mdollase = std::to_string(MD); mdollase+= "_MDpara" ; mdollase+= std::to_string(r); mdollase+= "_MDperc"; mdollase += std::to_string(f);

    
    double chTol= 2.0;
    string charge =  std::to_string(chTol);
    double tolerance =5;
    string angletol =  std::to_string(tolerance);
    int steps =100;  //this is the grid for microstructure generation
    int grainNo =50;  //this is the number of grains , if you change the number, you might need to manually change it at few places, do find and replace in that case for "50"
    string grains =  std::to_string(grainNo);
    int Iterations=100;
    string totalmic =  std::to_string(Iterations);
    int Iterationseclevel = 100;
    string iterationsinmic =  std::to_string(Iterationseclevel);
    string date = "12April2021";
    string commonpart;
    commonpart = totalmic; commonpart+= "_R" ; commonpart+= random;commonpart+= "_S3" ; commonpart+= sigma; commonpart+= "_MD" ; commonpart+= mdollase; commonpart+= "_charge"; commonpart+= charge; commonpart+= "_angtol"; commonpart+= angletol; commonpart+= date;
    
   
    
   

    
   srand(int(time(0)));
   
    double *meanpercmax = new(std::nothrow) double[Iterations];
    double *meanperctot = new(std::nothrow) double[Iterations];
    double *meanpercmaxstddev = new(std::nothrow) double[Iterations];
    double *meanperctotstddev = new(std::nothrow) double[Iterations];
    
    
    for(int iterstart=0; iterstart<Iterations; iterstart+=1)
    {
        double *xCoor = new(std::nothrow) double[steps*steps*steps];
        double *yCoor = new(std::nothrow) double[steps*steps*steps];
        double *zCoor = new(std::nothrow) double[steps*steps*steps];

        double *grainCX = new(std::nothrow) double[grainNo];
        double *grainCY = new(std::nothrow) double[grainNo];
        double *grainCZ = new(std::nothrow) double[grainNo];

        double *belong = new(std::nothrow) double[steps*steps*steps];

        double *grainCenterCx = new(std::nothrow) double[grainNo];
        double *grainCenterCy = new(std::nothrow) double[grainNo];
        double *grainCenterCz = new(std::nothrow) double[grainNo];
        double *pointbelongtograin = new(std::nothrow) double[grainNo];

        double *grainSize = new(std::nothrow) double[grainNo];
        double *grainRad  = new(std::nothrow) double[grainNo];


    //this is to ensure that no two grain centers are very close
    
        generateCoordinates(steps, xCoor, yCoor, zCoor);

        generateGrainCenters(grainNo, steps, grainCX, grainCY, grainCZ);

        assignPointsToGrains(
            steps, grainNo,
            xCoor, yCoor, zCoor,
            grainCX, grainCY, grainCZ,
            belong
        );

        computeGrainCentersFromPoints(
            steps, grainNo,
            xCoor, yCoor, zCoor,
            belong,
            grainCenterCx, grainCenterCy, grainCenterCz,
            pointbelongtograin
        );

        computeGrainSizeAndRadius(
            steps, grainNo,
            belong,
            grainSize, grainRad
        );

        double *gbX     = new(std::nothrow) double[steps*steps*steps];
        double *gbY     = new(std::nothrow) double[steps*steps*steps];
        double *gbZ     = new(std::nothrow) double[steps*steps*steps];
        double *gb_GID1 = new(std::nothrow) double[steps*steps*steps];
        double *gb_GID2 = new(std::nothrow) double[steps*steps*steps];
        double *gb_ID   = new(std::nothrow) double[steps*steps*steps];

        double *A = new(std::nothrow) double[steps*steps*steps];
        double *B = new(std::nothrow) double[steps*steps*steps];
        double *C = new(std::nothrow) double[steps*steps*steps];

        double *meangbX = new(std::nothrow) double[steps*steps*steps];
        double *meangbY = new(std::nothrow) double[steps*steps*steps];
        double *meangbZ = new(std::nothrow) double[steps*steps*steps];

    
        int gbPOINTno = findGrainBoundaries(
            steps,
            belong,
            xCoor, yCoor, zCoor,
            gbX, gbY, gbZ,
            gb_GID1, gb_GID2, gb_ID
        );

        int countGB = fitGrainBoundaryPlanes(
            gbPOINTno,
            gbX, gbY, gbZ,
            gb_GID1, gb_GID2, gb_ID,
            A, B, C,
            meangbX, meangbY, meangbZ
        );
        printf ("the number of grain boundaries is %d\n",countGB);

        double *GBx = new(std::nothrow) double[countGB];
        double *GBy = new(std::nothrow) double[countGB];
        double *GBz = new(std::nothrow) double[countGB];
        double *GBID1 = new(std::nothrow) double[countGB];
        double *GBID2 = new(std::nothrow) double[countGB];

        countGB = computeGrainBoundaryNormals(
            gbPOINTno,
            A, B, C,
            gb_GID1, gb_GID2,
            GBx, GBy, GBz,
            GBID1, GBID2
        );

  
    int pocketID1, pocketID2, pocketID3;  //for preferential picking up of Euler angles
    pocketID1 = rand() % 91;
    pocketID2 = rand() % 91 ;
    pocketID3= rand() % 91 ;
        
    
    ////allocating euler angles
    double *EA1 = new(std::nothrow) double[grainNo];                      //Arrays for the three euler angle; this one is for phi1/phi2/phi
    double *EA2= new(std::nothrow) double[grainNo];                       //this one is for phi1/phi2/phi
    double *EA3 = new(std::nothrow) double[grainNo];                      //this one is for phi1/phi2/phi
    
    //keeping the first grain the global axes////
        if(randomweigh == 1)
            assignRandomEulerAngles(grainNo, EA1, EA2, EA3);

        if(MD == 1)
            assignMDEulerAngles(
                grainNo, r, f,
                htex, ktex, ltex,
                maxno,
                EA1, EA2, EA3
            );

        if(s3 == 1)
            applySigma3Weighting(
                gbweightno,
                countGB,
                a1, a2, a3,
                theta2,
                GBID1, GBID2,
                EA1, EA2, EA3
            );

    ////The microstructure code ends here/////
        double meanvaluenet = 0;
        double meanvaluetotalnet = 0;
        double standerddeviationmax = 0;
        double standerddeviationtot = 0;

        runPercolationForMicrostructure(
            grainNo,
            Iterationseclevel,
            countGB,
            chTol,
            tolerance,
            GBx,
            GBy,
            GBz,
            GBID1,
            GBID2,
            EA1,
            EA2,
            EA3,
            meanvaluenet,
            meanvaluetotalnet,
            standerddeviationmax,
            standerddeviationtot
        );

       

    
    ///////FINAL RESULT//////////////
        ///
        meanpercmax[iterstart] = meanvaluenet;
        meanperctot[iterstart] = meanvaluetotalnet;


        
        delete[] belong;
        delete[] xCoor;
        delete[] yCoor;
        delete[] zCoor;
        delete[] grainCX;
        delete[] grainCZ;
        delete[] grainCY;
        delete[] grainSize;
        delete[] grainRad;
        delete[] gbX;
        delete[] gbY;
        delete[] gbZ;
        delete[] gb_GID1;
        delete[] gb_GID2;
        delete[] gb_ID;
        delete[] A;
        delete[] B;
        delete[] GBID2;
        delete[] GBID1;
        delete[] GBx;
        delete[] GBy;
        delete[] GBz;
        delete[] C;
        delete[] meangbX;
        delete[] meangbY;
        delete[] meangbZ;
        delete[] EA1;
        delete[] EA2;
        delete[] EA3;
        delete[] grainCenterCx;
        delete[] grainCenterCy;
        delete[] grainCenterCz;
        delete[] pointbelongtograin;
    }
    //calculation of average values for entire microstructure
    
    
    
    
    //saving all the files
    writeBinaryFile(
        "PL_" + commonpart,
        meanpercmax,
        Iterations
    );


    writeBinaryFile(
        "PL_tot_" + commonpart,
        meanperctot,
        Iterations
    );




//freeeing all the space

delete[] meanpercmax;
delete[] meanperctot;


return 0;
}
