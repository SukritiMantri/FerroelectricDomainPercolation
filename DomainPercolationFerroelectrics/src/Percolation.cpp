//
//  Percolation.cpp
//  DomainPercolationFerroelectrics
//
//  Created by Sukriti Mantri on 2/10/26.
//

#include "Percolation.h"
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include "Constants.h"

using namespace std;

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
)
{
    double sumN = 0;
    double sumX = 0;
    double sumX2 = 0;
    double sumNtot = 0;
    double sumXtot = 0;
    double sumX2tot = 0;

 
    for(int iter2 =0; iter2 < Iterationseclevel; iter2+=1)
    {
    
double *grainPerclength = new(std::nothrow) double[grainNo]; //this stores how many grains has the thread percolated before it crosses this grain
double *favoriteneigh = new(std::nothrow) double[grainNo]; //this stores the branching factor for each grain
double *grainpercolationlengthMax = new(std::nothrow) double[grainNo];   //this is the max end to end percolatio possible for that nucleating grain
double *grainpercolationlengthTot = new(std::nothrow) double[grainNo]; //this is total web size for that nucleating point

double ThreadGrainArray[1000][50] = {NAN};
    ///////FINAL RESULT//////////////

    double meanvalue=0, meanvaluetot=0;
//double modearray[50]={0};

for(int graincirc=0; graincirc<50; graincirc+=1)
{
    int grainNuc;   //this will end up in random selection of nucleating grain
    grainNuc = graincirc;
    double domainWall1[50],domainWall2[50],domainWall3[50], pold11[50],pold12[50], pold13[50], pold21[50], pold22[50], pold23[50];
    double grainThreadID[50]; //this stores the thread that croesses a grain
    double grainCheck[50];  //this stores whether this grain allowed percolation=1 or not =0 or was not checked =2
    double grainlist[50]; // contains the list of grains that need to checked
    double threadPercolation[1000];  //percolation length of a particular thread with the index as the thread number, that is why 1000 as I dont know what it is
    double threadfinalGrain[1000]; // this saves the final grain that thread is at the moment;
    for (int var = 0; var < 1000; var+=1){threadPercolation[var] = NAN; threadfinalGrain[var]=NAN;}
            
    int threadID =0; //this is the inital thread ID for the nucleated grain

    for(int i=0; i<grainNo; i+=1)  //for all the grains
        {
            domainWall1[i] = NAN;
            domainWall2[i] = NAN;
            domainWall3[i] = NAN;
            pold11[i] = NAN;
            pold12[i] = NAN;
            pold13[i] = NAN;
            pold21[i] = NAN;
            pold22[i] = NAN;
            pold23[i] = NAN;
            
            grainCheck[i] = 2;
           // favoriteneigh[i]=NAN;  //assuming the intial branching factor to be NAN fro all grains
            grainThreadID[i] = -1;  //
            grainlist[i] = -1;   //assumed initial value for all the grains, the value -1 says that havent been checked yet.
         //   grainPerclength[i] = 0;
        }
        

    grainPerclength[grainNuc] = 0;//why is this zero ###3// also is percoaltion length a grain boundary data ###?? I realized you dont need bf at all if you associate percolation length to gbs
    
    //now do you know the misorientation fo this grain, we dont know which grain is it?
    double phi = EA1[grainNuc];          //Euler angles of first grain
    double psi = EA3[grainNuc];
    double theta = EA2[grainNuc];
    double Matrix[3][3];              //rot matrix for the first grain


    Matrix[0][0] = cos(psi)*cos(phi)-cos(theta)*sin(phi)*sin(psi);
    Matrix[0][1] = cos(psi)*sin(phi)+cos(theta)*cos(phi)*sin(psi);
    Matrix[0][2] = sin(psi)*sin(theta);
    Matrix[1][0] = -sin(psi)*cos(phi)-cos(theta)*sin(phi)*cos(psi);
    Matrix[1][1] = -sin(psi)*sin(phi)+cos(theta)*cos(phi)*cos(psi);
    Matrix[1][2] = cos(psi)*sin(theta);
    Matrix[2][0] = sin(theta)*sin(phi);
    Matrix[2][1] = -sin(theta)*cos(phi);
    Matrix[2][2] = cos(theta);
    
        //pick randomly any one domain wall
    
    
    int domNo = rand() % 6; //this gives the randomly selected domain wall for first grain

    double planeNormalsGrain1[6][3] = {{1,1,0},{1,-1,0}, {1, 0,1},{-1, 0,1},{0, 1, 1},{0, -1, 1}};
    double planeNormalsGrain2[1][3]  = {0,0,0};
    
//this gives you the transformed normals in grain 2.
    
    for(int k=0;k<3;k+=1)
        {
            for(int j=0;j<3;j+=1){
                planeNormalsGrain2[0][k]+=planeNormalsGrain1[domNo][j]*Matrix[k][j];}
        }
    
    domainWall1[grainNuc] = planeNormalsGrain2[0][0];       //domain wall for first grain
    domainWall2[grainNuc] = planeNormalsGrain2[0][1];
    domainWall3[grainNuc] = planeNormalsGrain2[0][2];
    
    //randomly picking up a pol vector combination as well
    //findign pol vectors in grain 2
    //randomly picking up which of the four combinations you want : 0=> positive and rightb order, 1=> negative and right order, 2=> positive and opp order, 3=> negative and opp order
    int combNo = rand() % 4;
    
    double polVectors1[6][3] = {{1,0,0}, {-1,0,0}, {0,1,0}, {0,-1,0}, {0,0,1}, {0,0,-1}};
    double polVectors2[6][3] = {0} ;
    
    for(int i=0;i<6;i+=1)
        {
            for(int k=0;k<3;k+=1)
                {
                    for(int j=0;j<3;j+=1)
                    {
                        polVectors2[i][k]+=polVectors1[i][j]*Matrix[k][j];}
                }
        }
   
    
    //all 6 domain walls correspond to pol vec 1 and pol vec 2
    //double charge1[6][2][3];
    double charge2[6][2][3];
    
    for( int i=0; i<3; i+=1)
        {
            charge2[0][0][i] = polVectors2[0][i]; charge2[0][1][i] = polVectors2[2][i];
            charge2[1][0][i] = polVectors2[0][i]; charge2[1][1][i] = polVectors2[3][i];
            charge2[2][0][i] = polVectors2[0][i]; charge2[2][1][i] = polVectors2[4][i];
            charge2[3][0][i] = polVectors2[0][i]; charge2[3][1][i] = polVectors2[5][i];
            charge2[4][0][i] = polVectors2[2][i]; charge2[4][1][i] = polVectors2[4][i];
            charge2[5][0][i] = polVectors2[2][i]; charge2[5][1][i] = polVectors2[5][i];
        
        }
    
    if( combNo==0)
        {
            pold11[grainNuc] =  charge2[domNo][0][0];
            pold12[grainNuc] =  charge2[domNo][0][1];
            pold13[grainNuc] = charge2[domNo][0][2];
            pold21[grainNuc] = charge2[domNo][1][0];
            pold22[grainNuc] = charge2[domNo][1][1];
            pold23[grainNuc] = charge2[domNo][1][2];
        }
    else if( combNo==1)
        {
            pold11[grainNuc] =  -charge2[domNo][0][0];
            pold12[grainNuc] =  -charge2[domNo][0][1];
            pold13[grainNuc] = -charge2[domNo][0][2];
            pold21[grainNuc] = -charge2[domNo][1][0];
            pold22[grainNuc] = -charge2[domNo][1][1];
            pold23[grainNuc] = -charge2[domNo][1][2];
        }
    else if( combNo==2)
        {
            pold11[grainNuc] =  charge2[domNo][1][0];
            pold12[grainNuc] =  charge2[domNo][1][1];
            pold13[grainNuc] = charge2[domNo][1][2];
            pold21[grainNuc] = charge2[domNo][0][0];
            pold22[grainNuc] = charge2[domNo][0][1];
            pold23[grainNuc] = charge2[domNo][0][2];
        }
    else
        {
            pold11[grainNuc] =  -charge2[domNo][1][0];
            pold12[grainNuc] =  -charge2[domNo][1][1];
            pold13[grainNuc] = -charge2[domNo][1][2];
            pold21[grainNuc] = -charge2[domNo][0][0];
            pold22[grainNuc] = -charge2[domNo][0][1];
            pold23[grainNuc] = -charge2[domNo][0][2];}
    
//we now have all the information for the nucleating grain

//now we will have to run through all the neighbours of the nucleating grain...
    int graintestID = 0;  //this is me saving the grains that I want to test
    grainlist[0] = grainNuc;  //first in the grain list is the nucleating grain that I want to check
    int countThread=0;

    grainThreadID[grainNuc] =0;  //the first grain i.e., the nucleating grain has zero thread ID
    
    for (int graintest = 0;  grainlist[graintest]!=-1 ; graintest= graintest+1) //looping over the grains I have saved
        {
            //cout << "the testing grain is = " << graintest << "\n";
            //cout << "the grain test ID at eh momemnt is = " << graintestID << "\n";
            //cout << "the grain that is being testes = " << grainlist[graintest] << "\n";

            int gr1ID;  //its first value will be the nucleating grain
            //  if(grainlist[graintest]==-1)// this means this grain has not been tested therfoer it cant be gr1ID
            //     break;
            gr1ID= grainlist[graintest] ;  //if GRAINTEST=0 this is nucleating grain
            grainCheck[gr1ID] = 1;   //why are you assuming that this will continue, shouldnt this be zero to start with
            int gr2ID;     //next grain's ID
 

            //Finding all the neighbours of the nucleating grain///
            int countNeigh=0; //count for the number of neighbours just for that grain
            int neighID[100];    //stores the grain boundary id for all those grain neighbours
            for(int i=0; i<countGB; i+=1)   //testibg all the grain boundaries to see if either side of the grain is gr1ID
            {
                neighID[i]=NAN;
                if(GBID1[i] == gr1ID)                  //i gives you the id of the grain boundary
                {
                    neighID[countNeigh]=i;   //saving the grain boundary ID
                    countNeigh = countNeigh+1;
                
                }
            
                else if(GBID2[i] == gr1ID)
                {
                    neighID[countNeigh]=i;
                    countNeigh = countNeigh+1;
                }
            
              //checking the next grain boundary
            
            }
    
            //cout << "the number of neighbours  = " << countNeigh <<  "\n";    //CHECK 1
        
            int favneigh = -1;   //this is your initial branching factor for this grain
            for(int selNeigh = 0; selNeigh<countNeigh; selNeigh+=1)   //covering all the neighbours of this grain
                {
          
                    int id=0;  //this is to check whihc side was the nucleating grain so that we can see whihc grain we need to check percolation for
        
               //the id of this grain is given by neighID[selNeigh]////
                    {
                        if(GBID1[neighID[selNeigh]]==gr1ID)     ////something wrong here
                            id=2;//something is wrond HERE#####
                        else if(GBID2[neighID[selNeigh]]==gr1ID)
                            id=1;
                        else
                            printf("soemhting is wrong here\n");   //CHECK
                    }
        
        
              //grID is the id of the selected neighbour for its euler angles
                    {
                        if(id==1)
                            gr2ID = GBID1[neighID[selNeigh]];
                        else if(id==2)
                            gr2ID = GBID2[neighID[selNeigh]];
                        else
                            printf("something is wrong here\n");   //CHECK
                    }
        
                    //inital assumptions for energy
                    double outputCharge = 1000;    //default value
                    double outputAngle = 90;   //default value
                    double output = 0;    //default value
                    double grainboundaryNormal[3] ;   //grainboundary space variable
                    double normGrainboundarynormal;
        
                    //Finding the gb normal between our first two grains
                    grainboundaryNormal[0] = GBx[neighID[selNeigh]]; //////this time you dont need ANGLES
                    grainboundaryNormal[1] = GBy[neighID[selNeigh]];
                    grainboundaryNormal[2] = GBz[neighID[selNeigh]];
                    normGrainboundarynormal = sqrt(pow(grainboundaryNormal[0],2)+pow(grainboundaryNormal[1],2)+ pow(grainboundaryNormal[2],2));
                    grainboundaryNormal[0] =  grainboundaryNormal[0]/normGrainboundarynormal;
                    grainboundaryNormal[1] = grainboundaryNormal[1]/normGrainboundarynormal;
                    grainboundaryNormal[2] = grainboundaryNormal[2]/normGrainboundarynormal;
       
            
                    if(isnan(domainWall1[gr2ID]))    //checking if that grain already has a domain wall or not
                    {
            
                        phi = EA1[gr2ID];    //FOR THE GRSAIN WITH GRAIN  ID next to ith grain as we already know the resukt for second grain with i=1
                        psi = EA3[gr2ID];
                        theta = EA2[gr2ID];
            
                        double Matrix[3][3];
            
                        Matrix[0][0] = cos(psi)*cos(phi)-cos(theta)*sin(phi)*sin(psi);
                        Matrix[0][1] = cos(psi)*sin(phi)+cos(theta)*cos(phi)*sin(psi);
                        Matrix[0][2] = sin(psi)*sin(theta);
                        Matrix[1][0] = -sin(psi)*cos(phi)-cos(theta)*sin(phi)*cos(psi);
                        Matrix[1][1] = -sin(psi)*sin(phi)+cos(theta)*cos(phi)*cos(psi);
                        Matrix[1][2] = cos(psi)*sin(theta);
                        Matrix[2][0] = sin(theta)*sin(phi);
                        Matrix[2][1] = -sin(theta)*cos(phi);
                        Matrix[2][2] = cos(theta);
            
                        //Set all the (110) plane normals for [100] oriented grain and calculate plane normals in grain 2
                        // creating my own function for matrix multiplication
                        //finding domain walls in grain 2//
                        double planeNormalsGrain1[6][3] = {{1,1,0},{1,-1,0}, {1, 0,1},{-1, 0,1},{0, 1, 1},{0, -1, 1}};
                        double planeNormalsGrain2g[6][3]  = {0};
            
                        //this gives you the transformed normals in grain 2.
                        for(int i=0;i<6;i+=1)
                        {
                            for(int k=0;k<3;k+=1)
                            {
                                for(int j=0;j<3;j+=1)
                                planeNormalsGrain2g[i][k]+=planeNormalsGrain1[i][j]*Matrix[k][j];
                            }
                        }
    
                        //pol vecs of the otehr grain//
                        //findign pol vectors in grain 2
                        double polVectors1[6][3] = {{1,0,0}, {-1,0,0}, {0,1,0}, {0,-1,0}, {0,0,1}, {0,0,-1}};
                        double polVectors2g[6][3] = {0} ;
            
                        for(int i=0;i<6;i+=1)
                        {
                            for(int k=0;k<3;k+=1)
                            {
                                for(int j=0;j<3;j+=1)
                                polVectors2g[i][k]+=polVectors1[i][j]*Matrix[k][j];
                            }
                        }
            

                        double charge2g[6][2][3];
            
                        for( int i=0; i<3; i+=1)
                        {
                
                
                            charge2g[0][0][i] = polVectors2g[0][i]; charge2g[0][1][i] = polVectors2g[2][i];
                            charge2g[1][0][i] = polVectors2g[0][i]; charge2g[1][1][i] = polVectors2g[3][i];
                            charge2g[2][0][i] = polVectors2g[0][i]; charge2g[2][1][i] = polVectors2g[4][i];
                            charge2g[3][0][i] = polVectors2g[0][i]; charge2g[3][1][i] = polVectors2g[5][i];
                            charge2g[4][0][i] = polVectors2g[2][i]; charge2g[4][1][i] = polVectors2g[4][i];
                            charge2g[5][0][i] = polVectors2g[2][i]; charge2g[5][1][i] = polVectors2g[5][i];
                
                        }
                        //NOW WE HAVE ALL THE INFORMATION ABOUT THE SECOND GRAIN AS WELL
            
            
                   
            
            
                    ////////////////////////////////////////angle calc./////////////////////////////
            
                        double lineInt1[3];  //line of intersection of domain wall plane in grain A and grain boundary plane
                        double lineIntNorm1;
                        double lineInt2[3]; //line of intersection of domain wall plane in grain B and grain boundary plane, rem there are 6 possibilities for that
                        double lineIntNorm2;
                        double Angle[6]={90};
            
                        double charge[6]={1000};
                        int DWident=0;   //this represents the finally selected domain wall to continue
                        int orderparam=0;    //this represents which pol vecto combination type ie, +ve/-ve or ordered /reverse
                        for(int l=0; l<6; l+=1)
                        {
                            double dotproduct=0;
                            lineInt1[0]= domainWall2[gr1ID]*grainboundaryNormal[2]-grainboundaryNormal[1]*domainWall3[gr1ID];
                            lineInt1[1]=domainWall3[gr1ID]*grainboundaryNormal[0]-domainWall1[gr1ID]*grainboundaryNormal[2];
                            lineInt1[2]=domainWall1[gr1ID]*grainboundaryNormal[1]-domainWall2[gr1ID]*grainboundaryNormal[0];
                            lineIntNorm1 = sqrt((pow(lineInt1[0],2))+(pow(lineInt1[1],2))+(pow(lineInt1[2],2)));
                            //calculting line of intersection 2 for all 6 possible domain walls in grain 1
                            lineInt2[0]=planeNormalsGrain2g[l][1]*grainboundaryNormal[2]-grainboundaryNormal[1]*planeNormalsGrain2g[l][2];
                            lineInt2[1]=planeNormalsGrain2g[l][2]*grainboundaryNormal   [0]-planeNormalsGrain2g[l][0]*grainboundaryNormal[2];
                            lineInt2[2]=planeNormalsGrain2g[l][0]*grainboundaryNormal[1]-planeNormalsGrain2g[l][1]*grainboundaryNormal[0];
                            lineIntNorm2 = sqrt((pow(lineInt2[0],2))+(pow(lineInt2[1],2))+(pow(lineInt2[2],2)));
                    
                            //calculting the angle between the two line of intersections
                    
                            dotproduct= ((lineInt1[0]*lineInt2[0]) + (lineInt1[1]*lineInt2[1]) + (lineInt1[2]*lineInt2[2]))/(lineIntNorm1*lineIntNorm2);
                            //calculting the angle between the two line of intersections
                    
                    
                            if(dotproduct <0)
                                dotproduct = - dotproduct;// we need the modulus for comparison
                    
                            if(dotproduct <=1.00005 && dotproduct >=1)
                                dotproduct=1;
                    
                            dotproduct = 180*acos(dotproduct)/pi ;
                
                    
                            Angle[l] = dotproduct;   //saving the angle for each domain wall of grain 2
                    
                    
                    
                            if(Angle[l] < tolerance)
                        
                            {
                                //cout << "angle less thna tolerance" << "\n";
                                output = output+1;
                                double Mincharge;
                                double polG1D1[3], polG2D1[3], polG1D2[3], polG2D2[3];
                        
                                double dotres[4];
                                double chargeres[4];
                        
                        
                                polG1D1[0] = pold11[gr1ID];
                                polG1D1[1] =pold12[gr1ID];
                                polG1D1[2] =pold13[gr1ID];
                                polG1D2[0] = pold21[gr1ID];
                                polG1D2[1] =pold22[gr1ID];
                                polG1D2[2] =pold23[gr1ID];
                        
                                for( int loop=0; loop<3; loop+=1)
                                {
                                    polG2D1[loop] = charge2g[l][0][loop]; polG2D2[loop] = charge2g[l][1][loop];
                            
                                }
                                double normp1 = sqrt(polG1D1[0]*polG1D1[0]+polG1D1[1]*polG1D1[1]+polG1D1[2]*polG1D1[2]);
                                double normp2 = sqrt(polG1D2[0]*polG1D2[0]+polG1D2[1]*polG1D2[1]+polG1D2[2]*polG1D2[2]);
                                double normp3 = sqrt(polG2D1[0]*polG2D1[0]+polG2D1[1]*polG2D1[1]+polG2D1[2]*polG2D1[2]);
                                double normp4 = sqrt(polG2D2[0]*polG2D2[0]+polG2D2[1]*polG2D2[1]+polG2D2[2]*polG2D2[2]);
                        
                                dotres[0] = (polG1D1[0]*grainboundaryNormal[0] + polG1D1[1]*grainboundaryNormal[1] + polG1D1[2]*grainboundaryNormal[2])/normp1;
                                dotres[1] = (polG1D2[0]*grainboundaryNormal[0] + polG1D2[1]*grainboundaryNormal[1] + polG1D2[2]*grainboundaryNormal[2])/normp2;
                                dotres[2] = (polG2D1[0]*grainboundaryNormal[0] + polG2D1[1]*grainboundaryNormal[1] + polG2D1[2]*grainboundaryNormal[2])/normp3;
                                dotres[3] = (polG2D2[0]*grainboundaryNormal[0] + polG2D2[1]*grainboundaryNormal[1] + polG2D2[2]*grainboundaryNormal[2])/normp4;
                        
                                chargeres[0] = fabs(dotres[0]-dotres[2])+fabs(dotres[1]-dotres[3]);  //order +ve
                                chargeres[1] = fabs(dotres[0]+dotres[2])+fabs(dotres[1]+dotres[3]);  //order -ve
                                chargeres[2] = fabs(dotres[0]-dotres[3])+fabs(dotres[1]-dotres[2]);  //reverse order +ve
                                chargeres[3] = fabs(dotres[0]+dotres[3])+fabs(dotres[1]+dotres[2]); //reverse order -ve
                        
                        
                                int countorder;
                                Mincharge = 1000;
                                for( int loop=0; loop<4; loop+=1)
                                {
                                    if(Mincharge > chargeres[loop])
                                    {
                                        Mincharge = chargeres[loop];
                                        countorder = loop;
                                    }
                                }
                        
                                charge[l] = Mincharge;   //saving the charge for that domain wall
                        
                                if(outputCharge> charge[l])
                                {
                                    outputAngle = Angle[l];
                                    outputCharge=charge[l];
                                    DWident = l;
                                    orderparam = countorder;
                           // outputStrain[m]=strain[k][l];
                                }
                            }//if ends
                    //KEEPING THE DOMAIN WALL KNOWN//
                        }
                            domainWall1[gr2ID] = planeNormalsGrain2g[DWident][0];
                            domainWall2[gr2ID] = planeNormalsGrain2g[DWident][1];
                            domainWall3[gr2ID] = planeNormalsGrain2g[DWident][2];
                            outputAngle = Angle[DWident];
                    
                            if(orderparam==0)
                            {
                        
                                pold11[gr2ID] = charge2g[DWident][0][0];
                                pold12[gr2ID] = charge2g[DWident][0][1];
                                pold13[gr2ID] = charge2g[DWident][0][2];
                                pold21[gr2ID] = charge2g[DWident][1][0];
                                pold22[gr2ID] = charge2g[DWident][1][1];
                                pold23[gr2ID] = charge2g[DWident][1][2];
                            }
                            else if( orderparam==1)
                            {
                        
                                pold11[gr2ID] = -charge2g[DWident][0][0];
                                pold12[gr2ID] = -charge2g[DWident][0][1];
                                pold13[gr2ID] = -charge2g[DWident][0][2];
                                pold21[gr2ID] = -charge2g[DWident][1][0];
                                pold22[gr2ID] = -charge2g[DWident][1][1];
                                pold23[gr2ID] = -charge2g[DWident][1][2];
                            }
                            else if(orderparam==2)
                            {
                        
                                pold11[gr2ID] = charge2g[DWident][1][0];
                                pold12[gr2ID] = charge2g[DWident][1][1];
                                pold13[gr2ID] = charge2g[DWident][1][2];
                                pold21[gr2ID] = charge2g[DWident][0][0];
                                pold22[gr2ID] = charge2g[DWident][0][1];
                                pold23[gr2ID] = charge2g[DWident][0][2];
                            }
                            else
                            {
                       
                                pold11[gr2ID] = -charge2g[DWident][1][0];
                                pold12[gr2ID] = -charge2g[DWident][1][1];
                                pold13[gr2ID] = -charge2g[DWident][1][2];
                                pold21[gr2ID] = -charge2g[DWident][0][0];
                                pold22[gr2ID] = -charge2g[DWident][0][1];
                                pold23[gr2ID] = -charge2g[DWident][0][2];
                            }
                        }
                    
            
                    if((output > 0) &&(outputAngle<tolerance) && (outputCharge<chTol)) //this grain is good to go
                    {
                    
                        if(favneigh==-1){favneigh=0;} //this is more to do with the nucleatinf grain
                        else{favneigh = favneigh+1; } //the first favneigth =0, its default was -1
                        //cout << "the brancing factor has increased to = " << favneigh << "\n";
                        graintestID = graintestID+1;   //increasing the number of grains we want to test
                        grainCheck[gr2ID] =1;    //making this grain the one where percolation was posible and has been tested
                        grainlist[graintestID] = gr2ID;    //this is adding thr grain ID of the next grain to be tested
                        //cout << "the grain that has enetered the list = " << grainlist[graintestID] << "\n";
                        grainPerclength[gr2ID] = grainPerclength[gr1ID]+1;  // the percolation length of the second grain is first grain +1
                        // threadID = grainThreadID[gr1ID];//WHY IS THIS LIKE THIS
                        if(favneigh==0)   //if this is the first neighbour beign tested than same thread will continue
                            {grainThreadID[gr2ID] = grainThreadID[gr1ID]; threadID = grainThreadID[gr2ID];}
                        else  //if new thread starts
                            {grainThreadID[gr2ID] = countThread+favneigh; threadID = grainThreadID[gr2ID];}  //THIS IS WHERE THE PROBLEM IS 1 probs here
                        threadfinalGrain[threadID] = gr2ID;
                        threadPercolation[threadID] = grainPerclength[gr2ID];
                        //what is it is emamanting from a grain //
                        ThreadGrainArray[threadID][int(threadPercolation[threadID])] = gr2ID;
                   
                        int indmax = grainPerclength[gr2ID];
                        for(int ind2=0; ind2<indmax; ind2=ind2+1)
                        {
                            ThreadGrainArray[threadID][ind2] = ThreadGrainArray[int(grainThreadID[gr1ID])][ind2];
                        }
                        for(int ind2=0; ind2<=indmax; ind2=ind2+1)
                        {
                            //cout << "["<< threadID << "]"<<"["<< ind2 << "]"<< "=" << ThreadGrainArray[threadID][ind2] << "\n";;
                        }

                        //cout << "the threadGrainArray for "<< threadID << " and percolation length " << threadPercolation[threadID] << "is " << gr2ID << "\n";
                        //cout << "the thread's previous name for the grain" << gr1ID << "is = " << grainThreadID[gr1ID] <<"\n";
                        //cout << "the thread's new name for the grain" << gr2ID<<" is = " << threadID << "\n";
                    }
                    else // if the percolation is not happening
                        grainCheck[gr2ID] =0;
            
        
                }  //move to the next neighbour
            //this is the branching factor for that grain

                if(favneigh !=-1) //that is if some conitnuation atleast   happened
                {
                    favoriteneigh[gr1ID] = (favneigh+1);
                    countThread = countThread + favneigh;
                }
       
        
        }  //you have covered all the neighbours and now move to the next neighbour to study it for other of its own neighbours

    //cout << "the total grains tested = " << (graintestID) << "\n";
 
    //cout << "the count thread is " << countThread << "\n";
    double maxpercolation = 0;
    double roadtomaxper=0;
    double maxpercthread[2];
    //finding the maximum percolation length
    for(int i=0; i<countThread; i+=1)
    {
        for (int j=i+1; j<=countThread; j+=1)
        {
            roadtomaxper = threadPercolation[j] + threadPercolation[i]; //added the total thread lenghts initially
            for(int thread1=1;thread1<=threadPercolation[i];thread1+=1)
            {
                for(int thread2=1;thread2<=threadPercolation[j];thread2+=1)
                {
                    if(ThreadGrainArray[i][thread1]==ThreadGrainArray[j][thread2]){roadtomaxper = roadtomaxper - 2;} //subtacting one by one if there is a common member in the two threads thisseems like ERROR,shouldT I be subtarcting it two times and nto one
                }
            }
            if(roadtomaxper>maxpercolation){maxpercolation = roadtomaxper; maxpercthread[0]=i; maxpercthread[1]=j;}
        }
    
    }

    //cout << "maximum percolation lengthn is " << maxpercolation <<  "and the two threads are (" << maxpercthread[0]<<","<< maxpercthread[1]<< ")\n";
    //may be grainpercolation length would be the best judge

    grainpercolationlengthMax[graincirc] = maxpercolation;
    grainpercolationlengthTot[graincirc] = graintestID;
    meanvalue = meanvalue + grainpercolationlengthMax[graincirc];
    meanvaluetot = meanvaluetot + grainpercolationlengthTot[graincirc];
    
}
    meanvalue = round(meanvalue/grainNo); meanvaluetot = round(meanvaluetot/grainNo);
        double standdevnorm = 0;
        double standdevtot = 0;
        for(int i=0; i<grainNo; i+=1)
        {
            standdevnorm = standdevnorm +  pow((grainpercolationlengthMax[i] - meanvalue),2);
            standdevtot = standdevtot + pow((grainpercolationlengthTot[i] - meanvaluetot),2);
        }
        standdevnorm = sqrt(standdevnorm/grainNo);
        standdevtot = sqrt(standdevtot/grainNo);
            sumN = sumN + grainNo;
            sumX = sumX + meanvalue*grainNo;
            sumX2 = sumX2+ pow(standdevnorm,2)*(grainNo-1)+ (pow((meanvalue*grainNo),2)/grainNo);
            
            
            sumNtot = sumNtot + grainNo;
            sumXtot = sumXtot + meanvaluetot*grainNo;
            sumX2tot = sumX2tot+ pow(standdevtot,2)*(grainNo-1)+ (pow((meanvaluetot*grainNo),2)/grainNo);
            
            
            
            delete[] grainpercolationlengthMax;
            delete[] grainpercolationlengthTot;
            delete[] grainPerclength;
            delete[] favoriteneigh;
        }   //iterations ove rthe same microstructure ends
    
    
    meanvaluenet = sumX / sumN;
    standerddeviationmax = sqrt((sumX2-pow(sumX,2)/sumN)/(sumN-1)); //Combine SD = sqrt((txx-tx2/tn) / (tn-1))
    meanvaluetotalnet = sumXtot/sumNtot;
    standerddeviationtot = sqrt((sumX2tot-pow(sumXtot,2)/sumNtot)/(sumNtot-1));
  
}
