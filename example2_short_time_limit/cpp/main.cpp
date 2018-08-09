//
//  main.cpp
//  diffusion_myelin_exchange
//
//  Created by magda on 2/24/17.
//  Copyright (c) 2017 Honghsi. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <iomanip>
#include <time.h>
#include <cstdlib>
#include <algorithm>
#include <string>
#include <complex>
#include <string>
#include <vector>
#include "RNG.h"
#include "diffusion_lib.h"
    
    using namespace std;
    
#define Pi 3.14159265
    
    //********** Define tissue parameters **********
#define dt 3.0e5 //7.5e5        // Time step in ps
    const int TN=5e2;//16e5;    // Number of time steps
    const double TT=TN*dt;      // Total time
    const long NPar=5*1000;         // Number of particles
    const int timepoints=5e2;//1e3; // Number of time points to record
    const int Tstep=TN/timepoints;
    
#define Din 2.0e-9          // Diffusion coefficient inside the cylinders in m^2/s
#define Dout 2.0e-9         // Diffusion coefficient outside the cylinders in m^2/s
#define Dmc 0.5e-9          // Diffusion coefficient in myelin, circumferential, in m^2/s
#define kappa 0.0           // Permeability of a lipid bi-layer in um/ms
#define n 5.0             // Number density of lipid bi-layer in 1/um
    const double Dmr=0.2e-9;//kappa/n*1e-9;        // Diffusivity in myelin, radial, in m^2/s
    const double stepIN=sqrt(4*dt*Din);     // Step size in IAS in micrometers
    const double stepOUT=sqrt(4*dt*Dout);   // Step size in EAS in micrometers
    const double stepMC=sqrt(4*dt*Dmc);     // Step size in myelin in um, circumferential
    const double stepMR=sqrt(4*dt*Dmr);     // Step size in myelin in um, radial

    int main(int argc, char *argv[]) {
        
        clock_t begin=clock();
        clock_t end=clock();
        
        // Define index number
        int i=0, j=0, k=0;
        
        // Initialize RNG
        RNGvalues v;
        v = init_KISS();
        
        //********** Load mictostructure **********
        unsigned int NPix=0, NAx=0, **APix;
        double res=0.0, *xCir, *yCir, *rCir;
        
        // resolution
        ifstream myfile1 ("/Users/magda/Desktop/diffusion_simulation/cpp_practice/diffusion_myelin_exchange_ellipse_propagator_speed_up_v5/diffusion_myelin_exchange/phantom_res.txt", ios::in);
        myfile1>>res;
        myfile1.close();
        
        // Pixel # along each side
        ifstream myfile2 ("/Users/magda/Desktop/diffusion_simulation/cpp_practice/diffusion_myelin_exchange_ellipse_propagator_speed_up_v5/diffusion_myelin_exchange/phantom_NPix.txt", ios::in);
        myfile2>>NPix;
        myfile2.close();
        
        // Pixelized matrix A indicating axon labels
        APix=new unsigned int *[NPix];
        for (i=0; i<NPix; i++){
            APix[i]=new unsigned int[NPix];
            for (j = 0;j<NPix;j++){
                APix[i][j]=0;
            }
        }
        ifstream myfile3 ("/Users/magda/Desktop/diffusion_simulation/cpp_practice/diffusion_myelin_exchange_ellipse_propagator_speed_up_v5/diffusion_myelin_exchange/phantom_APix.txt", ios::in);
        for (i=0; i<NPix; i++){
            for (j=0; j<NPix; j++){
                myfile3>>APix[i][j];
            }
        }
        myfile3.close();
        
        // Number of axons
        ifstream myfile4 ("/Users/magda/Desktop/diffusion_simulation/cpp_practice/diffusion_myelin_exchange_ellipse_propagator_speed_up_v5/diffusion_myelin_exchange/phantom_NAx.txt", ios::in);
        myfile4>>NAx;
        myfile4.close();
        
        // Circle center of x-coordinate
        ifstream myfile5 ("/Users/magda/Desktop/diffusion_simulation/cpp_practice/diffusion_myelin_exchange_ellipse_propagator_speed_up_v5/diffusion_myelin_exchange/phantom_xCir.txt", ios::in);
        i=0; xCir=new double[NAx];
        while (!myfile5.eof()&(i<NAx)){
            myfile5>>xCir[i];
            i++;
        }
        myfile5.close();
        
        // Circle center of y-coordinate
        ifstream myfile6 ("/Users/magda/Desktop/diffusion_simulation/cpp_practice/diffusion_myelin_exchange_ellipse_propagator_speed_up_v5/diffusion_myelin_exchange/phantom_yCir.txt", ios::in);
        i=0; yCir=new double[NAx];
        while (!myfile6.eof()&(i<NAx)){
            myfile6>>yCir[i];
            i++;
        }
        myfile6.close();
        
        // Circle radius
        ifstream myfile7 ("/Users/magda/Desktop/diffusion_simulation/cpp_practice/diffusion_myelin_exchange_ellipse_propagator_speed_up_v5/diffusion_myelin_exchange/phantom_rCir.txt", ios::in);
        i=0; rCir= new double[NAx];
        while (!myfile7.eof()&(i<NAx)){
            myfile7>>rCir[i];
            i++;
        }
        myfile7.close();
        
        // The smallest number, which is > NAx, in the base 10
        unsigned int Nmax=0;
        ifstream myfile8 ("/Users/magda/Desktop/diffusion_simulation/cpp_practice/diffusion_myelin_exchange_ellipse_propagator_speed_up_v5/diffusion_myelin_exchange/phantom_Nmax.txt", ios::in);
        myfile8>>Nmax;
        myfile8.close();
        
        // The smallest number, which is > NAx, in the base 10
        double gratio=0;
        ifstream myfile9 ("/Users/magda/Desktop/diffusion_simulation/cpp_practice/diffusion_myelin_exchange_ellipse_propagator_speed_up_v5/diffusion_myelin_exchange/phantom_gratio.txt", ios::in);
        myfile9>>gratio;
        myfile9.close();
        
        /*
        vector<double> x(2,0), xc(2,0), xt(2,0), vi(2,0);
        double rc=0, g=0.75, dx=0.25;
        xc[0]=0;    xc[1]=0; rc=0.05;
        x[0]=0.1/sqrt(2);   x[1]=0.1/sqrt(2);
        //vi[0]=-1/sqrt(2);    vi[1]=-1/sqrt(2);
        //xt=elasticCollisionEAS(x, vi, dx, xc, rc);
        //cout<<xt[0]<<'\t'<<xt[1]<<endl;
        //cout<<0.25/sqrt(2)<<endl;
        
        double vC=1, phi=Pi/8, dxC=phi*vecDistance(x, xc);
        double vR=-1, dxR=0.3;
        xt=diffuseRadial(x, vR, dxR, xc);
        cout<<xt[0]<<'\t'<<xt[1]<<endl;
        cout<<-0.2/sqrt(2)<<endl;
         
        //cout<<0.1*cos(Pi*3/8)<<'\t'<<0.1*sin(Pi*3/8)<<endl;
        */
        
        //cout<<stepOUT<<'\t'<<res<<'\t'<<stepOUT/res<<endl;
        //cout<<NPix<<endl;
        
        //begin=clock();
        
        //********** Initialize Particle Positions in EAS *********
        const double probOUT=Pi/4*res*kappa/Dout*1e-9;// Probability constant for a particle from EAS to myelin
        const double probIN=Pi/4*res*kappa/Din*1e-9; // Probability constant for a particle from IAS to myelin
        const double probMY=Pi/4*res*kappa/Dmr*1e-9; // Probability constant for a particle from myelin to IAS/EAS
        double xiPar[NPar], yiPar[NPar];        // Initial positions of particles
        double xRand=0, yRand=0;                // Random numbers
        unsigned int a=0;                       // Element of APix matrix
        int ax1=0, ax2=0;                       // Labels of axons close to the particle
        vector<double> xyPar(2,0), xyCir(2,0);  // Particle position and Circle center
        vector<int> xyParG(2,0);                // Particle position on a grid
        bool instruction1=false, instruction2=false;     // Flags for a particle inside an axon
        k=0;
        while(k<NPar){
            //xiPar[k]=0.5; yiPar[k]=0.5;
            //k++;
            
            
            v=JKISS(v);
            xRand=v.r;
            v=JKISS(v);
            yRand=v.r;
            xyPar[0]=xRand; xyPar[1]=yRand;
            
            // Determine axons around the particle
            xyParG=pixPosition(xyPar,NPix);
            a=APix[xyParG[0]][xyParG[1]];
            ax1=a%Nmax; ax2=a/Nmax;
            
            // If the particle is in IAS, take the position
            instruction1=false; instruction2=false;
            if ( ax1 ){
                xyCir[0]=xCir[ax1-1]; xyCir[1]=yCir[ax1-1];
                instruction1 = inIAS(xyPar, xyCir, rCir[ax1-1],gratio);
            }
            
            if ( ax2 ){
                xyCir[0]=xCir[ax2-1]; xyCir[1]=yCir[ax2-1];
                instruction2 = inIAS(xyPar, xyCir, rCir[ax2-1],gratio);
            }
            
            if ( instruction1==true | instruction2==true ){
                for (i=0; i<1000; i++) {
                    xiPar[k+i]=xRand; yiPar[k+i]=yRand;
                    
                }
                k=k+1000;
                //xiPar[k+i]=xRand; yiPar[k+i]=yRand
                //k++;
            }
            
            
            
        }
        
        // ********** Simulate diffusion **********
        
        //double MOM2[timepoints]={0}, MOM4[timepoints]={0};  // Phase induced by diffusion
        //double myelinTime[timepoints][NPar]={0};            // Time spent in myelin
        int xjmp[NPar]={0}, yjmp[NPar]={0};                 // Jump numbers across boundaries
        //double dX[timepoints]={0}, dY[timepoints]={0};      // Displacement
        vector<double> xyi(2,0), xyt(2,0), xyTmp(2,0);      // Position before/after diffusion
        vector<double> xyCollision(2,0);                    // Position after collision
        vector<int> xytG(2,0), xyTmpG(2,0);                 // Position on grid after diffusion
        xyCir[0]=0; xyCir[1]=0;                             // Position of the axon center
        a=0; unsigned int aTmp=0;                           // Element of APix
        int ax[4]={0}; bool instruction[4]={0};                  // Axon label
        double phi=0, probPerm=0;                           // Diffusion direction, probability to permeate the myelin sheath
        vector<double> vp(2,0);                             // Nomalized diffusion velocity
        double axDis[4]={0}, axDisMin=0;                    // Distance to axon surface, its min
        int ax_hit=0;                                       // The axon label hit by the particle
        double probPermMY=stepMR/res*probMY;    // Probability of permeating from myelin to EAS or IAS
        double probPermIN=stepIN/res*probIN;    // Probability of permeating from IAS to Myelin
        double probPermOUT=stepOUT/res*probOUT; // Probability of permeating from EAS to Myelin
        //cout<<probPermMY<<endl<<probPermIN<<endl<<probPermOUT<<endl;
        
        const int Nbval=100;
        double bval[Nbval]={0}, qt=0;                              // b-value (ms/um^2), q=\gamma * g * \delta (1/um)
        double dx=0, dy=0, xPhase=0, yPhase=0;
        
        for (i=0; i<Nbval; i++) {
            bval[i]=0.1*(i+1);
        }
        //double (x/y)t_total[timepoints][NPar]={0};
        
        
        double **xt_total, **yt_total;
        xt_total=new double *[timepoints];
        yt_total=new double *[timepoints];
        for (i=0; i<timepoints; i++){
            xt_total[i]=new double[NPar];
            yt_total[i]=new double[NPar];
            for (j = 0;j<NPar;j++){
                xt_total[i][j]=0;
                yt_total[i][j]=0;
            }
        }
        
        double dxSquare[timepoints]={0}, dySquare[timepoints]={0}, dxQuartic[timepoints]={0}, dyQuartic[timepoints]={0}, gxSignalRe[timepoints][Nbval]={0}, gxSignalIm[timepoints][Nbval]={0}, gySignalRe[timepoints][Nbval]={0}, gySignalIm[timepoints][Nbval]={0};
        
        begin=clock();
        
        for (k=0; k<NPar; k++){         // Loop over initial positions
            
            xyi[0]=xiPar[k]; xyi[1]=yiPar[k];               // Position before diffusion
            xyt=xyi;                                        // Position after diffusion
            xytG=pixPosition(xyt, NPix);                    // Position on grid after diffusion
            
            //cout<<"Number"<<k<<endl;
            
            for (i=0; i<TN; i++){       // Loop over time points
                
                // The axons close to the previous step
                a=APix[xytG[0]][xytG[1]];
                ax1=a%Nmax, ax2=a/Nmax;
                
                // Check if the particle is in EAS
                instruction1=false; instruction2=false;
                if ( ax1 ) {
                    xyCir[0]=xCir[ax1-1]; xyCir[1]=yCir[ax1-1];
                    instruction1=inAxon(xyt, xyCir, rCir[ax1-1]);
                }
                
                if ( ax2 ) {
                    xyCir[0]=xCir[ax2-1]; xyCir[1]=yCir[ax2-1];
                    instruction2=inAxon(xyt, xyCir, rCir[ax2-1]);
                }
                
                if ( instruction1==false & instruction2==false ) { // In EAS
                    //cout<<"EAS"<<endl;
                    ax[0]=ax1; ax[1]=ax2;
                    
                    //Diffusion in EAS
                    
                    // Initialize velocity
                    v=JKISS(v);
                    phi=2*Pi*v.r;                               // Diffusion direction
                    vp[0]=cos(phi); vp[1]=sin(phi);             // Normalized diffusion velocity
                    
                    // Primitive position after diffusion
                    xyTmp[0]=xyt[0]+stepOUT/res*vp[0];
                    xyTmp[1]=xyt[1]+stepOUT/res*vp[1];
                    xyTmpG=pixPosition(xyTmp, NPix);
                    aTmp=APix[xyTmpG[0]][xyTmpG[1]];
                    ax[2]=aTmp%Nmax, ax[3]=aTmp/Nmax;
                    
                    // Determine if the segment(xyt,xyTmp) overlaps with the axon
                    for (j=0; j<4; j++) {
                        instruction[j]=false;
                        if ( ax[j] ) {
                            xyCir[0]=xCir[ax[j]-1];
                            xyCir[1]=yCir[ax[j]-1];
                            instruction[j]=stepEAS2Axon(xyt, xyTmp, xyCir, rCir[ax[j]-1]);
                        }
                    }
                    //cout<<instruction[0]<<'\t'<<instruction[1]<<'\t'<<instruction[2]<<'\t'<<instruction[3]<<endl;
                    //cout<<ax[0]<<'\t'<<ax[1]<<'\t'<<ax[2]<<'\t'<<ax[3]<<endl<<endl;
                    if (instruction[0]==false & instruction[1]==false & instruction[2]==false & instruction[3]==false ) {
                        // Diffusion without hitting an axon
                        xyt=xyTmp;
                        //cout<<"no hit\t"<<j<<endl;
                    }
                    else {
                        // Determine the axon
                        for (j=0; j<4; j++) {
                            axDis[j]=100;
                            if ( instruction[j] ) {
                                xyCir[0]=xCir[ax[j]-1];
                                xyCir[1]=yCir[ax[j]-1];
                                axDis[j]=disOutSheath(xyTmp, xyCir, rCir[ax[j]-1]);
                            }
                        }
                        axDisMin=*min_element(axDis, axDis+4);
                        for (j=0; j<4; j++) {
                            if ( axDis[j]==axDisMin ) {
                                ax_hit=ax[j];
                            }
                        }
                        xyCir[0]=xCir[ax_hit-1];
                        xyCir[1]=yCir[ax_hit-1];
                        
                        v=JKISS(v);
                        //probPerm=probEAS2Myelin(xyt, vp, stepOUT/res, xyCir, rCir[ax_hit-1], probOUT); //cout<<ax_hit<<endl<<probPerm<<endl;
                        if (v.r<probPermOUT) {
                            // Inelastic collision from EAS to myelin outer sheath
                            //xyt=xyTmp;
                            xyt=inelasticCollisionEAS(xyt, vp, stepOUT/res, xyCir, rCir[ax_hit-1], stepMR/res);
                        }
                        else {
                            // Elastic collision in EAS
                            xyCollision=elasticCollisionEAS(xyt, vp, stepOUT/res, xyCir, rCir[ax_hit-1]);
                            
                            // Use xyTmp to save the present position
                            xyTmp=xyt;
                            
                            // Renew the next step after elastic position
                            xyt=xyCollision;
                            
                            // Cancel this step if bounce twice
                            xyTmpG=pixPosition(xyCollision, NPix);
                            aTmp=APix[xyTmpG[0]][xyTmpG[1]];
                            ax[2]=aTmp%Nmax, ax[3]=aTmp/Nmax;
                            
                            xyCir[0]=xCir[ax[2]-1];
                            xyCir[1]=yCir[ax[2]-1];
                            if ( inAxon(xyCollision, xyCir, rCir[ax[2]-1])) {
                                xyt=xyTmp;
                            }
                            
                            xyCir[0]=xCir[ax[3]-1];
                            xyCir[1]=yCir[ax[3]-1];
                            if ( inAxon(xyCollision, xyCir, rCir[ax[3]-1])) {
                                xyt=xyTmp;
                            }
                        }
                    }
                }
                else {
                    // Check whether the particle is in IAS
                    instruction1=false; instruction2=false;
                    if ( ax1 ) {
                        xyCir[0]=xCir[ax1-1]; xyCir[1]=yCir[ax1-1];
                        instruction1=inMyelin(xyt, xyCir, rCir[ax1-1],gratio);
                    }
                    if ( ax2 ) {
                        xyCir[0]=xCir[ax2-1]; xyCir[1]=yCir[ax2-1];
                        instruction2=inMyelin(xyt, xyCir, rCir[ax2-1],gratio);
                    }
                    
                    if ( instruction1==false & instruction2==false ) { // In IAS
                        // Determine the axon
                        instruction1=false; instruction2=false;
                        if ( ax1 ) {
                            xyCir[0]=xCir[ax1-1]; xyCir[1]=yCir[ax1-1];
                            instruction1=inIAS(xyt, xyCir, rCir[ax1-1],gratio);
                        }
                        if ( ax2 ) {
                            xyCir[0]=xCir[ax2-1]; xyCir[1]=yCir[ax2-1];
                            instruction2=inIAS(xyt, xyCir, rCir[ax2-1],gratio);
                        }
                        
                        if ( instruction1 ) {
                            ax[0]=ax1;
                        }
                        else if ( instruction2 ){
                            ax[0]=ax2;
                        }
                        else {
                            ax[0]=0;
                            cout<<"error: particle in IAS has no axon label."<<endl;
                        }
                        
                        // Diffusion in IAS
                        
                        // Initialize velocity
                        v=JKISS(v);
                        phi=2*Pi*v.r;                               // Diffusion direction
                        vp[0]=cos(phi); vp[1]=sin(phi);             // Normalized diffusion velocity
                        
                        // Primitive position after diffusion
                        xyTmp[0]=xyt[0]+stepIN/res*vp[0];
                        xyTmp[1]=xyt[1]+stepIN/res*vp[1];
                        
                        // Determine whether the segment(xyt,xyTmp) overlaps with the myelin
                        xyCir[0]=xCir[ax[0]-1];
                        xyCir[1]=yCir[ax[0]-1];
                        instruction[0]=!inIAS(xyTmp, xyCir, rCir[ax[0]-1], gratio);
                        
                        if (instruction[0]==false) {
                            // Diffusion without hitting the myelin
                            xyt=xyTmp;
                        }
                        else {
                            ax_hit=ax[0];
                            xyCir[0]=xCir[ax_hit-1];
                            xyCir[1]=yCir[ax_hit-1];
                            
                            v=JKISS(v);
                            //probPerm=probIAS2Myelin(xyt, vp, stepIN/res, xyCir, rCir[ax_hit-1], gratio, probIN);
                            
                            if (v.r<probPermIN) {
                                // Inelastic collision from IAS to myelin inner sheath
                                //xyt=xyTmp;
                                xyt=inelasticCollisionIAS(xyt, vp, stepIN/res, xyCir, rCir[ax_hit-1], gratio, stepMR/res);
                            }
                            else {
                                // Elastic collision in IAS
                                xyCollision=elasticCollisionIAS(xyt, vp, stepIN/res, xyCir, rCir[ax_hit-1], gratio);
                                
                                // Cancel this step if bounce twice
                                if ( inIAS(xyCollision, xyCir, rCir[ax_hit-1], gratio) ) {
                                    xyt=xyCollision;
                                }
                            }
                        }
                    }
                    else {  // In myelin
                        //cout<<"Myelin"<<endl;
                        //cout<<"instruction1="<<instruction1<<'\t'<<"instruction2="<<instruction2<<endl;
                        //cout<<"ax1="<<ax1<<'\t'<<"ax2="<<ax2<<endl;
                        
                        // Determine the axon
                        instruction1=false; instruction2=false;
                        if ( ax1 ) {
                            xyCir[0]=xCir[ax1-1]; xyCir[1]=yCir[ax1-1];
                            instruction1=inMyelin(xyt, xyCir, rCir[ax1-1],gratio);
                        }
                        if ( ax2 ) {
                            xyCir[0]=xCir[ax2-1]; xyCir[1]=yCir[ax2-1];
                            instruction2=inMyelin(xyt, xyCir, rCir[ax2-1],gratio);
                        }
                        
                        if ( instruction1 ) {
                            ax[0]=ax1;
                        }
                        else if ( instruction2 ){
                            ax[0]=ax2;
                        }
                        else {
                            ax[0]=0;
                            cout<<"error: particle in Myelin has no axon label."<<endl;
                        }
                        xyCir[0]=xCir[ax[0]-1];
                        xyCir[1]=yCir[ax[0]-1];
                        
                        // Diffusion in myelin
                        
                        // Initialize velocity, vp[0] = phi, vp[1] = randR
                        v=JKISS(v);
                        vp[0]=2*Pi*v.r;
                        
                        v=JKISS(v);
                        vp[1]=v.r;
                        
                        // Diffuse inside myelin with a compressed ellipse propagator
                        xyt=diffuseMyelinEllipseCircum(xyt, stepMC/res, vp[0], xyCir);
                        
                        xyTmp=diffuseMyelinEllipseRadial(xyt, stepMR/res, vp[0], xyCir, vp[1]);
                        
                        if ( inMyelin(xyTmp, xyCir, rCir[ax[0]-1], gratio) ) { // No hitting to myelin sheath
                            xyt=xyTmp;
                        }
                        else if ( inIAS(xyTmp, xyCir, rCir[ax[0]-1], gratio) ){ // Hitting the inner myelin sheath
                            v=JKISS(v);
                            //probPerm=probMyelin2IAS(xyt, stepMR/res, xyCir, rCir[ax[0]-1], gratio, probMY); //cout<<ax[0]<<endl<<probPerm<<endl;
                            if (v.r<probPermMY) {     // Permeate through myelin inner sheath
                                xyt=permeateMyelin2IAS(xyt, stepMR/res*fabs(sin(vp[0])), xyCir, rCir[ax[0]-1], gratio, stepIN/res);
                            }
                            else {                  // Elastic collision from myelin to the inner sheath
                                xyCollision=elasticCollisionMyelin2IAS(xyt, stepMR/res*fabs(sin(vp[0])), xyCir, rCir[ax[0]-1], gratio);
                                // Cancel this step if bounce the outer sheath
                                if ( inMyelin(xyCollision, xyCir, rCir[ax[0]-1], gratio) ) {
                                    xyt=xyCollision;
                                }
                            }
                        }
                        else {  // Hitting the outer myelin sheath
                            v=JKISS(v);
                            //probPerm=probMyelin2EAS(xyt, stepMR/res, xyCir, rCir[ax[0]-1], probMY); //cout<<probPerm<<endl;
                            if (v.r<probPermMY) {     // Permeate through myelin outer sheath
                                xyTmp=xyt;
                                
                                xyt=permeateMyelin2EAS(xyt, stepMR/res*fabs(sin(vp[0])), xyCir, rCir[ax[0]-1], stepOUT/res);
                                
                                xytG=pixPosition(xyt, NPix);
                                a=APix[xytG[0]][xytG[1]];
                                ax1=a%Nmax, ax2=a/Nmax;
                                
                                xyCir[0]=xCir[ax1-1];
                                xyCir[1]=yCir[ax1-1];
                                if ( inAxon(xyt, xyCir, rCir[ax1-1])) {
                                    xyt=xyTmp;
                                }
                                
                                xyCir[0]=xCir[ax2-1];
                                xyCir[1]=yCir[ax2-1];
                                if ( inAxon(xyt, xyCir, rCir[ax2-1])) {
                                    xyt=xyTmp;
                                }
                            }
                            else {                  // Elastic collision from myelin to the outer sheath
                                xyCollision=elasticCollisionMyelin2EAS(xyt, stepMR/res*fabs(sin(vp[0])), xyCir, rCir[ax[0]-1]);
                                // Cancel this step if bounce the inner sheath
                                if ( inMyelin(xyCollision, xyCir, rCir[ax[0]-1], gratio) ) {
                                    xyt=xyCollision;
                                }
                            }
                        }
                        
                    }
                }
                
                // Comply to the periodic boundary condition
                if (xyt[0]>1) {
                    xyt[0]-=1;
                    xjmp[k]+=1;
                }
                else if (xyt[0]<0) {
                    xyt[0]+=1;
                    xjmp[k]-=1;
                }
                
                if (xyt[1]>1) {
                    xyt[1]-=1;
                    yjmp[k]+=1;
                }
                else if (xyt[1]<0) {
                    xyt[1]+=1;
                    yjmp[k]-=1;
                }
                
                if ( (i%Tstep) ==0 ) {
                    dx=(xyt[0]+xjmp[k]-xyi[0])*res;
                    dy=(xyt[1]+yjmp[k]-xyi[1])*res;
                    dxSquare[i/Tstep]+=dx*dx;
                    dySquare[i/Tstep]+=dy*dy;
                    
                    dxQuartic[i/Tstep]+=dx*dx*dx*dx;
                    dyQuartic[i/Tstep]+=dy*dy*dy*dy;
                    
                    for (j=0; j<Nbval; j++) {
                    qt=sqrt(bval[j]/(i*dt*1e-9));
                    xPhase=-qt*dx;
                    yPhase=-qt*dy;
                    
                    gxSignalRe[i/Tstep][j]+=cos(xPhase);
                    gxSignalIm[i/Tstep][j]+=sin(xPhase);
                    gySignalRe[i/Tstep][j]+=cos(yPhase);
                    gySignalIm[i/Tstep][j]+=sin(yPhase);
                    }
                    
                    xt_total[i/Tstep][k]=xyt[0];
                    yt_total[i/Tstep][k]=xyt[1];
                    
                }
                
                                
                xytG=pixPosition(xyt, NPix);                // Position on grid after diffusion
            }
            
        }
        
        end=clock();
        cout << "Done! Elpased time "<<double(diffclock(end,begin)) << " s"<< endl;
        
        
        ofstream fxout("/Users/magda/Desktop/diffusion_simulation/cpp_practice/diffusion_myelin_exchange_ellipse_propagator_speed_up_v5/diffusion_myelin_exchange/x_diffusion");
        ofstream fyout("/Users/magda/Desktop/diffusion_simulation/cpp_practice/diffusion_myelin_exchange_ellipse_propagator_speed_up_v5/diffusion_myelin_exchange/y_diffusion");
        for (i=0; i<timepoints; i++) {
            for (k=0; k<NPar; k++) {
                if (k==NPar-1) {
                    //fxout<<xyt_total[0][i][k]<<endl;
                    //fyout<<xyt_total[1][i][k]<<endl;
                    fxout<<xt_total[i][k]<<endl;
                    fyout<<yt_total[i][k]<<endl;
                }
                else {
                    //fxout<<xyt_total[0][i][k]<<"\t";
                    //fyout<<xyt_total[1][i][k]<<"\t";
                    fxout<<xt_total[i][k]<<"\t";
                    fyout<<yt_total[i][k]<<"\t";
                }
            }
        }
        fxout.close();
        fyout.close();
        
        
        ofstream fdx2out("/Users/magda/Desktop/diffusion_simulation/cpp_practice/diffusion_myelin_exchange_ellipse_propagator_speed_up_v5/diffusion_myelin_exchange/dx2_diffusion");
        ofstream fdy2out("/Users/magda/Desktop/diffusion_simulation/cpp_practice/diffusion_myelin_exchange_ellipse_propagator_speed_up_v5/diffusion_myelin_exchange/dy2_diffusion");
        
        ofstream fdx4out("/Users/magda/Desktop/diffusion_simulation/cpp_practice/diffusion_myelin_exchange_ellipse_propagator_speed_up_v5/diffusion_myelin_exchange/dx4_diffusion");
        ofstream fdy4out("/Users/magda/Desktop/diffusion_simulation/cpp_practice/diffusion_myelin_exchange_ellipse_propagator_speed_up_v5/diffusion_myelin_exchange/dy4_diffusion");
        
        ofstream fgxSRout("/Users/magda/Desktop/diffusion_simulation/cpp_practice/diffusion_myelin_exchange_ellipse_propagator_speed_up_v5/diffusion_myelin_exchange/gxSigRe_diffusion");
        ofstream fgxSIout("/Users/magda/Desktop/diffusion_simulation/cpp_practice/diffusion_myelin_exchange_ellipse_propagator_speed_up_v5/diffusion_myelin_exchange/gxSigIm_diffusion");
        
        ofstream fgySRout("/Users/magda/Desktop/diffusion_simulation/cpp_practice/diffusion_myelin_exchange_ellipse_propagator_speed_up_v5/diffusion_myelin_exchange/gySigRe_diffusion");
        ofstream fgySIout("/Users/magda/Desktop/diffusion_simulation/cpp_practice/diffusion_myelin_exchange_ellipse_propagator_speed_up_v5/diffusion_myelin_exchange/gySigIm_diffusion");
        
        for (i=0; i<timepoints; i++) {
            fdx2out<<dxSquare[i]/NPar<<endl;
            fdy2out<<dySquare[i]/NPar<<endl;
            
            fdx4out<<dxQuartic[i]/NPar<<endl;
            fdy4out<<dyQuartic[i]/NPar<<endl;
            
            for (j=0; j<Nbval; j++) {
                if (j==Nbval-1) {
                    fgxSRout<<gxSignalRe[i][j]/NPar<<endl;
                    fgxSIout<<gxSignalIm[i][j]/NPar<<endl;
                    
                    fgySRout<<gySignalRe[i][j]/NPar<<endl;
                    fgySIout<<gySignalIm[i][j]/NPar<<endl;
                }
                else {
                    fgxSRout<<gxSignalRe[i][j]/NPar<<"\t";
                    fgxSIout<<gxSignalIm[i][j]/NPar<<"\t";
                    
                    fgySRout<<gySignalRe[i][j]/NPar<<"\t";
                    fgySIout<<gySignalIm[i][j]/NPar<<"\t";
                }
            }
        }
        fdx2out.close();
        fdy2out.close();
        fdx4out.close();
        fdy4out.close();
        
        fgxSRout.close();
        fgxSIout.close();
        
        fgySRout.close();
        fgySIout.close();
        
        ofstream paraout ("/Users/magda/Desktop/diffusion_simulation/cpp_practice/diffusion_myelin_exchange_ellipse_propagator_speed_up_v5/diffusion_myelin_exchange/sim_para");
        paraout<<dt*1e-9<<endl<<TN<<endl<<NPar<<endl<<timepoints<<endl;
        paraout<<Din*1e9<<endl<<Dout*1e9<<endl<<kappa<<endl<<Dmr*1e9<<endl<<Dmc*1e9<<endl;
        for (i=0; i<Nbval; i++) {
            paraout<<bval[i]<<"\t";
        }
        paraout.close();


        
}

