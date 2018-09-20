//
//  main.cpp
//  diffusion_myelin_exchange
//
//  Update Journal:
//  -- 6/26/2017: massive job version
//  -- 7/20/2017: do not divide DW signal & cumulants by b0 signal, record particle number in each compartment
//
//  Created by Hong-Hsi Lee in February, 2017.
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

    int main(int argc, char *argv[]) {
        
        clock_t begin=clock();
        clock_t end=clock();
        
        // Define index number
        int i=0, j=0, k=0;
        
        // Initialize RNG
        RNGvalues v;
        v = init_KISS();
        
        //********** Load mictostructure **********
        
        double dt=0.0;              // Time step in ms
        int TN=0;                   // Number of time steps
        const int NPar=5e4;         // Number of particles
        const int timepoints=1e3;   // Number of time points to record
        const int Nbval=10;         // Number of b-values for DKI
        
        double Din=0;               // Diffusion coefficient inside the cylinders in um^2/ms
        double Dout=0;              // Diffusion coefficient outside the cylinders in um^2/ms
        double Dmc=0;               // Diffusion coefficient in myelin, circumferential, in um^2/ms
        double Dmz=0;               // Diffusion coefficient in myelin, parallel to axon, in um^2/ms
        double kappa=0;             // Permeability of a lipid bi-layer in um/ms
        //#define n 5.0             // Number density of lipid bi-layer in 1/um
        double Dmr=0;//kappa/n;     // Diffusivity in myelin, radial, in um^2/ms
        
        string root, target;        // root for loading packing, target for saving files
        unsigned int NPix=0, NAx=0, **APix;
        double res=0.0, *xCir, *yCir, *rCir;
        
        // root
        ifstream myfile_root ("root.txt", ios::in);
        myfile_root>>root;
        
        // simulation parameter
        ifstream myfile0 ("simParamInput.txt", ios::in);
        myfile0>>dt; myfile0>>TN;
        myfile0>>Din; myfile0>>Dout; myfile0>>kappa; myfile0>>Dmr; myfile0>>Dmc; myfile0>>Dmz;        
        
        int Tstep=TN/timepoints;
        double stepIN=sqrt(4.0*dt*Din);     // Step size in IAS in µm
        double stepINZ=sqrt(2.0*dt*Din);    // Step size in IAS along z direction in µm
        double stepOUT=sqrt(4.0*dt*Dout);   // Step size in EAS in µm
        double stepOUTZ=sqrt(2.0*dt*Dout);  // Step size in EAS along z direction in µm
        double stepMC=sqrt(4.0*dt*Dmc);     // Step size in myelin in µm, circumferential
        double stepMR=sqrt(4.0*dt*Dmr);     // Step size in myelin in µm, radial
        double stepMZ=sqrt(2.0*dt*Dmz);     // Step size in myelin in µm, axial
        
        // resolution
        ifstream myfile1 (root+"/phantom_res.txt", ios::in);
        myfile1>>res;
        myfile1.close();
        
        // Pixel # along each side
        ifstream myfile2 (root+"/phantom_NPix.txt", ios::in);
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
        ifstream myfile3 (root+"/phantom_APix.txt", ios::in);
        for (i=0; i<NPix; i++){
            for (j=0; j<NPix; j++){
                myfile3>>APix[i][j];
            }
        }
        myfile3.close();
        
        // Number of axons
        ifstream myfile4 (root+"/phantom_NAx.txt", ios::in);
        myfile4>>NAx;
        myfile4.close();
        
        // Circle center of x-coordinate
        ifstream myfile5 (root+"/phantom_xCir.txt", ios::in);
        i=0; xCir=new double[NAx];
        while (!myfile5.eof()&(i<NAx)){
            myfile5>>xCir[i];
            i++;
        }
        myfile5.close();
        
        // Circle center of y-coordinate
        ifstream myfile6 (root+"/phantom_yCir.txt", ios::in);
        i=0; yCir=new double[NAx];
        while (!myfile6.eof()&(i<NAx)){
            myfile6>>yCir[i];
            i++;
        }
        myfile6.close();
        
        // Circle radius
        ifstream myfile7 (root+"/phantom_rCir.txt", ios::in);
        i=0; rCir= new double[NAx];
        while (!myfile7.eof()&(i<NAx)){
            myfile7>>rCir[i];
            i++;
        }
        myfile7.close();
        
        // The smallest number, which is > NAx, in the base 10
        unsigned int Nmax=0;
        ifstream myfile8 (root+"/phantom_Nmax.txt", ios::in);
        myfile8>>Nmax;
        myfile8.close();
        
        // g-ratio
        double gratio=0;
        ifstream myfile9 (root+"/phantom_gratio.txt", ios::in);
        myfile9>>gratio;
        myfile9.close();
        
        //********** Initialize Particle Positions at circle center *********
        const double probOUT=Pi/4*res*kappa/Dout;       // Probability constant from EAS to myelin
        const double probIN=Pi/4*res*kappa/Din;         // Probability constant from IAS to myelin
        const double probMY=Pi/4*res*kappa/Dmr;         // Probability constant from myelin to IAS/EAS
        double xiPar[NPar]={0}, yiPar[NPar]={0}, ziPar[NPar]={0}; // Initial position
        double xRand=0, yRand=0, zRand=0;               // Random number
        unsigned int a=0;                               // Element of APix matrix
        int ax1=0, ax2=0;                               // Label of axons close to the particle
        vector<double> xyPar(2,0), xyCir(2,0);          // Particle position and Circle center
        vector<int> xyParG(2,0);                        // Particle position on a grid
        bool instruction1=false, instruction2=false;    // true: inside an axon, false: outside an axon
        k=0;
        double NParIAS[timepoints]={0}, NParEAS[timepoints]={0}, NParMY[timepoints]={0};
        while(k<NPar){
            xiPar[k]=0.5; yiPar[k]=0.5;
            k++;
            /*
            v=JKISS(v);
            xRand=v.r;
            v=JKISS(v);
            yRand=v.r;
            v=JKISS(v);
            zRand=v.r;
            xyPar[0]=xRand; xyPar[1]=yRand;
            
            // Whether the particle is close to axons
            xyParG=pixPosition(xyPar,NPix);
            a=APix[xyParG[0]][xyParG[1]];
            ax1=a%Nmax; ax2=a/Nmax;
            
            // If the particle is in IAS, take the initial position
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
            //if ( instruction1==false & instruction2==false ){
                xiPar[k]=xRand; yiPar[k]=yRand; ziPar[k]=zRand;
                k++;
            }*/
        }
        
        // ********** Simulate diffusion **********
        
        int xjmp[NPar]={0}, yjmp[NPar]={0}, zjmp[NPar]={0}; // Number of corssing the boundary box
        vector<double> xyi(3,0), xyt(3,0), xyTmp(3,0);      // Position before/after diffusion
        vector<double> xyCollision(3,0);                    // Position after collision
        vector<int> xytG(2,0), xyTmpG(2,0);                 // Position on grid after diffusion
        xyCir[0]=0; xyCir[1]=0;                             // Position of the axon center
        a=0; unsigned int aTmp=0;                           // Element of APix
        int ax[4]={0}; bool instruction[4]={0};             // Axon label
        double phi=0, sinPhi=0;                             // Diffusion direction and its sine
        vector<double> vp(3,0);                             // Nomalized diffusion velocity
        double axDis[4]={0}, axDisMin=0;                    // Distance to axon surface, and its min
        int ax_hit=0;                                       // Label of the axon hit by the particle
        double probPermMY=stepMR/res*probMY;                // Probability of permeation from myelin to EAS or IAS
        double probPermIN=stepIN/res*probIN;                // Probability of permeation from IAS to Myelin
        double probPermOUT=stepOUT/res*probOUT;             // Probability of permeation from EAS to Myelin
        // Whether the particle is in Myelin or not while recording the signal 
        bool notInMyelin = false, atIAS=false, atEAS=false, atMyelin=false;               
        
        double bval[Nbval]={0}, qt=0;                       // b-value (ms/µm^2), q=\gamma * g * \delta (1/µm)
        double dx=0, dy=0, dz=0, xPhase=0, yPhase=0, zPhase=0;
        
        for (i=0; i<Nbval; i++) {
            bval[i]=0.1*(i+1);
        }
        
        double **xt_total, **yt_total, **zt_total;
        xt_total=new double *[timepoints];
        yt_total=new double *[timepoints];
        zt_total=new double *[timepoints];
        for (i=0; i<timepoints; i++){
            xt_total[i]=new double[NPar];
            yt_total[i]=new double[NPar];
            zt_total[i]=new double[NPar];
            for (j = 0;j<NPar;j++){
                xt_total[i][j]=0;
                yt_total[i][j]=0;
                zt_total[i][j]=0;
            }
        }
        
        double dxSquare[timepoints]={0}, dySquare[timepoints]={0}, dzSquare[timepoints]={0}, dxQuartic[timepoints]={0}, dyQuartic[timepoints]={0}, dzQuartic[timepoints]={0}, gxSignalRe[timepoints][Nbval]={0}, gxSignalIm[timepoints][Nbval]={0}, gySignalRe[timepoints][Nbval]={0}, gySignalIm[timepoints][Nbval]={0}, gzSignalRe[timepoints][Nbval]={0}, gzSignalIm[timepoints][Nbval]={0};
        double NParAS[timepoints]={0};
        
        begin=clock();
        
        for (k=0; k<NPar; k++){         // Loop over each particle from its initial position
            
            xyi[0]=xiPar[k]; xyi[1]=yiPar[k]; xyi[2]=ziPar[k];      // Position before diffusion
            xyt=xyi;                                                // Position after diffusion
            xytG=pixPosition(xyt, NPix);                            // Position on grid after diffusion
            
            atIAS=false; atEAS=false; atMyelin=false;
            
            for (i=0; i<TN; i++){       // Loop over time points
                
                // Axons close to the particle in  previous step
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
                    ax[0]=ax1; ax[1]=ax2;
                    
                    //Diffusion in EAS
                    
                    // Initialize velocity
                    v=JKISS(v);
                    phi=2*Pi*v.r;                               // Diffusion direction
                    vp[0]=cos(phi); vp[1]=sin(phi);             // Normalized diffusion velocity
                    
                    v=JKISS(v);
                    if (v.r<0.5) { vp[2]=1.0; }
                    else { vp[2]=-1.0; }
                    
                    // Primitive position after diffusion
                    xyTmp[0]=xyt[0]+stepOUT/res*vp[0];
                    xyTmp[1]=xyt[1]+stepOUT/res*vp[1];
                    xyTmp[2]=xyt[2]+stepOUTZ/res*vp[2];
                    xyTmpG=pixPosition(xyTmp, NPix);
                    aTmp=APix[xyTmpG[0]][xyTmpG[1]];
                    ax[2]=aTmp%Nmax, ax[3]=aTmp/Nmax;
                    
                    // Decide whether the segment(xyt,xyTmp) overlaps with the axon
                    for (j=0; j<4; j++) {
                        instruction[j]=false;
                        if ( ax[j] ) {
                            xyCir[0]=xCir[ax[j]-1];
                            xyCir[1]=yCir[ax[j]-1];
                            instruction[j]=stepEAS2Axon(xyt, xyTmp, xyCir, rCir[ax[j]-1]);
                        }
                    }
                    if (instruction[0]==false & instruction[1]==false & instruction[2]==false & instruction[3]==false ) {
                        // Diffusion without hitting any axon
                        xyt=xyTmp;
                        atEAS=true; atIAS=false; atMyelin=false;
                    }
                    else {
                        // Decide the axon hit by the particle
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
                        if (v.r<probPermOUT) {
                            // Inelastic collision from EAS to myelin outer sheath
                            xyt=inelasticCollisionEAS(xyt, vp, stepOUT/res, xyCir, rCir[ax_hit-1], stepMR/res, stepOUTZ/res, stepMZ/res);
                            atEAS=false; atIAS=false; atMyelin=true;
                        }
                        else {
                            atEAS=true; atIAS=false; atMyelin=false;
                            // Elastic collision in EAS
                            xyCollision=elasticCollisionEAS(xyt, vp, stepOUT/res, xyCir, rCir[ax_hit-1], stepOUTZ/res);
                            
                            // Use xyTmp to save the present position
                            xyTmp=xyt;
                            
                            // Renew the next step after elastic position
                            xyt=xyCollision;
                            
                            // Cancel this step if bouncing twice
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
                        // Decide the axon where the particle stays in
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
                        
                        v=JKISS(v);
                        if (v.r<0.5) { vp[2]=1; }
                        else { vp[2]=-1; }
                        
                        // Primitive position after diffusion
                        xyTmp[0]=xyt[0]+stepIN/res*vp[0];
                        xyTmp[1]=xyt[1]+stepIN/res*vp[1];
                        xyTmp[2]=xyt[2]+stepINZ/res*vp[2];
                        
                        // Determine whether the segment(xyt,xyTmp) overlaps with the myelin
                        xyCir[0]=xCir[ax[0]-1];
                        xyCir[1]=yCir[ax[0]-1];
                        instruction[0]=!inIAS(xyTmp, xyCir, rCir[ax[0]-1], gratio);
                        
                        if (instruction[0]==false) {
                            // Diffusion in IAS without hitting the myelin
                            xyt=xyTmp;
                            atEAS=false; atIAS=true; atMyelin=false;
                        }
                        else {
                            ax_hit=ax[0];
                            xyCir[0]=xCir[ax_hit-1];
                            xyCir[1]=yCir[ax_hit-1];
                            
                            v=JKISS(v);
                            
                            if (v.r<probPermIN) {
                                // Inelastic collision from IAS to myelin inner sheath
                                xyt=inelasticCollisionIAS(xyt, vp, stepIN/res, xyCir, rCir[ax_hit-1], gratio, stepMR/res, stepINZ/res, stepMZ/res);
                                atEAS=false; atIAS=false; atMyelin=true;
                            }
                            else {
                                atEAS=false; atIAS=true; atMyelin=false;
                                // Elastic collision in IAS
                                xyCollision=elasticCollisionIAS(xyt, vp, stepIN/res, xyCir, rCir[ax_hit-1], gratio, stepINZ/res);
                                
                                // Cancel this step if bouncing twice
                                if ( inIAS(xyCollision, xyCir, rCir[ax_hit-1], gratio) ) {
                                    xyt=xyCollision;
                                }
                            }
                        }
                    }
                    else {  // In myelin
                        // Decide the axon's myelin sheath where the axon stays in
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
                        
                        // Initialize velocity, vp[0] = phi, vp[1] = rand, vp[2] = ±1
                        v=JKISS(v);
                        vp[0]=2*Pi*v.r;
                        
                        v=JKISS(v);
                        vp[1]=v.r;
                        
                        v=JKISS(v);
                        if (v.r<0.5) { vp[2]=1; }
                        else { vp[2]=-1; }
                        
                        // Diffuse inside myelin with a compressed ellipse propagator
                        xyt=diffuseMyelinEllipseCircum(xyt, stepMC/res, vp[0], xyCir);
                        
                        xyTmp=diffuseMyelinEllipseRadial(xyt, stepMR/res, vp[0], xyCir, vp[1], stepMZ/res*vp[2]);
                        
                        if ( inMyelin(xyTmp, xyCir, rCir[ax[0]-1], gratio) ) { // Not hitting the myelin sheath
                            xyt=xyTmp;
                            atEAS=false; atIAS=false; atMyelin=true;
                        }
                        else if ( inIAS(xyTmp, xyCir, rCir[ax[0]-1], gratio) ){ // Hitting the inner myelin sheath
                            v=JKISS(v);
                            if (v.r<probPermMY) {     // Permeate through inner myelin sheath
                                sinPhi = fabs(sin(vp[0]));
                                xyt=permeateMyelin2IAS(xyt, stepMR/res*sinPhi, xyCir, rCir[ax[0]-1], gratio, stepIN/res*sinPhi, stepINZ/res*vp[2], stepMZ/res*vp[2]);
                                atEAS=false; atIAS=true; atMyelin=false;
                            }
                            else {                    // Elastic collision from myelin to inner myelin sheath
                                atEAS=false; atIAS=false; atMyelin=true;
                                sinPhi = fabs(sin(vp[0]));
                                xyCollision=elasticCollisionMyelin2IAS(xyt, stepMR/res*sinPhi, xyCir, rCir[ax[0]-1], gratio, stepMZ/res*vp[2]);
                                // Cancel this step if bouncing the outer sheath
                                if ( inMyelin(xyCollision, xyCir, rCir[ax[0]-1], gratio) ) {
                                    xyt=xyCollision;
                                }
                            }
                        }
                        else {  // Hitting the outer myelin sheath
                            v=JKISS(v);
                            if (v.r<probPermMY) {     // Permeate through outer myelin sheath
                                atEAS=true; atIAS=false; atMyelin=false;
                                xyTmp=xyt;
                                
                                sinPhi = fabs(sin(vp[0]));
                                xyt=permeateMyelin2EAS(xyt, stepMR/res*sinPhi, xyCir, rCir[ax[0]-1], stepOUT/res*sinPhi, stepOUTZ/res*vp[2], stepMZ/res*vp[2]);
                                
                                // Cancel this step if hitting another axon
                                xytG=pixPosition(xyt, NPix);
                                a=APix[xytG[0]][xytG[1]];
                                ax1=a%Nmax, ax2=a/Nmax;
                                
                                xyCir[0]=xCir[ax1-1];
                                xyCir[1]=yCir[ax1-1];
                                if ( inAxon(xyt, xyCir, rCir[ax1-1])) {
                                    xyt=xyTmp;
                                    atEAS=false; atIAS=false; atMyelin=true;
                                }
                                
                                xyCir[0]=xCir[ax2-1];
                                xyCir[1]=yCir[ax2-1];
                                if ( inAxon(xyt, xyCir, rCir[ax2-1])) {
                                    xyt=xyTmp;
                                    atEAS=false; atIAS=false; atMyelin=true;
                                }
                            }
                            else {                  // Elastic collision from myelin to outer myelin sheath
                                atEAS=false; atIAS=false; atMyelin=true;
                                sinPhi = fabs(sin(vp[0]));
                                xyCollision=elasticCollisionMyelin2EAS(xyt, stepMR/res*sinPhi, xyCir, rCir[ax[0]-1], stepMZ/res*vp[2]);
                                // Cancel this step if hitting the inner myelin sheath
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
                
                if (xyt[2]>1) {
                    xyt[2]-=1;
                    zjmp[k]+=1;
                }
                else if (xyt[2]<0) {
                    xyt[2]+=1;
                    zjmp[k]-=1;
                }
                
                if ( (i%Tstep) ==0 ) {
                    NParIAS[i/Tstep]+=atIAS;
                    NParEAS[i/Tstep]+=atEAS;
                    NParMY[i/Tstep]+=atMyelin;
                    
                    notInMyelin=atEAS | atIAS;
                    NParAS[i/Tstep]+=notInMyelin;
                    
                    dx=(xyt[0]+xjmp[k]-xyi[0])*res;
                    dy=(xyt[1]+yjmp[k]-xyi[1])*res;
                    dz=(xyt[2]+zjmp[k]-xyi[2])*res;
                    dxSquare[i/Tstep]+=dx*dx*notInMyelin;
                    dySquare[i/Tstep]+=dy*dy*notInMyelin;
                    dzSquare[i/Tstep]+=dz*dz*notInMyelin;
                    
                    dxQuartic[i/Tstep]+=dx*dx*dx*dx*notInMyelin;
                    dyQuartic[i/Tstep]+=dy*dy*dy*dy*notInMyelin;
                    dzQuartic[i/Tstep]+=dz*dz*dz*dz*notInMyelin;
                    
                    for (j=0; j<Nbval; j++) {
                        qt=sqrt(bval[j]/(i*dt));
                        xPhase=-qt*dx;
                        yPhase=-qt*dy;
                        zPhase=-qt*dz;
                        
                        gxSignalRe[i/Tstep][j]+=cos(xPhase)*notInMyelin;
                        gxSignalIm[i/Tstep][j]+=sin(xPhase)*notInMyelin;
                        gySignalRe[i/Tstep][j]+=cos(yPhase)*notInMyelin;
                        gySignalIm[i/Tstep][j]+=sin(yPhase)*notInMyelin;
                        gzSignalRe[i/Tstep][j]+=cos(zPhase)*notInMyelin;
                        gzSignalIm[i/Tstep][j]+=sin(zPhase)*notInMyelin;

                    }
                    
                    xt_total[i/Tstep][k]=xyt[0]*res;
                    yt_total[i/Tstep][k]=xyt[1]*res;
                    zt_total[i/Tstep][k]=xyt[2]*res;
                    
                }
                xytG=pixPosition(xyt, NPix);                // Position on grid after diffusion
            }
            
        }
        
        end=clock();
        cout << "Done! Elpased time "<<double(diffclock(end,begin)) << " s"<< endl;
        
        
        ofstream fxout("x_diffusion");
        ofstream fyout("y_diffusion");
        ofstream fzout("z_diffusion");
        for (i=0; i<timepoints; i++) {
            for (k=0; k<NPar; k++) {
                if (k==NPar-1) {
                    fxout<<xt_total[i][k]<<endl;
                    fyout<<yt_total[i][k]<<endl;
                    fzout<<zt_total[i][k]<<endl;
                }
                else {
                    fxout<<xt_total[i][k]<<"\t";
                    fyout<<yt_total[i][k]<<"\t";
                    fzout<<zt_total[i][k]<<"\t";
                }
            }
        }
        fxout.close();
        fyout.close();
        fzout.close();
        
        
        ofstream fdx2out("dx2_diffusion");
        ofstream fdy2out("dy2_diffusion");
        ofstream fdz2out("dz2_diffusion");
        
        ofstream fdx4out("dx4_diffusion");
        ofstream fdy4out("dy4_diffusion");
        ofstream fdz4out("dz4_diffusion");
        
        ofstream fgxSRout("gxSigRe_diffusion");
        ofstream fgxSIout("gxSigIm_diffusion");
        
        ofstream fgySRout("gySigRe_diffusion");
        ofstream fgySIout("gySigIm_diffusion");
        
        ofstream fgzSRout("gzSigRe_diffusion");
        ofstream fgzSIout("gzSigIm_diffusion");
        
        ofstream fNparASout("NParAS");
        
        ofstream fNparIASout("NParIAS");
        ofstream fNparEASout("NParEAS");
        ofstream fNparMYout("NParMY");
        
        for (i=0; i<timepoints; i++) {
            fdx2out<<dxSquare[i]<<endl;
            fdy2out<<dySquare[i]<<endl;
            fdz2out<<dzSquare[i]<<endl;
            
            fdx4out<<dxQuartic[i]<<endl;
            fdy4out<<dyQuartic[i]<<endl;
            fdz4out<<dzQuartic[i]<<endl;
            
            fNparASout<<NParAS[i]<<endl;
            fNparIASout<<NParIAS[i]<<endl;
            fNparEASout<<NParEAS[i]<<endl;
            fNparMYout<<NParMY[i]<<endl;
            
            for (j=0; j<Nbval; j++) {
                if (j==Nbval-1) {
                    fgxSRout<<gxSignalRe[i][j]<<endl;
                    fgxSIout<<gxSignalIm[i][j]<<endl;
                    
                    fgySRout<<gySignalRe[i][j]<<endl;
                    fgySIout<<gySignalIm[i][j]<<endl;
                    
                    fgzSRout<<gzSignalRe[i][j]<<endl;
                    fgzSIout<<gzSignalIm[i][j]<<endl;
                }
                else {
                    fgxSRout<<gxSignalRe[i][j]<<"\t";
                    fgxSIout<<gxSignalIm[i][j]<<"\t";
                    
                    fgySRout<<gySignalRe[i][j]<<"\t";
                    fgySIout<<gySignalIm[i][j]<<"\t";
                    
                    fgzSRout<<gzSignalRe[i][j]<<"\t";
                    fgzSIout<<gzSignalIm[i][j]<<"\t";
                }
            }
        }
        fdx2out.close();
        fdy2out.close();
        fdz2out.close();
        fdx4out.close();
        fdy4out.close();
        fdz4out.close();
        
        fgxSRout.close();
        fgxSIout.close();
        
        fgySRout.close();
        fgySIout.close();
        
        fgzSRout.close();
        fgzSIout.close();
        
        fNparASout.close();
        fNparIASout.close();
        fNparEASout.close();
        fNparMYout.close();
        
        ofstream paraout ("sim_para");
        paraout<<dt<<endl<<TN<<endl<<NPar<<endl<<timepoints<<endl;
        paraout<<Din<<endl<<Dout<<endl<<kappa<<endl<<Dmr<<endl<<Dmc<<endl<<Dmz<<endl;
        for (i=0; i<Nbval; i++) {
            paraout<<bval[i]<<"\t";
        }
        paraout.close();
}

