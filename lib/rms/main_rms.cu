//  diffusion_3D_RMS
//
//  Update Journal:
//  -- 03/11/2019: equal step length random leap, one fiber, IAS, 3D, cuda version
//  -- 04/25/2019: implement mitochondria, high permeability, short T2, same diffusivity as IAS
//  -- 01/29/2020: implement generalized realistic microstructure simulator (RMS): elastic reflection, water exchange (no permeability) and T2 relaxation
//  -- 07/13/2020: PGSE signal and thrust pointer
//  -- 01/10/2021: I/O .bin instead of .txt
//  -- 03/04/2021: get parameters from command line
//  -- 03/24/2021: create anisotropic free diffusivity
//  Created by Hong-Hsi Lee, Massachusetts General Hospital
//  Modified by Ricardo Coronado-Leija, New York University

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
#include <unistd.h> // small
#include <getopt.h> // large
#include <cstdint>  // uint8_t

#include <cuda.h>
#include <curand_kernel.h>
#include <cuda_runtime.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>

using namespace std;
    
#define Pi 3.14159265
// #define timepoints 1000
#define Ngrad_max 2000 
#define nite 4
#define Nc_max 3
#define DEB 0
#define Nbin 400

// ********** cuda kernel **********
__device__ double atomAdd(double* address, double val)
{
    unsigned long long int* address_as_ull =
    (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val +
                                             __longlong_as_double(assumed)));
        
        // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (assumed != old);
    
    return __longlong_as_double(old);
}

__global__ void setup_kernel(curandStatePhilox4_32_10_t *state, unsigned long seed){
    int idx = threadIdx.x+blockDim.x*blockIdx.x;
    curand_init(seed, idx, 0, &state[idx]);
}
__global__ void propagate(
    curandStatePhilox4_32_10_t *state, 
    double *sig0, 
    double *sig,   // pgse
    double *sigRe, // narrow pulse
    double *dx2, 
    double *dx4, 
    double *NPar_count,
    double *NPar_bin,
    const double *TD, 
    const int TN, 
    const double *stepDa, 
    const double *stepDr, 
    const double *T2, 
    const double *Pij, 
    const int Nc,  
    const int NPar, 
    const double res, 
    const double dt,
    const uint8_t *APix,
    const int NPix1, 
    const int NPix2, 
    const int NPix3, 
    const double *bvec,
    const double *bval, 
    const double *grad,
    const double *Delta,
    const double *delta,   
    const double *echotime,
    const int Ngrad,  
    const int timepoints,
    const bool pgse_flag,
    const int initFlag){

// cuda
int idx    = threadIdx.x + blockDim.x * blockIdx.x;
int stride = blockDim.x * gridDim.x;
curandStatePhilox4_32_10_t localstate=state[idx];    
// time step of the points that will be saved   
int Tstep  = TN/timepoints;

for(int k = idx; k < NPar; k += stride){

// ########################################################################################### //
    double step[Nc_max] = {0};
    for(int kk = 0; kk < Nc_max; kk++){
       step[kk] = max(stepDa[kk],stepDr[kk]);
       }

    #if DEB == 1
    for(int kk = 0; kk < Nc_max; kk++){
    printf("normalized axial  step size compartment %d = %.8f\n",kk,stepDa[kk]);
    printf("normalized radial step size compartment %d = %.8f\n",kk,stepDr[kk]);
    printf("normalized max    step size compartment %d = %.8f\n",kk,step[kk]);
    }  // kk  
    #endif

    // Random number
    double vRand = 0;    
    // Particle position on a grid
    int xParGi[3] = {0}, xParGj[3] = {0};
        
    // xi:   initial particle position
    // xt:   particle position at the i-th step
    // tt:   distance between particle and y-z, x-z, x-y box wall
    // vt:   a unit vector indicating the hopping direction
    // xTmp: a temporary variable to save the position
    double xi[3] = {0}, xt[3] = {0}, tt[3] = {0}, vt[3] = {0};
    // double xTmp[3] = {0};
        
    // for new direction    
    double cos_theta = 0, sin_theta = 0;
    int tidx = 0, nidx = 0;
     
    // size of the box
    int NPix[3] = {0};
    NPix[0] = NPix1; NPix[1] = NPix2; NPix[2] = NPix3;

    // Signal weighted by T2 relaxation
    double s0 = 0;

    // The flags of hitting the medium boundary
    bool flip0[3] = {false}, flip1[3] = {false};
    // bool flip0_tmp[3]={false}, flip1_tmp[3]={false};
    // bool flipx0=false, flipx1=false, flipy0=false, flipy1=false, flipz0=false, flipz1=false;

    // q = \gamma * g * \delta (1/µm)
    double qx = 0;                 // narrow pulse
    double phase = 0; // wide pulse
    double dx = 0, dy = 0, dz = 0;
    double x1[Ngrad_max*3] = {0};
    double pm[3] = {0}; pm[0] = 1; pm[1] = 1; pm[2] = 1;

    // fstep: remaining fraction of one step
    // tmp:   temporary variable
    // tmin:  the shortest distance between the particle and y-z, x-z, x-y plane (% of step)
    double fstep = 0, tmp = 0, tmin = 0; //, maxStep;

    // Elements of lookup table APix
    int ai = 0, aj = 0;
    unsigned long long int aidx;

    // The box wall hit by the particle. 1:y-z plane, 2:x-z plane, 3:x-y plane
    int ii_hit = 0;

    // The time staying in compartments
    double t[Nc_max] = {0};
//        printf("step size=%.4f\n",step[0]);
//        printf("step size=%.4f\n",step[1]);    
//        printf("before\n");
    // ********** Initialize Particle Positions inside compartments (Apix[index] > 0) ********** //
    while (1){
        if ( initFlag == 4 ) {
            xi[0] = 0.5*static_cast<double>(NPix1);
            xi[1] = 0.5*static_cast<double>(NPix2);
            xi[2] = 0.5*static_cast<double>(NPix3);
            break;
        } else {
        // random initial position
        xi[0] = curand_uniform_double(&localstate)*static_cast<double>(NPix1);
        xi[1] = curand_uniform_double(&localstate)*static_cast<double>(NPix2);
        xi[2] = curand_uniform_double(&localstate)*static_cast<double>(NPix3);
        
        // Whether the particle is inside compartments
        xParGi[0] = floor(xi[0]); 
        xParGi[1] = floor(xi[1]); 
        xParGi[2] = floor(xi[2]);
        aidx = (unsigned long long int)( (unsigned long long int) ( ((unsigned long long int)(NPix2) )*((unsigned long long int)(NPix1))*((unsigned long long int)(xParGi[2])) + ((unsigned long long int)(NPix1))*((unsigned long long int)(xParGi[1])) ) + (unsigned long long int)(xParGi[0]) );

        if ( initFlag == 3 ) {
            if ( APix[ (unsigned long long int)(aidx) ] !=0 ){ break; }
        } else {
            if ( APix[ (unsigned long long int)(aidx) ] == initFlag ){ break; } 
        }
        
        }
    } // while
    
    // printf("after\n");
    // ********** Simulate diffusion ********** //
    // updating particle position
    xt[0] = xi[0]; 
    xt[1] = xi[1]; 
    xt[2] = xi[2];
    // position of the particle on the voxelized geometry
    xParGi[0]=floor(xt[0]); 
    xParGi[1]=floor(xt[1]); 
    xParGi[2]=floor(xt[2]);
    // compartment
    aidx = (unsigned long long int)( (unsigned long long int) ( ((unsigned long long int)(NPix2) )*((unsigned long long int)(NPix1))*((unsigned long long int)(xParGi[2])) + ((unsigned long long int)(NPix1))*((unsigned long long int)(xParGi[1])) ) + (unsigned long long int)(xParGi[0]) );
    ai = (int)( APix[ (unsigned long long int)(aidx) ] );
        
    #if DEB == 1    
    printf("Particle Initial Position = %.2f, %.2f, %.2f\n",xt[0],xt[1],xt[2]);
    printf("Current Voxel = %d, %d, %d : %d\n",xParGi[0],xParGi[1],xParGi[2],ai);
    #endif 

    for(int i = 0; i < TN; i++){
    // ========================================================================================= //
    // ========================================================================================= //
    // ========================================================================================= //
        fstep = 1.0;
        // set boundary hitting flags to false (for over three directions)
        for(int jj = 0; jj < 3; jj++){
            flip0[jj] = false; 
            flip1[jj] = false;
            }
        // flipx0=false; flipx1=false; flipy0=false; flipy1=false; flipz0=false; flipz1=false;
        
        // Random elevation 
        vRand     = curand_uniform_double(&localstate);
        cos_theta = 1.0 - 2.0*vRand;
        sin_theta = 2.0*sqrt(vRand*(1.0 - vRand)); //sin(acos(cos_theta));
        // random azimuthal
        vRand = curand_uniform_double(&localstate);
        vt[0] = sin_theta*cos(2.0*Pi*vRand); // x
        vt[1] = sin_theta*sin(2.0*Pi*vRand); // y
        vt[2] = cos_theta;                   // z
       
        #if DEB == 1    
        printf("New direction = %.2f, %.2f, %.2f\n",vt[0],vt[1],vt[2]);
        #endif 

//         // rescaling new direction according to stepDa, stepDr
//         maxStep = max(stepDa[ai],stepDr[ai]);
//         vt[0]   = stepDr[ai]*vt[0]/maxStep;
//         vt[1]   = stepDr[ai]*vt[1]/maxStep;
//         vt[2]   = stepDa[ai]*vt[2]/maxStep;
// 
//         #if DEB == 1    
//         printf("Scaled direction = %.2f, %.2f, %.2f\n",vt[0],vt[1],vt[2]);
//         #endif 
     
            
        for(int j = 0; (j < nite) && (fstep > 0); j++){ // several attempts ? 
            
            // check hitting walls in the voxel    
            tmin   =  2.0; 
            ii_hit = -1.0;
            for(int ii = 0; ii < 3; ii++){
                if(vt[ii] > 0.0){ // positive movement   i       i+1
                    // distance particle to next voxel   |  *<--->|
                    tmp = static_cast<double>(xParGi[ii]) + 1.0 - xt[ii];
                    } 
                else if(vt[ii] < 0.0){ // negative movement      i       i+1 
                    // distance particle to current voxel (neg)  |<--->*  |         
                    tmp = static_cast<double>(xParGi[ii]) - xt[ii];
                    } 
                else{
                    tmp = 2.0; // vt[ii] == 0 ?
                    } // if-elseif-else     
            // tmp=static_cast<double>(xParGi[ii]) + fmax(0.0,static_cast<double>(vt[ii]>0.0)) - xt[ii];
                if(fabs(tmp) > fabs(vt[ii])){
                    // movement of the particle is still inside of the voxel
                    tt[ii] = 2.0;
                    } 
                else{
                    // movement could get particle outside of the voxel
                    tt[ii] = max(0.0,tmp/vt[ii]); // which percent (positive) of vt is tmp ?
                    } // if-else
                // particle crossed at least one voxel-wall (just save th shortest distance)
                // tmin is the percentage of the step needed to get to the respective wall    
                if(tt[ii] < tmin){
                    tmin   = tt[ii];
                    ii_hit = ii;
                    } // if

                #if DEB == 1    
                printf("ii = %d, vt = %.2f, tmp = %.2f, tt = %.2f, tmin = %.2f, ii_hit = %d \n",
                        ii,vt[ii],tmp,tt[ii],tmin,ii_hit);
                #endif  

                } // i
                
            // if ( ii_hit<0 ) {
            //    printf("Error1: walker does not encounter the box wall, tmin=%.4f, xGi=%i,%i,%i, xGj=%i,%i,%i\n",tmin,xParGi[0],xParGi[1],xParGi[2],xParGj[0],xParGj[1],xParGj[2]);
            //    break;
            // }
                
            // Should be impossible (tmin initial is zero, then it is assigned 2.0 or max(0,+num))
            if(tmin < 0.0){
                printf("Error: walker jumps into wrong direction, tmin = %.4f\n",tmin);
                break;
                }
            
            // Update Position 
            if( (fstep*step[ai]) >= tmin){ // just move until reaching the closest wall
                fstep = fstep - tmin/step[ai]; // remaining fraction of the step
                xt[0] = xt[0] + tmin*vt[0];
                xt[1] = xt[1] + tmin*vt[1];
                xt[2] = xt[2] + tmin*vt[2];
                t[ai] = t[ai] + tmin/step[ai];
                } 
            else{ // move the full step
                xt[0] = xt[0] + step[ai]*fstep*vt[0];
                xt[1] = xt[1] + step[ai]*fstep*vt[1];
                xt[2] = xt[2] + step[ai]*fstep*vt[2];
                t[ai] = t[ai] + fstep;

                #if DEB == 1    
                printf("Final Position = %.2f, %.2f, %.2f\n",xt[0],xt[1],xt[2]);
                #endif 

               // printf("Error3: walker does not encounter the box wall, tmin=%.4f, xGi=%i,%i,%i, xGj=%i,%i,%i, vt=%.4f, ii_hit=%i\n",tmin,xParGi[0],xParGi[1],xParGi[2],xParGj[0],xParGj[1],xParGj[2],vt[ii_hit],ii_hit);
                break;
                } // if-else
                
            // set position on the lookup table     
            xParGj[0] = xParGi[0]; 
            xParGj[1] = xParGi[1]; 
            xParGj[2] = xParGi[2];

            #if DEB == 1    
            printf("Position Closest Wall %d = %.2f, %.2f, %.2f : fstep=%f t=%f\n",i,xt[0],xt[1],xt[2],fstep,t[ai]);
            #endif 

           // if ( (~flipx0) && (~flipx1) && (~flipy0) && (~flipy1) && (~flipz0) && (~flipz1) ) {
           // if ( (~flip0[0]) && (~flip1[0]) && (~flip0[1]) && (~flip1[1]) && (~flip0[2]) && (~flip1[2]) ) {

            // I believe this goes here, otherwise we would be accessing a negative indice
            // but originally it was after the section on
            // behavior of the particle on the wall
            if ( ii_hit < 0 ){
                printf("Error: walker does not encounter the box wall, tmin=%.4f, xGi=%i,%i,%i, xGj=%i,%i,%i\n",tmin,xParGi[0],xParGi[1],xParGi[2],xParGj[0],xParGj[1],xParGj[2]);
                break;
                } // if

            // behavior of the particle on the wall
            if(vt[ii_hit] >= 0.0){ // positive movement
                xParGj[ii_hit] = xParGj[ii_hit] + 1; // go next voxel
                if(xParGj[ii_hit]  >= NPix[ii_hit]){ // if upper boundary of the medium reached 
                    xParGj[ii_hit]  = NPix[ii_hit] - 1;
                    flip1[ii_hit]   = true;          // mirror boundary conditions
                    vt[ii_hit]      = -vt[ii_hit];   // oposite direction 
                    continue;
                    } // if
                } // if 
            else{ // negative movement
                xParGj[ii_hit] = xParGj[ii_hit] - 1; // go previous voxel
                if(xParGj[ii_hit] < 0){              // if bottom boundary of the medium reached   
                    xParGj[ii_hit] = 0;
                    flip0[ii_hit]  = true;           // mirror boundary conditions
                    vt[ii_hit]     = -vt[ii_hit];    // oposite direction 
                    continue;
                    } // if
                } // else
                
            // get compartment at current voxel    
            aidx = (unsigned long long int)( (unsigned long long int) ( ((unsigned long long int)(NPix2) )*((unsigned long long int)(NPix1))*((unsigned long long int)(xParGj[2])) + ((unsigned long long int)(NPix1))*((unsigned long long int)(xParGj[1])) ) + (unsigned long long int)(xParGj[0]) );    

            aj = (int) ( APix[ (unsigned long long int)(aidx) ] );
               
            // Checking permeability
            if(Pij[ai*Nc + aj] > (1 - 1e-10) ){ // same compartment: full permeability
                ai = aj; // update voxel of hitted wall
                xParGi[ii_hit] = xParGj[ii_hit]; 
                } 
            else if(Pij[ai*Nc + aj] < 1e-10){ // different compartments: no permeability 
                vt[ii_hit] = -vt[ii_hit]; // change direction of diffusion
               // printf("ai=%i, aj=%i, xGj=%i,%i,%i\n",ai,aj,xParGj[0],xParGj[1],xParGj[2]);
               } 
            else{ // different compartment: some permeability
                // printf("ai=%i, aj=%i, Pij=%.4f\n",ai,aj,Pij[ai*Nc+aj]);
                vRand = curand_uniform_double(&localstate);
                if(vRand < Pij[ai*Nc+aj]){
                    ai = aj; // update voxel of hitted wall
                    xParGi[ii_hit] = xParGj[ii_hit];
                    } 
                else{
                    vt[ii_hit] = -vt[ii_hit]; // change direction of diffusion
                    }
                } // if-elseif-else

            #if DEB == 1    
            printf("New direction 2 = %.2f, %.2f, %.2f\n",vt[0],vt[1],vt[2]);
            printf("Current Voxel 2 = %d, %d, %d : %d\n",xParGi[0],xParGi[1],xParGi[2],ai);
            #endif 
                
            } // j - several attempts
            
        // Apply flipping due to mirroring boundary conditions (why outside j and xi ??? )
        for(int jj = 0; jj < 3; jj++){
            if(flip0[jj]){
                xi[jj] = -xi[jj];
                pm[jj] = -pm[jj];
                } // if
            if(flip1[jj]){
                xi[jj] = 2.0*static_cast<double>(NPix[jj])-xi[jj];
                pm[jj] = -pm[jj];
                } // if
            // reset flipping    
            flip0[jj] = false; 
            flip1[jj] = false;
            } // jj

        #if DEB == 1    
        printf("Current Position i %d = %.2f, %.2f, %.2f\n",i,xi[0],xi[1],xi[2]);
        printf("Current Position t %d = %.2f, %.2f, %.2f\n",i,xt[0],xt[1],xt[2]);
        #endif 

    if(pgse_flag){  // Gradient parameters of PGSE: Delta, delta, |g|, gx, gy, gz  (only saved if asked)

        // Displacement
        dx = (xt[0] - xi[0])*res*pm[0];
        dy = (xt[1] - xi[1])*res*pm[1];
        dz = (xt[2] - xi[2])*res*pm[2];

        // add phase 
        for(int j = 0; j < Ngrad; j++){
            // First Pulse
            if(static_cast<double>(i+1)*dt <= delta[j]){
                x1[3*j]+=dx; x1[3*j+1]+=dy; x1[3*j+2]+=dz;
            } // 1st
            else if ( (static_cast<double>(i+1)*dt > Delta[j]) & (static_cast<double>(i+1)*dt <= (Delta[j]+delta[j])) ){
                x1[3*j]-=dx; x1[3*j+1]-=dy; x1[3*j+2]-=dz;
            } // 2nd
        } // j

//         // add phase
//         for(int j = 0; j < Ngrad; j++){
//             // Second pulse
//             if( ( static_cast<double>(i + 1)*dt >= Delta[j] ) & ( static_cast<double>(i + 1)*dt < (Delta[j] + delta[j]) ) ){
//                 phase[j] -= grad[j] * (dx*bvec[j*3] + dy*bvec[j*3 + 1] + dz*bvec[j*3 + 2])*dt;
//                 } // 2nd 
//             } // j   
        
        // Readout (here it may be better to change t for echotime) and have one for each ngrad
        if ( i == (TN-1) ){ // just saves the last time
            s0 = 0.0;
            for(int j = 0; j < Nc; j++){
                s0 += (t[j]/T2[j]);
                } // j
            s0 = exp(-1.0*s0);
            // pgse signal
            for(int j = 0; j < Ngrad; j++){
                phase = grad[j] * ( x1[j*3]*bvec[j*3] + x1[j*3+1]*bvec[j*3+1] + x1[j*3+2]*bvec[j*3+2] )*dt;
                atomAdd(&sig[j],s0*cos(phase));
                } // j
        } // if
    } // pgse

    if ( (i%Tstep) == 0 ) { // Save moment tensor for dx^2 and dx^4, and signal for the b-table
        
        // T2 Relaxation
        s0 = 0.0;
        for(int j = 0; j < Nc; j++){
            s0 = s0 + (t[j]/T2[j]);
            } // j
        s0 = exp(-1.0*s0); // s0 = 1.0;
        
        tidx = i/Tstep;
        nidx = Nc*tidx + ai;
        
        atomAdd(&sig0[tidx],s0);       // Update s0 signal 
        atomAdd(&NPar_count[nidx],1);  // Update npar (per compartment) counter 
        
        if (initFlag==4){
        nidx = floor(sqrt((xt[0]-static_cast<double>(NPix1)/2.0)*(xt[0]-static_cast<double>(NPix1)/2.0) + (xt[1]-static_cast<double>(NPix2)/2.0)*(xt[1]-static_cast<double>(NPix2)/2.0) + (xt[2]-static_cast<double>(NPix3)/2.0)*(xt[2]-static_cast<double>(NPix3)/2.0))/(static_cast<double>(NPix1)/2.0/static_cast<double>(Nbin)));
        if (nidx<=Nbin){ atomAdd(&NPar_bin[Nbin*tidx+nidx],1); }
        }

        // Displacement
        dx = (xt[0] - xi[0])*res*pm[0];
        dy = (xt[1] - xi[1])*res*pm[1];
        dz = (xt[2] - xi[2])*res*pm[2];

        #if DEB == 1    
        printf("Current Position i %d = %.4f, %.4f, %.4f\n",i,xi[0]*res,xi[1]*res,xi[2]*res);
        printf("Current Position t %d = %.4f, %.4f, %.4f\n",i,xt[0]*res,xt[1]*res,xt[2]*res);
        printf("Current Displace t %d = %.4f, %.4f, %.4f\n",tidx,dx,dy,dz);
        #endif 
     
        // Second Order Moment
        atomAdd(&dx2[6*tidx+0],s0*dx*dx);
        atomAdd(&dx2[6*tidx+1],s0*dx*dy);
        atomAdd(&dx2[6*tidx+2],s0*dx*dz);
        atomAdd(&dx2[6*tidx+3],s0*dy*dy);
        atomAdd(&dx2[6*tidx+4],s0*dy*dz);
        atomAdd(&dx2[6*tidx+5],s0*dz*dz);
        
        // Fourth Order Moment
        atomAdd(&dx4[15*tidx+0] ,s0*dx*dx*dx*dx);
        atomAdd(&dx4[15*tidx+1] ,s0*dx*dx*dx*dy);
        atomAdd(&dx4[15*tidx+2] ,s0*dx*dx*dx*dz);
        atomAdd(&dx4[15*tidx+3] ,s0*dx*dx*dy*dy);
        atomAdd(&dx4[15*tidx+4] ,s0*dx*dx*dy*dz);
        atomAdd(&dx4[15*tidx+5] ,s0*dx*dx*dz*dz);
        atomAdd(&dx4[15*tidx+6] ,s0*dx*dy*dy*dy);
        atomAdd(&dx4[15*tidx+7] ,s0*dx*dy*dy*dz);
        atomAdd(&dx4[15*tidx+8] ,s0*dx*dy*dz*dz);
        atomAdd(&dx4[15*tidx+9] ,s0*dx*dz*dz*dz);
        atomAdd(&dx4[15*tidx+10],s0*dy*dy*dy*dy);
        atomAdd(&dx4[15*tidx+11],s0*dy*dy*dy*dz);
        atomAdd(&dx4[15*tidx+12],s0*dy*dy*dz*dz);
        atomAdd(&dx4[15*tidx+13],s0*dy*dz*dz*dz);
        atomAdd(&dx4[15*tidx+14],s0*dz*dz*dz*dz);
        
        // Diffusion Signal Narrow Pulse Limit (this will always be saved)
        // Acording to Rafael Patiño:
        // \delta Phase = a(t)*gamma*g(t)*z(t)*(\delta t)
        // Phase        = \sum_i \delta Phase_i
        // Signal       = \sum_i e^{ i * Phase _i}
        // Signal_{re}  = \sum_i cos( Phase_i )
        for(int j = 0; j < Ngrad; j++){ // bvalue = t*q^2 => q = sqrt( bvalue/t )
            qx = sqrt( bval[j] / TD[tidx] ) * (dx*bvec[j*3] + dy*bvec[j*3+1] + dz*bvec[j*3+2]);
            atomAdd(&sigRe[Ngrad*tidx+j] , s0*cos(qx));
            } // j

        } // if Tstep (save state)
    // ========================================================================================= //
    // ========================================================================================= //
    // ========================================================================================= //
    } // steps i 
// ########################################################################################### //    
} // particles k 
state[idx]=localstate;
} // function

// .......................................................................... //
// ########################################################################## //
// ########################################################################## //
// ########################################################################## //
// .......................................................................... //

int show_help(){
printf(
// === Non Optional Arguments === //
"rms [ options ] input btable output\n"
"\n"
"\tinput \n"
"\t\tname of the input (binary) file with the substrate/medium.\n"
"\tbtable\n"
"\t\tname of the btable (.txt) file to generate the diffusion signal.\n" 
"\t\tEach line of the file must be (with b in ms/um^2):\n"
"\t\t\t x1 y1 z1 b1\n"
"\t\t\t x2 y2 z2 b2\n"
"\t\t\t x3 y3 z3 b3\n"
"\t\t\t .  .  .  . \n"
"\t\t\t .  .  .  . \n"
"\t\t\t .  .  .  . \n"
"\t\t\t xn yn zn bn\n"
"\toutput\n"
"\t\tname of the output (binary) file with the results of the simulation.\n"     
"\n"
"Realistic Microstructure Simulator (RMS).\n"
"Gets D(t), K(t) and the time dependent diffusion signal for the given substrate.\n"
"\n"
// === Optional Arguments === //
"Options:\n"
"\n"
// === More Needed Optional Arguments === //
"\t -particles (p)\n"
"\t\tNumber of random walkers to simulate.\n" 
"\t\t(default: -particles 1e6).\n" 
"\t -time (t)\n"
"\t\ttotal time of the simulation, in ms.\n" 
"\t\t(default: -time 100).\n" 
"\t -dintra (i)\n"
"\t\tintra axonal space (IAS) diffusivity at time = 0, in um^2/ms.\n" 
"\t\t(default: -dintra 2).\n" 
"\t -dextra (e)\n"
"\t\textra axonal space (EAS) diffusivity at time = 0, in um^2/ms.\n" 
"\t\t(default: -dextra 2).\n" 
"\t -voxstep (v)\n"
"\t\tvoxstep takes values between 0 and 1.\n"
"\t\tthe length step is equal to voxstep * (voxel size).\n"    
"\t\tthe time step is then equal to (length step)^2/(6*D).\n"  
"\t\t(default: -voxstep 0.9).\n" 
"\t -mspoints (m)\n"
"\t\tNumber of points per ms to sample D(t), K(t) and diffusion signal.\n"
"\t\t(default: -mspoints 10).\n" 
"\t -space (s)\n"
"\t\tSelect in which space the simulation will be performed:\n"
"\t\t1 = IAS, 2 = EAS.\n"
"\t\t(default: -space 1).\n" 
"\t -pgse (g)\n"
"\t\tgenerates the diffusion signal for a PGSE sequence, saving it with the name specified, \n"
"\t\t\\Delta (D), \\delta (d) and TE should be indicated in the btable.\n"
"\t\tEach line of the file must be (with b in ms/um^2):\n"
"\t\t\t x1 y1 z1 b1 D1 d1 TE1\n"
"\t\t\t x2 y2 z2 b2 D2 d2 TE2\n"
"\t\t\t x3 y3 z3 b3 D3 d3 TE3\n"
"\t\t\t .  .  .  .  .  .   . \n"
"\t\t\t .  .  .  .  .  .   . \n"
"\t\t\t .  .  .  .  .  .   . \n"
"\t\t\t xn yn zn bn Dn dn TEn\n"
"\t -help (h)\n"
"\t\tshow this help\n"
"\n"
"Author:\n"
"Hong Hsi Lee (ORCID 0000-0002-3663-6559)\n"
"Ricardo Coronado-Leija\n"
"\n"
"References:\n"
"Lee, et al., Journal of Neuroscience Methods 2021 (doi:10.1016/j.jneumeth.2020.109018)\n"
"Fieremans & Lee., Neuroimage 2018 (doi:10.1016/j.neuroimage.2018.06.046)\n"
"\n");
return 0;
}    

int main(int argc, char *argv[]) {

// To check time
clock_t begin = clock();
clock_t end   = clock();

// Initialize seed for RNG
unsigned long seed = 0;
FILE *urandom;
urandom = fopen("/dev/random", "r");
fread(&seed,sizeof(seed),1,urandom);
fclose(urandom);

// .......................................................................... //
// ########################################################################## //
// .......................................................................... //

int option_index; // getopt_long_only stores the option index here.
int goloval;      // getopt_long_only returns the option value here

int    NPar       = 1000000;      // # particles
double dtime      = 100.0;        // total diffusion time 
double vstep      = 0.9;          // length_step = vstep*voxel_size;    
int    Nc         = 3;            // # compartments
double Dintra     = 2.0;          // ICS diffusion at time = 0 ms
double Dextra     = 2.0;          // ECS diffusion at time = 0 ms
double t2intra    = 1e20;         // ICS T2 time 
double t2extra    = 1e20;         // ECS T2 time
double kappa      = 0.0;          // permeability, um/ms
int    mspoints   = 10;           // points per mili-second
int    initFlag   = 3;            // 1 = ICS, 2 = ECS, 3 = ICS + ECS, 4 = center
bool   pgse_flag  = false;        // use PGSE sequence
bool   debug_flag = false;        // for debugging
bool   help_flag  = false;        // show help
char namePGSE[500];               // name output PGSE signal

// ======================== Reading Input Parameters ======================== //

// options
struct option long_options[] = {
// These options don’t set a flag. We distinguish them by their indices.
{"particles",       required_argument, 0, 'p'},
{"time",            required_argument, 0, 't'},
{"voxstep",         required_argument, 0, 'v'},
{"compartments",    required_argument, 0, 'c'},
{"dintra",          required_argument, 0, 'i'},
{"dextra",          required_argument, 0, 'e'},
{"t2intra",         required_argument, 0, 'r'},
{"t2extra",         required_argument, 0, 'a'},
{"permeability",    required_argument, 0, 'k'},
{"mspoints",        required_argument, 0, 'm'},
{"space",           required_argument, 0, 's'},
{"pgse",            required_argument, 0, 'g'},
{"debug",           no_argument      , 0, 'd'},
{"help",            no_argument      , 0, 'h'},
{0, 0, 0, 0}
};

while(true){
// Detecting the next option
goloval = getopt_long_only (argc, argv, "p:t:v:c:i:e:r:m:s:g:dh", long_options, &option_index);

// Detect the end of the options and break the while. 
if (goloval == -1) break;

switch (goloval){
  case 'p':
     NPar = atoi(optarg);
     break;
  case 't':
     dtime = atof(optarg);   
     break;
  case 'v':
     vstep = atof(optarg);   
     if(vstep < 0.0001 || 0.9999 < vstep){
     cout << "vstep = " << vstep << "is not a valid value,"
          << " it should have a value between 0 and 1, settting to 0.9" << endl;
     vstep = 0.9;   
     }
     break;
  case 'c':
     Nc = atoi(optarg);
     break;
  case 'i':
     Dintra = atof(optarg);     
  case 'e':
     Dextra = atof(optarg);          
     break;
  case 'r':
     t2intra = atof(optarg);     
     break;
  case 'a':
     t2extra = atof(optarg);     
     break;
  case 'k':
     kappa = atof(optarg);     
     break;
  case 'm':
     mspoints = atoi(optarg);
     if(mspoints < 1){
     cout << "mspoints = " << mspoints << "is not a valid value,"
          << " it should be an integer larger than 1, settting to 10" << endl; 
     mspoints = 10;
     }  
     break;
  case 's':
     initFlag = atoi(optarg);
     if(initFlag < 1 || initFlag > 4){ 
       cout << "Incorrect option, settting to IAS + EAS." << endl;
       initFlag = 1;
       } //
     break;     
  case 'g':
     pgse_flag = true;
     strcpy(namePGSE,optarg); 
     cout << "pgse_flag: " << pgse_flag << " optarg: " << namePGSE << endl;
     break;     
  case 'd':
     debug_flag = true;
     break;     
  case 'h':
     help_flag = true;
     break;
  case '?':
     //getopt_long already printed an error message.
     break;
  default:
     abort();
  } // switch
} // while (1)

// Showing the help
if(help_flag || argc == 1){
show_help();
exit(1);
}

// ======================================================================= //      
// Remaining command line arguments (not options).
if(argc - optind != 3){
printf("Expecting 3 arguments %d provided. Use option -help to show the help.\n",(argc - optind));
return 1;
}

// --- Non-Optional arguments --- //
char name_input[500], basename_output[500], name_btable[500];

sprintf(name_input     ,"%s",argv[optind++]);
sprintf(name_btable    ,"%s",argv[optind++]);
sprintf(basename_output,"%s",argv[optind++]);

// .......................................................................... //
// ########################################################################## //
// .......................................................................... //

// =========================== Load Mictostructure ========================== //

unsigned int i, j, k, Nbs, Nbvec, NPix1, NPix2, NPix3, vs;
double Nbvals, res, dl, dt, TN_TP; 
int TN, timetemp, timepoints;

// ===== Read b-table ===== //
cout << "Load b-table: " << name_btable << endl;
Nbs = 0; // counter elements in btable
ifstream myfile(name_btable,ios::in);
if(myfile.is_open()){
while(!myfile.eof()){
myfile >> Nbvals; Nbs++;
} // while
myfile.close();
} // if
else{
cout << "ERROR: Can't open file: " << name_btable << endl; 
exit(1); 
}
// decompose in columns (pgse 7 or narrow 4)
if(pgse_flag){
Nbvec = Nbs/7;
Nbs   = Nbvec*7;
if(Nbvec > Ngrad_max){
cout << "ERROR: Maxim number of ngrad should be: " << Ngrad_max << endl; 
} // if
} // pgse
else{
Nbvec = Nbs/4;
Nbs   = Nbvec*4;
} // narrow
// diffusion = mu m^2 / ms => bvalue = ms / mu m^2 => gyromagnetic ratio = 1 (s^(-1) T^(-1)) => 1/ (ms * mT)
// Saving btab directly on a host vector
// Narrow: [gx gy gz bval] 
// PGSE: gx, gy, gz, |g|, Delta, delta, TE 
thrust::host_vector<double> grad(Nbvec);   // gradient (mT/mu m)
thrust::host_vector<double> Delta(Nbvec);  // Delta    (ms)
thrust::host_vector<double> delta(Nbvec);  // delta    (ms)
thrust::host_vector<double> TE(Nbvec);     // TE       (ms)
thrust::host_vector<double> bval(Nbvec);   // bvector  (unitary)
thrust::host_vector<double> bvec(Nbvec*3); // bvalue   (ms / mu m^2)
ifstream myfile1(name_btable,ios::in);
// reading btable
if(pgse_flag){
cout << "Computing bvalues" << endl;
for(i = 0; i < Nbvec; i++) {
myfile1 >> bvec[i*3+0]; 
myfile1 >> bvec[i*3+1]; 
myfile1 >> bvec[i*3+2];
myfile1 >> bval[i];
myfile1 >> Delta[i];
myfile1 >> delta[i];
myfile1 >> TE[i]; 
grad[i] = sqrt(bval[i]/delta[i]/delta[i]/( Delta[i] - delta[i]/3.0 ));
} // i 
} // pgse
else{
for(i = 0; i < Nbvec; i++) {
myfile1 >> bvec[i*3+0]; 
myfile1 >> bvec[i*3+1]; 
myfile1 >> bvec[i*3+2];
myfile1 >> bval[i]; 
} // i  
} // narrow
myfile1.close();
// // debug
// for(i = 0; i < Nbvec; i++){
// cout << bvec[i*3+0] << " " << bvec[i*3+1] << " "  << bvec[i*3+2] << " " << bval[i] << endl;
// } // i 
cout << "N elements = " << Nbs << " Nbvec = " << Nbvec << endl;

// ===== Read medium substrate (.bin) ===== //
cout << "Reading Substrate: " << name_input << endl;
FILE *finput;   
finput = fopen(name_input, "rb");
if(!finput){
cout << "ERROR: Can't open file: " << name_input << endl; 
exit(1); 
} // if
fread(&NPix1,sizeof(unsigned int),1,finput); 
fread(&NPix2,sizeof(unsigned int),1,finput); 
fread(&NPix3,sizeof(unsigned int),1,finput); 
fread(&vs   ,sizeof(unsigned int),1,finput); 

uint8_t *AcharPix  = new uint8_t[NPix1*NPix2*NPix3];
fread(AcharPix,sizeof(uint8_t),NPix1*NPix2*NPix3,finput);
fclose(finput);

thrust::host_vector<uint8_t> AcharPixB(NPix1*NPix2*NPix3);
for(i = 0; i < NPix1*NPix2*NPix3; i++){
    AcharPixB[i] = AcharPix[i];
} // i

// Copying raw array directly to a device_vector APix
// thrust::device_vector<bool> d_APix(AcharPixB,AcharPixB+NPix1*NPix2*NPix3); 
thrust::device_vector<uint8_t> d_APix=AcharPixB;
res = ((double)vs)/1000.0;
cout << "nx: " << NPix1 << " ny: " << NPix2 << " nz: " << NPix3 << " vs: " << vs << endl;
delete[] AcharPix;

dl         = vstep*res;
// all compartments have equal D for now
dt         = (double) ( Dintra > Dextra ? (dl*dl)/(6.0*Dintra) :  (dl*dl)/(6.0*Dextra) ); 
TN         = ceil(dtime/dt); 
timetemp   = round( ((double)TN)/dtime);
timepoints = (mspoints <= timetemp) ? round(dtime*((double)mspoints) ) : TN;
TN_TP      = ((double)TN)/((double)timepoints);
TN         = ceil( ((double)TN)/((double)timepoints) )*timepoints;

char *initFlagS = (char*)(initFlag == 1 ? "ICS" : ( initFlag == 2 ? "ECS" : ( initFlag == 3 ? "ICS+ECS" : ( initFlag == 4 ? "Center" : "NOT VALID") ) ) );

cout << endl;
cout << "Simulation Parameters (RMS)"      << endl       << endl;
cout << "Random Walkers: "                 << NPar       << endl;
cout << "Voxel Size: "                     << res        << " (um)"    << endl;
cout << "Length Step: "                    << dl         << " (um)"    << endl;
cout << "Time Step: "                      << dt         << " (ms)"    << endl;
cout << "Total Time Simulation: "          << dtime      << " (ms)"    << endl;
cout << "Number of Time Steps: "           << TN         << endl;
cout << "Dintra(t=0): "                    << Dintra     << " um^2/ms" << endl;
cout << "Dextra(t=0): "                    << Dextra     << " um^2/ms" << endl;
cout << "Permeability: "                   << kappa      << " um/ms"   << endl;
cout << "T2intra: "                        << t2intra    << " ms"      << endl;
cout << "T2extra: "                        << t2extra    << " ms"      << endl;
cout << "Simulated Space: "                << initFlagS  << endl;
cout << "Points per mili-second: "         << mspoints   << endl;
cout << "Total Points to Save: "           << timepoints << endl;
cout << "Total Steps/Points: "             << TN_TP      << endl;
cout << endl;

// .......................................................................... //
// ########################################################################## //
// .......................................................................... //

thrust::host_vector<double> Da(Nc);
thrust::host_vector<double> Dr(Nc);
thrust::host_vector<double> T2(Nc);
thrust::host_vector<double> stepDa(Nc); 
thrust::host_vector<double> stepDr(Nc); 
thrust::host_vector<double> Pij(Nc*Nc);
        
Da[0] = 0.0; Da[1] = Dintra; Da[2] = Dextra;
Dr[0] = 0.0; Dr[1] = Dintra; Dr[2] = Dextra;

T2[0] = 1e20/dt;
T2[1] = t2intra/dt;
T2[2] = t2extra/dt;
    
// Step size in IAS in µm (normalized)
for(i = 0; i < Nc; i++){
stepDa[i] = sqrt(6.0*dt*Da[i])/res;
stepDr[i] = sqrt(6.0*dt*Dr[i])/res;
#if DEB == 1 
cout << "step size axial  = " << stepDa[i] << endl;
cout << "step size radial = " << stepDr[i] << endl;
#endif
} // i   

k = 0;
for (i=0; i<Nc; i++){
    for (j=0; j<Nc; j++){
        if ( (i==0) || (j==0) ) {
            Pij[k]=0.0;
        } else if ( i==j ) {
            Pij[k]=1.0;
        } else {
            Pij[k]= kappa * (stepDa[i]*res) / Da[i] * 2/3 / ( 1 + kappa/2*(stepDa[i]*res/Da[i] + stepDa[j]*res/Da[j]) * 2/3 );
//             Pij[k]=fmin(1.0,sqrt(D[j]/D[i]));
        }
        cout<<"permeation probability = i" << i << " j" << j << " Pij" << Pij[k]<<endl;
        k++;
    } // j
} // i
        
thrust::host_vector<double> TD(timepoints);
for (i = 0 ; i < timepoints; i++){
TD[i]=(i*(TN/timepoints)+1)*dt;
}

// ********** Simulate diffusion **********

// Initialize state of RNG
int blockSize = 64;
int numBlocks = (NPar + blockSize - 1) / blockSize;
cout<<numBlocks<<endl<<blockSize<<endl;

thrust::device_vector<curandStatePhilox4_32_10_t> devState(numBlocks*blockSize);
setup_kernel<<<numBlocks, blockSize>>>(devState.data().get(),seed);

// Initialize output
thrust::host_vector<double> sig0(timepoints);
thrust::host_vector<double> sig(Nbvec);
thrust::host_vector<double> sigRe(timepoints*Nbvec);
// thrust::host_vector<double> sigIm(timepoints*Nbvec);
thrust::host_vector<double> dx2(timepoints*6);
thrust::host_vector<double> dx4(timepoints*15);
thrust::host_vector<double> NPar_count(timepoints*Nc);
thrust::host_vector<double> NPar_bin(timepoints*Nbin);
for (i=0;i<timepoints;i++){ sig0[i]=0; }
for (i=0;i<Nbvec;i++){ sig[i]=0; }
for (i=0;i<timepoints*Nbvec;i++){ sigRe[i]=0; }
// for (i=0;i<timepoints*Nbvec;i++){ sigIm[i]=0; }    
for (i=0;i<timepoints*6;i++){ dx2[i]=0; }
for (i=0;i<timepoints*15;i++){ dx4[i]=0; }
for (i=0;i<timepoints*Nc;i++){ NPar_count[i]=0; }
for (i=0;i<timepoints*Nbin;i++){ NPar_bin[i]=0; }

// Move data from host to device
thrust::device_vector<double> d_sig0=sig0;
thrust::device_vector<double> d_sig=sig;
thrust::device_vector<double> d_sigRe=sigRe;
// thrust::device_vector<double> d_sigIm=sigIm;
thrust::device_vector<double> d_dx2=dx2;
thrust::device_vector<double> d_dx4=dx4;
thrust::device_vector<double> d_NPar_count = NPar_count;
thrust::device_vector<double> d_NPar_bin = NPar_bin;
thrust::device_vector<double> d_TD=TD;
thrust::device_vector<double> d_stepDa=stepDa;
thrust::device_vector<double> d_stepDr=stepDr;
thrust::device_vector<double> d_T2=T2;
thrust::device_vector<double> d_Pij=Pij;
thrust::device_vector<double> d_grad=grad;
thrust::device_vector<double> d_Delta=Delta;
thrust::device_vector<double> d_delta=delta;
thrust::device_vector<double> d_TE=TE;
thrust::device_vector<double> d_bval=bval;
thrust::device_vector<double> d_bvec=bvec;
// thrust::device_vector<int> d_APix=APix;

//        double *NPar_count; cudaMallocManaged(&NPar_count,sizeof(double)); NPar_count[0] = 0;
// Parallel computation
begin=clock();
propagate<<<numBlocks, blockSize>>>(devState.data().get(),
                                    d_sig0.data().get(), d_sig.data().get(), d_sigRe.data().get(),
                                    d_dx2.data().get(), d_dx4.data().get(),
                                    d_NPar_count.data().get(),
                                    d_NPar_bin.data().get(),
                                    d_TD.data().get(), TN,
                                    d_stepDa.data().get(), d_stepDr.data().get(), d_T2.data().get(), d_Pij.data().get(),
                                    Nc, NPar, res, dt,
                                    d_APix.data().get(), NPix1, NPix2, NPix3,
                                    d_bvec.data().get(), d_bval.data().get(), d_grad.data().get(),
                                    d_Delta.data().get(), d_delta.data().get(), d_TE.data().get(),
                                    Nbvec, timepoints, pgse_flag, initFlag);
cudaDeviceSynchronize();
end=clock();
cout << "Done! Elpased time "<<double((end-begin)/CLOCKS_PER_SEC) << " s"<< endl;


thrust::copy(d_sig0.begin(), d_sig0.end(), sig0.begin());
thrust::copy(d_sig.begin(), d_sig.end(), sig.begin());
thrust::copy(d_sigRe.begin(), d_sigRe.end(), sigRe.begin());
// thrust::copy(d_sigIm.begin(), d_sigIm.end(), sigIm.begin());
thrust::copy(d_dx2.begin(), d_dx2.end(), dx2.begin());
thrust::copy(d_dx4.begin(), d_dx4.end(), dx4.begin());
thrust::copy(d_NPar_count.begin(), d_NPar_count.end(), NPar_count.begin());
thrust::copy(d_NPar_bin.begin(), d_NPar_bin.end(), NPar_bin.begin());
        
// .......................................................................... //
// ########################################################################## //
// ########################################################################## //
// ########################################################################## //
// .......................................................................... //

if(debug_flag){
// Save number of particles per compartment
char nameDB[500];
sprintf(nameDB,"%s_NPar_count.txt",basename_output);
ofstream fNpout(nameDB);
for(i = 0; i < timepoints; i++){
for(j = 0; j < Nc; j++){
if(j == Nc-1){
fNpout << NPar_count[i*Nc+j] << endl;
}
else{
fNpout<<NPar_count[i*Nc+j]<<"\t";
} // if-else
} // j
} // i
fNpout.close();
} // if debug

// saving PGSE signal
if(pgse_flag){
// copy    
double *SignalPGSE = new double[Nbvec];
for(k = 0; k < Nbvec; k++){ 
SignalPGSE[k] = sig[k];   
// cout << sig[k] << endl ;
}
// saving
// for now it will be saved in a different file  
FILE * pgseFile;      
pgseFile = fopen (namePGSE, "wb");
fwrite (SignalPGSE,sizeof(double),Nbvec,pgseFile);
fclose (pgseFile);
delete[] SignalPGSE;
} // pgse

// Saving output results
// timepoints(1)+TDout(timepoints)+dx2(timepoints*6)+dx4(timepoints*15)+sig0(timepoints)+sig(timepoints)*Nbvec);
unsigned int h, y, x;
unsigned int nmat     = 1;
unsigned int nparams  = 15;
unsigned int ncolumns = 1 + 6 + 15 + 1 + Nbvec;
if ( initFlag == 4 ) { ncolumns += Nbin; }
unsigned int matelem  = ncolumns*timepoints;
unsigned int ntotal   = nmat*ncolumns*timepoints+(nmat*3)+nparams;
double *Out = new double[ntotal];

// Some Parameters (16) of the Simulation
y = 0;
Out[y] = nparams;    y++;  // 01. Number of initialization points (randomly will be 1)
Out[y] = nmat;       y++;  // 02. Number of initialization points (randomly will be 1)
Out[y] = ncolumns;   y++;  // 03. Number of columns of each ini matrix (= 1 x 6 x 15 x Nbvecs)
Out[y] = timepoints; y++;  // 04. Number of time points being saved
Out[y] = dt;         y++;  // 05. Time step in ms
Out[y] = TN;         y++;  // 06. Number of time steps
Out[y] = NPar;       y++;  // 07. Number of random walkers
Out[y] = Nbvec;      y++;  // 08. Number of Diffusion Encoding Gradients
Out[y] = Dintra;     y++;  // 09. Diffusion coefficient inside the axon in um^2/ms
Out[y] = Dextra;     y++;  // 10. Diffusion coefficient outside the axon in um^2/ms
Out[y] = kappa;      y++;  // 11. Permeability of a lipid bi-layer in um/ms
Out[y] = initFlag;   y++;  // 12. Ini pos: 1=ICS, 2=ECS, 3=ICS+ECS+myelin, 4=center
Out[y] = res;        y++;  // 13. Box size
Out[y] = t2intra;    y++;  // 14. T2 inside cell, ms
Out[y] = t2extra;    y++;  // 15. T2 outside cell, ms

// Initialization Points (y is the offset)
for(k = 0; k < nmat; k++){ 
Out[y] = -1; y++;  
Out[y] = -1; y++;  
Out[y] = -1; y++;  
// cout << Out[y-3] << " " << Out[y-2] << " " << Out[y-1] << " " << endl ;
}

// Set outputs to save on a .bin file
for(k = 0; k < nmat; k++){
for(i = 0; i < timepoints; i++){    
// Timepoints
j      = 0;     
x      = y + matelem*k + timepoints*j + i; // (timepoints*j = 0) ToFix: TD should be saved only once
Out[x] = TD[i];
j++;
// sig0   
x      = y + matelem*k + timepoints*j + i; // 
Out[x] = sig0[i];
j++;
// Moment Tensor for dx^2 
for(h = 0; h < 6; h++){
x      = y + matelem*k + timepoints*j + i;     
Out[x] = dx2[k*timepoints*6 + i*6 + h]; 
j++;
} // h
// Moment Tensor for dx^4
for(h = 0; h < 15; h++){
x      = y + matelem*k + timepoints*j + i;     
Out[x] = dx4[k*timepoints*15 + i*15 + h]; 
j++;
} // h
// Signal Re
for(h = 0; h < Nbvec; h++){
x      = y + matelem*k + timepoints*j + i;      
Out[x] = sigRe[k*timepoints*Nbvec + i*Nbvec + h]; 
j++;
} // h

if ( initFlag == 4 ) {
// Particle density in bin
for(h = 0; h < Nbin; h++){
x      = y + matelem*k + timepoints*j + i;
Out[x] = NPar_bin[k*timepoints*Nbin + i*Nbin + h];
j++;
} // h
} // if
        
} // i
} // k

//char nameFile[300];
FILE * outFile;
//sprintf(nameFile,"%s_results.bin",basename_output);
outFile = fopen (basename_output, "wb");
fwrite (Out,sizeof(double),ntotal,outFile);
fclose (outFile);

delete[] Out;

return 0;
} // main


