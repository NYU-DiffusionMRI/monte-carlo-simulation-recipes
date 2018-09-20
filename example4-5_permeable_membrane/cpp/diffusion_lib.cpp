//
//  diffusion_lib.cpp
//
//  Created by Hong-Hsi Lee in February, 2017.
//

#include "diffusion_lib.h"
#include "RNG.h"
#include <math.h>
#include <vector>

using namespace std;

double diffclock(clock_t clock1,clock_t clock2)
{
	double diffticks=clock1-clock2;
	double diffms=(diffticks)/CLOCKS_PER_SEC;
	return diffms;
}

void vecSubtract (const vector<double> &x, const vector<double> &y, vector<double> &z) {
    // Subtraction of two vectors x-y = z
    z[0]=x[0]-y[0];
    z[1]=x[1]-y[1];
}

inline double vecDistance (const vector<double> &x, const vector<double> &y) {
    return sqrt( pow(x[0]-y[0], 2) + pow(x[1]-y[1], 2) );
}

inline double vecInProduct (const vector<double> &x, const vector<double> &y) {
    // Inner product of two vector
    return ( x[0]*y[0] + x[1]*y[1] );
}

inline double vecNorm (const vector<double> &x) {
    return sqrt( x[0]*x[0] + x[1]*x[1] );
}

vector<int> pixPosition ( vector<double> x, unsigned int NPix ) {
    vector<int> xPix(2,0);
    
    if ( x[0]<0 ) { x[0]+=1; }
    if ( x[0]>1 ) { x[0]-=1; }
    
    if ( x[1]<0 ) { x[1]+=1; }
    if ( x[1]>1 ) { x[1]-=1; }
    
    xPix[0]=floor(x[0]*NPix);
    xPix[1]=floor(x[1]*NPix);
    return xPix;
}

void translateXc ( const vector<double> &x, vector<double> &xc) {
    // Translate circle center xc to make it as close to the position x as possible
    int ti=0, tj=0;
    double disSquare = (x[0]-xc[0])*(x[0]-xc[0])+(x[1]-xc[1])*(x[1]-xc[1]), disSquare_t=0;
    for (int i=-1; i<2; i++) {
        for (int j=-1; j<2; j++) {
            disSquare_t=(x[0]-xc[0]-i)*(x[0]-xc[0]-i)+(x[1]-xc[1]-j)*(x[1]-xc[1]-j);
            if (disSquare_t<disSquare) {
                disSquare=disSquare_t;
                ti=i;
                tj=j;
            }
        }
    }
    xc[0]+=ti;
    xc[1]+=tj;
}

bool inAxon ( const vector<double> &x, vector<double> xc, const double &rc ) {
    // If the point x is in the circle (xc,rc), return 1; if not, return 0.
    
    // Translate circle center xc to make it as close to the position xt as possible
    translateXc(x,xc);
    
    return ( ( (x[0]-xc[0])*(x[0]-xc[0]) + (x[1]-xc[1])*(x[1]-xc[1]) ) <= rc*rc );    
}

bool inMyelin ( const vector<double> &x, vector<double> xc, const double &rc, const double &g ) {
    // If the point x is in the myelin (xc,rc,rc*g), return 1; if not, return 0.
    
    // Translate circle center xc to make it as close to the position xt as possible
    translateXc(x,xc);
    
    double disSquare = (x[0]-xc[0])*(x[0]-xc[0]) + (x[1]-xc[1])*(x[1]-xc[1]);
    return ( ( disSquare <= rc*rc ) & ( disSquare >=rc*rc*g*g ) );    
}

bool inIAS ( const vector<double> &x, vector<double> xc, const double &rc, const double &g ) {
    // If the point x is in the Intra-axonal space (xc,rc*g), return 1; if not, return 0.
    
    // Translate circle center xc to make it as close to the position xt as possible
    translateXc(x,xc);
    return ( ( (x[0]-xc[0])*(x[0]-xc[0]) + (x[1]-xc[1])*(x[1]-xc[1]) ) < rc*rc*g*g ) ;
}

bool stepEAS2Axon (const vector<double> &xi, const vector<double> &xt, vector<double> xc, const double &rc ) {
    // If segment(xi,xt) overlaps with circle (xc,rc), return 1; if not, return 0.
    
    // Translate circle center xc to make it as close to the position xt as possible
    translateXc(xt,xc);
    
    // If xt is in axon, segment overlaps with circle
    if ( ( (xt[0]-xc[0])*(xt[0]-xc[0]) + (xt[1]-xc[1])*(xt[1]-xc[1]) ) <= rc*rc ) {
        return true;
    }
    
    // L: a*x[0]+b*x[1]+c=0, d=distance of xc to L
    double a=xi[1]-xt[1], b=xt[0]-xi[0], c=xi[0]*xt[1]-xi[1]*xt[0], sos_ab= a*a+b*b, d=fabs( a*xc[0]+b*xc[1]+c )/sqrt( sos_ab );
    if (d>rc) {
        return false;
    }
    
    // xl: a point on L closest to xc
    vector<double> xl(2,0);
    xl[0] = ( b*( b*xc[0]-a*xc[1] )-a*c )/( sos_ab );
    xl[1] = ( a*( -b*xc[0]+a*xc[1] )-b*c )/( sos_ab );
    
    // xl is in axon, but xi and xt are not in axon.
    return ( ( (xi[0]-xl[0])*(xt[0]-xl[0])+(xi[1]-xl[1])*(xt[1]-xl[1]) ) <= 0 );    
}

double disOutSheath ( const vector<double> &x, vector<double> xc, const double &rc ) {
    // The distance between the point x and surface of the circle (xc,rc)
    
    // Translate circle center xc to make it as close to the position x as possible
    translateXc(x,xc);
    
    return fabs( sqrt( pow(x[0]-xc[0],2) + pow(x[1]-xc[1],2) ) - rc );
}

double disInSheath ( const vector<double> &x, vector<double> xc, const double &rc, const double &g ) {
    // The distance between the point x and surface of the circle (xc,rc*g)
    
    // Translate circle center xc to make it as close to the position x as possible
    translateXc(x,xc);
    
    return fabs( sqrt( pow(x[0]-xc[0],2) + pow(x[1]-xc[1],2) ) - rc*g );
}

vector<double> elasticCollisionEAS (const vector<double> &x, const vector<double> &v, const double &dx, vector<double> xc, const double &rc, const double &dz) {
    // Elastic collision of a particle x in EAS onto an circle (xc,rc) with a velocity v and a diffusion length dx
    
    // Translate circle center xc to make it as close to the position (x + dx*v) as possible
    vector<double> xTmp(3,0);
    xTmp[0]=x[0]+dx*v[0];
    xTmp[1]=x[1]+dx*v[1];
    xTmp[2]=x[2]+dz*v[2];
    translateXc(xTmp,xc);
    
    // distance( x+t*v, xc )==rc, solve t
    double a=0,b=0,c=0,t1=0,t2=0,t=0;
    a=v[0]*v[0] + v[1]*v[1];
    b=2*(x[0]-xc[0])*v[0] + 2*(x[1]-xc[1])*v[1];
    c=pow(x[0]-xc[0],2) + pow(x[1]-xc[1],2) - rc*rc;
    
    // xt: final position, xm: contact point on myelin, n: surface unit vector, v1: reflecting velocity, discri: discriminant
    vector<double> xt(3,0), xm(2,0), n(2,0), v1(2,0);
    double discri=b*b-4*a*c;
    if (discri<0) {          // Does not hit the axon
        xt=xTmp;
    }
    else {
        discri=sqrt(discri);
        t1=0.5/a*( -b+discri );
        t2=0.5/a*( -b-discri );
        t=min(t1,t2);
        if ( (t>=dx) | (t<0) ) {            // Does not hit the axon
            xt=xTmp;
        }
        else {                  // Hit the axon
            xm[0]=x[0]+t*v[0];
            xm[1]=x[1]+t*v[1];
            
            // n parallel to (xm-xc_t)
            n[0]=(xm[0]-xc[0])/vecDistance(xm, xc);
            n[1]=(xm[1]-xc[1])/vecDistance(xm, xc);
            
            // v1 = v - 2(v \dot n)*n
            v1[0]=v[0]-2*vecInProduct(v,n)*n[0];
            v1[1]=v[1]-2*vecInProduct(v,n)*n[1];
            
            // xt = xm + (dx-t)*v1
            xt[0]=xm[0]+(dx-t)*v1[0];
            xt[1]=xm[1]+(dx-t)*v1[1];
            xt[2]=xTmp[2];
        }
    }
    return xt;
}

vector<double> inelasticCollisionEAS (const vector<double>  &x, const vector<double> &v, const double &dx, vector<double> xc, const double &rc, const double &dxR, const double &dz, const double &dzMY) {
    // Inelastic collision of a particle x in EAS onto an circle (xc,rc) with a velocity v and a diffusion length dx. The dxR is the step size in myelin
    
    // Translate circle center xc to make it as close to the position (x + dx*v) as possible
    vector<double> xTmp(3,0);
    xTmp[0]=x[0]+dx*v[0];
    xTmp[1]=x[1]+dx*v[1];
    xTmp[2]=x[2]+dz*v[2];
    translateXc(xTmp,xc);
    
    // distance( x+t*v, xc )==rc, solve t
    double a=0,b=0,c=0,t1=0,t2=0,t=0;
    a=v[0]*v[0] + v[1]*v[1];
    b=2*(x[0]-xc[0])*v[0] + 2*(x[1]-xc[1])*v[1];
    c=pow(x[0]-xc[0],2) + pow(x[1]-xc[1],2) - rc*rc;
    
    // xt: final position, xm: contact point on myelin, discri: discriminant
    vector<double> xt(3,0), xm(2,0), n(2,0);
    double discri=b*b-4*a*c;
    if (discri<0) {          // Does not hit the axon
        xt=xTmp;
    }
    else {
        discri=sqrt(discri);
        t1=0.5/a*( -b+discri );
        t2=0.5/a*( -b-discri );
        t=min(t1,t2);
        if ( (t>=dx) | (t<0) ) {            // Does not hit the axon
            xt=xTmp;
        }
        else {                  // Hit the axon
            // xm = x + t*v;
            xm[0]=x[0]+t*v[0];
            xm[1]=x[1]+t*v[1];
            
            // n parallel to (xm-xc_t)
            n[0]=(xm[0]-xc[0])/vecDistance(xm, xc);
            n[1]=(xm[1]-xc[1])/vecDistance(xm, xc);
            
            // xt = xc + (xm-xc)*rt/rc;
            double rt = rc - (1-t/dx)*dxR*fabs(vecInProduct(n, v));      // Distance between xt and xc
            xt[0] = xc[0] + ( xm[0]-xc[0] )*rt/rc;
            xt[1] = xc[1] + ( xm[1]-xc[1] )*rt/rc;
            xt[2] = x[2]+t/dx*dz*v[2] + (1-t/dx)*dzMY*v[2];
        }
    }
    return xt;
}

vector<double> elasticCollisionIAS (const vector<double> &x, const vector<double> &v, const double &dx, vector<double> xc, const double &rc, const double &g, const double &dz) {
    // Elastic collision of a particle x in IAS onto an circle (xc,rc*g) with a velocity v and a diffusion length dx
    
    // Translate circle center xc to make it as close to the position x as possible
    translateXc(x,xc);
    
    // distance( x+t*v, xc )==rc*g, solve t
    double a=0,b=0,c=0,t1=0,t2=0,t=0;
    a=v[0]*v[0] + v[1]*v[1];
    b=2*(x[0]-xc[0])*v[0] + 2*(x[1]-xc[1])*v[1];
    c=pow(x[0]-xc[0],2) + pow(x[1]-xc[1],2) - rc*rc*g*g;
    
    // xt: final position, xm: contact point on myelin, n: surface unit vector, v1: reflecting velocity
    vector<double> xt(3,0), xm(2,0), n(2,0), v1(2,0);
    double discri=b*b-4*a*c;
    if (discri<0) {          // Does not hit the axon
        xt[0]=x[0]+dx*v[0];
        xt[1]=x[1]+dx*v[1];
        xt[2]=x[2]+dz*v[2];
    }
    else {
        discri=sqrt(discri);
        t1=0.5/a*( -b+discri );
        t2=0.5/a*( -b-discri );
        t=max(t1,t2);   // Choose the t>0 assuring that particle goes to the correct direction
        if (t>=dx) {            // Does not hit the axon
            xt[0]=x[0]+dx*v[0];
            xt[1]=x[1]+dx*v[1];
            xt[2]=x[2]+dz*v[2];
        }
        else {                  // Hit the axon
            // xm = x + t*v;
            xm[0]=x[0]+t*v[0];
            xm[1]=x[1]+t*v[1];
            
            // n parallel to (xm-xc_t)
            n[0]=(xm[0]-xc[0])/vecDistance(xm, xc);
            n[1]=(xm[1]-xc[1])/vecDistance(xm, xc);
            
            // v1 = v - 2(v \dot n)*n
            v1[0]=v[0]-2*vecInProduct(v,n)*n[0];
            v1[1]=v[1]-2*vecInProduct(v,n)*n[1];
            
            // xt = xm + (dx-t)*v1
            xt[0]=xm[0]+(dx-t)*v1[0];
            xt[1]=xm[1]+(dx-t)*v1[1];
            xt[2]=x[2]+dz*v[2];
        }
    }
    return xt;
}

vector<double> inelasticCollisionIAS (const vector<double> &x, const vector<double> &v, const double &dx, vector<double> xc, const double &rc, const double &g, const double &dxR, const double &dz, const double &dzMY) {
    // Inelastic collision of a particle x in IAS onto an circle (xc,rc*g) with a velocity v and a diffusion length dx
    
    // Translate circle center xc to make it as close to the position x as possible
    translateXc(x,xc);
    
    // distance( x+t*v, xc )==rc*g, solve t
    double a=0,b=0,c=0,t1=0,t2=0,t=0;
    a=v[0]*v[0] + v[1]*v[1];
    b=2*(x[0]-xc[0])*v[0] + 2*(x[1]-xc[1])*v[1];
    c=pow(x[0]-xc[0],2) + pow(x[1]-xc[1],2) - rc*rc*g*g;
    
    // xt: final position, xm: contact point on myelin, discri: discriminant
    vector<double> xt(3,0), xm(2,0), n(2,0);
    double discri=b*b-4*a*c;
    if (discri<0) {          // Does not hit the axon
        xt[0]=x[0]+dx*v[0];
        xt[1]=x[1]+dx*v[1];
        xt[2]=x[2]+dz*v[2];
    }
    else {
        discri=sqrt(discri);
        t1=0.5/a*( -b+discri );
        t2=0.5/a*( -b-discri );
        t=max(t1,t2);   // Choose the t>0 assuring that particle goes to the correct direction
        if (t>=dx) {            // Does not hit the axon
            xt[0]=x[0]+dx*v[0];
            xt[1]=x[1]+dx*v[1];
            xt[2]=x[2]+dz*v[2];
        }
        else {                  // Hit the axon
            // xm = x + t*v;
            xm[0]=x[0]+t*v[0];
            xm[1]=x[1]+t*v[1];
            
            // n parallel to (xm-xc_t)
            n[0]=(xm[0]-xc[0])/vecDistance(xm, xc);
            n[1]=(xm[1]-xc[1])/vecDistance(xm, xc);
            
            // xt = xc + (xm-xc)*rt/rc/g;
            double rt = rc*g + (1-t/dx)*dxR*fabs(vecInProduct(n, v));      // Distance between xt and xc
            xt[0] = xc[0] + ( xm[0]-xc[0] )*rt/(rc*g);
            xt[1] = xc[1] + ( xm[1]-xc[1] )*rt/(rc*g);
            xt[2] = x[2]+t/dx*dz*v[2] + (1-t/dx)*dzMY*v[2];
        }
    }
    return xt;
}

vector<double> diffuseMyelinEllipseCircum (const vector<double> &x, const double &dxC, const double &phi, vector<double> xc) {
    // Diffuse with an compressed ellipse propagator, particle at x, step size dxC and dxR in circumferential and radial, parameter phi, circle center xc, random number randR to determine the radial direction
    
    // Translate circle center xc to make it as close to the position x as possible
    translateXc(x,xc);
    
    double r0=vecDistance(x, xc), theta=dxC*cos(phi)/r0, cos_theta=cos(theta), sin_theta=sin(theta);
    
    vector<double> xTmp(3,0);   // xt: final position
    vecSubtract(x, xc, xTmp);          // xTmp = x-xc_t
    // Diffuse circumferentially, rotation around xc_t with theta
    xc[0]+=xTmp[0]*cos_theta-xTmp[1]*sin_theta;
    xc[1]+=xTmp[0]*sin_theta+xTmp[1]*cos_theta;
    
    xTmp[0]=xc[0]; xTmp[1]=xc[1]; xTmp[2]=x[2];
    
    return xTmp;
}

vector<double> diffuseMyelinEllipseRadial (const vector<double> &x, const double &dxR, const double &phi, vector<double> xc, const double &randR, const double &dzMY) {
    // Diffuse with an compressed ellipse propagator, particle at x, step size dxC and dxR in circumferential and radial, parameter phi, circle center xc, random number randR to determine the radial direction
    
    // Translate circle center xc to make it as close to the position x as possible
    translateXc(x,xc);
    
    double r0=vecDistance(x, xc), dr=dxR*fabs(sin(phi));
    
    if ( randR<(0.5-0.25*dr/r0) ) {
        dr=-dr;
    }
    
    vector<double> xTmp(3,0);   // xt: final position
    
    vecSubtract(x, xc, xTmp);         // xTmp = xt-xc_t
    double tmp=vecNorm(xTmp);
    // Normalize xTmp: unit vector along xt-xc_t
    // Diffuse radially with a step size dr
    xTmp[0]=x[0]+xTmp[0]/tmp*dr;
    xTmp[1]=x[1]+xTmp[1]/tmp*dr;
    xTmp[2]=x[2]+dzMY;
    
    return xTmp;
}

vector<double> elasticCollisionMyelin2IAS (const vector<double> &x, const double &dr, vector<double> xc, const double &rc, const double &g, const double &dzMY) {
    // Elastic collision of a particle at position x in myelin, hitting the inner sheath (xc,rc*g), with a step size dr
    
    // Translate circle center xc to make it as close to the position x as possible
    //vector<double> xc_t(2,0);
    translateXc(x,xc);
    
    //vector<double> xt(2,0);         // Final position
    double rt=0, r_in=rc*g;         // Distance between xt and xc, r_in = inner radius
    double r0=vecDistance(x, xc);   // Distance between x and xc
    vector<double> xt(3,0);
    if ( dr<(r0-r_in) ) {      // If no hitting to the inner sheath
        rt = r0-dr;
    }
    else {                      // Hit to the inner sheath
        rt = r_in + dr - ( r0 - r_in );
    }
    xt[0] = xc[0] + ( x[0]-xc[0] )*rt/r0;
    xt[1] = xc[1] + ( x[1]-xc[1] )*rt/r0;
    xt[2] = x[2] + dzMY;
    
    return xt;
}

vector<double> permeateMyelin2IAS (const vector<double> &x, const double &dr, vector<double> xc, const double &rc, const double &g, const double &dxIAS, const double &dz, const double &dzMY) {
    // A particle at position x in myelin, permeating through the inner sheath (xc,rc*g), with a step size dr, and dxIAS is the step size in IAS
    
    // Translate circle center xc to make it as close to the position x as possible
    translateXc(x,xc);
    
    double rt=0, r_in=rc*g;         // Distance between xt and xc, r_in = inner radius
    double r0=vecDistance(x, xc);   // Distance between x and xc
    vector<double> xt(3,0);
    if ( dr<(r0-r_in) ) {      // If no hitting to the inner sheath
        rt = r0-dr;
        xt[2] = x[2] + dzMY;
    }
    else {                      // Hit to the inner sheath
        rt = r_in - ( 1 - (r0-r_in)/dr )*dxIAS;
        xt[2] = x[2] + (r0-r_in)/dr*dzMY + ( 1 - (r0-r_in)/dr )*dz;
    }
    xt[0] = xc[0] + ( x[0]-xc[0] )*rt/r0;
    xt[1] = xc[1] + ( x[1]-xc[1] )*rt/r0;
    
    return xt;
}

vector<double> elasticCollisionMyelin2EAS (const vector<double> &x, const double &dr, vector<double> xc, const double &rc, const double &dzMY) {
    // Elastic collision of a particle at position x in myelin, hitting the outer sheath (xc,rc), with a step size dr
    
    // Translate circle center xc to make it as close to the position x as possible
    translateXc(x,xc);
    
    //vector<double> xt(2,0);         // Final position
    double rt=0;                    // Distance between xt and xc
    double r0=vecDistance(x, xc);   // Distance between x and xc
    vector<double> xt(3,0);
    if ( dr<(rc-r0) ) {        // If no hitting to the outer sheath
        rt = r0+dr;
    }
    else {                      // Hit to the outer sheath
        rt = rc - ( dr - ( rc - r0 ) );
    }
    xt[0] = xc[0] + ( x[0]-xc[0] )*rt/r0;
    xt[1] = xc[1] + ( x[1]-xc[1] )*rt/r0;
    xt[2] = x[2] + dzMY;
    
    return xt;
}

vector<double> permeateMyelin2EAS (const vector<double> &x, const double &dr, vector<double> xc, const double &rc, const double &dxEAS, const double &dz, const double &dzMY) {
    // A particle at position x in myelin, permeating through the outer sheath (xc,rc), with a step size dr, and dxEAS is the step size in EAS
    
    // Translate circle center xc to make it as close to the position x as possible
    translateXc(x,xc);
    
    double rt=0;         // Distance between xt and xc
    double r0=vecDistance(x, xc);   // Distance between x and xc
    vector<double> xt(3,0);
    if ( dr<(rc-r0) ) {      // If no hitting to the inner sheath
        rt = r0+dr;
        xt[2] = x[2] + dzMY;
    }
    else {                      // Hit to the inner sheath
        rt = rc + ( 1 - (rc-r0)/dr )*dxEAS;
        xt[2] = x[2] + (rc-r0)/dr*dzMY + ( 1 - (rc-r0)/dr )*dz;
    }
    xt[0] = xc[0] + ( x[0]-xc[0] )*rt/r0;
    xt[1] = xc[1] + ( x[1]-xc[1] )*rt/r0;
    
    return xt;
}



























