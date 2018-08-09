//
//  diffusion_lib.cpp
//  practice
//
//  Created by magda on 2/16/17.
//  Copyright (c) 2017 Honghsi. All rights reserved.
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

/*
vector<double> vecSum (vector<double> x, vector<double> y) {
    // Summation of two vectors
    
    //vector<double> sum(x.size(),0);
    //for (int i=0; i<x.size(); i++) {
    //    sum[i]=x[i]+y[i];
    //}
    //return sum;
    
    x[0]+=y[0];
    x[1]+=y[1];
    return x;
}
 */

void vecSubtract (const vector<double> &x, const vector<double> &y, vector<double> &z) {
    // Subtraction of two vectors x-y = z
    /*
    vector<double> sub(x.size(),0);
    for (int i=0; i<x.size(); i++) {
        sub[i]=x[i]-y[i];
    }
    return sub;
     */
    z[0]=x[0]-y[0];
    z[1]=x[1]-y[1];
    //x[0]-=y[0];
    //x[1]-=y[1];
    //return x;
}

inline double vecDistance (const vector<double> &x, const vector<double> &y) {
    /*
    double distance=0;
    for (int i=0; i<x.size(); i++) {
        distance+=pow(x[i]-y[i],2);
    }
    distance = sqrt(distance);
    return distance;
     */
    //return sqrt( (x[0]-y[0])*(x[0]-y[0]) + (x[1]-y[1])*(x[1]-y[1]) );
    return sqrt( pow(x[0]-y[0], 2) + pow(x[1]-y[1], 2) );
}

inline double vecInProduct (const vector<double> &x, const vector<double> &y) {
    // Inner product of two vector
    /*
    double inner_product=0;
    for (int i=0; i<x.size(); i++) {
        inner_product+=x[i]*y[i];
    }
    return inner_product;
     */
    return ( x[0]*y[0] + x[1]*y[1] );
}

inline double vecNorm (const vector<double> &x) {
    /*
    double norm=0;
    for (int i=0; i<x.size(); i++) {
        norm+=x[i]*x[i];
    }
    norm=sqrt(norm);
    return norm;
     */
    return sqrt( x[0]*x[0] + x[1]*x[1] );
}

vector<int> pixPosition ( vector<double> x, unsigned int NPix ) {
    vector<int> xPix(2,0);
    
    if ( x[0]<0 ) { x[0]+=1; }
    if ( x[0]>1 ) { x[0]-=1; }
    
    if ( x[1]<0 ) { x[1]+=1; }
    if ( x[1]>1 ) { x[1]-=1; }
    
    xPix[0]=floor(x[0]*NPix); //xPix[0]%=NPix;
    xPix[1]=floor(x[1]*NPix); //xPix[1]%=NPix;
    return xPix;
}

/*
inline bool needTranslateXc (vector<double> xc, double rc) {
    return ( (xc[0]+rc) > 1 ) | ( xc[0] < rc ) | ( (xc[1]+rc) > 1 ) | ( xc[1] < rc );
}
 */

void translateXc ( const vector<double> &x, vector<double> &xc) {
    // Translate circle center xc to make it as close to the position x as possible
    //vector<double> xc_t(2,0);
    int ti=0, tj=0;
    //double dis=vecDistance(x, xc), dis_t=dis;
    double disSquare = (x[0]-xc[0])*(x[0]-xc[0])+(x[1]-xc[1])*(x[1]-xc[1]), disSquare_t=0;
    for (int i=-1; i<2; i++) {
        for (int j=-1; j<2; j++) {
            //xc_t[0]=xc[0]+i;
            //xc_t[1]=xc[1]+j;
            disSquare_t=(x[0]-xc[0]-i)*(x[0]-xc[0]-i)+(x[1]-xc[1]-j)*(x[1]-xc[1]-j);
            //disSquare_t=(x[0]-xc_t[0])*(x[0]-xc_t[0])+(x[1]-xc_t[1])*(x[1]-xc_t[1]);
            if (disSquare_t<disSquare) {
                disSquare=disSquare_t;
                //xc_t[0]=xc[0]+i;
                //xc_t[1]=xc[1]+j;
                ti=i;
                tj=j;
            }
        }
    }
    xc[0]+=ti;
    xc[1]+=tj;
    
    //return xc;
}

bool inAxon ( const vector<double> &x, vector<double> xc, const double &rc ) {
    // If the point x is in the circle (xc,rc), return 1; if not, return 0.
    
    // Translate circle center xc to make it as close to the position xt as possible
    //vector<double> xc_t(2,0);
    translateXc(x,xc);
    
    //if ( sqrt( pow(x[0]-xc_t[0], 2) + pow(x[1]-xc_t[1], 2) )<=rc) {
    return ( ( (x[0]-xc[0])*(x[0]-xc[0]) + (x[1]-xc[1])*(x[1]-xc[1]) ) <= rc*rc );
    //return ( vecDistance(x, xc) <= rc );
    
}

bool inMyelin ( const vector<double> &x, vector<double> xc, const double &rc, const double &g ) {
    // If the point x is in the myelin (xc,rc,rc*g), return 1; if not, return 0.
    
    // Translate circle center xc to make it as close to the position xt as possible
    //vector<double> xc_t(2,0);
    translateXc(x,xc);
    
    //if ( (vecDistance(x, xc_t)<=rc) & (vecDistance(x, xc_t)>=rc*g) ) {
    double disSquare = (x[0]-xc[0])*(x[0]-xc[0]) + (x[1]-xc[1])*(x[1]-xc[1]);
    return ( ( disSquare <= rc*rc ) & ( disSquare >=rc*rc*g*g ) );
    //double dis=vecDistance(x, xc);
    //return ( ( dis <= rc) & ( dis >= rc*g ) );
    
}

bool inIAS ( const vector<double> &x, vector<double> xc, const double &rc, const double &g ) {
    // If the point x is in the Intra-axonal space (xc,rc*g), return 1; if not, return 0.
    
    // Translate circle center xc to make it as close to the position xt as possible
    //vector<double> xc_t(2,0);
    translateXc(x,xc);
    return ( ( (x[0]-xc[0])*(x[0]-xc[0]) + (x[1]-xc[1])*(x[1]-xc[1]) ) < rc*rc*g*g ) ;
    
    //return ( vecDistance(x, xc)<rc*g );
    
}

bool stepEAS2Axon (const vector<double> &xi, const vector<double> &xt, vector<double> xc, const double &rc ) {
    // If segment(xi,xt) overlaps with circle (xc,rc), return 1; if not, return 0.
    
    // Translate circle center xc to make it as close to the position xt as possible
    //vector<double> xc_t(2,0);
    translateXc(xt,xc);
    
    // If xt is in axon, segment overlaps with circle
    //if ( inAxon(xt,xc,rc)) {
    //if ( vecDistance(xt, xc) <= rc ) {

    if ( ( (xt[0]-xc[0])*(xt[0]-xc[0]) + (xt[1]-xc[1])*(xt[1]-xc[1]) ) <= rc*rc ) {
        return true;
    }
    
    // L: a*x[0]+b*x[1]+c=0, d=distance of xc to L
    double a=xi[1]-xt[1], b=xt[0]-xi[0], c=xi[0]*xt[1]-xi[1]*xt[0], sos_ab= a*a+b*b, d=fabs( a*xc[0]+b*xc[1]+c )/sqrt( sos_ab );
    //cout<<"distance to L ="<<d<<endl;
    //cout<<"a="<<a<<",b="<<b<<",c="<<c<<endl;
    if (d>rc) {
        return false;
    }
    
    // xl: a point on L closest to xc
    vector<double> xl(2,0);
    xl[0] = ( b*( b*xc[0]-a*xc[1] )-a*c )/( sos_ab );
    xl[1] = ( a*( -b*xc[0]+a*xc[1] )-b*c )/( sos_ab );
    //cout<<"closest point xl="<<xl[0]<<','<<xl[1]<<endl;
    
    // xl is in axon, but xi and xt are not in axon.
    return ( ( (xi[0]-xl[0])*(xt[0]-xl[0])+(xi[1]-xl[1])*(xt[1]-xl[1]) ) <= 0 );
    
    //vector<double> r1(2,0), r2(2,0);
    //r1=vecSubtract(xi, xl); r2=vecSubtract(xt, xl);
    //return ( vecInProduct(r1, r2)<=0 );
    
}

/*
bool stepIAS2nonIAS (vector<double> xt, vector<double> xc, double rc, double g ) {
    // If segment(xi,xt) overlaps with non-IAS = !circle (xc,rc*g), return 1; if not, return 0.
    
 
    // Translate circle center xc to make it as close to the position xt as possible
    //vector<double> xc_t(2,0);
    //xc_t=translateXc(xt,xc);
 
    
    // If xt is in non-IAS, segment overlaps with non-IAS
    return !inIAS(xt,xc,rc,g);
    
}


bool stepMyelin2EAS (vector<double> xt, vector<double> xc, double rc)  {
    // If segment(xi,xt) overlaps with EAS = !circle (xc,rc), return 1; if not, return 0.
    
 
    // Translate circle center xc to make it as close to the position xt as possible
    //vector<double> xc_t(2,0);
    //xc_t=translateXc(xt,xc);
 
    
    // If xt is in EAS, segment overlaps with EAS
    return !inAxon(xt,xc,rc);
    
}

bool stepMyelin2IAS (vector<double> xt, vector<double> xc, double rc, double g )  {
    // If segment(xi,xt) overlaps with IAS = circle (xc,rc), return 1; if not, return 0.
    
 
    // Translate circle center xc to make it as close to the position xt as possible
    //vector<double> xc_t(2,0);
    //xc_t=translateXc(xt,xc);
 
    
    // If xt is in IAS, segment overlaps with IAS
    return inIAS(xt,xc,rc,g);
    
}
 */

double disOutSheath ( const vector<double> &x, vector<double> xc, const double &rc ) {
    // The distance between the point x and surface of the circle (xc,rc)
    
    // Translate circle center xc to make it as close to the position x as possible
    //vector<double> xc_t(2,0);
    translateXc(x,xc);
    
    //return fabs( sqrt( (x[0]-xc[0])*(x[0]-xc[0]) + (x[1]-xc[1])*(x[1]-xc[1]) ) - rc );
    return fabs( sqrt( pow(x[0]-xc[0],2) + pow(x[1]-xc[1],2) ) - rc );
}

double disInSheath ( const vector<double> &x, vector<double> xc, const double &rc, const double &g ) {
    // The distance between the point x and surface of the circle (xc,rc*g)
    
    // Translate circle center xc to make it as close to the position x as possible
    //vector<double> xc_t(2,0);
    translateXc(x,xc);
    
    return fabs( sqrt( pow(x[0]-xc[0],2) + pow(x[1]-xc[1],2) ) - rc*g );
}

/*
double disCenter ( vector<double> x, vector<double> xc) {
    // The distance between the point x and the circle center xc
    
    // Translate circle center xc to make it as close to the position x as possible
    //vector<double> xc_t(2,0);
    translateXc(x,xc);
    
    return sqrt( pow(x[0]-xc[0],2) + pow(x[1]-xc[1],2) ) ;
}
 */

vector<double> elasticCollisionEAS (const vector<double> &x, const vector<double> &v, const double &dx, vector<double> xc, const double &rc) {
    // Elastic collision of a particle x in EAS onto an circle (xc,rc) with a velocity v and a diffusion length dx
    
    // Translate circle center xc to make it as close to the position (x + dx*v) as possible
    vector<double> xTmp(2,0);
    xTmp[0]=x[0]+dx*v[0];
    xTmp[1]=x[1]+dx*v[1];
    translateXc(xTmp,xc);
    
    // distance( x+t*v, xc )==rc, solve t
    double a=0,b=0,c=0,t1=0,t2=0,t=0;
    a=v[0]*v[0] + v[1]*v[1];
    b=2*(x[0]-xc[0])*v[0] + 2*(x[1]-xc[1])*v[1];
    c=pow(x[0]-xc[0],2) + pow(x[1]-xc[1],2) - rc*rc;
    //cout<<"a,b,c\t"<<a<<'\t'<<b<<'\t'<<c<<endl;
    
    // xt: final position, xm: contact point on myelin, n: surface unit vector, v1: reflecting velocity, discri: discriminant
    vector<double> xt(2,0), xm(2,0), n(2,0), v1(2,0);
    double discri=b*b-4*a*c;
    if (discri<0) {          // Does not hit the axon
        //xt[0]=x[0]+dx*v[0];
        //xt[1]=x[1]+dx*v[1];
        xt=xTmp;
        //cout<<"case 1, do not hit,\t"<<xt[0]<<'\t'<<xt[1]<<endl;
    }
    else {
        discri=sqrt(discri);
        t1=0.5/a*( -b+discri );
        t2=0.5/a*( -b-discri );
        t=min(t1,t2);
        //cout<<"t\t"<<t1<<'\t'<<t2<<'\t'<<t<<'\t'<<endl;
        if ( (t>=dx) | (t<0) ) {            // Does not hit the axon
            //xt[0]=x[0]+dx*v[0];
            //xt[1]=x[1]+dx*v[1];
            xt=xTmp;
            //cout<<"case 2, do not hit,\t"<<xt[0]<<'\t'<<xt[1]<<endl;
        }
        else {                  // Hit the axon
            //cout<<"case 3, hit"<<endl;
            // xm = x + t*v;
            xm[0]=x[0]+t*v[0];
            xm[1]=x[1]+t*v[1];
            //cout<<"xm\t"<<xm[0]<<'\t'<<xm[1]<<endl;
            
            // n parallel to (xm-xc_t)
            n[0]=(xm[0]-xc[0])/vecDistance(xm, xc);
            n[1]=(xm[1]-xc[1])/vecDistance(xm, xc);
            //cout<<"n\t"<<n[0]<<'\t'<<n[1]<<endl;
            
            // v1 = v - 2(v \dot n)*n
            v1[0]=v[0]-2*vecInProduct(v,n)*n[0];
            v1[1]=v[1]-2*vecInProduct(v,n)*n[1];
            //cout<<"v1\t"<<v1[0]<<'\t'<<v1[1]<<endl;
            
            // xt = xm + (dx-t)*v1
            xt[0]=xm[0]+(dx-t)*v1[0];
            xt[1]=xm[1]+(dx-t)*v1[1];
            //cout<<"xt\t"<<xt[0]<<'\t'<<xt[1]<<endl;
        }
    }
    return xt;
}

vector<double> inelasticCollisionEAS (const vector<double>  &x, const vector<double> &v, const double &dx, vector<double> xc, const double &rc, const double &dxR) {
    // Inelastic collision of a particle x in EAS onto an circle (xc,rc) with a velocity v and a diffusion length dx. The dxR is the step size in myelin
    
    // Translate circle center xc to make it as close to the position (x + dx*v) as possible
    vector<double> xTmp(2,0);
    xTmp[0]=x[0]+dx*v[0];
    xTmp[1]=x[1]+dx*v[1];
    translateXc(xTmp,xc);
    
    // distance( x+t*v, xc )==rc, solve t
    double a=0,b=0,c=0,t1=0,t2=0,t=0;
    a=v[0]*v[0] + v[1]*v[1];
    b=2*(x[0]-xc[0])*v[0] + 2*(x[1]-xc[1])*v[1];
    c=pow(x[0]-xc[0],2) + pow(x[1]-xc[1],2) - rc*rc;
    //cout<<"a,b,c\t"<<a<<'\t'<<b<<'\t'<<c<<endl;
    
    // xt: final position, xm: contact point on myelin, discri: discriminant
    vector<double> xt(2,0), xm(2,0), n(2,0);
    double discri=b*b-4*a*c;
    if (discri<0) {          // Does not hit the axon
        //xt[0]=x[0]+dx*v[0];
        //xt[1]=x[1]+dx*v[1];
        xt=xTmp;
        //cout<<"case 1, do not hit,\t"<<xt[0]<<'\t'<<xt[1]<<endl;
    }
    else {
        discri=sqrt(discri);
        t1=0.5/a*( -b+discri );
        t2=0.5/a*( -b-discri );
        t=min(t1,t2);
        //cout<<"t\t"<<t1<<'\t'<<t2<<'\t'<<t<<'\t'<<endl;
        if ( (t>=dx) | (t<0) ) {            // Does not hit the axon
            //xt[0]=x[0]+dx*v[0];
            //xt[1]=x[1]+dx*v[1];
            xt=xTmp;
            //cout<<"case 2, do not hit,\t"<<xt[0]<<'\t'<<xt[1]<<endl;
        }
        else {                  // Hit the axon
            //cout<<"case 3, hit"<<endl;
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
            //xt=xm;
        }
    }
    return xt;
}

vector<double> elasticCollisionIAS (const vector<double> &x, const vector<double> &v, const double &dx, vector<double> xc, const double &rc, const double &g) {
    // Elastic collision of a particle x in IAS onto an circle (xc,rc*g) with a velocity v and a diffusion length dx
    
    // Translate circle center xc to make it as close to the position x as possible
    //vector<double> xc_t(2,0);
    translateXc(x,xc);
    
    // distance( x+t*v, xc )==rc*g, solve t
    double a=0,b=0,c=0,t1=0,t2=0,t=0;
    a=v[0]*v[0] + v[1]*v[1];
    b=2*(x[0]-xc[0])*v[0] + 2*(x[1]-xc[1])*v[1];
    c=pow(x[0]-xc[0],2) + pow(x[1]-xc[1],2) - rc*rc*g*g;
    //cout<<"a,b,c\t"<<a<<'\t'<<b<<'\t'<<c<<endl;
    
    // xt: final position, xm: contact point on myelin, n: surface unit vector, v1: reflecting velocity
    vector<double> xt(2,0), xm(2,0), n(2,0), v1(2,0);
    double discri=b*b-4*a*c;
    if (discri<0) {          // Does not hit the axon
        xt[0]=x[0]+dx*v[0];
        xt[1]=x[1]+dx*v[1];
        //cout<<"case 1, do not hit,\t"<<xt[0]<<'\t'<<xt[1]<<endl;
    }
    else {
        discri=sqrt(discri);
        t1=0.5/a*( -b+discri );
        t2=0.5/a*( -b-discri );
        t=max(t1,t2);   // Choose the t>0 assuring that particle goes to the correct direction
        //cout<<"t\t"<<t1<<'\t'<<t2<<'\t'<<t<<'\t'<<endl;
        if (t>=dx) {            // Does not hit the axon
            xt[0]=x[0]+dx*v[0];
            xt[1]=x[1]+dx*v[1];
            //cout<<"case 2, do not hit,\t"<<xt[0]<<'\t'<<xt[1]<<endl;
        }
        else {                  // Hit the axon
            //cout<<"case 3, hit"<<endl;
            // xm = x + t*v;
            xm[0]=x[0]+t*v[0];
            xm[1]=x[1]+t*v[1];
            //cout<<"xm\t"<<xm[0]<<'\t'<<xm[1]<<endl;
            
            // n parallel to (xm-xc_t)
            n[0]=(xm[0]-xc[0])/vecDistance(xm, xc);
            n[1]=(xm[1]-xc[1])/vecDistance(xm, xc);
            //cout<<"n\t"<<n[0]<<'\t'<<n[1]<<endl;
            
            // v1 = v - 2(v \dot n)*n
            v1[0]=v[0]-2*vecInProduct(v,n)*n[0];
            v1[1]=v[1]-2*vecInProduct(v,n)*n[1];
            //cout<<"v1\t"<<v1[0]<<'\t'<<v1[1]<<endl;
            
            // xt = xm + (dx-t)*v1
            xt[0]=xm[0]+(dx-t)*v1[0];
            xt[1]=xm[1]+(dx-t)*v1[1];
            //cout<<"xt\t"<<xt[0]<<'\t'<<xt[1]<<endl;
        }
    }
    return xt;
}

vector<double> inelasticCollisionIAS (const vector<double> &x, const vector<double> &v, const double &dx, vector<double> xc, const double &rc, const double &g, const double &dxR) {
    // Inelastic collision of a particle x in IAS onto an circle (xc,rc*g) with a velocity v and a diffusion length dx
    
    // Translate circle center xc to make it as close to the position x as possible
    //vector<double> xc_t(2,0);
    translateXc(x,xc);
    
    // distance( x+t*v, xc )==rc*g, solve t
    double a=0,b=0,c=0,t1=0,t2=0,t=0;
    a=v[0]*v[0] + v[1]*v[1];
    b=2*(x[0]-xc[0])*v[0] + 2*(x[1]-xc[1])*v[1];
    c=pow(x[0]-xc[0],2) + pow(x[1]-xc[1],2) - rc*rc*g*g;
    //cout<<"a,b,c\t"<<a<<'\t'<<b<<'\t'<<c<<endl;
    
    // xt: final position, xm: contact point on myelin, discri: discriminant
    vector<double> xt(2,0), xm(2,0), n(2,0);
    double discri=b*b-4*a*c;
    if (discri<0) {          // Does not hit the axon
        xt[0]=x[0]+dx*v[0];
        xt[1]=x[1]+dx*v[1];
        //cout<<"case 1, do not hit,\t"<<xt[0]<<'\t'<<xt[1]<<endl;
    }
    else {
        discri=sqrt(discri);
        t1=0.5/a*( -b+discri );
        t2=0.5/a*( -b-discri );
        t=max(t1,t2);   // Choose the t>0 assuring that particle goes to the correct direction
        //cout<<"t\t"<<t1<<'\t'<<t2<<'\t'<<t<<'\t'<<endl;
        if (t>=dx) {            // Does not hit the axon
            xt[0]=x[0]+dx*v[0];
            xt[1]=x[1]+dx*v[1];
            //cout<<"case 2, do not hit,\t"<<xt[0]<<'\t'<<xt[1]<<endl;
        }
        else {                  // Hit the axon
            //cout<<"case 3, hit"<<endl;
            // xm = x + t*v;
            xm[0]=x[0]+t*v[0];
            xm[1]=x[1]+t*v[1];
            //cout<<"xm\t"<<xm[0]<<'\t'<<xm[1]<<endl;
            
            // n parallel to (xm-xc_t)
            n[0]=(xm[0]-xc[0])/vecDistance(xm, xc);
            n[1]=(xm[1]-xc[1])/vecDistance(xm, xc);
            
            // xt = xc + (xm-xc)*rt/rc/g;
            double rt = rc*g + (1-t/dx)*dxR*fabs(vecInProduct(n, v));      // Distance between xt and xc
            xt[0] = xc[0] + ( xm[0]-xc[0] )*rt/(rc*g);
            xt[1] = xc[1] + ( xm[1]-xc[1] )*rt/(rc*g);
        }
    }
    return xt;
}
/*
vector<double> diffuseCircum (vector<double> x, double vC, double dxC, vector<double> xc) {
    // Diffuse circumferentially around a center xc from x with a step dxC and a direction vC (1: counterclockwise, -1: clockwise)
    
    // Translate circle center xc to make it as close to the position x as possible
    //vector<double> xc_t(2,0);
    translateXc(x,xc);
    
    vector<double> xt(2,0),xTmp(2,0);   // xt: final position, xTmp = x-xc_t;
    double phi=0;                       // Azimuthal angle
    vecSubtract(x, xc, xTmp);          // Translation
    phi = dxC/vecNorm(xTmp)*vC;         // phi = s/r * (\pm 1)
    double cos_phi=cos(phi), sin_phi=sin(phi);
    // Rotation around xc_t with phi
    xt[0]=xc[0]+xTmp[0]*cos_phi-xTmp[1]*sin_phi;
    xt[1]=xc[1]+xTmp[0]*sin_phi+xTmp[1]*cos_phi;
    
    return xt;
}

vector<double> diffuseRadial (vector<double> x, double vR, double dxR, vector<double> xc) {
    // Diffuse away from or toward a center xc from x with a step dxR and a direction vR (1: away from the center, -1: toward the center)
    
    // Translate circle center xc to make it as close to the position x as possible
    //vector<double> xc_t(2,0);
    translateXc(x,xc);
    
    vector<double> xt(2,0),n(2,0);   // xt: final position, n: unit vector parallel to x-xc_t
    double tmp=0;
    vecSubtract(x, xc, n); tmp=vecNorm(n);
    n[0]/=tmp;
    n[1]/=tmp;
    xt[0]=x[0]+n[0]*dxR*vR;
    xt[1]=x[1]+n[1]*dxR*vR;
    
    return xt;
}
 */

vector<double> diffuseMyelinEllipseCircum (const vector<double> &x, const double &dxC, const double &phi, vector<double> xc) {
    // Diffuse with an compressed ellipse propagator, particle at x, step size dxC and dxR in circumferential and radial, parameter phi, circle center xc, random number randR to determine the radial direction
    
    // Translate circle center xc to make it as close to the position x as possible
    //vector<double> xc_t(2,0);
    translateXc(x,xc);
    
    double r0=vecDistance(x, xc), theta=dxC*cos(phi)/r0, cos_theta=cos(theta), sin_theta=sin(theta);
    /*
    if ( randR<(0.5-0.25*dr/r0) ) {
        dr=-dr;
    }
     */
    
    vector<double> xTmp(2,0);   // xt: final position
    vecSubtract(x, xc, xTmp);          // xTmp = x-xc_t
    // Diffuse circumferentially, rotation around xc_t with theta
    xc[0]+=xTmp[0]*cos_theta-xTmp[1]*sin_theta;
    xc[1]+=xTmp[0]*sin_theta+xTmp[1]*cos_theta;
    
    /*
    vecSubtract(xt, xc, xTmp);         // xTmp = xt-xc_t
    double tmp=vecNorm(xTmp);
    // Normalize xTmp: unit vector along xt-xc_t
    //xTmp[0]/=tmp;
    //xTmp[1]/=tmp;
    // Diffuse radially with a step size dr
    xt[0]+=xTmp[0]/tmp*dr;
    xt[1]+=xTmp[1]/tmp*dr;
     */
    
    return xc;
}

vector<double> diffuseMyelinEllipseRadial (const vector<double> &x, const double &dxR, const double &phi, vector<double> xc, const double &randR) {
    // Diffuse with an compressed ellipse propagator, particle at x, step size dxC and dxR in circumferential and radial, parameter phi, circle center xc, random number randR to determine the radial direction
    
    // Translate circle center xc to make it as close to the position x as possible
    //vector<double> xc_t(2,0);
    translateXc(x,xc);
    
    double r0=vecDistance(x, xc), dr=dxR*fabs(sin(phi));
    
    if ( randR<(0.5-0.25*dr/r0) ) {
        dr=-dr;
    }
    
    vector<double> xTmp(2,0);   // xt: final position
    
    vecSubtract(x, xc, xTmp);         // xTmp = xt-xc_t
    double tmp=vecNorm(xTmp);
    // Normalize xTmp: unit vector along xt-xc_t
    //xTmp[0]/=tmp;
    //xTmp[1]/=tmp;
    // Diffuse radially with a step size dr
    xc[0]=x[0]+xTmp[0]/tmp*dr;
    xc[1]=x[1]+xTmp[1]/tmp*dr;
    
    return xc;
}

vector<double> elasticCollisionMyelin2IAS (const vector<double> &x, const double &dr, vector<double> xc, const double &rc, const double &g) {
    // Elastic collision of a particle at position x in myelin, hitting the inner sheath (xc,rc*g), with a step size dr
    
    // Translate circle center xc to make it as close to the position x as possible
    //vector<double> xc_t(2,0);
    translateXc(x,xc);
    
    //vector<double> xt(2,0);         // Final position
    double rt=0, r_in=rc*g;         // Distance between xt and xc, r_in = inner radius
    double r0=vecDistance(x, xc);   // Distance between x and xc
    if ( dr<(r0-r_in) ) {      // If no hitting to the inner sheath
        rt = r0-dr;
    }
    else {                      // Hit to the inner sheath
        rt = r_in + dr - ( r0 - r_in );
    }
    xc[0] = xc[0] + ( x[0]-xc[0] )*rt/r0;
    xc[1] = xc[1] + ( x[1]-xc[1] )*rt/r0;
    
    return xc;
    
    //xt[0] = xc[0] + ( x[0]-xc[0] )*rt/r0;
    //xt[1] = xc[1] + ( x[1]-xc[1] )*rt/r0;
    
    //return xt;
}

vector<double> permeateMyelin2IAS (const vector<double> &x, const double &dr, vector<double> xc, const double &rc, const double &g, const double &dxIAS) {
    // A particle at position x in myelin, permeating through the inner sheath (xc,rc*g), with a step size dr, and dxIAS is the step size in IAS
    
    // Translate circle center xc to make it as close to the position x as possible
    translateXc(x,xc);
    
    double rt=0, r_in=rc*g;         // Distance between xt and xc, r_in = inner radius
    double r0=vecDistance(x, xc);   // Distance between x and xc
    if ( dr<(r0-r_in) ) {      // If no hitting to the inner sheath
        rt = r0-dr;
    }
    else {                      // Hit to the inner sheath
        rt = r_in - ( 1 - (r0-r_in)/dr )*dxIAS;
    }
    xc[0] = xc[0] + ( x[0]-xc[0] )*rt/r0;
    xc[1] = xc[1] + ( x[1]-xc[1] )*rt/r0;
    
    return xc;
}

vector<double> elasticCollisionMyelin2EAS (const vector<double> &x, const double &dr, vector<double> xc, const double &rc) {
    // Elastic collision of a particle at position x in myelin, hitting the outer sheath (xc,rc), with a step size dr
    
    // Translate circle center xc to make it as close to the position x as possible
    //vector<double> xc_t(2,0);
    translateXc(x,xc);
    
    //vector<double> xt(2,0);         // Final position
    double rt=0;                    // Distance between xt and xc
    double r0=vecDistance(x, xc);   // Distance between x and xc
    if ( dr<(rc-r0) ) {        // If no hitting to the outer sheath
        rt = r0+dr;
    }
    else {                      // Hit to the outer sheath
        rt = rc - ( dr - ( rc - r0 ) );
    }
    xc[0] = xc[0] + ( x[0]-xc[0] )*rt/r0;
    xc[1] = xc[1] + ( x[1]-xc[1] )*rt/r0;
    
    return xc;
    
    //xt[0] = xc[0] + ( x[0]-xc[0] )*rt/r0;
    //xt[1] = xc[1] + ( x[1]-xc[1] )*rt/r0;
    
    //return xt;
}

vector<double> permeateMyelin2EAS (const vector<double> &x, const double &dr, vector<double> xc, const double &rc, const double &dxEAS) {
    // A particle at position x in myelin, permeating through the outer sheath (xc,rc), with a step size dr, and dxEAS is the step size in EAS
    
    // Translate circle center xc to make it as close to the position x as possible
    translateXc(x,xc);
    
    double rt=0;         // Distance between xt and xc
    double r0=vecDistance(x, xc);   // Distance between x and xc
    if ( dr<(rc-r0) ) {      // If no hitting to the inner sheath
        rt = r0+dr;
    }
    else {                      // Hit to the inner sheath
        rt = rc + ( 1 - (rc-r0)/dr )*dxEAS;
    }
    xc[0] = xc[0] + ( x[0]-xc[0] )*rt/r0;
    xc[1] = xc[1] + ( x[1]-xc[1] )*rt/r0;
    
    return xc;
}
/*
double probEAS2Myelin (const vector<double> &x, const vector<double> &v, const double &dx, vector<double> xc, const double &rc, const double &probUnit) {
    // Probability for a particle at x from EAS to myelin, circle (xc,rc), velocity v, step size dx, res*kappa/Dout = probUnit
    
    double prob=0;
    
    // Translate circle center xc to make it as close to the position (x + dx*v) as possible
    vector<double> xTmp(2,0);
    xTmp[0]=x[0]+dx*v[0];
    xTmp[1]=x[1]+dx*v[1];
    translateXc(xTmp,xc);
    
    // distance( x+t*v, xc )==rc, solve t
    double a=0,b=0,c=0,t1=0,t2=0,t=0;
    a=v[0]*v[0] + v[1]*v[1];
    b=2*(x[0]-xc[0])*v[0] + 2*(x[1]-xc[1])*v[1];
    c=pow(x[0]-xc[0],2) + pow(x[1]-xc[1],2) - rc*rc;
    //cout<<"a,b,c\t"<<a<<'\t'<<b<<'\t'<<c<<endl;
    
    // xt: final position, xm: contact point on myelin, n: surface unit vector, v1: reflecting velocity
    vector<double> xt(2,0), xm(2,0), n(2,0);
    double discri=b*b-4*a*c;
    if (discri<0) {          // Does not hit the axon
        prob=0;
    }
    else {
        discri=sqrt(discri);
        t1=0.5/a*( -b+discri );
        t2=0.5/a*( -b-discri );
        t=min(t1,t2);
        //cout<<"t\t"<<t1<<'\t'<<t2<<'\t'<<t<<'\t'<<endl;
        if (t>=dx | t<0) {            // Does not hit the axon
            prob=0;
        }
        else {                  // Hit the axon
            
            // xm = x + t*v;
            xm[0]=x[0]+t*v[0];
            xm[1]=x[1]+t*v[1];
            //cout<<"xm\t"<<xm[0]<<'\t'<<xm[1]<<endl;
            
            // n parallel to (xm-xc_t)
            n[0]=(xm[0]-xc[0])/vecDistance(xm, xc);
            n[1]=(xm[1]-xc[1])/vecDistance(xm, xc);
            //cout<<"n\t"<<n[0]<<'\t'<<n[1]<<endl;
            
            
            prob=dx*probUnit*fabs(vecInProduct(n, v));      // Approximated solution, and assume that <2*t> ~ dx
            prob=prob/(1+prob);                             // Exact solution
        }
    }
    return prob;
}

double probIAS2Myelin (const vector<double> &x, const vector<double> &v, const double &dx, vector<double> xc, const double &rc, const double &g, const double &probUnit) {
    // Probability for a particle at x from IAS to myelin, circle (xc,rc), velocity v, step size dx, gratio g, res*kappa/Din = probUnit
    
    double prob=0;
    
    // Translate circle center x to make it as close to the position (x + dx*v) as possible
    //vector<double> xc_t(2,0);
    translateXc(x,xc);
    
    // distance( x+t*v, xc )==rc, solve t
    double a=0,b=0,c=0,t1=0,t2=0,t=0;
    a=v[0]*v[0] + v[1]*v[1];
    b=2*(x[0]-xc[0])*v[0] + 2*(x[1]-xc[1])*v[1];
    c=pow(x[0]-xc[0],2) + pow(x[1]-xc[1],2) - rc*rc*g*g;
    //cout<<"a,b,c\t"<<a<<'\t'<<b<<'\t'<<c<<endl;
    
    // xt: final position, xm: contact point on myelin, n: surface unit vector, v1: reflecting velocity
    vector<double> xt(2,0), xm(2,0), n(2,0);
    double discri=b*b-4*a*c;
    if (discri<0) {          // Does not hit the axon
        prob=0;
    }
    else {
        discri=sqrt(discri);
        t1=0.5/a*( -b+discri );
        t2=0.5/a*( -b-discri );
        t=max(t1,t2);   // Choose the t>0 assuring that particle goes to the correct direction
        //cout<<"t\t"<<t1<<'\t'<<t2<<'\t'<<t<<'\t'<<endl;
        if (t>=dx) {            // Does not hit the axon
            prob=0;
        }
        else {                  // Hit the axon
            
            // xm = x + t*v;
            xm[0]=x[0]+t*v[0];
            xm[1]=x[1]+t*v[1];
            //cout<<"xm\t"<<xm[0]<<'\t'<<xm[1]<<endl;
            
            // n parallel to (xm-xc_t)
            n[0]=(xm[0]-xc[0])/vecDistance(xm, xc);
            n[1]=(xm[1]-xc[1])/vecDistance(xm, xc);
            //cout<<"n\t"<<n[0]<<'\t'<<n[1]<<endl;
            
            
            prob=dx*probUnit*fabs(vecInProduct(n, v));      // Approximated solution, and assume that <2*t> ~ dx
            prob=prob/(1+prob);                             // Exact solution
        }
    }
    return prob;
}

double probMyelin2EAS (vector<double> x, double dxR, vector<double> xc, double rc, double probUnit) {
    // Probability for a particle at x from Myelin to EAS, circle (xc,rc), step size dxR, res*kappa/Dmr = probUnit
    
    // Translate circle center x to make it as close to the position x as possible
    //vector<double> xc_t(2,0);
    //xc=translateXc(x,xc);
    
    double t=disOutSheath(x, xc, rc);
    
    if (t>dxR) {
        return 0;
    }
    else {
        return dxR*probUnit/(1+dxR*probUnit);                                  // Exact solution, and assume that <2*t> ~ dxR
    }
    
}

double probMyelin2IAS (vector<double> x, double dxR, vector<double> xc, double rc, double g, double probUnit) {
    // Probability for a particle at x from Myelin to IAS, circle (xc,rc), step size dxR, res*kappa/Dmr = probUnit
    
    // Translate circle center x to make it as close to the position x as possible
    //vector<double> xc_t(2,0);
    //xc=translateXc(x,xc);
    
    double t=disInSheath(x, xc, rc, g);
    
    if (t>dxR) {
        return 0;
    }
    else {
        return dxR*probUnit/(1+dxR*probUnit);                                  // Exact solution, and assume that <2*t> ~ dxR
    }
    
}

double probMyelinOutward (const vector<double> &x, const double &dxR, vector<double> xc) {
    // Probability for a particle at x in the myelin diffuse outwardly, with a step size dxR, toward the circle center xc
    
    // Translate circle center x to make it as close to the position x as possible
    //vector<double> xc_t(2,0);
    translateXc(x,xc);
    
    
    return 0.5 + 0.25*dxR/vecDistance(x, xc);
    //double r0=vecDistance(x, xc);
    
    //return 0.5+0.25*dxR/r0;
    
}
 */




























