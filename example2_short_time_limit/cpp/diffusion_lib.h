//
//  diffusion_lib.h
//  practice
//
//  Created by magda on 2/16/17.
//  Copyright (c) 2017 Honghsi. All rights reserved.
//

#ifndef __practice__diffusion_lib__
#define __practice__diffusion_lib__

#include <iostream>
#include <vector>

using namespace std;

double diffclock(clock_t ,clock_t );

//vector<double> vecSum (vector<double>, vector<double>);

void vecSubtract (const vector<double> &, const vector<double> &, vector<double> &);

inline double vecDistance (const vector<double> &, const vector<double> &);

inline double vecInProduct (const vector<double> &, const vector<double> &);

inline double vecNorm (const vector<double> &);

vector<int> pixPosition (vector<double>, unsigned int);

//inline bool needTranslateXc (vector<double>, double);

void translateXc (const vector<double> &, vector<double> &);

bool inAxon (const vector<double> &, vector<double>, const double &);

bool inMyelin (const vector<double> &, vector<double>, const double &, const double &);

bool inIAS ( const vector<double> &, vector<double>, const double &, const double &);

bool stepEAS2Axon (const vector<double> &, const vector<double> &, vector<double>, const double &);

//bool stepIAS2nonIAS (vector<double>, vector<double>, double, double);

//bool stepMyelin2EAS (vector<double>, vector<double>, double);

//bool stepMyelin2IAS (vector<double>, vector<double>, double, double);

double disOutSheath ( const vector<double> &, vector<double>, const double &);

double disInSheath ( const vector<double> &, vector<double>, const double &, const double &);

//double disCenter ( vector<double>, vector<double>);

vector<double> elasticCollisionEAS (const vector<double> &, const vector<double> &, const double &, vector<double>, const double &);

vector<double> inelasticCollisionEAS (const vector<double> &, const vector<double> &, const double &, vector<double>, const double &, const double &);

vector<double> elasticCollisionIAS (const vector<double> &, const vector<double> &, const double &, vector<double>, const double &, const double &);

vector<double> inelasticCollisionIAS (const vector<double> &, const vector<double> &, const double &, vector<double>, const double &, const double &, const double &);

//vector<double> diffuseCircum (vector<double>, double, double, vector<double>);

//vector<double> diffuseRadial (vector<double>, double, double, vector<double>);

vector<double> diffuseMyelinEllipseCircum (const vector<double> &, const double &, const double &, vector<double>);

vector<double> diffuseMyelinEllipseRadial (const vector<double> &, const double &, const double &, vector<double>, const double &);

vector<double> elasticCollisionMyelin2IAS (const vector<double> &, const double &, vector<double>, const double &, const double &);

vector<double> permeateMyelin2IAS (const vector<double> &, const double &, vector<double>, const double &, const double &, const double &);

vector<double> elasticCollisionMyelin2EAS (const vector<double> &, const double &, vector<double>, const double &);

vector<double> permeateMyelin2EAS (const vector<double> &, const double &, vector<double>, const double &, const double &);

/*
double probEAS2Myelin (const vector<double> &, const vector<double> &, const double &, vector<double>, const double &, const double &);

double probIAS2Myelin (const vector<double> &, const vector<double> &, const double &, vector<double>, const double &, const double &, const double &);

double probMyelin2EAS (vector<double>, double, vector<double>, double, double);

double probMyelin2IAS (vector<double>, double, vector<double>, double, double, double);

double probMyelinOutward (const vector<double> &, const double &, vector<double>);
 */


#endif /* defined(__practice__diffusion_lib__) */
