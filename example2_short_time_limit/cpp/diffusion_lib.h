//
//  diffusion_lib.h
//  practice
//
//  Created by Hong-Hsi Lee on 2/16/17.
//  Copyright (c) 2017 Hong-Hsi Lee. All rights reserved.
//

#ifndef __practice__diffusion_lib__
#define __practice__diffusion_lib__

#include <iostream>
#include <vector>

using namespace std;

double diffclock(clock_t ,clock_t );

void vecSubtract (const vector<double> &, const vector<double> &, vector<double> &);

inline double vecDistance (const vector<double> &, const vector<double> &);

inline double vecInProduct (const vector<double> &, const vector<double> &);

inline double vecNorm (const vector<double> &);

vector<int> pixPosition (vector<double>, unsigned int);

void translateXc (const vector<double> &, vector<double> &);

bool inAxon (const vector<double> &, vector<double>, const double &);

bool inMyelin (const vector<double> &, vector<double>, const double &, const double &);

bool inIAS ( const vector<double> &, vector<double>, const double &, const double &);

bool stepEAS2Axon (const vector<double> &, const vector<double> &, vector<double>, const double &);

double disOutSheath ( const vector<double> &, vector<double>, const double &);

double disInSheath ( const vector<double> &, vector<double>, const double &, const double &);

vector<double> elasticCollisionEAS (const vector<double> &, const vector<double> &, const double &, vector<double>, const double &);

vector<double> inelasticCollisionEAS (const vector<double> &, const vector<double> &, const double &, vector<double>, const double &, const double &);

vector<double> elasticCollisionIAS (const vector<double> &, const vector<double> &, const double &, vector<double>, const double &, const double &);

vector<double> inelasticCollisionIAS (const vector<double> &, const vector<double> &, const double &, vector<double>, const double &, const double &, const double &);

vector<double> diffuseMyelinEllipseCircum (const vector<double> &, const double &, const double &, vector<double>);

vector<double> diffuseMyelinEllipseRadial (const vector<double> &, const double &, const double &, vector<double>, const double &);

vector<double> elasticCollisionMyelin2IAS (const vector<double> &, const double &, vector<double>, const double &, const double &);

vector<double> permeateMyelin2IAS (const vector<double> &, const double &, vector<double>, const double &, const double &, const double &);

vector<double> elasticCollisionMyelin2EAS (const vector<double> &, const double &, vector<double>, const double &);

vector<double> permeateMyelin2EAS (const vector<double> &, const double &, vector<double>, const double &, const double &);

#endif /* defined(__practice__diffusion_lib__) */
