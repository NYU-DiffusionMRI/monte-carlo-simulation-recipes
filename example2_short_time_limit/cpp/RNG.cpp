//
//  RNG.cpp
//  practice
//
//  Created by magda on 2/16/17.
//  Copyright (c) 2017 Honghsi. All rights reserved.
//

#include "RNG.h"
#include <iostream>
#include <stdio.h>
//#include <random>

using namespace std;

#define Pi 3.14159265



/***********************RANDOM NUMBER GENERATOR (KISS) BEGIN****************************/
//RNG JKISS from http://www0.cs.ucl.ac.uk/staff/d.jones/GoodPracticeRNG.pdf

//Seed initialization using entropy pool in /dev/urandom.  This should only be called once.
unsigned int devrand()
{
	
    unsigned int seed=0;
	FILE *urandom;
    
	urandom = fopen ("/dev/random", "r");
	fread (&seed, sizeof (seed), 1, urandom);
	fclose(urandom);
	cout<<"seed is " << seed<<endl;
    
	return seed;
}

//Initialize RNG
RNGvalues init_KISS()
{
	
	RNGvalues v;
    
	unsigned int x,y,z,c;
	//unsigned long long t;
    
	x=devrand();
	while (!(y=devrand())){} /* y must not be zero! */ z = devrand();
	/* We don’t really need to set c as well but let's anyway... */
	/* NOTE: offset c by 1 to avoid z=c=0 */
	c = devrand() % 698769068 + 1; /* Should be less than 698769069 */
    
	v.x = x;
	v.y = y;
	v.z = z;
	v.c = c;
	return v;
}

//RNG script

RNGvalues JKISS(RNGvalues v){
    
    unsigned long long t;
    
    v.x = 314527869 * v.x + 1234567;
    v.y ^= v.y << 5; v.y ^= v.y >> 7; v.y ^= v.y << 22;
    t = 4294584393ULL * v.z + v.c;
    v.c = t >> 32;
    v.z = static_cast<unsigned int>(t);
    
    v.r = (v.x+v.y+v.z)/ 4294967295.0;

    return v;
}

/***********************RANDOM NUMBER GENERATOR (KISS) END ******************************/


/********** Mersenne Twister RNG *********/
/*
RNGvalues init_randMT () {
    RNGvalues v;
    v.x=0; v.y=0; v.z=0; v.c=0;
    v.r=0;
    return v;
}

RNGvalues randMT (RNGvalues v) {
    mt19937 generator (v.x);
    uniform_real_distribution<double> dis(0.0, 1.0);
    v.r=dis(generator);
    return v;
}
*/





