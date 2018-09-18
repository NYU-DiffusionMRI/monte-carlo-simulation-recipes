//
//  RNG.h
//
//  Created by Hong-Hsi Lee in February, 2017.
//

#ifndef __practice__RNG__
#define __practice__RNG__

#include <iostream>
//#include <random>

//Define structure used in RNG
struct RNGvalues
{
	unsigned int x;	   //seed
	unsigned int y;    //seed
	unsigned int z;    //seed
	unsigned int c;    //seed
	double r;   //random number output from RNG
};

unsigned int devrand();
RNGvalues init_KISS();
RNGvalues JKISS(RNGvalues);


#endif /* defined(__practice__RNG__) */
