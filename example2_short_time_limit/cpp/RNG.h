//
//  RNG.h
//  practice
//
//  Created by Hong-Hsi Lee on 2/16/17.
//  Copyright (c) 2017 Hong-Hsi Lee. All rights reserved.
//

#ifndef __practice__RNG__
#define __practice__RNG__

#include <iostream>

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
