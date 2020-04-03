// =================================================================
// Licensed Materials - Property of IBM
//
// "Restricted Materials of IBM"
//
// BCM-YKT-07-18-2017
//
// (C) Copyright IBM Corp. 2005-2017  All rights reserved
//
// US Government Users Restricted Rights -
// Use, duplication or disclosure restricted by
// GSA ADP Schedule Contract with IBM Corp.
//
// =================================================================

#ifdef HAVE_MPI
#include <mpi.h>
#endif
#include "NumInt.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>

void numInt::initializeIterator(int numvars_, double dT_)
{
  dT = dT_; 
  x.increaseSizeTo(numvars_);
}

void RK4Phased::initializeIterator(int numvars_, double dT_)
{
     numInt::initializeIterator(numvars_,dT_);
     dT2=dT_/2.0; dT6 = dT_/6.0;
     dx1.increaseSizeTo(numvars_);
     dx2.increaseSizeTo(numvars_);
     dx3.increaseSizeTo(numvars_);
     x1.increaseSizeTo(numvars_);  
}

void RK4Phased::callIteratePhase1()
{
  derivs(x,dx1); 

  for (i1=x1.begin(),i2=x.begin(),i3=dx1.begin();i1!=x1.end();i1++,i2++,i3++) 
    *i1 = *i2 + dT2**i3;
}

void RK4Phased::callIteratePhase2()
{
  derivs(x1,dx2);

  for (i1=x1.begin(),i2=x.begin(),i3=dx2.begin();i1!=x1.end();i1++,i2++,i3++) 
    *i1 = *i2 + dT2**i3;
}

void RK4Phased::callIteratePhase3()
{
  derivs(x1,dx3);

  for (i1=x1.begin(),i2=x.begin(),i3=dx3.begin(),i4=dx2.begin();i1!=x1.end();i1++,i2++,i3++,i4++)
  { 
    *i1 = *i2 + dT**i3;
    *i3 += *i4;
  }
}

void RK4Phased::callIteratePhase4()
{
  derivs(x1,dx2);
  for (i1=dx1.begin(),i2=x.begin(),i3=dx3.begin(),i4=dx2.begin();i1!=dx1.end();i1++,i2++,i3++,i4++)
   *i2 = *i2 + dT6*(*i1 + *i4 + 2.0**i3);

}

