#ifndef NUMINT_H
#define NUMINT_H
#include "Copyright.h"

#include "rndm.h"
#include "ShallowArray.h"

#define DPREC 10.0e-12

class numInt
{
 protected :

  virtual void initializeIterator(int, double) = 0;
  virtual void callIteratePhase1() = 0;
  virtual void callIteratePhase2() = 0;
  virtual void callIteratePhase3() = 0;
  virtual void callIteratePhase4() = 0;
  virtual void prePhase1() = 0;
  virtual void prePhase2() = 0;
  virtual void prePhase3() = 0;
  virtual void prePhase4() = 0;
  virtual void derivs(const ShallowArray< double > &, ShallowArray< double > &) = 0;
  virtual void flushVars(const ShallowArray< double > &) = 0;

 protected:
  ShallowArray< double > x;
  double dT;
};


class RK4Phased  : public numInt
{
  protected:
  
  void initializeIterator(int, double);
  void callIteratePhase1();
  void callIteratePhase2();
  void callIteratePhase3();
  void callIteratePhase4();
  void prePhase1(){flushVars(x);}
  void prePhase2(){flushVars(x1);}
  void prePhase3(){flushVars(x1);}
  void prePhase4(){flushVars(x1);}

 private:

  ShallowArray< double > dx1, dx2, dx3, x1;
  ShallowArray< double >::iterator i1, i2, i3, i4;
  double dT2, dT6;
};