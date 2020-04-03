#include "Lens.h"
#include "MahonUnit.h"
#include "CG_MahonUnit.h"
#include "GridLayerData.h"
#include "rndm.h"
#include "NumInt.h"

#define SHD getSharedMembers()
#define ITER getSimulation().getIteration()
#define TIME ITER*getSharedMembers().deltaT
#define RANK getSimulation().getRank()
#define NUMVARS 12

template <typename T>
bool Threshold(T val, T thresh)
{
  return val>thresh;
}

template <typename T>
T sigmoid(const T & V, const T Vb, const T k) 
{
return 1.0/(1.0 + exp(-1.0*(V - Vb)/k));
}

template <typename T>
T pow(const T & val, const int & q)
{
  int i = 0;
  T retval = 1.0;
  while(i<q) {retval*=val;i++;}
  return val;
}


template <typename T>
T IonChannel(const T & V, const T & m, const T & h, const T g, const T Vb, 
	     const int p, const int q) 
{
  return g*pow(m,p)*pow(h,q)*(V-Vb);
}

template <typename T>
T IonChannel31(const T & V, const T & m, const T & h, const T g, const T Vb) 
{
  return g*m*m*m*h*(V-Vb);
}

template <typename T>
T IonChannel4(const T & V, const T & m, const T g, const T Vb) 
{
  return g*m*m*m*m*(V-Vb);
}

template <typename T>
T IonChannel(const T & V, const T & m, const T & h, const T g, const T Vb) 
{
  return g*m*h*(V-Vb);
}

template <typename T>
T IonChannel(const T & V, const T & m, const T g, const T Vb) 
{
  return g*m*(V-Vb);
}


template <typename T>
T IonChannel(const T & V, const T g, const T Vb) 
{
  return g*(V-Vb);
}


template <typename T>
T ratefcn(const T & x, const T xb, const T t) 
{
  return (xb - x)/t;
}


template <typename T>
T taufcn(const T & x, const T tau1, const T phi, const T sig0) 
{
  const T val1 = (x - phi)/sig0;
  return tau1/(exp(val1) + exp(val1*-1.0));
}


template <typename T>
T Ashtaufcn(const T & x)
{
  const T val1 = (x + 38.2)/28.0;
  return 1790.0 + 2930.0*exp(val1*val1*-1.0)*val1;
}



template <typename T> //wang buzaki 96 V shifted 7mv
T Kmalpha(const T & x)
{
  const T val1 = x + 27.0;
  return -0.01*val1/(exp(-0.1*val1) - 1.0);
}

template <typename T> //wang buzaki 96  V shifted 7mv
T Kmbeta(const T & x)
{
  return 0.125*exp((x+37.0)/-80.0);
}

template <typename T> //wang buzaki 96 V shifted 7mv
T Namalpha(const T & x)
{
  const T val1 = x + 28.0;
  return -0.1*val1/(exp(-0.1*val1) - 1.0);
}

template <typename T> //wang buzaki 96  V shifted 7mv
T Nambeta(const T & x)
{
  return 4.0*exp((x+53.0)/-18.0);
}

template <typename T> //wang buzaki 96 V shifted 7mv
T Nahalpha(const T & x)
{
  return 0.07*exp(-(x+51.0)/20.0);
}

template <typename T> //wang buzaki 96  V shifted 7mv
T Nahbeta(const T & x)
{
  return 1.0/(exp(-0.1*(x+21.0)) + 1.0);
}


template <typename T> //wang buzaki 96
T gatefcn(const T & x, const T alpha, const T beta, const T scale)
{
  return scale * (alpha*(1.0-x) - beta*x);
}


template <typename T> //wang buzaki 96
T gatefcnInstant(const T alpha, const T beta)
{
  return alpha / (alpha + beta);
}

template <typename T>
T Tadj(const T q10, const T cels, const T Etemp)
{
  return pow(q10, (cels-Etemp)/10.0);
}



template <typename T>
T TadjAdj(const T x, const T tadj)
{
  return x/tadj;
}

#define q10 2.5
#define CELSIUS 37
#define ETEMPKAF 22 
#define ETEMPNAS 21 
#define TADJPKAF 3.953
#define TADJPNAS 4.332


#define GLEAK 0.075
#define VLEAK -90.0 
#define GNAS 0.11
#define VNAS 40.0
#define GNAP 0.02
#define VNAP 45.0
#define GKRP 0.42
#define VKRP -77.5
#define GAS 0.32
#define VAS -85.0
#define GAF 0.09
#define VAF -73.0
#define GKIR 0.15
#define VKIR -90.0
#define GKCHAN 6.0
#define VKCHAN -90.0
#define GNACHAN 35.0
#define VNACHAN 55.0

#define NASMTHE -16.0 
#define NASMK 9.4 
#define NASMTAU1 637.8
#define NASMPHI -33.5
#define NASMSIG0 26.3

#define NAPMTHE -47.8 
#define NAPMK 3.1 
#define NAPMTAU1 1.0

#define KRPMTHE -13.4 
#define KRPMK 12.1 
#define KRPMTAU1 206.2 
#define KRPMPHI -53.9 
#define KRPMSIG0 26.5

#define KRPHTHE -55.0 
#define KRPHK -19.0 

#define ASHTHE -78.8
#define ASHK -10.4

#define XASMTHE -25.6
#define XASMK 13.3
#define XASMTAU1 131.4
#define XASMPHI -37.4
#define XASMSIG0 27.3

#define AFMTHE -33.1 
#define AFMK 7.5
#define AFMTAU1 1.0

#define AFHTHE -70.4 
#define AFHK -7.6
#define AFHTAU1 25.0

#define KIRMTHE -100.0
#define KIRMK -10.0
#define KIRMTAU1 0.01

#define VCL -80.0
#define SYNA 2
#define SYNB 0.1

#define CAPAC 1.0

void MahonUnit::derivs(const ShallowArray< double > & x, ShallowArray< double > & dx)
{

  const double & V = x[0];
  double & dV = dx[0];
  const double & Nasm = x[1];
  double & dNasm = dx[1];
  const double & Napm = x[2];
  double & dNapm = dx[2];
  const double & Krpm = x[3];
  double & dKrpm = dx[3];
  const double & Krph = x[4];
  double & dKrph = dx[4];
  const double & Asm = x[5];
  double & dAsm = dx[5];
  const double & Ash = x[6];
  double & dAsh = dx[6];
  const double & Afm = x[7];
  double & dAfm = dx[7];
  const double & Afh = x[8];
  double & dAfh = dx[8];
  const double & Km = x[9];
  double & dKm = dx[9];
  const double & Nah = x[10];
  double & dNah = dx[10];

  const double & g = x[11];
  double & dg = dx[11];


  ShallowArray<Input>::iterator iter, end=MSNNetInps.end();
  double drive=0;
  for (iter=MSNNetInps.begin(); iter!=end; ++iter) {
    drive += *(iter->input)*iter->weight;
  }
  
  dg = SYNA*((double) (V>0.0))*(1.0 - g) - synb*g;

  dV = drive*(V - VCL);

  dV -= IonChannel<double>(V,GLEAK,VLEAK);

  dV -= IonChannel<double>(V,Nasm,GNAS,VNAS);
  dV -= IonChannel<double>(V,Napm,GNAP,VNAP);
  dV -= IonChannel<double>(V,Krpm,Krph,GKRP,VKRP);
  dV -= IonChannel<double>(V,Asm,Ash,GAS,VAS);
  dV -= IonChannel<double>(V,Afm,Afh,GAF,VAF);

  dV -= IonChannel<double>(V,sigmoid<double>(V,KIRMTHE,KIRMK),GKIR,VKIR);

  dV -= IonChannel4<double>(V,Km,GKCHAN,VKCHAN);
  dV -= IonChannel31<double>
    (V,gatefcnInstant<double>(Namalpha<double>(V),Nambeta<double>(V)),Nah,GNACHAN,VNACHAN);

  dV -= (*drivinp)*V;
  dV /=CAPAC;

  dNasm = ratefcn<double>(Nasm,sigmoid<double>(V,NASMTHE,NASMK),
			  TadjAdj<double>(taufcn<double>(V,NASMTAU1,NASMPHI,NASMSIG0),TADJPNAS));

  dNapm = ratefcn<double>(Napm,sigmoid<double>(V,NAPMTHE,NAPMK),TadjAdj<double>(NAPMTAU1,TADJPKAF));

  dKrpm = ratefcn<double>(Krpm,sigmoid<double>(V,KRPMTHE,KRPMK),
			  TadjAdj<double>(taufcn<double>(V,KRPMTAU1,KRPMPHI,KRPMSIG0),TADJPKAF));

  const double val1 = Ashtaufcn(V);

  dKrph = ratefcn<double>(Krph,sigmoid<double>(V,KRPHTHE,KRPHK),TadjAdj<double>(val1*3.0,TADJPKAF));

  dAsh = ratefcn<double>(Ash,sigmoid<double>(V,ASHTHE,ASHK),TadjAdj<double>(val1,TADJPKAF));

  dAsm = ratefcn<double>(Asm,sigmoid<double>(V,XASMTHE,XASMK),
			 TadjAdj<double>(taufcn<double>(V,XASMTAU1,XASMPHI,XASMSIG0),TADJPKAF));

  dAfm = ratefcn<double>(Afm,sigmoid<double>(V,AFMTHE,AFMK),TadjAdj<double>(AFMTAU1,TADJPKAF));

  dAfh = ratefcn<double>(Afh,sigmoid<double>(V,AFHTHE,AFHK),TadjAdj<double>(AFHTAU1,TADJPKAF));

  dKm = gatefcn<double>(Km,Kmalpha<double>(V),Kmbeta<double>(V),5.0);

  dNah = gatefcn<double>(Nah,Nahalpha<double>(V),Nahbeta<double>(V),5.0);


}
 

void MahonUnit::initialize(RNG& rng) 
{
  initializeIterator(NUMVARS,SHD.deltaT);
  
  double & V = x[0];
  double & Nasm = x[1];
  double & Napm = x[2];
  double & Krpm = x[3];
  double & Krph = x[4];
  double & Asm = x[5];
  double & Ash = x[6];
  double & Afm = x[7];
  double & Afh = x[8];
  
  double & Km = x[9];  
  double & Nah = x[10];

  V=V_init;
  Nasm = sigmoid<double>(V,NASMTHE,NASMK);
  Napm = sigmoid<double>(V,NAPMTHE,NAPMK);
  Krpm = sigmoid<double>(V,KRPMTHE,KRPMK);
  Krph = sigmoid<double>(V,KRPHTHE,KRPHK);
  Ash = sigmoid<double>(V,ASHTHE,ASHK);
  Asm = sigmoid<double>(V,XASMTHE,XASMK);
  Afm = sigmoid<double>(V,AFMTHE,AFMK);
  Afh = sigmoid<double>(V,AFHTHE,AFHK);
  Km = gatefcnInstant<double>(Kmalpha<double>(V),Kmbeta<double>(V));
  Nah = gatefcnInstant<double>(Nahalpha<double>(V),Nahbeta<double>(V));
    
  std::ofstream ofs;
  ofs.open ("weights.dat", std::ofstream::out | std::ofstream::app);
  outputWeights(ofs);
  ofs.close();
  
}

void MahonUnit::update1(RNG& rng) 
{
  callIteratePhase1();
}

void MahonUnit::update2(RNG& rng) 
{
  callIteratePhase2();
}

void MahonUnit::update3(RNG& rng) 
{
  callIteratePhase3();
}

void MahonUnit::update4(RNG& rng) 
{
  callIteratePhase4();
}

void MahonUnit::flushVars(const ShallowArray< double > & x)
{
  g_out = x[11];
}

void MahonUnit::flushVars1(RNG& rng) 
{
  prePhase1();
}

void MahonUnit::flushVars2(RNG& rng) 
{
  prePhase2();
}

void MahonUnit::flushVars3(RNG& rng) 
{
  prePhase3();
}

void MahonUnit::flushVars4(RNG& rng) 
{
  prePhase4();
}


void MahonUnit::updateOutputs(RNG& rng) 
{
  if (x[0]>= SHD.spikethresh && var1 < SHD.spikethresh) spike = true;
  else spike = false;
  var1 = x[0];
  var2 = x[5];

  ShallowArray<Input>::iterator iter, end=MSNNetInps.end();
  double drive=0;
  for (iter=MSNNetInps.begin(); iter!=end; ++iter) {
    drive += *(iter->input)*iter->weight;
  }

  var3 = drive; 
}

bool MahonUnit::ConnectP1(const String& CG_direction, const String& CG_component, NodeDescriptor* CG_node, Edge* CG_edge, VariableDescriptor* CG_variable, Constant* CG_constant, CG_MahonUnitInAttrPSet* CG_inAttrPset, CG_MahonUnitOutAttrPSet* CG_outAttrPset) 
{ 
   RNG rng;
     unsigned preIdx=CG_node->getNode()->getIndex();
  unsigned postIdx=getNode()->getIndex();
  unsigned inWeightIdx=(postIdx*89423217)+(preIdx*37);
  rng.reSeedShared(connectionSeed+inWeightIdx);
  double in = drandom(rng);    
  if (in<CG_inAttrPset->connectionProb) return true;
  return false;
}

MahonUnit::~MahonUnit() 
{
}

void MahonUnit::setIndices(const String& CG_direction, const String& CG_component, NodeDescriptor* CG_node, Edge* CG_edge, VariableDescriptor* CG_variable, Constant* CG_constant, CG_MahonUnitInAttrPSet* CG_inAttrPset, CG_MahonUnitOutAttrPSet* CG_outAttrPset) 
{
  MSNNetInps[MSNNetInps.size()-1].row = getGlobalIndex()+1; // +1 is for Matlab
  MSNNetInps[MSNNetInps.size()-1].col = CG_node->getGlobalIndex()+1; 

}


void MahonUnit::outputWeights(std::ofstream& fs)
{
  ShallowArray<Input>::iterator iter, end=MSNNetInps.end();

  for (iter=MSNNetInps.begin(); iter!=end; ++iter)
    fs<<iter->row<<" "<<iter->col<<" "<<iter->weight<<std::endl;
}
