#include <math.h>
#include <float.h>
#include "libiad.h"
#include "mathlink.h"

static void error(char *s)
{
	MLEvaluate(stdlink,s);
	MLNextPacket(stdlink);
	MLNewPacket(stdlink);
  	MLPutFunction(stdlink, "List", 4);
	MLPutSymbol(stdlink,"$Failed");
	MLPutSymbol(stdlink,"$Failed");
	MLPutSymbol(stdlink,"$Failed");
	MLPutSymbol(stdlink,"$Failed");
}


:Begin:
:Function:       AD
:Pattern:        AD[n_Integer, a_Real, b_Real, g_Real, nslab_Real, ntop_Real, nbottom_Real]
:Arguments:      {n, a, b, g, nslab, ntop, nbottom}
:ArgumentTypes:  {Integer, Real, Real, Real, Real, Real, Real}
:ReturnType:     Manual
:End:
:Evaluate:       AD[n_Integer,a_,b_,g_,nslab_,ntop_,nbottom_] := AD[n,N[a],N[b],N[g],N[nslab],N[ntop],N[nbottom]]
:Evaluate:       AD::badQuadNum    = "The number of quadrature angles must be a multiple of 8."
:Evaluate:       AD::badAnisotropy = "Anisotropy must be between -1 and 1."
:Evaluate:       AD::badThickness  = "Optical thickness must be non-negative."
:Evaluate:       AD::badAlbedo     = "Albedo must be between zero and one."
:Evaluate:       AD::badSlabIndex  = "Index of refraction of slab must be between 0.1 and 10."
:Evaluate:       AD::badTopIndex   = "Index of refraction of top slide must be between 0.1 and 10."
:Evaluate:       AD::badBotIndex   = "Index of refraction of bottom slide must be between 0.1 and 10."

void AD(int n, double a, double b, double g, double nslab, double ntop, double nbottom)
{
  double UR1, UT1, URU, UTU;
  double rr[4];
  
  if ((a<0) || (a>1))              {error("Message[AD::badAlbedo]");     return;}
  if (b<0)                         {error("Message[AD::badThickness]");  return;}
  if ((g<=-1) || (g>=1))           {error("Message[AD::badAnisotropy]"); return;}
  if ((n/8)*8 != n)                {error("Message[AD::badQuadNum]");    return;}
  if ((nslab<0.1)  ||(nslab>10))   {error("Message[AD::badSlabIndex]");  return;}
  if ((ntop<0.1)   ||(ntop>10))    {error("Message[AD::badTopIndex]");   return;}
  if ((nbottom<0.1)||(nbottom>10)) {error("Message[AD::badBotIndex]");   return;}

  ez_RT(n,nslab,ntop,nbottom, a, b, g, &UR1,&UT1,&URU,&UTU);
  rr[0] = UR1;
  rr[1] = UT1;
  rr[2] = URU;
  rr[3] = UTU;
  
  MLPutRealList( stdlink, rr, 4);
}

:Begin:
:Function: 		 ADThick
:Pattern: 		 ADThick[n_Integer, a_Real, g_Real, nslab_Real, ntop_Real]
:Arguments: 	 {n, a, g, nslab, ntop}
:ArgumentTypes:  {Integer, Real, Real, Real, Real}
:ReturnType: 	 Manual
:End:
:Evaluate:       ADThick[n_Integer,a_,g_,nslab_,ntop_] := ADThick[n,N[a],N[g],N[nslab],N[ntop]]

void ADThick(int n, double a, double g, double nslab, double ntop)
{
  double UR1, UT1, URU, UTU;
  double rr[2];
  
  if ((a<0) || (a>1))              {error("Message[AD::badAlbedo]");     return;}
  if ((g<=-1) || (g>=1))           {error("Message[AD::badAnisotropy]"); return;}
  if ((n/8)*8 != n)                {error("Message[AD::badQuadNum]");    return;}
  if ((nslab<0.1)  ||(nslab>10))   {error("Message[AD::badSlabIndex]");  return;}
  if ((ntop<0.1)   ||(ntop>10))    {error("Message[AD::badTopIndex]");   return;}

  ez_RT(n,nslab,ntop,nslab, a, HUGE_VAL, g, &UR1,&UT1,&URU,&UTU);

  rr[0] = UR1;
  rr[1] = URU;
  
  MLPutRealList( stdlink, rr, 2);
}

:Begin:
:Function:       IAD
:Pattern:        IAD[UR1_Real, UT1_Real, Tc_Real, nslab_Real, nslide_Real]
:Arguments:      {UR1, UT1, Tc, nslab, nslide}
:ArgumentTypes:  {Real, Real, Real, Real, Real}
:ReturnType:     Manual
:End:
:Evaluate:       IAD[Rt_, Tt_, Tc_, nslab_, nslide_] := IAD[N[Rt], N[Tt], N[Tc], N[nslab], N[nslide]]

void IAD(double UR1, double UT1, double Tc, double nslab, double nslide)
{
  double a, b, g;
  double rr[4];
  int error;
  
  ez_Inverse_RT(nslab,nslide,UR1,UT1,Tc,&a,&b,&g,&error);
  rr[0] = a;
  rr[1] = b;
  rr[2] = g;
  rr[3] = error;
  
  MLPutRealList( stdlink, rr, 4);
}

:Begin:
:Function:       ADCone
:Pattern:        ADCone[n_Integer, a_Real, b_Real, g_Real, mu_Real, nslab_Real, ntop_Real, nbottom_Real]
:Arguments:      {n, a, b, g, mu, nslab, ntop, nbottom}
:ArgumentTypes:  {Integer, Real, Real, Real, Real, Real, Real, Real}
:ReturnType:     Manual
:End:
:Evaluate:       ADCone[n_Integer,a_,b_,g_,mu_,nslab_,ntop_,nbottom_] := ADCone[n,N[a],N[b],N[g],N[mu],N[nslab],N[ntop],N[nbottom]]

void ADCone(int n, double a, double b, double g, double mu, double nslab, double ntop, double nbottom)
{
  double UR1, UT1, URU, UTU;
  double rr[4];
  
  if ((a<0) || (a>1))              {error("Message[AD::badAlbedo]");     return;}
  if (b<0)                         {error("Message[AD::badThickness]");  return;}
  if ((g<=-1) || (g>=1))           {error("Message[AD::badAnisotropy]"); return;}
  if ((n/8)*8 != n)                {error("Message[AD::badQuadNum]");    return;}
  if ((nslab<0.1)  ||(nslab>10))   {error("Message[AD::badSlabIndex]");  return;}
  if ((ntop<0.1)   ||(ntop>10))    {error("Message[AD::badTopIndex]");   return;}
  if ((nbottom<0.1)||(nbottom>10)) {error("Message[AD::badBotIndex]");   return;}

  ez_RT_Cone(n,nslab,ntop,nbottom, a, b, g, mu, &UR1,&UT1,&URU,&UTU);
  rr[0] = UR1;
  rr[1] = UT1;
  rr[2] = URU;
  rr[3] = UTU;
  
  MLPutRealList( stdlink, rr, 4);
}

:Begin:
:Function: 		 ADLayers
:Pattern: 		 ADLayers[n_Integer, a_List, b_List, g_List, nslab_Real, ntop_Real, nbottom_Real]
:Arguments: 	 {n, a, b, g, nslab, ntop, nbottom}
:ArgumentTypes:  {Integer, RealList, RealList, RealList, Real, Real, Real}
:ReturnType: 	 Manual
:End:
:Evaluate:       ADLayers[n_Integer,a_List,b_List,g_List,nslab_,ntop_,nbottom_] := ADLayers[n,N[a],N[b],N[g],N[nslab],N[ntop],N[nbottom]]

void ADLayers(int n, double *a, long alen, double *b, long blen, double *g, long glen, double nslab, double ntop, double nbottom)
{
  double UR1, UT1, URU, UTU;
  double rr[4];

  RT_Layers(n,nslab,ntop,nbottom,alen,a,b,g, &UR1, &UT1, &URU, &UTU);

  rr[0] = UR1;
  rr[1] = UT1;
  rr[2] = URU;
  rr[3] = UTU;
  
  MLPutRealList( stdlink, rr, 4);
}

:Begin:
:Function:       IADSphere
:Pattern:        IADSphere[setup_List, analysis_List, rsphere_List, tsphere_List, meas_List]
:Arguments:      {setup, analysis, rsphere, tsphere, meas}
:ArgumentTypes:  {RealList, RealList, RealList, RealList, RealList}
:ReturnType:     Manual
:End:
:Evaluate:       IADSphere[setup_List, analysis_List, rsphere_List, tsphere_List, meas_List] := IADSphere[N[setup], N[analysis], N[rsphere], N[tsphere], N[meas]] 

void IADSphere(double *setup,      long len_setup, 
	             double *analysis, long len_analysis,
	             double *sphere_r,   long len_sphere_r,
	             double *sphere_t,   long len_sphere_t,
	             double *meas,       long len_meas)
{
  double results[4];
  
  Spheres_Inverse_RT(setup, analysis,sphere_r,sphere_t,meas,results);
    
  MLPutRealList( stdlink, results, 4);
}

int main(int argc, char* argv[])
{
	return MLMain(argc, argv);
}

