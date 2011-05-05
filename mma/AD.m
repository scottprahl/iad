(* :Name:Optics`AD*)

(* :Title:Adding Doubling Calculations*)

(* :Author:Scott Prahl*)

(* :Summary:This package finds reflection and transmission for slabs.*)

(* :Context:Optics`AD`*)

(* :Package Version:1.4*)

(* :Copyright:Copyright 2010, Oregon Medical Laser Center*)

(* :History:Originally written by Scott Prahl, November 2002 *)

(* :Discussion: A bunch of handy Adding-Doubling scattering routines *)

BeginPackage["Optics`AD`"]

AD::usage = "AD[N,a,b,g,nslab,ntop,nbottom] returns {UR1,UT1,URU,UTU} for a slab with albedo a, optical thickness b, scattering anisotropy g, an index of refraction nslab, a top slide with index ntop, and a bottom slide with index nbottom.  N is the number of quadrature angles."

ADCone::usage = "ADCone[a,mu], ADCone[a,b,mu], ADCone[a,b,g,mu], ADCone[a,b,g,mu,nslab], ADCone[a,b,g,nslab,ntop,nbottom], and ADCone[N,a,b,g,nslab,ntop,nbottom] returns {UR1,UT1,URU,UTU} for a slab with albedo a, optical thickness b, scattering anisotropy g, an index of refraction nslab, a top slide with index ntop, and a bottom slide with index nbottom.  N is the number of quadrature angles."

ADThick::usage = "ADThick[a], ADThick[a,g], ADThick[a,g,nslab],ADThick[a,g,nslab,ntop,nbottom], and ADThick[N,a,g,nslab,ntop,nbottom] returns {UR1,URU} for a semi-infinite medium with albedo a, optical thickness b, scattering anisotropy g, an index of refraction nslab, a top slide with index ntop, and a bottom slide with index nbottom.  N is the number of quadrature angles."

ADLayers::usage = "ADLayers[{a,...},{b,...},{g,...}], ADLayers[{a,...},{b,...},{g,...},nslab], ADLayers[{a,...},{b,...},{g,...},nslab,ntop,nbottom], and ADLayers[N,{a,...},{b,...},{g,...},nslab,ntop,nbottom] returns {UR1,UT1,URU,UTU} for a layered slab with albedos {a,...}, optical thickness {b,...}, scattering anisotropy {g,...}, an index of refraction nslab, a top slide with index ntop, and a bottom slide with index nbottom.  N is the number of quadrature angles."

RTLayers::usage = "RTLayers[{mua,...},{mus,...},{g,...},{d,...}], RTLayers[{mua,...},{mus,...},{g,...},{d,...},nslab], RTLayers[{mua,...},{mus,...},{g,...},{d,...},nslab,ntop,nbottom], and RTLayers[N,{mua,...},{mus,...},{g,...},{d,...},nslab,ntop,nbottom] returns {UR1,UT1,URU,UTU} for a layered slab with absorption coefficents {mua,...}, scattering coefficients {mus,...}, anisotropies {g,...}, and thicknesses {d,...}.  The index of refraction of each layer in the slab is nslab, a top slide with index ntop, and a bottom slide with index nbottom.  N is the number of quadrature angles."

UR1Layers::usage = "UR1Layers[{mua,...},{mus,...},{g,...},{d,...}], UR1Layers[{mua,...},{mus,...},{g,...},{d,...},nslab], UR1Layers[{mua,...},{mus,...},{g,...},{d,...},nslab,ntop,nbottom], and UR1Layers[N,{mua,...},{mus,...},{g,...},{d,...},nslab,ntop,nbottom] returns UR1 for a layered slab with absorption coefficents {mua,...}, scattering coefficients {mus,...}, anisotropies {g,...}, and thicknesses {d,...}.  The index of refraction of each layer in the slab is nslab, a top slide with index ntop, and a bottom slide with index nbottom.  N is the number of quadrature angles."

UT1Layers::usage = "UT1Layers[{mua,...},{mus,...},{g,...},{d,...}], UT1Layers[{mua,...},{mus,...},{g,...},{d,...},nslab], UT1Layers[{mua,...},{mus,...},{g,...},{d,...},nslab,ntop,nbottom], and UT1Layers[N,{mua,...},{mus,...},{g,...},{d,...},nslab,ntop,nbottom] returns UT1 for a layered slab with absorption coefficents {mua,...}, scattering coefficients {mus,...}, anisotropies {g,...}, and thicknesses {d,...}.  The index of refraction of each layer in the slab is nslab, a top slide with index ntop, and a bottom slide with index nbottom.  N is the number of quadrature angles."

UR1::usage = "UR1[a], UR1[a,b], UR1[a,b,g], UR1[a,b,g,nslab], UR1[a,b,g,nslab,ntop,nbottom], and UR1[N, a,b,g,nslab,ntop,nbottom] returns the total reflection for normal irradiance on a slab with albedo a, optical thickness b, scattering anisotropy g, an index of refraction nslab, a top slide with index ntop, and a bottom slide with index nbottom.  N is the number of quadrature angles."

UT1::usage = "UT1[a,b], UT1[a,b,g], UT1[a,b,g,nslab], UT1[a,b,g,nslab,ntop,nbottom], UT1[N,a,b,g,nslab,ntop,nbottom] returns the total transmission for normal irradiance on a slab with albedo a, optical thickness b, scattering anisotropy g, an index of refraction nslab, a top slide with index ntop, and a bottom slide with index nbottom.  N is the number of quadrature angles."

URU::usage = "URU[a], URU[a,b], URU[a,b,g], URU[a,b,g,nslab], URU[a,b,g,nslab,ntop,nbottom], URU[N,a,b,g,nslab,ntop,nbottom] returns the total reflection for diffuse irradiance on a slab with albedo a, optical thickness b, scattering anisotropy g, an index of refraction nslab, a top slide with index ntop, and a bottom slide with index nbottom.  N is the number of quadrature angles."

UTU::usage = "UTU[a,b], UTU[a,b,g], and UTU[a,b,g,nslab], UTU[a,b,g,nslab,ntop,nbottom], UTU[N,a,b,g,nslab,ntop,nbottom] returns the total transmission for diffuse irradiance on a slab with albedo a, optical thickness b, scattering anisotropy g, an index of refraction nslab, a top slide with index ntop, and a bottom slide with index nbottom.  N is the number of quadrature angles."

ADDefaultQuadraturePoints::usage = "ADDefaultQuadraturePoints is typically 16 and is the number of quadrature points used to numerically approximate integrals during adding-doubling calculations.  This number should be a multiple of 8 and at least 16."

IAD::usage = "IAD[Rt, Tt, Tc, nslab, nslide] calculates {a,b,g} for a slab with total reflectance Rt, total transmission Tt, unscattered transmission Tc. The index of refraction of the slab is nslab, the index of refraction of the top and bottom slides is nslide."

IADSphere::usage = "IADSphere[setup_List, analysis_List, rsphere_List, tsphere_List, meas_List, results_List] calculates the optical properties for a particular experiment."

Begin["`Private`"]

Install["Optics/External/AD"]

ADDefaultQuadraturePoints = 16;

ADCone[a_,b_,g_,mu_,nslab_,ntop_,nbottom_] := ADCone[ADDefaultQuadraturePoints,N[a],N[b],N[g],N[mu],N[nslab],N[ntop],N[nbottom]]
ADCone[a_,b_,g_,mu_,nslab_] := ADCone[ADDefaultQuadraturePoints,N[a],N[b],N[g],N[mu],N[nslab],1.0,1.0]
ADCone[a_,b_,g_,mu_] := ADCone[ADDefaultQuadraturePoints,N[a],N[b],N[g],N[mu],1.0,1.0,1.0]
ADCone[a_,b_,mu_] := ADCone[ADDefaultQuadraturePoints,N[a],N[b],0.0,N[mu],1.0,1.0,1.0]
ADCone[a_,mu_] := ADCone[ADDefaultQuadraturePoints,N[a],10000.0,0.0,N[mu],1.0,1.0,1.0]

ADThick[a_, g_:0.0, nslab_:1.0, ntop_:1.0, n_Integer:ADDefaultQuadraturePoints] := ADThick[n,N[a],N[g],N[nslab],N[ntop]]

RTLayers[mua_List, mus_List, g_List, d_List, nslab_:1.0, ntop_:1.0, nbottom_:1.0, n_Integer:ADDefaultQuadraturePoints] := 
    ADLayers[n,N[mus/(mua+mus)],N[d(mus+mua)],N[g],N[nslab],N[ntop],N[nbottom]]

UR1Layers[mua_List, mus_List, g_List, d_List, nslab_:1.0, ntop_:1.0, nbottom_:1.0, n_Integer:ADDefaultQuadraturePoints] := 
    ADLayers[n,N[mus/(mua+mus)],N[d(mus+mua)],N[g],N[nslab],N[ntop],N[nbottom]][[1]]

UT1Layers[mua_List, mus_List, g_List, d_List, nslab_:1.0, ntop_:1.0, nbottom_:1.0, n_Integer:ADDefaultQuadraturePoints] := 
    ADLayers[n,N[mus/(mua+mus)],N[d(mus+mua)],N[g],N[nslab],N[ntop],N[nbottom]][[2]]

UR1[a_, b_, g_:0.0, nslab_:1.0, ntop_:1.0, nbottom_:1.0, n_Integer:ADDefaultQuadraturePoints] := 
     AD[ADDefaultQuadraturePoints,a,b,g,nslab,ntop,nbottom][[1]]
UT1[a_, b_, g_:0.0, nslab_:1.0, ntop_:1.0, nbottom_:1.0, n_Integer:ADDefaultQuadraturePoints] := 
     AD[ADDefaultQuadraturePoints,a,b,g,nslab,ntop,nbottom][[2]]
URU[a_, b_, g_:0.0, nslab_:1.0, ntop_:1.0, nbottom_:1.0, n_Integer:ADDefaultQuadraturePoints] := 
     AD[ADDefaultQuadraturePoints,a,b,g,nslab,ntop,nbottom][[3]]
UTU[a_, b_, g_:0.0, nslab_:1.0, ntop_:1.0, nbottom_:1.0, n_Integer:ADDefaultQuadraturePoints] := 
     AD[ADDefaultQuadraturePoints,a,b,g,nslab,ntop,nbottom][[4]]
     
UR1[a_] := ADThick[ADDefaultQuadraturePoints,a,0.0,1.0,1.0][[1]]
URU[a_] := ADThick[ADDefaultQuadraturePoints,a,0.0,1.0,1.0][[2]]
	
IAD[Rt_, Tt_, Tc_, nslab_] := IAD[N[Rt], N[Tt], N[Tc], N[nslab], 1.0]
IAD[Rt_, Tt_, Tc_] := IAD[N[Rt], N[Tt], N[Tc], 1.0, 1.0]
IAD[Rt_, Tt_] := IAD[N[Rt], N[Tt], 0.0, 1.0, 1.0]
IAD[Rt_] := IAD[N[Rt], 0.0, 0.0, 1.0, 1.0]

      
End [ ]

EndPackage[ ]