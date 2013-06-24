(* ::Package:: *)

BeginPackage["WSScattering`"]


f::usage="f[k, \[Theta], m:938.272046, \
V0:-39.42463943157022-7.884927886314045*I, R:3, a:0.54] gives \
the scattering amplitude calculated via a partial wave \
expansion."

d\[Sigma]d\[CapitalOmega]::usage="d\[Sigma]d\[CapitalOmega][k, points:10, m:938.272046, \
V0:-39.42463943157022-7.884927886314045*I, R:3, a:0.54] gives a \
table with values of the differential scattering cross-section \
calculated from the scattering amplitude. The density of table \
entries increases with increasing value for points."

Eif0::usage="f[k, \[Theta], m:938.272046, \
V0:-39.42463943157022-7.884927886314045*I, R:3, a:0.54] gives \
the scattering amplitude calculated via the zeroth order eikonal \
expansion."

eik0d\[Sigma]d\[CapitalOmega]::usage="eik0d\[Sigma]d\[CapitalOmega][k, points:10, m:938.272046, \
V0:-39.42463943157022-7.884927886314045*I, R:3, a:0.54] gives a \
table with values of the differential scattering cross-section \
calculated from the zeroth order eikonal scattering amplitude. \
The density of table entries increases with increasing value for \
points."

Eif1::usage="f[k, \[Theta], m:938.272046, \
V0:-39.42463943157022-7.884927886314045*I, R:3, a:0.54] gives \
the scattering amplitude calculated via the first order eikonal \
expansion."

eik1d\[Sigma]d\[CapitalOmega]::usage="eik1d\[Sigma]d\[CapitalOmega][k, points:10, m:938.272046, \
V0:-39.42463943157022-7.884927886314045*I, R:3, a:0.54] gives a \
table with values of the differential scattering cross-section \
calculated from the first order eikonal scattering amplitude. \
The density of table entries increases with increasing value for \
points."

Eif2::usage="f[k, \[Theta], m:938.272046, \
V0:-39.42463943157022-7.884927886314045*I, R:3, a:0.54] gives \
the scattering amplitude calculated via the second order eikonal \
expansion."

eik2d\[Sigma]d\[CapitalOmega]::usage="eik2d\[Sigma]d\[CapitalOmega][k, points:10, m:938.272046, \
V0:-39.42463943157022-7.884927886314045*I, R:3, a:0.54] gives a \
table with values of the differential scattering cross-section \
calculated from the second order eikonal scattering amplitude. \
The density of table entries increases with increasing value for \
points."

Eif3::usage="f[k, \[Theta], m:938.272046, \
V0:-39.42463943157022-7.884927886314045*I, R:3, a:0.54] gives \
the scattering amplitude calculated via the third order eikonal \
expansion."

eik3d\[Sigma]d\[CapitalOmega]::usage="eik3d\[Sigma]d\[CapitalOmega][k, points:10, m:938.272046, \
V0:-39.42463943157022-7.884927886314045*I, R:3, a:0.54] gives a \
table with values of the differential scattering cross-section \
calculated from the third order eikonal scattering amplitude. \
The density of table entries increases with increasing value for \
points."


Begin["`Private`"]

\[HBar]c=197.32697178
V[V0_,R_,a_,r_]:=V0/(1+E^((r-R)/a))

s[m_,V0_,R_,a_,L_,k_]:=NDSolve[{u''[r]+(k^2-2*m/(\[HBar]c)^2*V[V0,R,a,r]-L*(L+1)/r^2)*u[r]==0,u[(1000*k)^(-1)]==0,u'[(1000*k)^(-1)]==8},u,{r,(1000*k)^(-1),5*R},MaxSteps->20000]
u[m_,V0_,R_,a_,L_,k_,r_]:=u[r]/.s[m,V0,R,a,L,k][[1]]

\[Beta][m_,V0_,R_,a_,L_,k_]:=5*R/Evaluate[u[m,V0,R,a,L,k,5*R]]*(D[u[m,V0,R,a,L,k,r],r]/.r->(5*R))-1
tan\[Delta][m_,V0_,R_,a_,L_,k_]:=N[(k*5*R*Derivative[0,1][SphericalBesselJ][L,k*5*R]-\[Beta][m,V0,R,a,L,k]*SphericalBesselJ[L,k*5*R])/(k*5*R*Derivative[0,1][SphericalBesselY][L,k*5*R]-\[Beta][m,V0,R,a,L,k]*SphericalBesselY[L,k*5*R])]

f[k_,\[Theta]_,m_:938.272046,V0_:-39.42463943157022-7.884927886314045*I,R_:3,a_:0.54]:=1/k*Sum[(2*L+1)*tan\[Delta][m,V0,R,a,L,k]/(1-I*tan\[Delta][m,V0,R,a,L,k])*LegendreP[L,Cos[\[Theta]]],{L,0,2*k*R}]

d\[Sigma]d\[CapitalOmega][k_,points_:10,m_:938.272046,V0_:-39.42463943157022-7.884927886314045*I,R_:3,a_:0.54]:=Table[{\[Theta],Abs[f[k,\[Theta],m,V0,R,a]]^2},{\[Theta],0,Pi,N[1/points]}]


\[Chi]0[m_,V0_,R_,a_,k_,b_]=Quiet[-2*m/(k*(\[HBar]c)^2)*NIntegrate[V[V0,R,a,Sqrt[b^2+z^2]],{z,0,15}]];
T0[m_,V0_,R_,a_,k_,b_]:=Exp[I*\[Chi]0[m,V0,R,a,k,b]]-1
Eif0[k_,\[Theta]_,m_:938.272046,V0_:-39.42463943157022-7.884927886314045*I,R_:3,a_:0.54]:=-I*k*NIntegrate[b*BesselJ[0,2*k*b*Sin[\[Theta]/2]]*T0[m,V0,R,a,k,b],{b,0,15}]

eik0d\[Sigma]d\[CapitalOmega][k_,points_:10,m_:938.272046,V0_:-39.42463943157022-7.884927886314045*I,R_:3,a_:0.54]:=Table[{N[\[Theta]/points],Abs[Eif0[k,\[Theta]/points,m,V0,R,a]]^2},{\[Theta],0,points*Pi,1}]


En[m_,k_]:=(k*\[HBar]c)^2/(2*m)
U[V0_,R_,a_,r_]:=V[V0,R,a,r]/V0
\[Epsilon][m_,V0_,k_]:=V0/(2*En[m,k])


\[Tau]1[m_,V0_,R_,a_,k_,b_]=Quiet[-k*\[Epsilon][m,V0,k]^2*NIntegrate[Evaluate[(1#+b*D[#,b])&[U[V0,R,a,Sqrt[z^2+b^2]]^2]],{z,0,15}]];
T1[m_,V0_,R_,a_,k_,b_]:=Exp[I*(\[Chi]0[m,V0,R,a,k,b]+\[Tau]1[m,V0,R,a,k,b])]-1
Eif1[k_,\[Theta]_,m_:938.272046,V0_:-39.42463943157022-7.884927886314045*I,R_:3,a_:0.54]:=-I*k*NIntegrate[b*BesselJ[0,2*k*b*Sin[\[Theta]/2]]*T1[m,V0,R,a,k,b],{b,0,15}]

eik1d\[Sigma]d\[CapitalOmega][k_,points_:10,m_:938.272046,V0_:-39.42463943157022-7.884927886314045*I,R_:3,a_:0.54]:=Table[{N[\[Theta]/points],Abs[Eif1[k,\[Theta]/points,m,V0,R,a]]^2},{\[Theta],0,points*Pi,1}]


\[Tau]2[m_,V0_,R_,a_,k_,b_]=Quiet[-k*\[Epsilon][m,V0,k]^3*NIntegrate[Evaluate[(1#+5/3*b*D[#,b]+1/3*b^2*D[#,{b,2}])&[U[V0,R,a,Sqrt[z^2+b^2]]^3]],{z,0,15}]-b*D[\[Chi]0[m,V0,R,a,k,b],b]^3/(24*k^2)];
\[Omega]2[m_,V0_,R_,a_,k_,b_]=Quiet[D[\[Chi]0[m,V0,R,a,k,b],b]*D[b*D[\[Chi]0[m,V0,R,a,k,b],b],b]/(8*k^2)];

T2[m_,V0_,R_,a_,k_,b_]:=Exp[I*(\[Chi]0[m,V0,R,a,k,b]+\[Tau]1[m,V0,R,a,k,b]+\[Tau]2[m,V0,R,a,k,b])]*Exp[-\[Omega]2[m,V0,R,a,k,b]]-1
Eif2[k_,\[Theta]_,m_:938.272046,V0_:-39.42463943157022-7.884927886314045*I,R_:3,a_:0.54]:=-I*k*NIntegrate[b*BesselJ[0,2*k*b*Sin[\[Theta]/2]]*T2[m,V0,R,a,k,b],{b,0,15}]

eik2d\[Sigma]d\[CapitalOmega][k_,points_:10,m_:938.272046,V0_:-39.42463943157022-7.884927886314045*I,R_:3,a_:0.54]:=Table[{N[\[Theta]/points],Abs[Eif2[k,\[Theta]/points,m,V0,R,a]]^2},{\[Theta],0,points*Pi,1}]


\[Tau]3[m_,V0_,R_,a_,k_,b_]=Quiet[-k*\[Epsilon][m,V0,k]^4*NIntegrate[Evaluate[(5/4#+11/4*b*D[#,b]+b^2*D[#,{b,2}]+1/12*b^3*D[#,{b,3}])&[U[V0,R,a,Sqrt[z^2+b^2]]^4]],{z,0,15}]-b*D[\[Tau]1[m,V0,R,a,k,b],b]*D[\[Chi]0[m,V0,R,a,k,b],b]^2/(8*k^2)];
\[Phi]3[m_,V0_,R_,a_,k_,b_]=Quiet[-k*\[Epsilon][m,V0,k]^2*NIntegrate[Evaluate[(1#+5/3*b*D[#,b]+1/3*b^2*D[#,{b,2}])&[(1/(2*k)*D[U[V0,R,a,Sqrt[z^2+b^2]],b])^2]],{z,0,15}]];
\[Omega]3[m_,V0_,R_,a_,k_,b_]=Quiet[(D[\[Chi]0[m,V0,R,a,k,b],b]*D[b*D[\[Tau]1[m,V0,R,a,k,b],b],b]+D[\[Tau]1[m,V0,R,a,k,b],b]*D[b*D[\[Chi]0[m,V0,R,a,k,b],b],b])/(8*k^2)];

T3[m_,V0_,R_,a_,k_,b_]:=Exp[I*(\[Chi]0[m,V0,R,a,k,b]+\[Tau]1[m,V0,R,a,k,b]+\[Tau]2[m,V0,R,a,k,b]+\[Tau]3[m,V0,R,a,k,b]+\[Phi]3[m,V0,R,a,k,b])]*Exp[-(\[Omega]2[m,V0,R,a,k,b]+\[Omega]3[m,V0,R,a,k,b])]-1
Eif3[k_,\[Theta]_,m_:938.272046,V0_:-39.42463943157022-7.884927886314045*I,R_:3,a_:0.54]:=-I*k*NIntegrate[b*BesselJ[0,2*k*b*Sin[\[Theta]/2]]*T3[k,b],{b,0,15}]

eik3d\[Sigma]d\[CapitalOmega][k_,points_:10,m_:938.272046,V0_:-39.42463943157022-7.884927886314045*I,R_:3,a_:0.54]:=Table[{N[\[Theta]/points],Abs[Eif3[k,\[Theta]/points,m,V0,R,a]]^2},{\[Theta],0,points*Pi,1}]

End[]


EndPackage[]
