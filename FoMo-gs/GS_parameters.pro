pro DL_parameters,n_0=n_0,T=T,B=B,n_b=n_b,angle=angle,length=length,Nf=Nf,Itot=Itot,Vtot=Vtot,f=f

; INPUT:
; n_0[length] = plasma number desity in cm-3 
; T[length] = plasma temperature in K
; B[length] = magnetic field strength in G
; n_0[length] = nonthermal number desity in cm-3
; angle[length] = B-LOS in deg  
; length = number of pixels alonf the LOS, defined by the cylinder depth 
; Nf = number of frequensies in GS spectrum
;
; OUTPUT:
; Itot[Nf] = total intensity (Stokes I) in SFU
; Vtot[Nf] = total polarization (Stokes V) in SFU
; f[Nf] = emission frequencies in GHz

;EXTERNAL LIBRARY DIRECTORY (main GS code):
libname='E:\Belgium\Wiki\FoMo-GS-code\MW_Transfer_64.dll'; for Win64-bit

;libname='/.../MWTransfer.so'; for Linux 64-bit

;Firstly, we create one-dimensional array describing the properties of a single volume element
 ParmIn=fltarr(29)
 ;Input parameters for GS intensity calculations
 ParmIn[0] = 25.e5^2.  ;Area, cm^2
 ParmIn[1] = 25.e5     ;the length of the pixel  along the line of sight [cm]
 ;ParmIn[2] = 1.e6      ;T_0, K !!!Will be Changed letter!!!!!!
 ;ParmIn[3] =0.05   ;\eps (not used in this casee)
 ;ParmIn[4] =4.0    ;\kappa (not used in this case)
 ParmIn[5] =16      ;number of integration nodes (the integration over energy using the trapezoidal method)
 ParmIn[6] =0.01     ;E_min, MeV
 ParmIn[7] =10.     ;E_max, MeV
 ParmIn[8]=0.5     ;E_break, MeV (not used in this case)
 ParmIn[9] =1.5     ;\delta_1 the power-law (PL) index in the single-power-law (SPL) distrib-n or low-energy PL index
 ParmIn[10]=3.0    ;\delta_2 the high-energy PL index in the DPL distribution 
 ;ParmIn[11]=1.6e9   ;n_0 - thermal electron density, cm^{-3} !!!Will be Changed letter
 ParmIn[12]=2.5e6    ;n_b - nonthermal electron density, cm^{-3}
 ;ParmIn[13]=10     ;B - magnetic field, G
 ;ParmIn[14]= angle  ;theta - the viewing angle relative to B, degrees, (must be < 90) 
 ParmIn[15]=1.e9    ;f0 - starting frequency to calculate spectrum in Hz--------------------------------!!!!!!!!!!!!!!!!!!!!!!!!!
 ParmIn[16]=0.2     ;d - logarithmic step in frequensy, f1=f0^d in Hz -----------------------------------!!!!!!!!!!!!!!!!!!!!!!!!!
 ParmIn[17]=4       ;distribution over energy (2-thermal; 3-Single-Power-Law, 4-DPL, 5-thermal/nonthermal, 6-Kappa)
 ParmIn[18]=Nf      ;number of frequencies (specified above)
 ParmIn[19]=1.       ;distribution over pitch-angle (0/1-Isotropic; 2-exp-loss-cone; 3-Gaussian-loss-cone; 4-Gauss; 5-Super-Gauss)
 ;ParmIn[20]=90      ;loss-cone boundary, degrees
 ;ParmIn[21]=90     ;beam direction (degrees) in GAU and SGA
 ;ParmIn[22]=0.4     ;\Delta\mu
 ;ParmIn[23]=1       ;a_4 in SGA (not used in this example)
 ParmIn[25]=0        ;f_cr/f_b (boundary freq. in the gyrofreq.,0 means that the spectrum is purely continuous)
 ParmIn[26]=ParmIn[25] ;f^WH_cr (not used in case of f^C_cr=0)
 ;ParmIn[27]=1        ;controls renormalisation in the hybrid code
 ParmIn[28]=2        ;Q-optimization on (with 2 bisection steps before including the term \ln Q)

 NSteps=length 
 Parms=fltarr(29, NSteps) ;the array of input parameters
 for i=0, NSteps-1 do begin
  Parms[*, i] =ParmIn  ;The parameters of all volume elements (pixels) are initially the same
;Parameters which are different for different pixels on the line of sight:
  Parms[2, i]=T(i)  
  Parms[11,i]=n_0(i) 
  Parms[12,i]=n_b(i)
  Parms[13,i]=B(i)
  Parms[14,i]=angle(i)
 endfor

RL=fltarr(7, Nf) ;the array for output parameters
res=call_external(libname, 'GET_MW', NSteps, Parms, RL, /f_value) ;calling the fast GS code

;Extracting the data from the output array (the exact coupling model is considered):
 f=RL[0, *]   ;emission frequency, GHz
 I_L=reform(RL[5, *]) ;left-polarized emission intensity (as observed from the Earth), sfu
 I_R=reform(RL[6, *]) ;right-polarized emission intensity (as observed from the Earth), sfu
 Itot=I_L+I_R ;total intensity (Stokes I)
 Vtot=I_L-I_R ;total polarization (Stokes V)

end
