;pro GS_emission, dimt=dimt, dimx=dimx, dimz=dimz, Nf=Nf, View_Angle=View_Angle, Radius=Radius

dimt=1 & dimx=204 & dimz=114 & Nf=11 & View_Angle=85. & Radius=42.

;Calculats GS intensity and polarisation spectra integrated along the LOS
;through inhomogeneous cylinder for each slice (X-position), Z-position, and time.
;Calls procedure 'GS_param.pro' where parameters for GS-emission calculation are defined.
;
; INPUT:
; dimt = number of points in time dimension
; dimx = dimension in x direction
; dimz = dimension in z direction
; Nf = number of frequensies in GS spectrum
; View_Angle = view angle in degree
; Radius = radius in pixels
; OUPUT:
; Produces IDL save files of GS-emission on the output of each slice
; Output sav file contains 3 variables (1 sav file per 1 slice):
; f[Nf] = emission frequencies
; Int[ray_num,Nf,dimt] = GS Intensity integrated along each LOS (ray_num) for each frequency (f), and time step (dimt)
; Pol[ray_num,Nf,dimt] = GS Polarisation integrated along each LOS for each frequency, and time step

; INPUT DIRECTORY:
	slicedir='E:\Belgium\Wiki\FoMo-GS-code\slices\'
	savedir='E:\Belgium\Wiki\FoMo-GS-code\'

window, 1, title='Total intensity', retain=2, xsize=500,ysize=600
Depth = findgen(Radius+1)
L = findgen(Radius+1)

;Calculate a depth of the cylinder along the LOS for different slices (at different X-positions)
Pat_angle=90.-View_Angle
csa=cos(Pat_angle*!dtor) & sna=sin(Pat_angle*!dtor) & tna=tan(Pat_angle*!dtor)
for X=0, Radius do begin
  h=X+1-Radius
  Depth(X) = 2*sqrt(Radius^2-h^2)/csa
  L(X)=fix(sqrt(Radius^2-h^2)+0.5)
  if (abs(Depth(X)) gt (sqrt((2*Radius)^2+dimz^2))) then Depth(X) = sqrt((2*Radius)^2+dimz^2)
endfor
start = long(Radius*2*sna)*2
finish = long(dimz*csa +0.5)

;number of rays along the slice
ray_num = Start + Finish + 1

Int=fltarr(ray_num,Nf,dimt) ;output intensity
Pol=fltarr(ray_num,Nf,dimt) ;output polarisation

;------------------------vectors LOS & vector B
WA=View_Angle*!dtor
Lx = sin(WA) & Ly = 0. & Lz = cos(WA)
Los=[Lx,Ly,Lz]
vectorB=fltarr(3,114)
An=fltarr(dimx,114)
;-------------------------

  for j=dimx/2-Radius,dimx/2 do begin
    for i=0, dimt-1 do begin
filename = 'slice_'+string(i,format='(i3.3)')+'t_'+string(j,format='(i3.3)')+'x'+'.sav'

;Restore slices with pertubed values in [CGS] unit system:
restore,slicedir+filename

Np1=reform(Np(0,*,*))      ;Np[1,dimy,dimz]- plasma number desity in cm-3
te1=reform(te(0,*,*))      ;Te[1,dimy,dimz]- plasma temperature in K
btot1=reform(btot(0,*,*))  ;Btot[1,dimy,dimz]- m.f. magnitute in G
Nb1=reform(Nb(0,*,*))
Br1=reform(br(0,*,*))
Bz1=reform(bz(0,*,*))

for k=0,dimx-1 do begin
fi=0.          
Y=j-dimx/2 & X=k-dimx/2 
;------------------------azimuth (fi) definition
if (Y eq 0) and (X gt 0) then fi=0.*!dtor 
if (Y eq 0) and (X le 0) then fi=180.*!dtor 
if (X eq 0) and (Y gt 0) then fi=90.*!dtor 
if (X eq 0) and (Y lt 0) then fi=270.*!dtor 
if (X ne 0) and (Y ne 0) then fi=atan(Y,X)

Bx=Br1(k,*)*cos(fi)
By=Br1(k,*)*sin(fi)
Bz=Bz1(k,*)

;------------------------B-LOS angle (An) definition 
vectorB=[[Bx(*)],[By(*)],[Bz(*)]] 
b=sqrt(Bx^2.+By^2.+Bz^2.) 
for zcoor=0,dimz-1 do begin  
 vecB=vectorB(zcoor,*)
 dotproduct = vecB # LOS
 An(k,zcoor) = !radeg*acos(dotproduct/b(zcoor))
endfor
endfor
;------------------------cut only pixels inside the cylinder in Y-direction 
l1=dimx/2-Radius & l2=dimx/2+Radius
np2=Np1(l1:l2,*) ;[l2-l1+1,114]
te2=te1(l1:l2,*)
Btot2=btot1(l1:l2,*)
nb2=nb1(l1:l2,*)
angle2=An(l1:l2,*)

length=sqrt((2*Radius)^2+dimz^2)  ;length along LOS in pix
print,'length=',length
y=indgen(length)*0 & z=indgen(length)*0
y1=indgen(length)*0 & z1=indgen(length)*0
;wset,1
;plot,z1,y1,xran=[-5,350],yran=[-5,90],/nodata,/xst,/yst,xtit='Z',ytit='Y'

i3=0
for m=-start,finish do begin
i2=0
For i1=0, length-1 do begin
z(i1)=long(m*csa+i1*sna+0.5) & y(i1)=long(-m*sna+i1*csa+0.5)
if (z(i1) ge 0) and (z(i1) le dimz) and (y(i1) ge 0) and (y(i1) lt 2*Radius) then begin
		z1(i2)=z(i1) & y1(i2)=y(i1)
			i2=i2+1
		endif
endfor
;print,'ray number',i3
if (i2 gt 0) then begin
; wset,1
; oplot,z1(0:i2-1),y1(0:i2-1);,xran=[-5,120],yran=[-5,90],/xst,/yst,xtit='Z',ytit='Y'
; wshow
length1=i2
Np3=fltarr(length1) & Te3=fltarr(length1) & Btot3=fltarr(length1) & Nb3=fltarr(length1)& angle3=fltarr(length1)
	Np3=Np2(y1[0:i2-1],z1[0:i2-1]) & Te3=te2(y1[0:i2-1],z1[0:i2-1]) 
	Btot3=Btot2(y1[0:i2-1],z1[0:i2-1]) & Nb3=Nb2(y1[0:i2-1],z1[0:i2-1])
  angle3=angle2(y1[0:i2-1],z1[0:i2-1])

;-------------sending parameters from selected pixels on the LOS to fast GS code
	DL_parameters,n_0=np3,T=Te3,B=Btot3,n_b=Nb3,angle=angle3,length=length1,Nf=Nf,Itot=Itot,Vtot=Vtot,f=f
	Int(i3,*,i)=Itot(*) ;[ray_num,Nf,dimt]
	Pol(i3,*,i)=Vtot(*) ;[ray_num,Nf,dimt]
	
wset,1 ;-------------plotting SG spectrum
plot, f, Int(i3,*,i), xtitle='Frequency, GHz', ytitle='Intensity, sfu' $
 ,tit=+'slice '+string(j,format='(i3.3)')+' time='+string(i,format='(i3.3)')+' z='+string(i3,format='(i3.3)') $
 ,psym=4,syms=1.2,/xlog,/ylog
 oplot, f, Int(i3,*,i)
endif
i3=i3+1
endfor
print,' time='+string(i,format='(i3.3)')
  endfor
 save,Int,Pol,f,filename=savedir+'GS_emission'+string(j,format='(i3.3)')+'x'+'.sav'
 print,'SLICE '+string(j,format='(i3.3)')+'  FINISHED'
    endfor
end

