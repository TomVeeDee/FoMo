pro visualeigf

; Time step, choice of constant variable and resolution of vector field relative to grid resolution

t=6
planecte = 'z'
coarsefactorvf = 3       ; (specify x/y/zlevel if plane x/y/z = cte in if-then clause)
varnb = 3                ; 3=density, 4=temperature, 5=vx, 6=vy, 7=vz

;--------------------------------------------------------------------------
;--------------------------------------------------------------------------

eigfcube = '/users/cpa/sgijsen/fomo/version_patrick_nov13/examples/data/datacubes_cpp/largeba5104completecyl/loc_data_cube_kink'+string(t,format="(i3.3)")+'.dat'
nx = 204
ny = 204
nz = 102

openr,lun,eigfcube,/get_lun
header=lonarr(4)
readf,lun,header
data=fltarr(header[0]+header[3],4244832)
readf,lun,data
close,lun
free_lun,lun

if varnb eq 3 then ti = 'Density'
if varnb eq 4 then ti = 'Temperature'
if varnb eq 5 then ti = 'vx'
if varnb eq 6 then ti = 'vy'
if varnb eq 7 then ti = 'vz'

;--------------------------------------------------------------------------
;--------------------------------------------------------------------------

if planecte eq 'z' then begin
  zlevel = 95.

zvalue = data[2,zlevel]
datazlevel=where(data[2,*] eq zvalue)
nd = size(datazlevel)
lim = nd[1]

var=fltarr(nx,ny)
vx=fltarr(nx,ny)
vy=fltarr(nx,ny)
vxrd=fltarr(nx/coarsefactorvf, ny/coarsefactorvf)
vyrd=fltarr(nx/coarsefactorvf, ny/coarsefactorvf)

for i=0,lim-1 do begin
   var[(i-(i mod ny))/ny,i mod ny] = data[varnb,datazlevel[i]]
   vx[(i-(i mod ny))/ny,i mod ny] = data[5,datazlevel[i]]
   vy[(i-(i mod ny))/ny,i mod ny] = data[6,datazlevel[i]]
endfor

for i = 1, nx/coarsefactorvf do begin
  for j = 1, ny/coarsefactorvf do begin
    vxrd[i-1,j-1]=vx[coarsefactorvf*i-1,coarsefactorvf*j-1]
    vyrd[i-1,j-1]=vy[coarsefactorvf*i-1,coarsefactorvf*j-1]
  endfor
endfor

x = findgen(nx/coarsefactorvf)*coarsefactorvf
y = findgen(ny/coarsefactorvf)*coarsefactorvf

endif   

;-----------------------------------------------------------
;-----------------------------------------------------------

if planecte eq 'y' then begin
  ylevel = 71.
  
yvalue = data[1,ylevel*nz]
dataylevel=where(data[1,*] eq yvalue)
nd = size(dataylevel)
lim = nd[1]

var=fltarr(nx,nz)
vx=fltarr(nx,nz)
vz=fltarr(nx,nz)
vxrd=fltarr(nx/coarsefactorvf, nz/coarsefactorvf)
vzrd=fltarr(nx/coarsefactorvf, nz/coarsefactorvf)

for i=0,lim-1 do begin
   var[(i-(i mod nz))/nz,i mod nz] = data[varnb,dataylevel[i]]
   vx[(i-(i mod nz))/nz,i mod nz] = data[5,dataylevel[i]]
   vz[(i-(i mod nz))/nz,i mod nz] = data[7,dataylevel[i]]
endfor

for i = 1, nx/coarsefactorvf do begin
  for j = 1, nz/coarsefactorvf do begin
    vxrd[i-1,j-1]=vx[coarsefactorvf*i-1,coarsefactorvf*j-1]
    vzrd[i-1,j-1]=vz[coarsefactorvf*i-1,coarsefactorvf*j-1]
  endfor
endfor

z = findgen(nz/coarsefactorvf)*coarsefactorvf
x = findgen(nx/coarsefactorvf)*coarsefactorvf
  
endif

;---------------------------------------------------------------
;---------------------------------------------------------------
  
if planecte eq 'x' then begin
  xlevel = 45.

xvalue = data[0,xlevel*ny*nz]
dataxlevel=where(data[0,*] eq xvalue)
nd = size(dataxlevel)
lim = nd[1]

var=fltarr(ny,nz)
vx=fltarr(ny,nz)
vy=fltarr(ny,nz)
vxrd=fltarr(ny/coarsefactorvf, nz/coarsefactorvf)
vyrd=fltarr(ny/coarsefactorvf, nz/coarsefactorvf)

for i=0,lim-1 do begin
   var[(i-(i mod nz))/nz,i mod nz] = data[varnb,dataxlevel[i]]
   vx[(i-(i mod nz))/nz,i mod nz] = data[5,dataxlevel[i]]
   vy[(i-(i mod nz))/nz,i mod nz] = data[6,dataxlevel[i]]
endfor

for i = 1, ny/coarsefactorvf do begin
  for j = 1, nz/coarsefactorvf do begin
    vxrd[i-1,j-1]=vx[coarsefactorvf*i-1,coarsefactorvf*j-1]
    vyrd[i-1,j-1]=vy[coarsefactorvf*i-1,coarsefactorvf*j-1]
  endfor
endfor

z = findgen(nz/coarsefactorvf)*coarsefactorvf
y = findgen(ny/coarsefactorvf)*coarsefactorvf
  
  endif

; -----------------------------------------------------
; -------  Making contour plot and vector field -------
; -----------------------------------------------------

window, 0
minValue = Min(var)
maxValue = Max(var)
nLevels = 12
position =  [0.125, 0.125, 0.9, 0.800]
cbposition = [0.125, 0.865, 0.9, 0.895]

cgDisplay, 600, 500, Title= ti ;+ ' in horizontal plane z=' + strtrim(string(zlevel),1)
cgLoadCT, 33, NColors=nLevels, Bottom=1, CLIP=[30,255]

contourLevels = cgConLevels(var, NLevels=nLevels, MinValue=minValue)
   
cgContour, var, /Fill, Levels=contourLevels, C_Colors=Bindgen(nLevels)+1B, $
   /OutLine, Position=position, XTitle=xtitle, YTitle=ytitle
 
cgColorbar, NColors=nlevels, Bottom=1, Position=cbposition, $
   Range=[Float(Round(MinValue*1000)/1000.), Float(Round(MaxValue*1000)/1000.)], Divisions=nLevels, $
   Title=cbTitle, TLocation='Top'
   
window, 1
contour, var,xtitle='x',ytitle='y',nlevels=30,/fill,/xstyle,/ystyle,title=ti

if planecte eq 'z' then begin
  v = VECTOR(vxrd, vyrd, x, y)
endif
if planecte eq 'y' then begin
  v = VECTOR(vxrd, vzrd, x, z)
endif
if planecte eq 'x' then begin
  v = VECTOR(vxrd, vyrd, y,z)
endif

end