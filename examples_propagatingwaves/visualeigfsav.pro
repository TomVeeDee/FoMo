pro visualeigfsav

; Specify time step in simulation, choice of plane and resolution of vector field relative to grid resolution

t=9
planecte = 'x'           ; Choice between slices of the form x=x0 ('x'), y=y0 ('y') and z=z0 ('z') (Determine x0,y0,z0 below)
coarsefactorvf = 6
varnb = 8            ; 3=density, 4=temperature, 5=bx, 6=by, 7=bz, 8=br, 9=vr, 10=vy, 11=vz, 12=vr

restore, '/users/cpa/sgijsen/FoMo/examples_propagatingwaves/cubes_B_sausage_all_ka2.24_'+string(t,format="(i3.3)")+'.sav'
restore, '/users/cpa/sgijsen/FoMo/examples_propagatingwaves/variables_propagating_sausage_0.012beta.sav'
restore, '/users/cpa/sgijsen/FoMo/examples_propagatingwaves/params_ka2.24.sav'

;--------------------------------------------------------------------------
;--------------------------------------------------------------------------

nx = dimx
ny = dimy
nz = dimz

; Select data to be plotted --- v_theta=0 for sausage mode

bx_cube=br_cube*(x_cube-r0)/dist_cube
by_cube=br_cube*(y_cube-r0)/dist_cube
vx_cube=vr_cube*(x_cube-r0)/dist_cube
vy_cube=vr_cube*(y_cube-r0)/dist_cube

if varnb eq 3 then begin 
  ti = 'Density'
  data = rh_cube
endif

if varnb eq 4 then begin 
  ti = 'Temperature'
  data = te_cube
endif

if varnb eq 5 then begin 
  ti = 'bx'
  data = bx_cube   ; -vt_cube*y_cube  (for kink mode)
endif

if varnb eq 6 then begin 
  ti = 'by'
  data = by_cube    ; +vt_cube*x_cube
endif

if varnb eq 7 then begin 
  ti = 'bz'
  data = bz_cube
endif

if varnb eq 8 then begin
  ti = 'br'
  data = br_cube
endif

if varnb eq 9 then begin 
  ti = 'vx'
  data = vx_cube   ; -vt_cube*y_cube  (for kink mode)
endif

if varnb eq 10 then begin 
  ti = 'vy'
  data = vy_cube    ; +vt_cube*x_cube
endif

if varnb eq 11 then begin 
  ti = 'vz'
  data = vz_cube
endif

if varnb eq 12 then begin
  ti = 'vr'
  data = vr_cube
endif

;--------------------------------------------------------------------------
;------------------------Make contour plot and vector field ---------------
;--------------------------------------------------------------------------


if planecte eq 'z' then begin
  zlevel = 66.
  data = reform(data[*,*,zlevel])
  
  bx=reform(vx_cube[*,*,zlevel])
  by=reform(by_cube[*,*,zlevel])
  bxrd=fltarr(nx/coarsefactorvf, ny/coarsefactorvf)
  byrd=fltarr(nx/coarsefactorvf, ny/coarsefactorvf)

  for i = 1, nx/coarsefactorvf do begin
    for j = 1, ny/coarsefactorvf do begin
      bxrd[i-1,j-1]=bx[coarsefactorvf*i-1,coarsefactorvf*j-1]
      byrd[i-1,j-1]=by[coarsefactorvf*i-1,coarsefactorvf*j-1]
    endfor
  endfor

  x = findgen(nx/coarsefactorvf)*coarsefactorvf
  y = findgen(ny/coarsefactorvf)*coarsefactorvf
  
  
  window, 0
  minValue = Min(data)
  maxValue = Max(data)
  nLevels = 10
  position =  [0.125, 0.125, 0.9, 0.800]
  cbposition = [0.125, 0.865, 0.9, 0.895]

  cgDisplay, 600, 500, Title= ti + ' in horizontal plane z=' + strtrim(string(zlevel),1)
  cgLoadCT, 33, NColors=nLevels, Bottom=1, CLIP=[30,255]

  contourLevels = cgConLevels(data, NLevels=nLevels, MinValue=minValue)
   
  cgContour, data, /Fill, Levels=contourLevels, C_Colors=Bindgen(nLevels)+1B, $
   /OutLine, Position=position, XTitle=xtitle, YTitle=ytitle
 
  cgColorbar, NColors=nlevels, Bottom=1, Position=cbposition, $
   Range=[Float(Round(MinValue*1000)/1000.), Float(Round(MaxValue*10)/10.)], Divisions=nLevels, $
   Title=cbTitle, TLocation='Top'
  
  ;window,1
  ;contour, data[*,*,zlevel],xtitle='x',ytitle='y',nlevels=30,/fill,/xstyle,/ystyle,title=ti
  v = VECTOR(bxrd, byrd, x, y)
  
  print, minValue
  
endif   


if planecte eq 'y' then begin
  ylevel = 71.
  data=reform(data[*,ylevel,*])
  
  bx=reform(bx_cube[*,ylevel,*])
  bz=reform(bz_cube[*,ylevel,*])
  bxrd=fltarr(nx/coarsefactorvf, nz/coarsefactorvf)
  bzrd=fltarr(nx/coarsefactorvf, nz/coarsefactorvf)

  for i = 1, nx/coarsefactorvf do begin
    for j = 1, nz/coarsefactorvf do begin
      vxrd[i-1,j-1]=vx[coarsefactorvf*i-1,coarsefactorvf*j-1]
      vzrd[i-1,j-1]=vz[coarsefactorvf*i-1,coarsefactorvf*j-1]
    endfor
  endfor

  x = findgen(nx/coarsefactorvf)*coarsefactorvf
  z = findgen(nz/coarsefactorvf)*coarsefactorvf
  
  ;contour, reform(data[*,ylevel,*]),xtitle='x',ytitle='z',nlevels=30,/fill,/xstyle,/ystyle,title=ti
  
  window, 0
  minValue = Min(data)
  maxValue = Max(data)
  nLevels = 10
  position =  [0.125, 0.125, 0.9, 0.800]
  cbposition = [0.125, 0.865, 0.9, 0.895]

  cgDisplay, 600, 500, Title= ti + ' in plane y=' + strtrim(string(ylevel),1)
  cgLoadCT, 33, NColors=nLevels, Bottom=1, CLIP=[30,255]

  contourLevels = cgConLevels(data, NLevels=nLevels, MinValue=minValue)
   
  cgContour, data, /Fill, Levels=contourLevels, C_Colors=Bindgen(nLevels)+1B, $
   /OutLine, Position=position, XTitle=xtitle, YTitle=ytitle
 
  cgColorbar, NColors=nlevels, Bottom=1, Position=cbposition, $
   Range=[Float(Round(MinValue*1000)/1000.), Float(Round(MaxValue*1000)/1000.)], Divisions=nLevels, $
   Title=cbTitle, TLocation='Top'
    
  v = VECTOR(bxrd, bzrd, x,z)
endif

  
if planecte eq 'x' then begin
  xlevel = 102.
  data=reform(data[xlevel,*,*])
  
  by=reform(by_cube[xlevel,*,*])
  bz=reform(bz_cube[xlevel,*,*])
  byrd=fltarr(ny/coarsefactorvf, nz/coarsefactorvf)
  bzrd=fltarr(ny/coarsefactorvf, nz/coarsefactorvf)

  for i = 1, ny/coarsefactorvf do begin
    for j = 1, nz/coarsefactorvf do begin
      byrd[i-1,j-1]=by[coarsefactorvf*i-1,coarsefactorvf*j-1]
      bzrd[i-1,j-1]=bz[coarsefactorvf*i-1,coarsefactorvf*j-1]
    endfor
  endfor

  y = findgen(ny/coarsefactorvf)*coarsefactorvf
  z = findgen(nz/coarsefactorvf)*coarsefactorvf
  
  ;contour, data[xlevel,*,*],xtitle='y',ytitle='z',nlevels=30,/fill,/xstyle,/ystyle,title=ti
  window, 0
  minValue = Min(data)
  maxValue = Max(data)
  nLevels = 10
  position =  [0.125, 0.125, 0.9, 0.800]
  cbposition = [0.125, 0.865, 0.9, 0.895]

  cgDisplay, 600, 500, Title= ti + ' in plane x=' + strtrim(string(xlevel),1)
  cgLoadCT, 33, NColors=nLevels, Bottom=1, CLIP=[30,255]

  contourLevels = cgConLevels(data, NLevels=nLevels, MinValue=minValue)
   
  cgContour, data, /Fill, Levels=contourLevels, C_Colors=Bindgen(nLevels)+1B, $
   /OutLine, Position=position, XTitle=xtitle, YTitle=ytitle
 
  cgColorbar, NColors=nlevels, Bottom=1, Position=cbposition, $
   Range=[Float(Round(MinValue*1000000.)/1000000.), Float(Round(MaxValue*1)/1.)], Divisions=nLevels, $
   Title=cbTitle, TLocation='Top'
  
  v = VECTOR(byrd, bzrd, y,z)
endif

end
