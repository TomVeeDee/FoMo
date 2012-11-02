
pro rotmat, dimx, dimy, dimz, ro, re, va, vae, co, ce, waka_ini, no_init_guess = no_init, high_t_ext = high_t

if n_params(0) lt 1 then begin
   print,'rotmat, dimx, dimy, dimz, ro, re, va, vae, co, ce, [no_init_guess = no_guess, high_ext_t = high_t]'
   return
endif

; Reads inital condition values from simulation run by Gruszecki et
; al. 2008, returns the parameters below and calculates an initial
; guess for solving the dispersion relation of the system.

; RETURNED VALUES:
; dimx = dimension along x axis of structure
; dimy = dimension along y axis of structure
; dimz = dimension along z axis of structure
; ro = internal density: rho[dimx/2,dimy/2,dimz/2]
; re = external density: rho[0,0,0]
; va = internal Alfven velocity: calculated at (dimx/2,dimy/2,dimz/2)
; vae = external Alfven velocity: calculated at (0,0,0)
; co = internal sound speed: calculated at (dimx/2,dimy/2,dimz/2)
; ce = external sound speed: calculated at (0,0,0)
; waka_ini = initial guess of w/k for solving dispersion relation
; OPTIONAL:
; if keyword 'no_init' is set then calculation of initial guess is avoided
; if keyword 'high_t' is set then external temperature is increased to
; 8*10^5 K

if keyword_set(no_init) then initial = 0 else initial = 1
if keyword_set(high_t) then val = 1 else val = 0

RESTORE,'/users/cpa/tomvd/data/forwardmodelling/sausage_marcin2/imagesave0000integratey.sav'
;RESTORE,'./imagesave0000integratey.sav' 
RESTORE,'/users/cpa/tomvd/data/forwardmodelling/sausage_marcin2/forpatrick.sav'
;RESTORE,'./forpatrick.sav'
;RESTORE,'/users/cpa/pantolin/Modeling/Bfield.sav'
RESTORE,'/users/cpa/pantolin/Modeling/Bfield_s.sav'

rho = data.rho
rho = reform(rho[0:201:20,*,0:100:2,0])
temperature = data.temperature
temperature = reform(temperature[0:201:20,*,0:100:2,0])
vx = vx[0:201:20,*,0:100:2]
vy = vy[0:201:20,*,0:100:2]
vz = vz[0:201:20,*,0:100:2]
proton=1.67262158*10^(-27.) ; kg
kboltz = 1.380658*10^(-23.) ; m^2 kg s^-2
;mu = 1.27 
;pe = rho/(mu*proton)*kboltz*temperature
pe = 2*rho/(proton)*kboltz*temperature
;bx = B[0].bx
;by = B[0].by
;bz = B[0].bz
;bx = bx[0:201:20,*,0:100:2]
;by = by[0:201:20,*,0:100:2]
;bz = bz[0:201:20,*,0:100:2]

sizes0=size(emission)

nx0 = sizes0[1]
ny0 = sizes0[2]
nz0 = sizes0[3]
nw0 = sizes0[4]
remission=fltarr(nz0,ny0,nx0,nw0)
rrho = fltarr(nz0,ny0,nx0)
rtemperature = fltarr(nz0,ny0,nx0)
rpe = fltarr(nz0,ny0,nx0)

for l=0,nw0-1 do begin
   for i=0,nx0-1 do remission[*,*,i,l]=transpose(reform(emission[i,*,*,l]))
endfor

sizes=size(remission)

nx = sizes[1]
ny = sizes[2]
nz = sizes[3]
nw = sizes[4]

sizv = size(vx)

nxv = sizv[3]
nyv = sizv[2]
nzv = sizv[1]
rvx = fltarr(nxv,nyv,nzv)
rvy = fltarr(nxv,nyv,nzv)
rvz = fltarr(nxv,nyv,nzv)
rbx = fltarr(nxv,nyv,nzv)
rby = fltarr(nxv,nyv,nzv)
rbz = fltarr(nxv,nyv,nzv)

for i=0,nzv-1 do begin
   rvx[*,*,i]=transpose(reform(vx[i,*,*]))
   rvy[*,*,i]=transpose(reform(vy[i,*,*]))
   rvz[*,*,i]=transpose(reform(vz[i,*,*]))
   rrho[*,*,i]=transpose(reform(rho[i,*,*]))
   rtemperature[*,*,i]=transpose(reform(temperature[i,*,*]))
   rpe[*,*,i]=transpose(reform(pe[i,*,*]))
   rbx[*,*,i]=transpose(reform(bx[i,*,*]))
   rby[*,*,i]=transpose(reform(by[i,*,*]))
   rbz[*,*,i]=transpose(reform(bz[i,*,*]))
endfor

gridx = findgen(nx)
gridy = findgen(ny)
gridz = findgen(nz)

if initial eq 1 then begin

   dim = 10000
   ka = 4.
   norm = 1.e5
   wa = findgen(dim+1)/float(dim)*3.*ka*1.e5/norm+6.*ka*1.e5/norm
   L = fltarr(dim+1)
   moa2 = fltarr(dim+1)
   mea2 = fltarr(dim+1)
   bslij = fltarr(dim+1)
   bslky = fltarr(dim+1)
   dbslij = fltarr(dim+1)
   dbslky = fltarr(dim+1)
   fracij = fltarr(dim+1)
   fracky = fltarr(dim+1)
   sigi = intarr(dim+1)
   sigk = intarr(dim+1)
   n = 0
   mup = 1.25663706*1.e-6
   gamma = 5./3.

   dimx = nx
   dimy = ny
   dimz = nz
   xm = round(dimx/2)
   ym = round(dimy/2)
   zm = round(dimz/2)

   if val eq 0 then begin
      rho_out = rrho[0,0,0]
      te_out = rtemperature[0,0,0]
      rho_in = rrho[xm,ym,zm]
      te_in = rtemperature[xm,ym,zm]
   endif else begin
      rho_out = rrho[0,0,0]
      te_out = double(8.e5)
      te_in = rtemperature[xm,ym,zm]
      rho_in = proton/(2.*kboltz)*(2*kboltz/proton*rho_out*te_out+((rBx[0,0,0]^2+rBy[0,0,0]^2+rBz[0,0,0]^2)-(rBx[xm,ym,zm]^2+rBy[xm,ym,zm]^2+rBz[xm,ym,zm]^2))/(2*mup))/te_in
   endelse

   ro = rho_in*norm^2
   re = rho_out*norm^2
   va = sqrt(rBx[xm,ym,zm]^2+rBy[xm,ym,zm]^2+rBz[xm,ym,zm]^2)/sqrt(mup*ro)
   vae = sqrt(rBx[0,0,0]^2+rBy[0,0,0]^2+rBz[0,0,0]^2)/sqrt(mup*re)
   co = sqrt(2*gamma*kboltz/proton*te_in)/1.e5
   ce = sqrt(2*gamma*kboltz/proton*te_out)/1.e5
   ct = co*va/sqrt(co^2+va^2)
   cte = ce*vae/sqrt(ce^2+vae^2)

   nummoa2 = (ka^2*co^2-wa^2)*(ka^2*va^2-wa^2)
   dnummoa2 = (ka^2*ct^2-wa^2)*(co^2+va^2)
   moa2 = nummoa2/dnummoa2
   nummea2 = (ka^2*ce^2-wa^2)*(ka^2*vae^2-wa^2)
   dnummea2 = (ka^2*cte^2-wa^2)*(ce^2+vae^2)
   mea2 = nummea2/dnummea2

   reg1 = where(moa2 gt 0. and mea2 gt 0.)
   reg2 = where(moa2 gt 0. and mea2 lt 0.)
   reg3 = where(moa2 lt 0. and mea2 gt 0.)
   reg4 = where(moa2 lt 0. and mea2 lt 0.)

   if reg1[0] ne -1 then begin
      bslij[reg1] = beseli(sqrt(abs(moa2[reg1])),n,/double)
      bslky[reg1] = beselk(sqrt(abs(mea2[reg1])),n,/double)
      dbslij[reg1] = dbeseli(sqrt(abs(moa2[reg1])),n)
      dbslky[reg1] = dbeselk(sqrt(abs(mea2[reg1])),n)
      fracij[reg1] = dbslij[reg1]/bslij[reg1]
      fracky[reg1] = dbslky[reg1]/bslky[reg1]
      sigi[reg1] = 1.
      sigk[reg1] = 1.
   endif
   if reg2[0] ne -1 then begin
      bslij[reg2] = beseli(sqrt(abs(moa2[reg2])),n,/double)
      bslky[reg2] = besely(sqrt(abs(mea2[reg2])),n,/double)
      dbslij[reg2] = dbeseli(sqrt(abs(moa2[reg2])),n)
      dbslky[reg2] = dbesely(sqrt(abs(mea2[reg2])),n)
      fracij[reg2] = dbslij[reg2]/bslij[reg2]
      fracky[reg2] = dbslky[reg2]/bslky[reg2]
      sigi[reg2] = 1.
      sigk[reg2] = -1.
   endif
   if reg3[0] ne -1 then begin
      bslij[reg3] = beselj(sqrt(abs(moa2[reg3])),n,/double)
      bslky[reg3] = beselk(sqrt(abs(mea2[reg3])),n,/double)
      dbslij[reg3] = dbeselj(sqrt(abs(moa2[reg3])),n)
      dbslky[reg3] = dbeselk(sqrt(abs(mea2[reg3])),n)
      fracij[reg3] = dbslij[reg3]/bslij[reg3]
      fracky[reg3] = dbslky[reg3]/bslky[reg3]
      sigi[reg3] = -1.
      sigk[reg3] = 1.
   endif
   if reg4[0] ne -1 then begin
      bslij[reg4] = beselj(sqrt(abs(moa2[reg4])),n,/double)
      bslky[reg4] = beselk(sqrt(abs(mea2[reg4])),n,/double)
;   bslky[reg4] = besely(sqrt(abs(mea2[reg4])),n,/double)
      dbslij[reg4] = dbeselj(sqrt(abs(moa2[reg4])),n)
      dbslky[reg4] = dbeselk(sqrt(abs(mea2[reg4])),n)
;   dbslky[reg4] = dbesely(sqrt(abs(mea2[reg4])),n)
      fracij[reg4] = dbslij[reg4]/bslij[reg4]
      fracky[reg4] = dbslky[reg4]/bslky[reg4]
      sigi[reg4] = -1.
      sigk[reg4] = -1.
   endif

   if reg1[0] ne -1 then L[reg1] = wa[reg1]^2*(sigi[reg1]*re*sqrt(abs(moa2[reg1]))*dbslij[reg1]*bslky[reg1]-ro*sigk[reg1]*sqrt(abs(mea2[reg1]))*dbslky[reg1]*bslij[reg1])-$
    ka^2*(re*vae^2*sigi[reg1]*sqrt(abs(moa2[reg1]))*dbslij[reg1]*bslky[reg1]-ro*va^2*sigk[reg1]*sqrt(abs(mea2[reg1]))*dbslky[reg1]*bslij[reg1])
   if reg2[0] ne -1 then L[reg2] = wa[reg2]^2*(sigi[reg2]*re*sqrt(abs(moa2[reg2]))*dbslij[reg2]*bslky[reg2]-sigk[reg2]*ro*sqrt(abs(mea2[reg2]))*dbslky[reg2]*bslij[reg2])-$
    ka^2*(re*vae^2*sigi[reg2]*sqrt(abs(moa2[reg2]))*dbslij[reg2]*bslky[reg2]-ro*va^2*sigk[reg2]*sqrt(abs(mea2[reg2]))*dbslky[reg2]*bslij[reg2])
   if reg3[0] ne -1 then L[reg3] = wa[reg3]^2*(re*sigi[reg3]*sqrt(abs(moa2[reg3]))*dbslij[reg3]*bslky[reg3]-ro*sigk[reg3]*sqrt(abs(mea2[reg3]))*dbslky[reg3]*bslij[reg3])-$
    ka^2*(re*vae^2*sigi[reg3]*sqrt(abs(moa2[reg3]))*dbslij[reg3]*bslky[reg3]-ro*va^2*sigk[reg3]*sqrt(abs(mea2[reg3]))*dbslky[reg3]*bslij[reg3])
   if reg4[0] ne -1 then L[reg4] = wa[reg4]^2*(re*sigi[reg4]*sqrt(abs(moa2[reg4]))*dbslij[reg4]*bslky[reg4]-ro*sigk[reg4]*sqrt(abs(mea2[reg4]))*dbslky[reg4]*bslij[reg4])-$
    ka^2*(re*vae^2*sigi[reg4]*sqrt(abs(moa2[reg4]))*dbslij[reg4]*bslky[reg4]-ro*va^2*sigk[reg4]*sqrt(abs(mea2[reg4]))*dbslky[reg4]*bslij[reg4])

   window,0
;plot,wa/ka,L/max(sqrt(abs(dnummoa2*dnummea2))),/xs,/ys,psym=3;,xr=[2,5],yr=[-1,1]
;locs = where(abs(L/max(sqrt(abs(dnummoa2*dnummea2)))) lt 1.e-5)

   plot,wa/ka,L,/xs,/ys,psym=3  ;,yr=[-0.01,0.01]
   locs = where(abs(L) lt 1.e-3)
   if locs[0] ne -1 then begin
      waka_ini = max(wa[locs]/ka[locs])
      print,'waka_ini = ', waka_ini 
   endif else begin
      print,'no solution'
   endelse

endif
;stop
end
