pro advectionttdiv

; Define equilibrium (units: SI with length in Mm)

   m=1.
   harmonic = 1.
   dimr = 260.
   dimth = 180.
   dimt = 48.
   rad = 6.6

   norm = 1.e6
   mup = 1.25663706*1.e-6/norm
   kboltz = 1.380658*10^(-23.)
   proton = 1.67262158*10^(-27.)
   gmma = 5./3.

   ne_out_cgs = 1.6667e10
   ne_in_cgs = 5.e10
   te_out = 1.e6
   te_in = 5.e6
   
   Bx_out = 0.
   By_out = 0.   
   Bz_out = 0.0182234     ; (Tesla)
   Bx_in = 0.
   By_in = 0.
   Bz_in = 0.018
      
   ro = ne_in_cgs * 1.e6 * proton *norm^3
   re = ne_out_cgs * 1.e6 * proton *norm^3
   vao = sqrt(Bx_in^2+By_in^2+Bz_in^2)/sqrt(mup*ro)        ; (Mm per second)
   vae = sqrt(Bx_out^2+By_out^2+Bz_out^2)/sqrt(mup*re)
   bo = vao*sqrt(mup*ro)
   be = vae*sqrt(mup*re)
   
   A = 0.05     ; A represents amplitude normalised by loop radius, i.e. A/R
   L = 5.      ; L represents L/R
   kz = harmonic*!pi/(L*rad)    ; kz is [1 / Mm] not normalised with respect to radius
   
   wao = vao * kz
   wae = vae * kz

; Kink frequency

  w = sqrt((ro*wao^2+re*wae^2)/(ro+re))
  
  ke = sqrt((wae^2-w^2)/(vae^2))
  ki = sqrt((-wao^2+w^2)/(vao^2))
  
  Co = 2.*ro*(w^2-wao^2)/ki*A*rad          ; Specifies the displacement of the tube boundary at t=0 (v(r=R,t=0)*dt)
  Ce = ki*ke*rad^2/2.*Co

; Define grids

  r = findgen(dimr)/(dimr-1.)*2.9*rad
  r[0] = r[1]/1.1
  th = findgen(dimth)/(dimth)*2*!pi
  z = L/2.*rad
  P = 2*!pi/w
  t = findgen(dimt)/(dimt-1.)*P
  dt = t[1]-t[0]
  
  rpos = fltarr(dimr, dimth)
  thpos = fltarr(dimr, dimth)
  loc = fltarr(dimr, dimth)

; Eigenfunctions with advection

filename = 'test_advectionTT/eigft'+string(0,format="(i3.3)")+'.dat'
 OPENW, lun, filename, /GET_LUN
 printf, lun, 2
 printf, lun, 8812
 printf, lun, long(dimr*dimth)
 printf, lun, 6

; Time zero

for i=0, dimr-1 do begin
  for j = 0, dimth-1 do begin
      rpos[i,j] = r[i]
      thpos[i,j] = th[j]
      xpos = rpos[i,j]*cos(thpos[i,j])
      ypos = rpos[i,j]*sin(thpos[i,j])     
        if rpos[i,j] lt rad then begin  ; specify loc, calculate eigenf, save position and eigenf, advect
                  
          loc[i,j] = 0
          vri_t = Co*ki*w/(2*ro*(w^2-wao^2))*cos(th[j])*sin(kz*z)*cos(w*t[0])
          vti_t = -Co*ki*w/(2*ro*(w^2-wao^2))*sin(th[j])*sin(kz*z)*cos(w*t[0])
          vxi_t = vri_t * cos(th[j]) - vti_t*sin(th[j])
          vyi_t = vri_t * sin(th[j]) + vti_t*cos(th[j])
          rhi_t = Co*ki*r[i]/(2*vao^2)*cos(th[j])*sin(kz*z)*sin(w*t[0])+ro
          ;pi_t = gmma*po/ro*rhi_t         ; po = 0 because we assumed pressureless equilibrium
          ;ti_t = pi_t/rhi_t ; Possibly correct for scaling factor
          bri_t = bo*kz*Co*ki/(2*ro*(w^2-wao^2))*cos(th[j])*cos(kz*z)*sin(w*t[0])
          bti_t = -bo*kz*Co*ki/(2*ro*(w^2-wao^2))*sin(th[j])*cos(kz*z)*sin(w*t[0])
          bzi_t = bo*Co*ki*r[i]/(2*ro*vao^2)*cos(th[j])*sin(kz*z)*sin(w*t[0])+bo  ; + int d{vz}/dz dt
          
          
          printf, lun, xpos, ypos, vxi_t, vyi_t, rhi_t, bri_t, bti_t, bzi_t
          
          rpos[i,j] = rpos[i,j] + vri_t * dt
          thpos[i,j] = thpos[i,j] + vti_t * dt
          if rpos[i,j] lt 0 then begin
                rpos[i,j] = -rpos[i,j]
                thpos[i,j] = thpos[i,j] + !pi
          endif          
          
        endif else begin
          
          loc[i,j] = 1         
          vre_t = Ce*ke*w/(re*(w^2-wae^2))*dbeselk(ke*r[i],m)*cos(th[j])*sin(kz*z)*cos(w*t[0])
          vte_t = -Ce*w*beselk(ke*r[i],m,/double)/(r[i]*re*(w^2-wae^2))*sin(th[j])*sin(kz*z)*cos(w*t[0])
          ;xre_t = Ce*ke/(re*(w^2-wae^2))*dbeselk(ke*r[i],m)*cos(th[j])*sin(kz*z)*sin(w*t[0])
          ;xte_t = -Ce*beselk(ke*r[i],m,/double)/(r[i]*re*(w^2-wae^2))*sin(th[j])*sin(kz*z)*sin(w*t[0])
          vxe_t = vre_t * cos(th[j]) - vte_t*sin(th[j])
          vye_t = vre_t * sin(th[j]) + vte_t*cos(th[j])
          rhe_t = Ce*beselk(ke*r[i],m,/double)/vae^2*cos(th[j])*sin(kz*z)*sin(w*t[0])+re
          ;pe_t = gmma*pe/re*rhe_t
          ;te_t = pe_t/rhe_t ; Possibly correction needed for scaling factor
          bre_t = be*kz*Ce*ke/(re*(w^2-wae^2))*dbeselk(ke*r[i],m)*cos(th[j])*cos(kz*z)*sin(w*t[0])
          bte_t = -be*kz*Ce*beselk(ke*r[i],m,/double)/(r[i]*re*(w^2-wae^2))*sin(th[j])*cos(kz*z)*sin(w*t[0])
          bze_t = be*Ce*beselk(ke*r[i],m,/double)/(re*vae^2)*cos(th[j])*sin(kz*z)*sin(w*t[0])+be
          
          printf, lun, xpos, ypos, vxe_t, vye_t, rhe_t, bre_t, bte_t, bze_t

          rpos[i,j] = rpos[i,j] + vre_t * dt
          thpos[i,j] = thpos[i,j] + vte_t * dt          
          if rpos[i,j] lt 0 then begin
                rpos[i,j] = -rpos[i,j]
                thpos[i,j] = thpos[i,j] + !pi
          endif       
          
     endelse
  endfor
endfor

CLOSE, lun
FREE_LUN, lun

; Time >> 0

for n=1, dimt-1 do begin
  filename = 'test_advectionTT/eigft'+string(n,format="(i3.3)")+'.dat'
 OPENW, lun, filename, /GET_LUN
 printf, lun, 2
 printf, lun, 8812+n
 printf, lun, long(dimr*dimth)
 printf, lun, 6
  
for i=0, dimr-1 do begin
  for j = 0, dimth-1 do begin
      xpos = rpos[i,j]*cos(thpos[i,j])
      ypos = rpos[i,j]*sin(thpos[i,j])     
        if loc[i,j] eq 0 then begin  ; specify loc, calculate eigenf, save position and eigenf, advect
                  
          vri_t = Co*ki*w/(2*ro*(w^2-wao^2))*cos(th[j])*sin(kz*z)*cos(w*t[n])
          vti_t = -Co*ki*w/(2*ro*(w^2-wao^2))*sin(th[j])*sin(kz*z)*cos(w*t[n])
          xri_t = Co*ki/(2*ro*(w^2-wao^2))*cos(th[j])*sin(kz*z)*sin(w*t[n])
          xti_t = -Co*ki/(2*ro*(w^2-wao^2))*sin(th[j])*sin(kz*z)*sin(w*t[n])
          vxi_t = vri_t * cos(th[j]) - vti_t*sin(th[j])
          vyi_t = vri_t * sin(th[j]) + vti_t*cos(th[j])
          rhi_t = Co*ki*r[i]/(2*vao^2)*cos(th[j])*sin(kz*z)*sin(w*t[n])+ro
          ;pi_t = gmma*po/ro*rhi_t         ; po = 0 because we assumed pressureless equilibrium
          ;ti_t = pi_t/rhi_t ; Possibly correct for scaling factor
          bri_t = bo*kz*Co*ki*w/(2*ro*(w^2-wao^2))*cos(th[j])*cos(kz*z)*cos(w*t[n])
          bti_t = -bo*kz*Co*ki*w/(2*ro*(w^2-wao^2))*sin(th[j])*cos(kz*z)*cos(w*t[n])
          bzi_t = Co*ki*r[i]*w/(2*ro*vao^2)*cos(th[j])*sin(kz*z)*cos(w*t[n])+bo  ; + d{vz}/dz

          printf, lun, xpos, ypos, vxi_t, vyi_t, rhi_t, bri_t, bti_t, bzi_t

          rpos[i,j] = rpos[i,j] + vri_t * dt
          thpos[i,j] = thpos[i,j] + vti_t * dt
          if rpos[i,j] lt 0 then begin
                rpos[i,j] = -rpos[i,j]
                thpos[i,j] = thpos[i,j] + !pi
          endif
                    
        endif else begin
          
          vre_t = -A*rad*w*(ke*rad)^2*dbeselk(ke*r[i],m)*cos(th[j])*sin(kz*z)*cos(w*t[n])
          vte_t = A*rad*w*(ke*rad)*beselk(ke*r[i],m,/double)/(r[i]/rad)*sin(th[j])*sin(kz*z)*cos(w*t[n])
          xre_t = Ce*ke/(re*(w^2-wae^2))*dbeselk(ke*r[i],m)*cos(th[j])*sin(kz*z)*sin(w*t[n])
          xte_t = -Ce*beselk(ke*r[i],m,/double)/(r[i]*re*(w^2-wae^2))*sin(th[j])*sin(kz*z)*sin(w*t[n])
          vxe_t = vre_t * cos(th[j]) - vte_t*sin(th[j])
          vye_t = vre_t * sin(th[j]) + vte_t*cos(th[j])
          rhe_t = Ce*beselk(ke*r[i],m,/double)/vae^2*cos(th[j])*sin(kz*z)*sin(w*t[n])+re
          ;pe_t = gmma*pe/re*rhe_t
          ;te_t = pe_t/rhe_t ; Possibly correction needed for scaling factor
          bre_t = be*kz*Ce*ke/(re*(w^2-wae^2))*dbeselk(ke*r[i],m)*cos(th[j])*cos(kz*z)*sin(w*t[n])
          bte_t = -be*kz*Ce*beselk(ke*r[i],m,/double)/(r[i]*re*(w^2-wae^2))*sin(th[j])*cos(kz*z)*sin(w*t[n])
          bze_t = be*Ce*beselk(ke*r[i],m,/double)/(re*vae^2)*cos(th[j])*sin(kz*z)*sin(w*t[n])+be
          
          printf, lun, xpos, ypos, vxe_t, vye_t, rhe_t, bre_t, bte_t, bze_t

          rpos[i,j] = rpos[i,j] + vre_t * dt
          thpos[i,j] = thpos[i,j] + vte_t * dt          
          if rpos[i,j] lt 0 then begin
                rpos[i,j] = -rpos[i,j]
                thpos[i,j] = thpos[i,j] + !pi
          endif       
          
      endelse
    endfor
  endfor       

print, 'time step '+strtrim(string(n),1)
CLOSE, lun
FREE_LUN, lun

endfor

end