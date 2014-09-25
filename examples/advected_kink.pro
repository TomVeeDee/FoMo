pro advected_kink

; Advected eigenfunctions, for processing in FoMo-C.
; Created by Stief Gijsen, Aug/Sep 2014.
; -----------------------------------------------------
; 
; This routine calculates the eigenfunctions (ne, T, vx, vy, vz, [bx, by, bz]) corresponding to
; standing MHD waves in straight cylinders. The azimuthal (m) and longitudinal (k_z = n * pi / L)
; wave numbers can be selected by the user. The density inside and outside the cylinder is constant,
; with the radial density profile discontinuous across the tube boundary. The eigenfunctions are calculated
; for a single wave period. The grid is advected with the solution: after each time step, 
; the new particle positions are calculated using a first-order integration.
; 
; The eigenfunctions are calculated by solving the dispersion relation for linear MHD waves
; in straight cylinders, specifically equation (8b) in Edwin & Roberts (1983).
; They are built up using div v = C * R(r) * cos(th)*sin(kz*z)*cos(w*t).
; The two integration constants alpha, beta are linked by the dispersion relation.
; The remaining free constant (in our case alpha) is chosen such that the displacement vector
; xi_r (r=R, th=0, z=L/2, t=P/4) = (A/R) * R  (a fixed percentage of the radius which can
; be specified by the user, f.e. 5%).
; 
; The eigenfunctions are written to a file, in a format supported for processing with FoMo-C.
; 
; -------------------------------------------------------------------------------------------------
; 
; 
; Define equilibrium (SI units in Megameter)

   m=1. ;azimuthal wave number
   harmonic = 1 ; used in longitudinal wave number (standing waves)
   dimr = 40.
   dimth = 20.   
   dimz = 20.
   dimt = 12./harmonic
   rad = 1. ; loop radius in Mm

   norm = 1.e6 ; Conversion factor between SI and lengths in Mm
   mup = 1.25663706*1.e-6/norm ;magnetic permeability rescaled to Si_Mm
   kboltz = 1.380658*10^(-23.)
   proton = 1.67262158*10^(-27.); proton mass
   gmma = 5./3.

   ne_out_cgs = 1.e9 ;electron number density outside loop
   ne_in_cgs = 3.e9 ;inside
   te_out = 3.e6 ;temperature outside
   te_in = 1.e6 ;inside
   
   Bx_out = 0.
   By_out = 0.
   Bz_out = 0.00228 ;Magnetic field with loop axis = z-axis of coordinate system
   Bx_in = 0.
   By_in = 0.
   Bz_in = 0.00228 ;Inside
   ; Note B, T, and ne must be chosen such that total pressure equilibrium is satisfied
      
   ro = ne_in_cgs * 1.e6 * proton *norm^3
   re = ne_out_cgs * 1.e6 * proton *norm^3
   vao = sqrt(Bx_in^2+By_in^2+Bz_in^2)/sqrt(mup*ro) ; Alfven speed inside
   vae = sqrt(Bx_out^2+By_out^2+Bz_out^2)/sqrt(mup*re); Alfven speed outside
   bo = vao*sqrt(mup*ro) ; magnitude of magnetic field inside
   be = vae*sqrt(mup*re) ; magnitude of magnetic field outside
   co = sqrt(2*gmma*kboltz/proton*te_in)/norm ; Sound speed inside: from ideal gas law with ne = electron density (hence factor two in numerator)
   ce = sqrt(2*gmma*kboltz/proton*te_out)/norm; soond speed outside
   ct = co*vao/sqrt(co^2+vao^2) ; tube speed inside
   cte = ce*vae/sqrt(ce^2+vae^2) ; tube speed outside
   
   po=ro*co^2/gmma ;pressure inside
   pe=re*ce^2/gmma ; pressure outside
   
   A = 0.15     ; A represents amplitude normalised by loop radius, i.e. A/R
   L = 200.      ; L represents L/R
   kz = harmonic*!pi/(L*rad)    ; kz is [1 / Mm] not normalised with respect to loop radius
   
   wao = vao * kz; Alfven frequency inside
   wae = vae * kz; Alfven frequency outside


; Calculate eigenfrequency (Solve transcendental equation numerically)

  dim=1000
  Q = fltarr(dim+1)
  
   vmin = vao
   va_mx = vae
   wa = findgen(dim+1)/float(dim)*(va_mx-vmin)*kz+vmin*kz ;Create array of omega-values in which to look for zeroes of dispersion function

   nummoa2 = (kz^2*co^2-wa^2)*(kz^2*vao^2-wa^2)
   dnummoa2 = (kz^2*ct^2-wa^2)*(co^2+vao^2)
   ki_2 = nummoa2/dnummoa2 ;radial wave number inside (m_0^2 in ER1983)
   nummea2 = (kz^2*ce^2-wa^2)*(kz^2*vae^2-wa^2)
   dnummea2 = (kz^2*cte^2-wa^2)*(ce^2+vae^2)
   ke_2 = nummea2/dnummea2 ;radial wave number outside (m_e^2 in ER1983)
   
   reg3 = where(ki_2 lt 0. and ke_2 gt 0.) ;Coronal solutions (finite for r->0 and vanishing for r-> \infty)
   
  ; Dispersion function: 
  Q = wa[reg3]^2*(re*sqrt(-ki_2[reg3])*dbeselj(sqrt(-ki_2[reg3])*rad,m)*beselk(sqrt(ke_2[reg3])*rad,m,/double)-ro*sqrt(ke_2[reg3])* $
    dbeselk(sqrt(ke_2[reg3])*rad,m)*beselj(sqrt(-ki_2[reg3])*rad,m,/double))-(re*wae^2*sqrt(-ki_2[reg3])*dbeselj(sqrt(-ki_2[reg3])*rad,m)* $
    beselk(sqrt(ke_2[reg3])*rad,m,/double)-ro*wao^2*sqrt(ke_2[reg3])*dbeselk(sqrt(ke_2[reg3])*rad,m)*beselj(sqrt(-ki_2[reg3])*rad,m,/double))
  
  ; Visualise dispersion function  
  ; Q2 = Q/1000.
  ; plot,wa,Q2,/xs,/ys,psym=3,yr=[-0.5,0.5]
   
  wi = wa(where(abs(Q) eq min(abs(Q)))) ; Find smallest value within wa-array
  w=wi[0] ; Sometimes, more values of wa correspond to wi, choose arbitrary one = eigenfrequency

    
  ke = sqrt((kz^2*ce^2-w^2)*(kz^2*vae^2-w^2)/((kz^2*cte^2-w^2)*(ce^2+vae^2))) ;radial wave number with eigenfrequency w
  ki = sqrt(-(kz^2*co^2-w^2)*(kz^2*vao^2-w^2)/((kz^2*ct^2-w^2)*(co^2+vao^2))) ;radial wave number inside
  
  alpha = (w^2-wao^2)*w^3/((w^2-kz^2*ct^2)*(co^2+vao^2)*dbeselj(ki*rad,m)*ki)*A*rad ; Specifies the displacement of the tube boundary at t=P/4, z=L/2, r=R, phi=0.
  beta = beselj(ki*rad,m,/double)/beselk(ke*rad,m,/double)*alpha ; Relation between alpha and beta

;---------------
; Define grids 
;---------------

  r = findgen(dimr)/(dimr-1.)*2.9*rad
  r[0] = r[1]/1.1
  th = findgen(dimth)/(dimth)*2*!pi
  z = findgen(dimz)/(dimz-1.)*L*rad
  P = 2*!pi/w
  t = findgen(dimt)/(dimt)*P
  dt = t[1]-t[0]
  
  
; Keep track of positions of each particle during advection:
  rpos = fltarr(dimr, dimth, dimz)
  thpos = fltarr(dimr, dimth, dimz)
  zpos = fltarr(dimr, dimth, dimz)
  loc = fltarr(dimr, dimth, dimz) 
; loc = boolean indicator: if loc=0, particle is inside loop, for loc=1 particle is outside loop.
;Particles that are inside loop at tme zero, remain inside after advection


; ---------------------------------------------------------
; Eigenfunctions with advection
; Time zero

filename = '/users/cpa/sgijsen/FoMo/eigft/eigft'+string(0,format="(i3.3)")+'.dat'
 OPENW, lun, filename, /GET_LUN
; First print general information about grid
  printf, lun, 3 ;number of dimensions
  printf, lun, 8812+0 ;indicator of 'eqtype' in FoMo-C (any integer larger than 8 suffices)
  printf, lun, long(dimr*dimth*dimz) ;number of grid points
  printf, lun, 5 ;number of variables

for k=0, dimz-1 do begin
for i=0, dimr-1 do begin
  for j = 0, dimth-1 do begin
      rpos[i,j,k] = r[i]
      thpos[i,j,k] = th[j]
      zpos[i,j,k] = z[k]
      xpos = rpos[i,j,k]*cos(thpos[i,j,k])
      ypos = rpos[i,j,k]*sin(thpos[i,j,k])     
        if rpos[i,j,k] lt rad then begin  ; Particle inside loop, or loc=0
                  
          loc[i,j,k] = 0
          vri_t = alpha*ki*(w^2-kz^2*ct^2)*(co^2+vao^2)*dbeselj(ki*r[i],m)/((w^2-wao^2)*w^2)*cos(th[j])*sin(kz*z[k])*cos(w*t[0])
          vti_t = -alpha*(w^2-kz^2*ct^2)*(co^2+vao^2)/((w^2-wao^2)*w^2*r[i])*beselj(ki*r[i],m,/double)*sin(th[j])*sin(kz*z[k])*cos(w*t[0])
          vxi_t = vri_t * cos(th[j]) - vti_t*sin(th[j]) ;Convert vr and v_th to vx and vy
          vyi_t = vri_t * sin(th[j]) + vti_t*cos(th[j])
          vzi_t = alpha*kz*co^2/w^2*beselj(ki*r[i],m,/double)*cos(th[j])*cos(kz*z[k])*cos(w*t[0])
          rhi_t = (alpha*ro/w*beselj(ki*r[i],m,/double)*cos(th[j])*sin(kz*z[k])*sin(w*t[0])+ro)
          pi_t = po*(1+gmma/ro*rhi_t)         
          ti_t = pi_t/rhi_t 
          rhi_t = rhi_t/(1.e6 * proton *norm^3) ; Rescaling of density to electrons / cm^3
          ti_t = ti_t*proton/kboltz*norm^2*3./16. ; Rescaling of temperature
          ;bri_t = bo*kz*alpha*ki/w*(w^2-kz^2*ct^2)*(co^2+vao^2)*dbeselj(ki*r[i],m)/((w^2-wao^2)*w^2)*cos(th[j])*cos(kz*z[k])*sin(w*t[0])
          ;bti_t = -bo*kz*alpha/w*(w^2-kz^2*ct^2)*(co^2+vao^2)/((w^2-wao^2)*w^2*r[i])*beselj(ki*r[i],m,/double)*sin(th[j])*cos(kz*z[k])*sin(w*t[0])
          ;bzi_t = bo*alpha/w*(1-kz^2*co^2/w^2)*beselj(ki*r[i],m,/double)*cos(th[j])*sin(kz*z[k])*sin(w*t[0])+bo
                   
          printf, lun, xpos, ypos, zpos[i,j,k], rhi_t, ti_t, vxi_t, vyi_t, vzi_t
          ;printf, lun, bri_t
          ;printf, lun, bti_t
          ;printf, lun, bzi_t
          
          ; Advection of particles: change of radial, azimuthal, and longitudinal position
          
          rpos[i,j,k] = rpos[i,j,k] + vri_t * dt
          thpos[i,j,k] = thpos[i,j,k] + vti_t * dt
          if rpos[i,j,k] lt 0 then begin
                rpos[i,j,k] = -rpos[i,j,k]
                thpos[i,j,k] = thpos[i,j,k] + !pi
          endif
          zpos[i,j,k] = zpos[i,j,k] + vzi_t * dt 
          
        endif else begin ; Particle outside loop, i.e. loc=1
          
          loc[i,j,k] = 1         
          vre_t = beta*ke*(w^2-kz^2*cte^2)*(ce^2+vae^2)*dbeselk(ke*r[i],m)/((w^2-wae^2)*w^2)*cos(th[j])*sin(kz*z[k])*cos(w*t[0])
          vte_t = -beta*(w^2-kz^2*cte^2)*(ce^2+vae^2)/((w^2-wae^2)*w^2*r[i])*beselk(ke*r[i],m,/double)*sin(th[j])*sin(kz*z[k])*cos(w*t[0])
          vxe_t = vre_t * cos(th[j]) - vte_t*sin(th[j])
          vye_t = vre_t * sin(th[j]) + vte_t*cos(th[j])
          vze_t = beta*kz*ce^2/w^2*beselk(ke*r[i],m,/double)*cos(th[j])*cos(kz*z[k])*cos(w*t[0])
          rhe_t = beta*re/w*beselk(ke*r[i],m,/double)*cos(th[j])*sin(kz*z[k])*sin(w*t[0])+re
          pe_t = pe*(1+gmma/re*rhe_t)         
          te_t = pe_t/rhe_t
          rhe_t = rhe_t/(1.e6 * proton *norm^3)
          te_t = te_t*proton/kboltz*norm^2*3./16.
          ;bre_t = be*kz*beta*ke/w*(w^2-kz^2*cte^2)*(ce^2+vae^2)*dbeselk(ke*r[i],m)/((w^2-wae^2)*w^2)*cos(th[j])*cos(kz*z[k])*sin(w*t[0])
          ;bte_t = -be*kz*beta/w*(w^2-kz^2*cte^2)*(ce^2+vae^2)/((w^2-wae^2)*w^2*r[i])*beselk(ke*r[i],m,/double)*sin(th[j])*cos(kz*z[k])*sin(w*t[0])
          ;bze_t = be*beta/w*(1-kz^2*ce^2/w^2)*beselk(ke*r[i],m,/double)*cos(th[j])*sin(kz*z[k])*sin(w*t[0])+be
          
          
          printf, lun, xpos, ypos, zpos[i,j,k], rhe_t, te_t, vxe_t, vye_t, vze_t
          ;printf, lun, bre_t
          ;printf, lun, bte_t
          ;printf, lun, bze_t

          rpos[i,j,k] = rpos[i,j,k] + vre_t * dt
          thpos[i,j,k] = thpos[i,j,k] + vte_t * dt          
          if rpos[i,j,k] lt 0 then begin
                rpos[i,j,k] = -rpos[i,j,k]
                thpos[i,j,k] = thpos[i,j,k] + !pi
          endif
          zpos[i,j,k] = zpos[i,j,k] + vze_t * dt
          
     endelse
  endfor
endfor
endfor

print, 'time step '+strtrim(string(0),1)
CLOSE, lun
FREE_LUN, lun

; -----------------------------------------
; Time > 0
; -----------------------------------------

for n=1, dimt-1 do begin
  filename = '/users/cpa/sgijsen/FoMo/eigft/eigft'+string(n,format="(i3.3)")+'.dat'
  OPENW, lun, filename, /GET_LUN
  printf, lun, 3
  printf, lun, 8812+n
  printf, lun, long(dimr*dimth*dimz)
  printf, lun, 5

  
  for k=0, dimz-1 do begin
  for i=0, dimr-1 do begin
  for j = 0, dimth-1 do begin
      xpos = rpos[i,j,k]*cos(thpos[i,j,k])
      ypos = rpos[i,j,k]*sin(thpos[i,j,k])     
        if loc[i,j,k] eq 0 then begin  ; specify loc, calculate eigenf, save position and eigenf, advect
                  
          vri_t = alpha*ki*(w^2-kz^2*ct^2)*(co^2+vao^2)*dbeselj(ki*r[i],m)/((w^2-wao^2)*w^2)*cos(th[j])*sin(kz*z[k])*cos(w*t[n])
          vti_t = -alpha*(w^2-kz^2*ct^2)*(co^2+vao^2)/((w^2-wao^2)*w^2*r[i])*beselj(ki*r[i],m,/double)*sin(th[j])*sin(kz*z[k])*cos(w*t[n])
          vxi_t = vri_t * cos(th[j]) - vti_t*sin(th[j])
          vyi_t = vri_t * sin(th[j]) + vti_t*cos(th[j])
          vzi_t = alpha*kz*co^2/w^2*beselj(ki*r[i],m,/double)*cos(th[j])*cos(kz*z[k])*cos(w*t[n])
          rhi_t = alpha*ro/w*beselj(ki*r[i],m,/double)*cos(th[j])*sin(kz*z[k])*sin(w*t[n])+ro
          pi_t = po*(1+gmma/ro*rhi_t)         
          ti_t = pi_t/rhi_t
          rhi_t = rhi_t/(1.e6 * proton *norm^3)
          ti_t = ti_t*proton/kboltz*norm^2*3./16.
          ;bri_t = bo*kz*ki*alpha/w*(w^2-kz^2*ct^2)*(co^2+vao^2)*dbeselj(ki*r[i],m)/((w^2-wao^2)*w^2)*cos(th[j])*cos(kz*z[k])*sin(w*t[n])
          ;bti_t = -bo*kz*alpha/w*(w^2-kz^2*ct^2)*(co^2+vao^2)/((w^2-wao^2)*w^2*r[i])*beselj(ki*r[i],m,/double)*sin(th[j])*cos(kz*z[k])*sin(w*t[n])
          ;bzi_t = bo*alpha/w*(1-kz^2*co^2/w^2)*beselj(ki*r[i],m,/double)*cos(th[j])*sin(kz*z[k])*sin(w*t[n])+bo
                   
          printf, lun, xpos, ypos, zpos[i,j,k], rhi_t, ti_t, vxi_t, vyi_t, vzi_t
          ;printf, lun, bri_t
          ;printf, lun, bti_t
          ;printf, lun, bzi_t


          rpos[i,j,k] = rpos[i,j,k] + vri_t * dt
          thpos[i,j,k] = thpos[i,j,k] + vti_t * dt
          if rpos[i,j,k] lt 0 then begin
                rpos[i,j,k] = -rpos[i,j,k]
                thpos[i,j,k] = thpos[i,j,k] + !pi
          endif
          zpos[i,j,k] = zpos[i,j,k] + vzi_t * dt 
          
                    
        endif else begin
          
          vre_t = beta*ke*(w^2-kz^2*cte^2)*(ce^2+vae^2)*dbeselk(ke*r[i],m)/((w^2-wae^2)*w^2)*cos(th[j])*sin(kz*z[k])*cos(w*t[n])
          vte_t = -beta*(w^2-kz^2*cte^2)*(ce^2+vae^2)/((w^2-wae^2)*w^2*r[i])*beselk(ke*r[i],m,/double)*sin(th[j])*sin(kz*z[k])*cos(w*t[n])
          vxe_t = vre_t * cos(th[j]) - vte_t*sin(th[j])
          vye_t = vre_t * sin(th[j]) + vte_t*cos(th[j])
          vze_t = beta*kz*ce^2/w^2*beselk(ke*r[i],m,/double)*cos(th[j])*cos(kz*z[k])*cos(w*t[n])
          rhe_t = beta*re/w*beselk(ke*r[i],m,/double)*cos(th[j])*sin(kz*z[k])*sin(w*t[n])+re
          pe_t = pe*(1+gmma/re*rhe_t)         
          te_t = pe_t/rhe_t
          rhe_t = rhe_t/(1.e6 * proton *norm^3)
          te_t = te_t*proton/kboltz*norm^2*3./16.
          ;bre_t = be*kz*ke*beta/w*(w^2-kz^2*cte^2)*(ce^2+vae^2)*dbeselk(ke*r[i],m)/((w^2-wae^2)*w^2)*cos(th[j])*cos(kz*z[k])*sin(w*t[n])
          ;bte_t = -be*kz*beta/w*(w^2-kz^2*cte^2)*(ce^2+vae^2)/((w^2-wae^2)*w^2*r[i])*beselk(ke*r[i],m,/double)*sin(th[j])*cos(kz*z[k])*sin(w*t[n])
          ;bze_t = be*beta/w*(1-kz^2*ce^2/w^2)*beselk(ke*r[i],m,/double)*cos(th[j])*sin(kz*z[k])*sin(w*t[n])+be
             
          printf, lun, xpos, ypos, zpos[i,j,k], rhe_t, te_t, vxe_t, vye_t, vze_t
          ;printf, lun, bre_t
          ;printf, lun, bte_t
          ;printf, lun, bze_t
          
          rpos[i,j,k] = rpos[i,j,k] + vre_t * dt
          thpos[i,j,k] = thpos[i,j,k] + vte_t * dt          
          if rpos[i,j,k] lt 0 then begin
                rpos[i,j,k] = -rpos[i,j,k]
                thpos[i,j,k] = thpos[i,j,k] + !pi
          endif
          zpos[i,j,k] = zpos[i,j,k] + vze_t * dt
          
      endelse
    endfor
  endfor       
  endfor 
CLOSE, lun
FREE_LUN, lun 
print, 'time step '+strtrim(string(n),1)
endfor

end