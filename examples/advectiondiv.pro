pro advectiondiv

; Eigenfunctions are built up from div v = C * R(r) * cos(th)*sin(kz*z)*cos(w*t)
; Constants C = Ci, Ce are linked by the dispersion relation
; The remaining free constant is chosen such that the displacement vector
; xi_r (r=R, th=0, z=L/2, t=P/4) = (A/R) * R  (a fixed percentage of the radius, f.e. 5%)
; 
; 
; Define equilibrium (units: SI with length in Mm)

   m=1.
   harmonic = 1
   dimr = 204.
   dimth = 204.
   dimt = 48./harmonic
   rad = 1.

   norm = 1.e6
   mup = 1.25663706*1.e-6/norm
   kboltz = 1.380658*10^(-23.)
   proton = 1.67262158*10^(-27.)
   gmma = 5./3.

   ne_out_cgs = 3.e9
   ne_in_cgs = 5.e10
   te_out = 2.e6
   te_in = 1.e7
   
   Bx_out = 0.
   By_out = 0.   
   Bz_out = 0.01158     ; (Tesla)
   Bx_in = 0.
   By_in = 0.
   Bz_in = 0.01000
      
   ro = ne_in_cgs * 1.e6 * proton *norm^3
   re = ne_out_cgs * 1.e6 * proton *norm^3
   vao = sqrt(Bx_in^2+By_in^2+Bz_in^2)/sqrt(mup*ro)        ; (Mm per second)
   vae = sqrt(Bx_out^2+By_out^2+Bz_out^2)/sqrt(mup*re)
   bo = vao*sqrt(mup*ro)
   be = vae*sqrt(mup*re)
   co = sqrt(2*gmma*kboltz/proton*te_in)/norm     ; Equals the usual definition, with factor 2 since we consider n as electron density
   ce = sqrt(2*gmma*kboltz/proton*te_out)/norm
   ct = co*vao/sqrt(co^2+vao^2)
   cte = ce*vae/sqrt(ce^2+vae^2)
   
   po=ro*co^2/gmma
   pe=re*ce^2/gmma
   
   A = 0.35     ; A represents amplitude normalised by loop radius, i.e. A/R
   L = 200.      ; L represents L/R
   kz = harmonic*!pi/(L*rad)    ; kz is [1 / Mm] not normalised with respect to radius
   
   wao = vao * kz
   wae = vae * kz

; Eigenfrequency (Patrick style)

  dim=10000
  Q = fltarr(dim+1)
  
   vmin = vao
   va_mx = vae
   wa = findgen(dim+1)/float(dim)*(va_mx-vmin)*kz+vmin*kz

   nummoa2 = (kz^2*co^2-wa^2)*(kz^2*vao^2-wa^2)
   dnummoa2 = (kz^2*ct^2-wa^2)*(co^2+vao^2)
   ki_2 = nummoa2/dnummoa2
   nummea2 = (kz^2*ce^2-wa^2)*(kz^2*vae^2-wa^2)
   dnummea2 = (kz^2*cte^2-wa^2)*(ce^2+vae^2)
   ke_2 = nummea2/dnummea2
   
   reg3 = where(ki_2 lt 0. and ke_2 gt 0.)
   
   
  Q = wa[reg3]^2*(re*sqrt(-ki_2[reg3])*dbeselj(sqrt(-ki_2[reg3])*rad,m)*beselk(sqrt(ke_2[reg3])*rad,m,/double)-ro*sqrt(ke_2[reg3])* $
    dbeselk(sqrt(ke_2[reg3])*rad,m)*beselj(sqrt(-ki_2[reg3])*rad,m,/double))-(re*wae^2*sqrt(-ki_2[reg3])*dbeselj(sqrt(-ki_2[reg3])*rad,m)* $
    beselk(sqrt(ke_2[reg3])*rad,m,/double)-ro*wao^2*sqrt(ke_2[reg3])*dbeselk(sqrt(ke_2[reg3])*rad,m)*beselj(sqrt(-ki_2[reg3])*rad,m,/double))
    
   Q2 = Q/norm
   plot,wa,Q2,/xs,/ys,psym=3,yr=[-0.1,0.1]
   
  wi = wa(where(abs(Q) eq min(abs(Q))))
  w=wi[0]

    
  ke = sqrt((kz^2*ce^2-w^2)*(kz^2*vae^2-w^2)/((kz^2*cte^2-w^2)*(ce^2+vae^2)))
  ki = sqrt(-(kz^2*co^2-w^2)*(kz^2*vao^2-w^2)/((kz^2*ct^2-w^2)*(co^2+vao^2)))
  
  alpha = (w^2-wao^2)*w^3/((w^2-kz^2*ct^2)*(co^2+vao^2)*dbeselj(ki*rad,m)*ki)*A*rad          ; Specifies the displacement of the tube boundary at t=0, z=L/2 (v(r=R,t=0)*dt)
  beta = beselj(ki*rad,m,/double)/beselk(ke*rad,m,/double)*alpha                        ; Extra factor omega to compensate xi_r ~ w^(-1) (Built up from div v!)

; Define grids

  r = findgen(dimr)/(dimr-1.)*2*rad
  r[0] = r[1]/1.1
  th = findgen(dimth)/(dimth)*2*!pi
  z = L/2.*rad
  P = 2*!pi/w
  t = findgen(dimt)/(dimt)*P
  dt = t[1]-t[0]
  
  rpos = fltarr(dimr, dimth)
  thpos = fltarr(dimr, dimth)
  loc = fltarr(dimr, dimth)

; Eigenfunctions with advection

filename = 'eigft'+string(0,format="(i3.3)")+'.dat'
 OPENW, lun, filename, /GET_LUN
 printf, lun, 2
 printf, lun, 8812
 printf, lun, long(dimr*dimth)
 printf, lun, 7

; Time zero

for i=0, dimr-1 do begin
  for j = 0, dimth-1 do begin
      rpos[i,j] = r[i]
      thpos[i,j] = th[j]
      xpos = rpos[i,j]*cos(thpos[i,j])
      ypos = rpos[i,j]*sin(thpos[i,j])     
        if rpos[i,j] lt rad then begin  ; specify loc, calculate eigenf, save position and eigenf, advect
                  
          loc[i,j] = 0
          vri_t = alpha*ki*(w^2-kz^2*ct^2)*(co^2+vao^2)*dbeselj(ki*r[i],m)/((w^2-wao^2)*w^2)*cos(th[j])*sin(kz*z)*cos(w*t[0])
          vti_t = -alpha*(w^2-kz^2*ct^2)*(co^2+vao^2)/((w^2-wao^2)*w^2*r[i])*beselj(ki*r[i],m,/double)*sin(th[j])*sin(kz*z)*cos(w*t[0])
          vxi_t = vri_t * cos(th[j]) - vti_t*sin(th[j])
          vyi_t = vri_t * sin(th[j]) + vti_t*cos(th[j])
          rhi_t = alpha*ro/w*beselj(ki*r[i],m,/double)*cos(th[j])*sin(kz*z)*sin(w*t[0])+ro
          pi_t = po*(1+gmma/ro*rhi_t)         
          ti_t = pi_t/rhi_t ; Possibly correct for scaling factor
          rhi_t = rhi_t/(1.e6 * proton *norm^3)
          ti_t = ti_t*proton/kboltz*norm^2*3./16.
          bri_t = bo*kz*alpha*ki/w*(w^2-kz^2*ct^2)*(co^2+vao^2)*dbeselj(ki*r[i],m)/((w^2-wao^2)*w^2)*cos(th[j])*cos(kz*z)*sin(w*t[0])
          bti_t = -bo*kz*alpha/w*(w^2-kz^2*ct^2)*(co^2+vao^2)/((w^2-wao^2)*w^2*r[i])*beselj(ki*r[i],m,/double)*sin(th[j])*cos(kz*z)*sin(w*t[0])
          bzi_t = bo*alpha/w*(1-kz^2*co^2/w^2)*beselj(ki*r[i],m,/double)*cos(th[j])*sin(kz*z)*sin(w*t[0])+bo
                   
          printf, lun, xpos, ypos, vxi_t, vyi_t, rhi_t, ti_t, bri_t, bti_t, bzi_t
          
          rpos[i,j] = rpos[i,j] + vri_t * dt
          thpos[i,j] = thpos[i,j] + vti_t * dt
          if rpos[i,j] lt 0 then begin
                rpos[i,j] = -rpos[i,j]
                thpos[i,j] = thpos[i,j] + !pi
          endif          
          
        endif else begin
          
          loc[i,j] = 1         
          vre_t = beta*ke*(w^2-kz^2*cte^2)*(ce^2+vae^2)*dbeselk(ke*r[i],m)/((w^2-wae^2)*w^2)*cos(th[j])*sin(kz*z)*cos(w*t[0])
          vte_t = -beta*(w^2-kz^2*cte^2)*(ce^2+vae^2)/((w^2-wae^2)*w^2*r[i])*beselk(ke*r[i],m,/double)*sin(th[j])*sin(kz*z)*cos(w*t[0])
          vxe_t = vre_t * cos(th[j]) - vte_t*sin(th[j])
          vye_t = vre_t * sin(th[j]) + vte_t*cos(th[j])
          rhe_t = beta*re/w*beselk(ke*r[i],m,/double)*cos(th[j])*sin(kz*z)*sin(w*t[0])+re
          pe_t = pe*(1+gmma/re*rhe_t)         
          te_t = pe_t/rhe_t ; Possibly correct for scaling factor
          rhe_t = rhe_t/(1.e6 * proton *norm^3)
          te_t = te_t*proton/kboltz*norm^2*3./16.
          bre_t = be*kz*beta*ke/w*(w^2-kz^2*cte^2)*(ce^2+vae^2)*dbeselk(ke*r[i],m)/((w^2-wae^2)*w^2)*cos(th[j])*cos(kz*z)*sin(w*t[0])
          bte_t = -be*kz*beta/w*(w^2-kz^2*cte^2)*(ce^2+vae^2)/((w^2-wae^2)*w^2*r[i])*beselk(ke*r[i],m,/double)*sin(th[j])*cos(kz*z)*sin(w*t[0])
          bze_t = be*beta/w*(1-kz^2*ce^2/w^2)*beselk(ke*r[i],m,/double)*cos(th[j])*sin(kz*z)*sin(w*t[0])+be
          
          printf, lun, xpos, ypos, vxe_t, vye_t, rhe_t, te_t, bre_t, bte_t, bze_t

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
  filename = 'eigft'+string(n,format="(i3.3)")+'.dat'
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
                  
          vri_t = alpha*ki*(w^2-kz^2*ct^2)*(co^2+vao^2)*dbeselj(ki*r[i],m)/((w^2-wao^2)*w^2)*cos(th[j])*sin(kz*z)*cos(w*t[n])
          vti_t = -alpha*(w^2-kz^2*ct^2)*(co^2+vao^2)/((w^2-wao^2)*w^2*r[i])*beselj(ki*r[i],m,/double)*sin(th[j])*sin(kz*z)*cos(w*t[n])
          vxi_t = vri_t * cos(th[j]) - vti_t*sin(th[j])
          vyi_t = vri_t * sin(th[j]) + vti_t*cos(th[j])
          rhi_t = alpha*ro/w*beselj(ki*r[i],m,/double)*cos(th[j])*sin(kz*z)*sin(w*t[n])+ro
          pi_t = po*(1+gmma/ro*rhi_t)         
          ti_t = pi_t/rhi_t ; Possibly correct for scaling factor
          rhi_t = rhi_t/(1.e6 * proton *norm^3)
          ti_t = ti_t*proton/kboltz*norm^2*3./16.
          bri_t = bo*kz*ki*alpha/w*(w^2-kz^2*ct^2)*(co^2+vao^2)*dbeselj(ki*r[i],m)/((w^2-wao^2)*w^2)*cos(th[j])*cos(kz*z)*sin(w*t[n])
          bti_t = -bo*kz*alpha/w*(w^2-kz^2*ct^2)*(co^2+vao^2)/((w^2-wao^2)*w^2*r[i])*beselj(ki*r[i],m,/double)*sin(th[j])*cos(kz*z)*sin(w*t[n])
          bzi_t = bo*alpha/w*(1-kz^2*co^2/w^2)*beselj(ki*r[i],m,/double)*cos(th[j])*sin(kz*z)*sin(w*t[n])+bo
                   
          printf, lun, xpos, ypos, vxi_t, vyi_t, rhi_t, ti_t, bri_t, bti_t, bzi_t

          rpos[i,j] = rpos[i,j] + vri_t * dt
          thpos[i,j] = thpos[i,j] + vti_t * dt
          if rpos[i,j] lt 0 then begin
                rpos[i,j] = -rpos[i,j]
                thpos[i,j] = thpos[i,j] + !pi
          endif
                    
        endif else begin
          
          vre_t = beta*ke*(w^2-kz^2*cte^2)*(ce^2+vae^2)*dbeselk(ke*r[i],m)/((w^2-wae^2)*w^2)*cos(th[j])*sin(kz*z)*cos(w*t[n])
          vte_t = -beta*(w^2-kz^2*cte^2)*(ce^2+vae^2)/((w^2-wae^2)*w^2*r[i])*beselk(ke*r[i],m,/double)*sin(th[j])*sin(kz*z)*cos(w*t[n])
          vxe_t = vre_t * cos(th[j]) - vte_t*sin(th[j])
          vye_t = vre_t * sin(th[j]) + vte_t*cos(th[j])
          rhe_t = beta*re/w*beselk(ke*r[i],m,/double)*cos(th[j])*sin(kz*z)*sin(w*t[n])+re
          pe_t = pe*(1+gmma/re*rhe_t)         
          te_t = pe_t/rhe_t ; Possibly correct for scaling factor
          rhe_t = rhe_t/(1.e6 * proton *norm^3)
          te_t = te_t*proton/kboltz*norm^2*3./16.
          bre_t = be*kz*ke*beta/w*(w^2-kz^2*cte^2)*(ce^2+vae^2)*dbeselk(ke*r[i],m)/((w^2-wae^2)*w^2)*cos(th[j])*cos(kz*z)*sin(w*t[n])
          bte_t = -be*kz*beta/w*(w^2-kz^2*cte^2)*(ce^2+vae^2)/((w^2-wae^2)*w^2*r[i])*beselk(ke*r[i],m,/double)*sin(th[j])*cos(kz*z)*sin(w*t[n])
          bze_t = be*beta/w*(1-kz^2*ce^2/w^2)*beselk(ke*r[i],m,/double)*cos(th[j])*sin(kz*z)*sin(w*t[n])+be
             
          printf, lun, xpos, ypos, vxe_t, vye_t, rhe_t, te_t, bre_t, bte_t, bze_t
          
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

print, 'time step '+strtrim(string(n),1)
endfor

end