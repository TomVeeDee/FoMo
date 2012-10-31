

pro smprof_cubes, dimz=dimz, dimt=dimt, gridx=gridx, rtot_t=rtot_t, te_t=te_t, rtot_sm_t, te_sm_t

if n_params(0) lt 1 then begin
   print,'smprof_cubes, dimz=dimz, dimt=dimt, gridx=gridx, rtot_t=rtot_t, te_t=te_t, rtot_sm_t, te_sm_t'
   return
endif

; Produces smooth density and temperature profiles based on sharp profiles
; INPUT:
; dimt = number of points in time dimension
; dimz = dimension in z direction (set equal to 1 wavelength)
; gridx = grid along x axis (radial direction)
; rtot_t = rtot(dimx,dimz,dimt) : total density
; te_t = te_t(dimx,dimy,dimz) : temperature
; OUTPUT:
; rtot_sm_t = rtot_sm_t(dimx,dimz,dimt) : total density with smooth profile
; te_sm_t = te_sm_t(dimx,dimz,dimt) : temperature with smooth profile

  te_i = te_t & te_sm_t = te_t
  rho_i = rtot_t & rtot_sm_t = rtot_t
  stp = 1.

  for j=0,dimt-1 do begin
     for i=0,dimz-1 do begin
        gr1 = ([max(te_t[1:-1,i,j]-te_t[0:-1,i,j]),!c])[1]
        ite1 = (te_t[gr1+5,i,j]-te_t[gr1-5,i,j])/10.*indgen(10)+te_t[gr1-5,i,j]
        iro1 = (rtot_t[gr1+5,i,j]-rtot_t[gr1-5,i,j])/10.*indgen(10)+rtot_t[gr1-5,i,j]
        gr2 = ([min(te_t[1:-1,i,j]-te_t[0:-1,i,j]),!c])[1]
        ite2 = (te_t[gr2+5,i,j]-te_t[gr2-5,i,j])/10.*indgen(10)+te_t[gr2-5,i,j]
        iro2 = (rtot_t[gr2+5,i,j]-rtot_t[gr2-5,i,j])/10.*indgen(10)+rtot_t[gr2-5,i,j]
        nrm_tan = ([max(atan((gridx-gridx[gr1])*stp)/!pi-atan((gridx-gridx[gr2])*stp)/!pi),!c])[0]
        tan_sm=(atan((gridx-gridx[gr1])*stp)/!pi-atan((gridx-gridx[gr2])*stp)/!pi)/nrm_tan
        te_i[gr1-5:gr1+4,i,j] = ite1   
        te_i[gr2-5:gr2+4,i,j] = ite2
        te_sm_t[*,i,j] = (te_i[*,i,j]-te_i[0,i,j])*tan_sm+te_i[0,i,j]
        rho_i[gr1-5:gr1+4,i,j] = iro1   
        rho_i[gr2-5:gr2+4,i,j] = iro2
        rtot_sm_t[*,i,j]=(rho_i[*,i,j]-rho_i[0,i,j])*tan_sm+rho_i[0,i,j]
     endfor
  endfor

end
