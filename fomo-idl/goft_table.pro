

pro goft_table, w0=w0, ion=ion

if keyword_set(ion) eq 0 then begin
   print,'goft_table, w0=w0, ion=ion'
   return
endif

; Given wavelength center of line (w0) and ion of interest (ion)
; produces the G(T,n) table for that wavelength transition

  ; determine wavelength transition:
  elements, w0=w0, ion=ion, logTm=logTm, enum=enum, inum=inum, ind=ind

  n_e_min = 1.e8
  n_e_max = 1.e11

;steplg = 0.005
steplg = 0.001

num = alog10(n_e_max/n_e_min)/steplg

n_e_lg = dindgen(num+1)/num*alog10(n_e_max/n_e_min)+alog10(n_e_min)

w0nm = string(round(w0),format='(i4.4)')
;openw,unit,'goft_table_frt_193.dat',/get_lun
;openw,unit,'goft_table_f2rt_171.dat',/get_lun
openw,unit,'goft_table_'+ion+'_'+w0nm+'.dat',/get_lun

for i=0,num do begin
   goft0=g_of_t(enum,inum,dens=n_e_lg[i],ioneq_file=concat_dir(concat_dir(!xuvtop,'ioneq'),'chianti.ioneq'),abund_file=concat_dir(concat_dir(!xuvtop,'abundance'),'sun_coronal.abund'),index=ind,/quiet)

   alogt0 = findgen(101)/20+4.
   wte = 0.75
   tmin = (logTm-wte)>4.
   tmax = (logTm+wte)<6.5
   pts = where(alogt0 ge tmin and alogt0 le tmax,npts)
   alogt1 = alogt0[pts]
   goft1 = goft0[pts]
   alogt2=findgen(201)/200*2*wte+tmin
   ion_interp,alogt1,goft1,alogt2,goft2
   npts2=n_elements(alogt2)

   printf,unit,n_e_lg[i],npts2
   printf,unit,alogt2
   printf,unit,goft2
endfor
free_lun,unit

end
