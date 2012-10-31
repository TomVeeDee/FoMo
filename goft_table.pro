

pro goft_table

; for Fe IX 171.073, set lin = 1
; for Fe XII 193.509, set lin = 2
lin = 1

  if lin eq 1 then begin
     enum = 26                  ; number for element
     inum = 9                   ; number of ion
     ind = 151                  ; index for wavelength of transition
  endif
  if lin eq 2 then begin
     enum = 26
     inum = 12
     ind = 844
  endif
n_e_min = 1.e8
n_e_max = 1.e10

;steplg = 0.005
steplg = 0.001

num = alog10(n_e_max/n_e_min)/steplg

n_e_lg = dindgen(num+1)/num*alog10(n_e_max/n_e_min)+alog10(n_e_min)

;openw,unit,'goft_table_frt_193.dat',/get_lun
openw,unit,'goft_table_f2rt_171.dat',/get_lun

for i=0,num do begin
   goft0=g_of_t(enum,inum,dens=n_e_lg[i],ioneq_file=concat_dir(concat_dir(!xuvtop,'ioneq'),'chianti.ioneq'),abund_file=concat_dir(concat_dir(!xuvtop,'abundance'),'sun_coronal.abund'),index=ind,/quiet)

   alogt0 = findgen(101)/20+4.
   pts = where(alogt0 ge 4.95 and alogt0 le 6.5,npts)
   alogt1 = alogt0[pts]
   goft1 = goft0[pts]
   alogt2=findgen(101)/100*1.55+4.95
   ion_interp,alogt1,goft1,alogt2,goft2
   npts2=n_elements(alogt2)

   printf,unit,n_e_lg[i],npts2
   printf,unit,alogt2
   printf,unit,goft2
endfor
free_lun,unit

end
