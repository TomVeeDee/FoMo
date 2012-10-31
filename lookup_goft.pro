
pro lookup_goft, ion=ion, n_e_lg=n_e_lg, logt=logt, goft_mat=goft_mat

; returns got_mat matrix 
; INPUT:
; 

  n_e_min = 1.e8
  n_e_max = 1.e10
;  steplg = 0.005
  steplg = 0.001
  num = alog10(n_e_max/n_e_min)/steplg

if ion eq 'fe_9' then  openr,unit,'goft_table_f2rt_171.dat',/get_lun
if ion eq 'fe_12' then  openr,unit,'goft_table_frt_193.dat',/get_lun
  readf,unit,n_e_0,pts
  goft_mat = fltarr(num+1,pts)
  n_e_lg = dblarr(num+1)
  logt = fltarr(pts)
  g_t = dblarr(pts)
  readf,unit,logt
  readf,unit,g_t
  goft_mat[0,*] = g_t
  n_e_lg[0] = n_e_0
  i = 1 

  while(eof(unit) eq 0) do begin
     readf,unit,n_e_0,pts
     readf,unit,logt
     readf,unit,g_t
     goft_mat[i,*] = g_t
     n_e_lg[i] = n_e_0
     i = i+1
  endwhile
  free_lun,unit
end
