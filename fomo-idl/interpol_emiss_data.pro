
pro interpol_emiss_data,n_e,te,ion=ion, w0=w0,emission_goft=emission_goft, goft=goft, g_logte=g_logte, g_logne=g_logne, tstep=tstep, sav=sav, sdir=sdir, filenm=filenm

emission_goft = n_e*0.d & interp_goft = n_e*0.d
siz = size(n_e)
dimx = siz[1] & dimy = siz[2]

lookup_goft, ion=ion, w0=w0, n_e_lg=g_logne, logt=g_logte, goft_mat=goft_mat, watom= watom, filenm=filenm

logne = alog10(n_e)
logte = alog10(te)
lgne_sort = sort(logne)
logne_sorted = logne[lgne_sort] 
logte_sorted = logte[lgne_sort] 
arr_ne = logne_sorted*0. & num_arrne = n_elements(arr_ne)
arr_gne = findgen(n_elements(g_logne))
arr_te = logte_sorted*0. & num_arrte = n_elements(arr_te)
arr_gte = findgen(n_elements(g_logte))
for i=0,num_arrne-1 do arr_ne[lgne_sort[i]] = interpol(arr_gne,g_logne,logne_sorted[i])
for i=0,num_arrte-1 do arr_te[lgne_sort[i]] = interpol(arr_gte,g_logte,logte_sorted[i])

for i=0,num_arrne-1 do begin
   interp_goft[lgne_sort[i]] = interpolate(goft_mat,arr_ne[lgne_sort[i]],arr_te[lgne_sort[i]],/grid)
   print,string(13b)+' % finished: ',float(i)*100./(num_arrne-1),format='(a,f4.0,$)'
endfor
emission_goft[lgne_sort] = interp_goft[lgne_sort] * n_e[lgne_sort]^2

if keyword_set(sav) eq 1 then begin
   ntstep = string(tstep,"(i4)")
   save,emission_goft,filename=sdir+'emission_goft'+filenm+'_'+ntstep
endif

end
