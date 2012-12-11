pro lineongrid_goft_tab, rh_s, te_s, wave=wave,nwave=nwave,minwave=minwave,maxwave=maxwave,ion=ion, w0=w0, emission_goft=emission_goft, goft=goft, wayemi=wayemi

if n_params(0) lt 1 then begin
   print,'Check input directories'
   print,'lineongrid_goft_tab, rh_s, te_s, wave=wave,nwave=nwave,minwave=minwave,maxwave=maxwave,ion=ion, w0=w0, emission_goft=emission_goft, goft=goft, wayemi=wayem'
   return
endif

; Calculates the emission at each point of a given numerical box by
; reading tabulated G(T) values produced by function goft_table.pro
; INPUT:
; rh_j = 2D or 3D array of mass density in kg/m^3, normalized by 1.e10
; te_j = 2D or 3D array of temperature in K, normalized by
; (protonmass/(2*kboltz)*1.e10)
; OPTIONAL:
; ion = ion for which to calculate emissivities (default = fe_9)
; w0 = wavelength center of line transition (default = 171.073)
; minwave = infimum of wavelength range (default = 171.0)
; maxwave = suppremum of wavelength range (default = 171.14)
; set keyword wayemi to 1 if emission calculated point by point, else
; set to 2 if calculated through binning of G(T)*ne^2 term (default = 2)
; OUTPUT:
; wave = wavelength array set to nwave pts, containing line transition	
; w0 = wavelength center
; emission_goft = array of emission (x, y, (z), lambda)

if keyword_set(wayemi) eq 0 then begin 
   way = 2 
endif else begin
   if wayemi eq 1 then way = 1
   if wayemi eq 2 then way = 2
endelse

; set ionization,abundance and DEM packages for CHIANTI:
ioneq_name = '/users/cpa/tomvd/ssw/packages/chianti/dbase/ioneq/chianti.ioneq'
abund_name = '/users/cpa/tomvd/ssw/packages/chianti/dbase/abundance/sun_coronal.abund'
dem_name = '/users/cpa/tomvd/ssw/packages/chianti/dbase/dem/active_region.dem'
;ioneq_name= !xuvtop+'/ioneq/chianti.ioneq'
;abund_name=!xuvtop+'/abundance/sun_photospheric.abund'
;abund_name=!xuvtop+'/abundance/sun_coronal.abund'
;dem_name=!xuvtop+'/dem/active_region.dem'

proton=1.67262158*10^(-27.)
kboltz = 1.380658*10^(-23.)
gamma = 5./3.
mu = 1.27
way = 2 ; calculate emission through binning of G(T)*ne^2 matrix

normro = 1.e10
;   const = 3.5991482e+13 ; = T[mid]*rho[mid]^(-gamma+1)
;   normte = const*normro^(-gamma+1)
normte = proton/(2*kboltz)*normro
rh = rh_s / normro
T = te_s * normte

n_e=rh/proton/1.e6 ; in cgs
logT = alog10(T>1.)
logne = alog10(n_e)

sizes=size(rh)
dims=sizes[0]
; if dims<2 exit
if dims eq 3 then begin
   nx=sizes[1]
   ny=sizes[2]
   nz=sizes[3]
endif else begin
   nx = sizes[1]
   nz = sizes[2]
endelse

;if (~(keyword_set(ion))) then begin
;   ion='fe_9'  
;   ion = 'fe_12'
;endif

if (~(keyword_set(nwave))) then begin
   nwave=100
endif
;if ion eq 'fe_9' then w0 = 171.073 ; wave center
;if ion eq 'fe_12' then w0 = 193.509
minwave = w0-0.07
maxwave = w0+0.07
wave=findgen(nwave)/(nwave-1)*(maxwave-minwave)+minwave

if (dims eq 2) then emission_goft=dblarr(nx,nz,nwave)
if (dims eq 3) then emission_goft=dblarr(nx,ny,nz,nwave)
if way eq 3 then begin
   lognen = fltarr(nx,ny,nz,nwave)
   for i=0,nwave-1 do lognen[*,*,*,i]=logne
endif

pe = n_e * kboltz * T * 1.e7 ; in cgs

ne_sort = sort(n_e)
n_e_sorted = n_e[ne_sort]
Tlg_sorted = logT[ne_sort]

; Read tabulated G(ne,T) values for given number density (n_e_lg) and
; temperature (t_lg) arrays
lookup_goft, ion=ion, w0=w0, n_e_lg=n_e_lg, logt=t_lg, goft_mat=goft_mat

; create a ch_synthetic structure called "singleline"
ch_synthetic,min(wave),max(wave),output=singleline,density=max(n_e),logt_isothermal=max(alog10(T)),ioneq_name=ioneq_name,sngl_ion=ion

goft=n_e*0.
; calculate emission point by point:
if way eq 1 then begin
   for i=0.,n_elements(n_e)-1 do begin
      lc_ne = ([min(abs(n_e_lg-alog10(n_e_sorted[i]))),!c])[1]
      lc_te = ([min(abs(t_lg-tlg_sorted[i])),!c])[1]
      coords=array_indices(n_e,ne_sort[i])
;      if dims eq 2 then goft[coords[0],coords[1]] = spline(t_lg,goft_mat[lc_ne,*],tlg_sorted[i])
;      if dims eq 3 then goft[coords[0],coords[1],coords[2]] = spline(t_lg,goft_mat[lc_ne,*],tlg_sorted[i])
      goft[ne_sort[i]] = spline(t_lg,goft_mat[lc_ne,*],tlg_sorted[i])
      singleline.logt_isothermal[0]=tlg_sorted[i]
      singleline.lines.int[0]=goft[ne_sort[i]]*n_e_sorted[i]^2
      make_chianti_spec, singleline, wave, out, abund_name=concat_dir(concat_dir(!xuvtop,'abundance'),'sun_coronal.abund')
      emission_goft[coords[0],coords[1],coords[2],*]=out.spectrum
      print,string(13b)+' % finished: ',float(i)*100./(n_elements(t)-1),format='(a,f4.0,$)'
   endfor
endif

; calculate emission through binning of G(T)*ne^2 matrix
if way eq 2 then begin
   emi=n_e*0.
   for i=0.,n_elements(n_e)-1 do begin
      lc_ne = ([min(abs(n_e_lg-alog10(n_e_sorted[i]))),!c])[1]
      lc_te = ([min(abs(t_lg-tlg_sorted[i])),!c])[1]
      goft[ne_sort[i]] = interpol(goft_mat[lc_ne,*],t_lg,tlg_sorted[i])
      emi[ne_sort[i]] = goft[ne_sort[i]]*n_e_sorted[i]^2
      print,string(13b)+' % finished: ',float(i)*100./(n_elements(n_e)-1),format='(a,f4.0,$)'
   endfor

   numemi = min([n_elements(uniq(emi)),10000]) ; 10000 pts are sufficient
   hist_emi = histogram(emi,nbins=numemi,locations=loc_emi,reverse_indices=Rem)
   nhemi = n_elements(hist_emi)
   if dims eq 3 then begin
      emin = fltarr(nx,ny,nz,nwave) 
      for i=0,nwave-1 do emin[*,*,*,i]=emi
   endif else begin
      emin = fltarr(nx,nz,nwave)
      for i=0,nwave-1 do emin[*,*,i]=emi
   endelse
   hist_emin = histogram(emin,nbins=numemi,locations=loc_emin,reverse_indices=Remn)
   for i=0.,nhemi-1 do begin
      if hist_emi[i] ne 0 then begin
         singleline.logt_isothermal[0] = max(logt[Rem[Rem[i]:Rem[i+1]-1]])
         singleline.lines.int[0] = loc_emi[i]
         make_chianti_spec, singleline, wave, out, abund_name=abund_name
         nRem = n_elements(Rem[Rem[i]:Rem[i+1]-1])
         for j=0.,nRem-1 do emission_goft[Remn[Remn[i]+j:Remn[i+1]-1:nRem]] = out.spectrum
      endif
      print,string(13b)+' % finished: ',float(i)*100./(nhemi-1),format='(a,f4.0,$)'
   endfor
endif

if way eq 3 then begin
   num = n_elements(n_e_lg)
   hist_ne=histogram(logne,nbins=401,locations=loc_ne,min=8,max=10,reverse_indices=Rne)
   hist_nen=histogram(lognen,nbins=401,locations=loc_nen,min=8,max=10,reverse_indices=Rnen)
   nhne = n_elements(hist_ne)
   for i=0.,nhne-1 do begin
      if hist_ne[i] ne 0 then begin
         singleline.logt_isothermal[0] = tlg_sorted[i]
         nRne = n_elements(Rne[Rne[i]:Rne[i+1]-1])
         stop
         for j=0.,nRne-1 do begin
            goft[Rne[Rne[i]+j]]=spline(t_lg,goft_mat[i,*],logt[Rne[Rne[i]+j]])
            singleline.lines.int[0]=goft[Rne[Rne[i]+j]]*(10.^loc_ne[i])^2
            make_chianti_spec, singleline, wave, out, abund_name=concat_dir(concat_dir(!xuvtop,'abundance'),'sun_coronal.abund')
            emission_goft[Rnen[Rnen[i]+j:Rnen[i+1]-1:nRne]] = out.spectrum
         endfor
      endif
      print,string(13b)+' % finished: ',float(i)*100./(nhne-1),format='(a,f4.0,$)'
   endfor
endif
        
end
