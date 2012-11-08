
pro readcubes,set=set,rh_cube_t,te_cube_t,vr_cube_t,vz_cube_t,br_cube_t,bz_cube_t,aa,gridx,gridy,gridz,gridz_ext,co,ce,ro,re,va,vae,dimt,tarr,r0,ka_rt,wk_rt,rh_cube_t_ext,wave,dim00d,dim30d,dim45d,dim60d,mag=mag

; Reads cubes produced by datacubes_wt.pro

if keyword_set(set) eq 0 then begin
 print,'readcubes,set=set,rh_cube_t,te_cube_t,vr_cube_t,vz_cube_t,br_cube_t,bz_cube_t,aa,gridx,gridy,gridz,gridz_ext,co,ce,ro,re,va,vae,dimt,tarr,r0,ka_rt,wk_rt,rh_cube_t_ext,wave,dim00d,dim30d,dim45d,dim60d [,mag=mag]'
   return
endif

if set eq 2 or set eq 23 then begin
   ka_n = 'ka2.24'
   dir = '/volume1/scratch/set2/'
   dir2 = dir
   ndz = 3
;   dir = '/users/cpa/pantolin/Modeling/cubes/set2/'
;   dir2 = '/volume1/scratch/set2/'
endif
if set eq 24 then begin
   ka_n = 'ka2.24_hgres'
   dir = '/volume1/scratch/set2/hgres/'
   ndz = 3
;   dir = '/users/cpa/pantolin/Modeling/cubes/set2/'
   dir2 = dir
;   dir2 = '/volume1/scratch/set2/'
endif
if set eq 25 then begin
   ka_n = 'ka2.24_hgres2d'
   dir = '/volume1/scratch/set2/hgres2d/'
   ndz = 3
;   dir = '/users/cpa/pantolin/Modeling/cubes/set2/'
   dir2 = dir
;   dir2 = '/volume1/scratch/set2/'
endif
if set eq 22 then begin
   ka_n = 'ka2.24_sm'
   dir = '/volume1/scratch/set2/sm/'
   dir2 = dir
   ndz = 3
;   dir = '/users/cpa/pantolin/Modeling/cubes/set2/'
;   dir2 = '/volume1/scratch/set2/'
endif
if set eq 21 then begin
   ka_n = 'ka2.24_highT'
   ndz = 3
   dir = '/volume1/scratch/set2/high_T/'
   dir2 = dir
;   dir = '/users/cpa/pantolin/Modeling/cubes/set2/'
endif
if set eq 3 or set eq 31 then begin
   ka_n = 'ka1.25'
   ndz = 2
   dir = '/volume1/scratch/set3/'
   dir2 = dir
;   dir = '/users/cpa/pantolin/Modeling/cubes/set3/'
endif
   
;dir =  '/Volumes/Karmeliet/Data/modeling/cubes/'
if set ne 24 then restore,dir+'cubes_'+ka_n+'_'+string(0,format="(i3.3)")+'.sav' else restore,dir+'slice_rh_'+ka_n+'_'+string(0,format="(i3.3)")+string(0,format='(i3.3)')+'x'+'.sav'
if set eq 22 then siz = size(rh_cube_sm) 
if set eq 2 or set eq 23 or set eq 21 or set eq 3 or set eq 25 or set eq 26 then siz = size(rh_cube)
if set eq 24 then siz = size(rho)
if set eq 25 or set eq 26 then ldimz = 1 else ldimz = 2
restore,dir+'params_'+ka_n+'.sav'

if set ne 25 and set ne 26 then files = file_search(dir2+'rslt_slice_*.sav',count=nfiles,/fully_qualify_path) 
if set eq 25 then files = file_search(dir2+'rslt_cubes_'+ka_n+'.fe9'+'*.sav',count=nfiles,/fully_qualify_path)
if set eq 26 then files = file_search(dir2+'rslt_cubes_'+ka_n+'.fe12'+'*.sav',count=nfiles,/fully_qualify_path)
restore,files[0]
dim00d = (size(image00d_ext))[ldimz]
dim30d = (size(image30d_ext))[ldimz]
dim45d = (size(image45d_ext))[ldimz]
dim60d = (size(image60d_ext))[ldimz]

if set ne 25 and set ne 26 then begin 
   rh_cube_t = fltarr(siz[2],siz[3],dimt)
   rh_cube_t_ext = fltarr(siz[2],siz[3]*ndz,dimt)
   te_cube_t = fltarr(siz[2],siz[3],dimt)
   if keyword_set(mag) then begin
      br_cube_t = fltarr(siz[2],siz[3],dimt)
      bz_cube_t = fltarr(siz[2],siz[3],dimt)
   endif
   if set ne 24 then begin
      vr_cube_t = fltarr(siz[2],siz[3],dimt)
      vz_cube_t = fltarr(siz[2],siz[3],dimt)
   endif
   gridz_ext = fltarr(siz[3]*ndz)
endif else begin
   rh_cube_t = fltarr(siz[1],siz[2],dimt)
   rh_cube_t_ext = fltarr(siz[1],siz[2]*ndz,dimt)
   te_cube_t = fltarr(siz[1],siz[2],dimt)
   vr_cube_t = fltarr(siz[1],siz[2],dimt)
   vz_cube_t = fltarr(siz[1],siz[2],dimt)
   if keyword_set(mag) then begin
      br_cube_t = fltarr(siz[1],siz[2],dimt)
      bz_cube_t = fltarr(siz[1],siz[2],dimt)
   endif
   gridz_ext = fltarr(siz[2]*ndz)
endelse

if set eq 2 or set eq 23 or set eq 21 or set eq 3 then begin
   for i=0,dimt-1 do begin
      restore,dir+'cubes_'+ka_n+'_'+string(i,format="(i3.3)")+'.sav'
      rh_cube_t[*,*,i]=rh_cube[dimx/2,*,*]
      te_cube_t[*,*,i]=te_cube[dimx/2,*,*]
      vr_cube_t[*,*,i]=vr_cube[dimx/2,*,*]
      vz_cube_t[*,*,i]=vz_cube[dimx/2,*,*]
      if keyword_set(mag) then begin
         br_cube_t[*,*,i] = br_cube[dimx/2,*,*]
         bz_cube_t[*,*,i] = bz_cube[dimx/2,*,*]
      endif
      print,string(13b)+' % finished: ',float(i)*100./(dimt-1),format='(a,f4.0,$)'
   endfor
endif
if set eq 24 then begin
   j=dimx/2
   for i=0,dimt-1 do begin
      restore,dir+'slice_rh_'+ka_n+'_'+string(i,format="(i3.3)")+string(j,format='(i3.3)')+'x'+'.sav'
      rh_cube_t[*,*,i]=rho[0,*,*]
      restore,dir+'slice_te_'+ka_n+'_'+string(i,format="(i3.3)")+string(j,format='(i3.3)')+'x'+'.sav'
      te_cube_t[*,*,i]=te[0,*,*]
;      restore,dir+'slice_vr_'+ka_n+'_'+string(i,format="(i3.3)")+string(j,format='(i3.3)')+'x'+'.sav'
;      vr_cube_t[*,*,i]=vr[0,*,*]
;      restore,dir+'slice_vz_'+ka_n+'_'+string(i,format="(i3.3)")+string(j,format='(i3.3)')+'x'+'.sav'
;      vz_cube_t[*,*,i]=vz[0,*,*]
      print,string(13b)+' % finished: ',float(i)*100./(dimt-1),format='(a,f4.0,$)'
   endfor
endif
if set eq 22 then begin
   for i=0,dimt-1 do begin
      restore,dir+'cubes_'+ka_n+'_'+string(i,format="(i3.3)")+'.sav'
      rh_cube_t[*,*,i]=rh_cube_sm[dimx/2,*,*]
      te_cube_t[*,*,i]=te_cube_sm[dimx/2,*,*]
      vr_cube_t[*,*,i]=vr_cube_sm[dimx/2,*,*]
      vz_cube_t[*,*,i]=vz_cube_sm[dimx/2,*,*]
      if keyword_set(mag) then begin
         br_cube_t[*,*,i] = br_cube_sm[dimx/2,*,*]
         bz_cube_t[*,*,i] = bz_cube_sm[dimx/2,*,*]
      endif      
      print,string(13b)+' % finished: ',float(i)*100./(dimt-1),format='(a,f4.0,$)'
   endfor
endif
if set eq 25 or set eq 26 then begin
   for i=0,dimt-1 do begin
      restore,dir+'cubes_'+ka_n+'_'+string(i,format="(i3.3)")+'.sav'
      rh_cube_t[*,*,i]=rh_cube
      te_cube_t[*,*,i]=te_cube
      vr_cube_t[*,*,i]=vr_cube
      vz_cube_t[*,*,i]=vz_cube
      if keyword_set(mag) then begin
         br_cube_t[*,*,i] = br_cube
         bz_cube_t[*,*,i] = bz_cube
      endif
      print,string(13b)+' % finished: ',float(i)*100./(dimt-1),format='(a,f4.0,$)'
   endfor
endif
if set eq 2 or set eq 21 or set eq 23 or set eq 22 or set eq 24 or set eq 25 or set eq 26 then begin
   gridz_ext[0:dimz-1] = gridz
   gridz_ext[dimz:2*dimz-1] = gridz[0:dimz-1]+gridz[dimz-1]+gridz[1]
   gridz_ext[dimz*2:*] = gridz_ext[dimz:dimz*2-1]+gridz[dimz-1]+gridz[1]
   rh_cube_t_ext[*,0:dimz-1,*] = rh_cube_t
   rh_cube_t_ext[*,dimz:2*dimz-1,*] = rh_cube_t
   rh_cube_t_ext[*,dimz*2:*,*] = rh_cube_t
endif
if set eq 3 then begin
   gridz_ext[0:dimz-1] = gridz
   gridz_ext[dimz:*] = gridz[0:dimz-1]+gridz[dimz-1]+gridz[1]   
   rh_cube_t_ext[*,0:dimz-1,*] = rh_cube_t
   rh_cube_t_ext[*,dimz:*,*] = rh_cube_t
endif

end
