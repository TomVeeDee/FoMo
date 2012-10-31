
pro rays, set=set, gridy=gridy, n_gridy=n_gridy, n_gridz=n_gridz, ngrid=ngrid, nray0, nray1, vr_cube_t=vr_cube_t, gridz=gridz

  if n_params(0) lt 1 then begin
     print,'do before:'
     print,'gridlos, gridx=gridz_ext, gridy=gridy, gridz=gridx, mua_d=mua_d, velx=0., vely=0., velz=0., n_gridz, n_gridy, ngrid, losvel_ext'
     print,'then:'
     print,'rays, set=set, gridy=gridy, n_gridy=n_gridy, n_gridz=n_gridz, ngrid=ngrid, ray00d_0, ray00d_1, vr_cube_t=vr_cube_t, gridz=gridz'
     return
  endif

  dimy = n_elements(gridy)
  if set eq 0 then begin
     z0 = 67. & y0 = 102.
     z1 = 89. & y1 = 102.
  endif
  if set eq 2 or set eq 21 or set eq 22 or set eq 23 or set eq 24 then begin
     siz = size(vr_cube_t)
     dimz = n_elements(gridz)
     vr_cube_t_ext = fltarr(siz[1],siz[2]*3,siz[3])
     vr_cube_t_ext[*,0:dimz-1,*] = vr_cube_t
     vr_cube_t_ext[*,dimz:2*dimz-1,*] = vr_cube_t
     vr_cube_t_ext[*,2*dimz:-1,*] = vr_cube_t
     z0 = ([min(abs(vr_cube_t_ext[siz[1]/2,round(siz[2]*5/4):siz[2]*7/4,0])),!c])[1]+round(siz[2]*5/4) & y0=dimy/2
     z1 = ([max(abs(vr_cube_t_ext[siz[1]/2,round(siz[2]*6/4):siz[2]*8/4,0])),!c])[1]+round(siz[2]*6/4) & y1=dimy/2
;     z0 = ([min(abs(vr_cube_t_ext[siz[1]/2,150:200,0])),!c])[1]+150 & y0=dimy/2
;     z1 = ([max(abs(vr_cube_t_ext[siz[1]/2,180:220,0])),!c])[1]+180 & y1=dimy/2
;     z0 = 170 & y0 = 101
;     z1 = 199 & y1 = 101
  endif
  if set eq 3 then begin
     siz = size(vr_cube_t)
     dimz = n_elements(gridz)
     vr_cube_t_ext = fltarr(siz[1],siz[2]*2,siz[3])
     vr_cube_t_ext[*,0:dimz-1,*] = vr_cube_t
     vr_cube_t_ext[*,dimz:-1,*] = vr_cube_t
     z0 = ([min(abs(vr_cube_t_ext[siz[1]/2,150:250,0])),!c])[1]+150 & y0 = dimy/2
     z1 = ([max(abs(vr_cube_t_ext[siz[1]/2,200:300,0])),!c])[1]+200 & y1 = dimy/2
;     z0 = 201 & y0 = 101
;     z1 = 252 & y1 = 101
  endif
  loczy0 = fltarr(n_elements(ngrid)-1,2)
  loczy1 = fltarr(n_elements(ngrid)-1,2)
  for i=0,n_elements(ngrid)-2 do begin
     loczy0[i,*]=[min(sqrt((n_gridy[ngrid[i]:ngrid[i+1]-1]-y0)^2+(n_gridz[ngrid[i]:ngrid[i+1]-1]-z0)^2)),!c]
     loczy1[i,*]=[min(sqrt((n_gridy[ngrid[i]:ngrid[i+1]-1]-y1)^2+(n_gridz[ngrid[i]:ngrid[i+1]-1]-z1)^2)),!c]
  endfor
  nray0 = ([min(loczy0[*,0]),!c])[1]
  nray1 = ([min(loczy1[*,0]),!c])[1]
  print,'for z0:',nray0
  print,'for z1:',nray1
end
