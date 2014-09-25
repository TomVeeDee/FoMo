pro visualadvection, t=t

varnb = 4                       ; 2=vx, 3=vy, 4=rh, 5=br, 6=bth, 7=bz
coarsefactorvf = 2
eigfcube = 'eigft'+string(t,format="(i3.3)")+'.dat'
nx = 204.
ny = 204.

openr,lun,eigfcube,/get_lun
header=lonarr(4)
readf,lun,header
data=fltarr(9,nx*ny)
readf,lun,data
close,lun
free_lun,lun

ti = 'Velocity field'

vxreg = fltarr(nx, ny)
vyreg = fltarr(nx, ny)
vxrd=fltarr(nx/coarsefactorvf, ny/coarsefactorvf)
vyrd=fltarr(nx/coarsefactorvf, ny/coarsefactorvf)

x=reform(data[0,*])
y=reform(data[1,*])
vx=reform(data[2,*])
vy=reform(data[3,*])
var=reform(data[varnb,*])

triangulate, x,y,tr,b
;contour, trigrid(x,y,vx,tr), xtitle='x',ytitle='y',nlevels=30,/fill,/xstyle,/ystyle,title=ti
contour, trigrid(x,y,var,tr), nlevels=30, /fill

vxreg = trigrid(x,y,vx, input=vxreg, tr)
vyreg = trigrid(x,y,vy, input=vyreg, tr)
x=findgen(nx)
y=findgen(ny)
;v = VECTOR(vxreg, vyreg, x, y, xtitle = 'x (*100km)', ytitle = 'y (*100km)', title='Velocity field at z=L/2, t=4, L/R=20')

for i = 1, nx/coarsefactorvf do begin
  for j = 1, ny/coarsefactorvf do begin
    vxrd[i-1,j-1]=vxreg[coarsefactorvf*i-1,coarsefactorvf*j-1]
    vyrd[i-1,j-1]=vyreg[coarsefactorvf*i-1,coarsefactorvf*j-1]
  endfor
endfor

x = findgen(nx/coarsefactorvf)*coarsefactorvf
y = findgen(ny/coarsefactorvf)*coarsefactorvf
wait, 0.2

if t mod 15 eq 0 then v = VECTOR(vxrd, vyrd, x, y, xrange=[0,50], yrange=[0,50], xtitle = 'x', ytitle = 'y', title='Velocity field at z=L/2, t=' + string(t))

end