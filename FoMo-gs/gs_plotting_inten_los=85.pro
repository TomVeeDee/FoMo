
; INPUT DIRECTORY:
  slicedir='E:\Belgium\Wiki\FoMo-GS-code\GS-slices\'
; OUTPUT DIRECTORY:
  savedir='E:\Belgium\Wiki\FoMo-GS-code\'

 dimy=42 & Radius = 42 & t=0
frequen = ['1','1.6','2.5','4','6.3','10','16','25','35','63','100']; GHz
n_f=[3,5,10]
gridy=findgen(dimy*2)/dimy
nf=n_elements(n_f)
sz=129 ;number of rays
Int_cube = fltarr(sz,nf,dimy)
Pol_cube = fltarr(sz,nf,dimy)
i=0
for j=61,102 do begin
  restore,slicedir+'GS_emission'+string(j,format='(i3.3)')+'x'+'.sav'
 
  It=Int[*,n_f]
  Pl=Pol[*,n_f]
  
  Int_cube[*,*,i] = It
  Pol_cube[*,*,i] = Pl
    i=i+1
endfor
gridz_ext=findgen(sz)/Radius 
ff=0
  cube_Int = fltarr(sz,dimy*2) & cube_Pol = fltarr(sz,dimy*2)
;glue array in y directions
  cube_Int [0:sz-1,0:dimy-1] = bytscl(Int_cube(*,ff,*),top=254)
  half_cylinder=fltarr(sz,dimy)
  half_cylinder=cube_Int(*,0:dimy-1)
  cube_Int [*,dimy:dimy*2-1] = reverse(half_cylinder,2)
cube_Pol [0:sz-1,0:dimy-1] = bytscl(Pol_cube(*,ff,*),top=254)
half_cylinder=cube_pol(*,0:dimy-1)
cube_Pol [*,dimy:dimy*2-1] = reverse(half_cylinder,2)

     set_plot,'ps'
     !P.thick = 5
     !x.thick = 3
     !y.thick = 3
device,filename=savedir+'GS_emission_LOS=85-angle.eps',/encapsulated,xsize=8,ysize=10,/inches,/portrait,set_font='Times',/tt_font,font_size=33,bits_per_pixel=8,/color
    loadct,3        
  cube_Int[0,0,0]=255
  cube_Pol[0,0,0]=255
titl_int =' f='+frequen(n_f(ff))+' GHz'  & tis_int = '[10!u-4!n SFU]'
;  titl_pol = 'Circular polarization degree' & tis_pol = '[10!u-4!n %]'
  
  !P.multi=[0,1,4,0,0]
  plot_image,cube_Int(*,*),gridz_ext,gridy,xtitle='z, Mm',ytitle='y, Mm',title=titl_int,charthick=2,$
  max=max(cube_int),min=min(cube_int), /iso,background=white, xticks=3, xtickname=['0','1','2','3'],yticks=2,ytickname=['0','1','2']
  chars = 0.7
  
   x1=0.53 & x2=.545 & y1=.825 & y2=.962
  color_bar, x1, x2, y1, y2, /normal, max=(max(Int_cube(*,ff,*))*1.e4),min=fix(min(Int_cube(*,ff,*))*1.e4),$
  ticks=3,color=1,chars=chars,charth=3,tit=tis_int,tickf='(i3)'
ff=1
titl_int =' f='+frequen(n_f(ff))+' GHz'  & tis_int = '[10!u-4!n SFU]'
  cube_Int = fltarr(sz,dimy*2) & cube_Pol = fltarr(sz,dimy*2)
  cube_Int [0:sz-1,0:dimy-1] = bytscl(Int_cube(*,ff,*));,top=254)
half_cylinder=fltarr(sz,dimy)
half_cylinder=cube_Int(*,0:dimy-1)
cube_Int [*,dimy:dimy*2-1] = reverse(half_cylinder,2)

plot_image,cube_Int(*,*),gridz_ext,gridy,xtitle='z, Mm',ytitle='y, Mm',title=titl_int,charthick=2,$
  max=max(cube_int),min=min(cube_int), /iso,background=white,  xticks=3, xtickname=['0','1','2','3'],yticks=2,ytickname=['0','1','2']
  chars = 0.7
y1=.575 & y2=.715

color_bar, x1, x2, y1, y2, /normal, max=(max(Int_cube(*,ff,*))*1.e4),min=fix(min(Int_cube(*,ff,*))*1.e4),$
  ticks=3,color=1,chars=chars,charth=3,tit=tis_int,tickf='(i3)
ff=2
titl_int =' f='+frequen(n_f(ff))+' GHz'  & tis_int = '[10!u-4!n SFU]'
  cube_Int = fltarr(sz,dimy*2) & cube_Pol = fltarr(sz,dimy*2)
  cube_Int [0:sz-1,0:dimy-1] = bytscl(Int_cube(*,ff,*))
  half_cylinder=fltarr(sz,dimy)
half_cylinder=cube_Int(*,0:dimy-1)
cube_Int [*,dimy:dimy*2-1] = reverse(half_cylinder,2)

plot_image,cube_Int(*,*),gridz_ext,gridy,xtitle='z, Mm',ytitle='y, Mm',title=titl_int,charthick=2,$
  max=max(cube_int),min=min(cube_int), /iso,background=white, xticks=3, xtickname=['0','1','2','3'],yticks=2,ytickname=['0','1','2']
  chars = 0.7
y1=.326 & y2=.464
  color_bar, x1, x2, y1, y2, /normal, max=(max(Int_cube(*,ff,*))*1.e4),min=fix(min(Int_cube(*,ff,*))*1.e5),$
  ticks=3,color=1,chars=chars,charth=3,tit=tis_int,tickf='(i3)'

     device,/close
     set_plot,'win'
!p.multi=0
end