function readgoftcube_dat,emissionsave,data,compress=compress

; reads the data in emissionsave files, not the triangulation
; invoke as data=readgoftcube('fileyouwanttoread')
; data will be a 4xn array, with n the number of grid points in x,y,\lambda

maxversionlength=40

openr,lun,emissionsave,/get_lun,compress=compress
character='a'
versionstring=''
count=0
while ((character ne '#') and (count lt maxversionlength)) do begin
	readu,lun,character
	versionstring+=character
	count++
endwhile
versionstring=strmid(versionstring,0,strlen(versionstring)-1)


ng=0L
dim=0L
nvars=0L
unitsize=0LL
chiantisize=0LL
abundsize=0LL
if (strpos(versionstring,'FoMo-v') eq -1) then begin
; this is the old version of the file
	point_lun,lun,0
	
	readu,lun,dim
	readu,lun,ng
	readu,lun,nvars
	readu,lun,chiantisize
	chiantibuffer=bytarr(chiantisize)
	readu,lun,chiantibuffer
	chiantifile=string(chiantibuffer)
	readu,lun,abundsize
	abundbuffer=bytarr(abundsize)
	readu,lun,abundbuffer
	abundfile=string(abundbuffer)
endif else begin
	readu,lun,dim
	readu,lun,ng
	readu,lun,nvars
	units=strarr(dim+nvars)
	for i=0,dim+nvars-1 do begin
		readu,lun,unitsize
		unitbuffer=bytarr(unitsize)
		readu,lun,unitbuffer
		units[i]=string(unitbuffer)
	endfor
	readu,lun,chiantisize
	chiantibuffer=bytarr(chiantisize)
	readu,lun,chiantibuffer
	chiantifile=string(chiantibuffer)
	readu,lun,abundsize
	abundbuffer=bytarr(abundsize)
	readu,lun,abundbuffer
	abundfile=string(abundbuffer)
endelse


data=fltarr(dim+nvars,ng)
temp=0.
for i=0,dim+nvars-1 do begin
	for j=0,ng-1 do begin
		readu,lun,temp
		data[i,j]=temp
	endfor
endfor

close,lun
free_lun,lun
return,data

end
