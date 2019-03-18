function readgoftcube_txt,emissionsave,data,compress=compress

; reads the data in emissionsave files, not the triangulation
; invoke as data=readgoftcube('fileyouwanttoread')
; data will be a 4xn array, with n the number of grid points in x,y,\lambda

openr,lun,emissionsave,/get_lun,compress=compress
versionheader=strarr(1)
readf,lun,versionheader

if (strpos(versionheader,'FoMo-v') eq -1) then begin
; this is the version of the file pre-3.4
	point_lun,lun,0
	header=lonarr(3)
	readf,lun,header
	stringheader=strarr(2)
	readf,lun,stringheader
endif else begin
; this is the new version of the file (v3.4)
	header=lonarr(3)
        readf,lun,header
	stringheader=strarr(3)
        readf,lun,stringheader
endelse
	
data=fltarr(header[0]+header[2],header[1])
readf,lun,data

close,lun
free_lun,lun
return,data

end
