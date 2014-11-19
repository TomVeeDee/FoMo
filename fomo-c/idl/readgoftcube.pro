function readgoftcube,emissionsave,data

; reads the data in emissionsave files, not the triangulation

openr,lun,emissionsave,/get_lun
header=lonarr(4)
readf,lun,header
data=fltarr(header[0]+header[3],header[2])
readf,lun,data

close,lun
free_lun,lun
return,data
end
