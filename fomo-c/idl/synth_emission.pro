; Consider more robust gaussfile function D.Y. 17 Oct 2014
function readgoftcube,infile,data

; reads the data in infile files, not the triangulation

openr,lun,infile,/get_lun
header=lonarr(4); ndim, datatype, ng,nvars
readf,lun,header
data=fltarr(header[0]+header[3],header[2])
readf,lun,data

close,lun
free_lun,lun
return,data
end


 pro synth_emission, infile,intensity,lvec,doppler,$
           emiss=emiss,sigma=sigma,chisq=chisq,lambda0=lambda0,aia=aia,verbose=verbose
;     datadir='~/'
;     infile=datadir+'propslow_out030t003'
;    aia=0
  if N_elements(infile) eq 0 then $
    message,'synth_emission,infile,intensity,lvec,doppler,emiss=emiss,sigma=sigma,chisq=chisq,aia=aia'
  if ~file_test(infile) then message,infile+' do not exist!'
  if keyword_set(aia) then aia=1 else aia=0 
  if keyword_set(verbose) then verbose=1 else verbose=0
  if aia and (N_elements(doppler) $
	    or N_elements(lvec) $
            or N_elements(emiss) $
            or N_elements(sigma) $
          or N_elements(chisq)) then message,'call synth_emission,infile,intensity,/aia'
  data=readgoftcube(infile); read in the emission file
  ;; reform data into emission profiles. 
  ; it is a 2D array, with 4 columns and N rows
  ; N = nx*ny*nl
  eps=1.0d-4; should be less than spectral resolution
  C0=299792.458d; speed of light km/s
  sizes=size(data)
  Nd=sizes[2]

  temp=where(abs(data[0,*]-data[0,0]) le eps,count); x-dimension
  nx=Nd/count

  temp=where(abs(data[1,*]-data[1,0]) le eps,count); y-dimension 
  ny=Nd/count

  temp=where(abs(data[2,*]-data[2,0]) le eps,count); lambda-dimension
  nl=Nd/count

  xvec=dblarr(nx) & yvec=dblarr(ny)
  lvec=dblarr(nl) & emiss=dblarr(nx,ny,nl)

  j=0 & k=0 & l=0

  for i=0,Nd-1 do begin
	  j=i/nl
	  k=i/nl/nx
	  j-=k*nx
	  l=i mod nl
	  ;if ((j mod nlambda) eq 0) then xvec[j]=data[0,i]
	  ;if ((k mod nlambda*nx) eq 0) then yvec[k]=data[1,i]
	  ;if ((j eq 0) and (k eq 0)) then lvec[l]=data[2,i]
	  xvec[j]=data[0,i]
	  yvec[k]=data[1,i]
	  lvec[l]=data[2,i]
	  emiss[j,k,l]=data[3,i]
  endfor
; fit data as gaussian shape
intensity=dblarr(nx,ny)
if aia then intensity=reform(emiss); simply assign values
if not aia then begin
    doppler=dblarr(nx,ny)
    sigma=dblarr(nx,ny)
    chisq=dblarr(nx,ny)

    if n_elements(lambda0) eq 0 then l0=mean(lvec) else l0=lambda0
    mem=max(emiss)
    errors=mean(emiss)
    ; 50,100, inside the loop
    ; 50 50, outside the loop 
    for i=0,nx-1 do begin
    ;for i=nx/2,nx/2 do begin
	    for j=0,ny-1 do begin
    ;	for j=ny/2,ny/2 do begin
		    localmem=max(emiss[i,j,*])
		    if (localmem ge mem/100000.) then begin
			    yfit=gaussfit(lvec,reform(emiss[i,j,*]),coeff,$
				nterms=3,chisq=ch,estimates=[localmem,l0,(max(lvec)-min(lvec))/4.],$
				measure_errors=replicate(errors,nl))
                          if verbose then begin 
                             print,i,j
                             plot,lvec,reform(emiss[i,j,*]),psym=4
                             oplot,lvec,yfit
                             wait,0.1
                          endif
			  ;  print,i,j,coeff
			    intensity[i,j]=coeff[0]
			    doppler[i,j]=(coeff[1]-l0)/l0*C0 ;  [km/s]
			    sigma[i,j]=abs(coeff[2])/l0*C0 ; km/s sigma/l0*c
			    chisq[i,j]=ch
                           ; print,ch
		    endif 
	    endfor
    endfor
endif
end


