pro update_tables,path=path,spectral_lines=spectral_lines,aia=aia,eit=eit

if (~keyword_set(path)) then path='.'

; find list of all files, and look for aia files (need special treatment) and normal chianti files
listoffiles=file_search(path+'/'+'goft_table_*.dat')
aiafiles=listoffiles[where(strmatch(listoffiles,'*aia*'))]
specfiles=listoffiles[where(~strmatch(listoffiles,'*aia*'))]

if keyword_set(spectral_lines) then begin
	; estimated time (with 4 cores), 3 days
	for i=0,n_elements(specfiles)-1 do begin
		; the creation of dummy files allows for parallel updating of the files
		; check if the file was processed already or claimed for working
		if ~(file_test(specfiles[i]+'.working') || file_test(specfiles[i]+'.updated')) then begin
			; claim the file to work on by creating a temporary file 
			spawn,'touch '+specfiles[i]+'.working'
			; initialise variables for reading
			ion=''
			w0=0.
			print,'Updating ',specfiles[i]
			; for each file, read the ion and rest wavelength
			openr,unit,specfiles[i],/get_lun
			readf,unit,ion
			readf,unit,w0
			free_lun,unit

			; even though all spectroscopic files should be done with 
			; sun_coronal_2012_schmelz.abund (because it is normalised in FoMo is other 
			; abundances are needed), let's check if the file name has abco. 
			if (strmatch(specfiles[i],'*abco*')) then begin
				goft_table,ion=ion,w0=w0,gotdir=path+'/',fact=100
				; the goft_table routine immediately overwrites the old file, and if this script fails, then we're left with an empty file
				; should we therefore store the result in a temp directory and then move it to this location?
			endif
			; create dummy file to indicate that the file has been updated, and remove dummy file
			spawn,'touch '+specfiles[i]+'.updated'
			spawn,'rm -f '+specfiles[i]+'.working'
		endif
	endfor
endif

if keyword_set(aia) then begin
	listofwvl=['304', '171', '193', '211', '335', '094', '131', '1600', '1700', '4500']
	listofwvl=['304', '171', '193', '211', '335', '094', '131']
	listofabund=['coronal','photospheric']
	listofabundext=['abco','abph']
	; estimated time, around 4 days
	for i=0,n_elements(listofabund)-1 do begin
		for j=0,n_elements(listofwvl)-1 do begin
			filename=path+'/'+'goft_table_aia'+listofwvl[j]+'_'+listofabundext[i]+'.dat'
			if ~(file_test(filename+'.working') || file_test(filename+'.updated')) then begin
				spawn,'touch '+filename+'.working'
				print,'Updating AIA ', listofwvl[j], ' for ',listofabund[i]
				make_aiaresponse, sngfilter=listofwvl[j],gotdir=path+'/',file_abund=listofabund[i]
				spawn,'touch '+filename+'.updated'
				spawn,'rm -f '+filename+'.working'
			endif
		endfor
	endfor
endif

if keyword_set(eit) then begin
	listofwvl=['304', '171', '195', '284']
	listofabund=['coronal','photospheric']
	listofabundext=['abco','abph']
	; estimated time, around 4 days
	for i=0,n_elements(listofabund)-1 do begin
		for j=0,n_elements(listofwvl)-1 do begin
			filename=path+'/'+'goft_table_eit'+listofwvl[j]+'_'+listofabundext[i]+'.dat'
			if ~(file_test(filename+'.working') || file_test(filename+'.updated')) then begin
				spawn,'touch '+filename+'.working'
				print,'Updating AIA ', listofwvl[j], ' for ',listofabund[i]
				make_aiaresponse, sngfilter=listofwvl[j],gotdir=path+'/',file_abund=listofabund[i]
				spawn,'touch '+filename+'.updated'
				spawn,'rm -f '+filename+'.working'
			endif
		endfor
	endfor
endif

end
