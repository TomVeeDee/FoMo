
; GENERATE DATACUBES OF A CYLINDER WITH A FAST SAUSAGE WAVE
pro main

rotmat, dimx, dimy, dimz, ro, re, va, vae, co, ce, bo, be, waka_ini

tubemodes, rho_int=ro, rho_ext=re, valfv_int=va, valv_ext=vae, cs_int=co, cs_ext=ce,  bo_int=bo, be_ext=be, waka_ini=waka_ini, waka_root, ka_root

; set dimensions of your model
dimx=204 & dimy=204

vel_modes, waka_root=waka_root, ka_root=ka_root, dimx=dimx, dimy=dimy, reg0, reg1, reg2, reg3, reg4, gridx, vr_md, vt_md, vz_md, pr_md, rr_md, br_md, bt_md, bz_md
dimt=30

; specify first the wavelength in velmod_wt.pro. Set keyword 'mag' for
; including magnetic field
velmod_wt, waka_root=waka_root, ka_root=ka_root, gridx=gridx, dimt=dimt, reg3=reg3, vr_md=vr_md, aa, wk_rt, ka_rt, kafix, dimz, tarr, gridz, vr_t, vz_t, rr_t, rtot_t, ptot_t, dispr, te_t, br_t, bz_t ;mag=mag

; set model = string with name of treated model:
;      model = 'base' corresponds to ka = 2.24, 
;      model = 'long' corresponds to ka = 1.25
;      model = 'high_res' corresponds to ka = 2.24, high spatial resolution
;      model = 'high_res2d' corresponds to ka = 2.24, high 2D spatial resolution
;      model = 'highT' corresponds to ka = 2.24, high external temperature

;!!! check output directory for following routines

; if smooth profile case then run smprof_cubes.pro:
;smprof_cubes, dimz, dimt, gridx, rtot_t, te_t, rtot_sm_t, te_sm_t
; specify first the wavelength in datacubes_wt.pro (add /smooth keyword for smooth profile)

; define whether to work with all variables (='all') or with single
; variables (='vr', 'vz', 'rh', or 'te')
sngcub = 'all'

; 2d or 3d with no smoothing. set keyword 'mag' for including magnetic field
datacubes_wt, rho_int=ro, rho_ext=re, valfv_int=va, valv_ext=vae, cs_int=co, cs_ext=ce, radius = aa, gridx=gridx, gridz=gridz, dimt=dimt, tarr=tarr, ka_rt=ka_rt, kafix=kafix, wk_rt=wk_rt, vr_t=vr_t, vz_t=vz_t, rtot_t=rtot_t, te_t=te_t, br_t=br_t, bz_t=bz_t, model=model, vr_cube, vz_cube, te_cube, rh_cube, br_cube, bz_cube, sngcub=sngcub, /save_cubes ;,mag=mag

; with smoothing:
;datacubes_wt, rho_int=ro, rho_ext=re, valfv_int=va, valv_ext=vae,
;cs_int=co, cs_ext=ce, radius = aa, gridx=gridx, gridz=gridz, dimt=dimt, tarr=tarr, ka_rt=ka_rt, kafix=kafix, wk_rt=wk_rt, vr_t=vr_t, vz_t=vz_t, rtot_t=rtot_sm_t, te_t=te_sm_t, model=model,vr_cube, vz_cube, te_cube, rh_cube

; for 3d datasets, divide in 2D (radial) slices:
divcubes,dimt=dimt,dimx=dimx,dimy=dimy,model=model,sngcub=sngcub;,/mag,/smooth

print,'for continuing to intensity calculation type .c'
stop

;----------------
; GENERATE INTENSITY CUBES FROM DATA

; set the model:
; set = case considered: 
;       set = 2 -> base model - 171: ka = 2.24 with 171 intensity
;       set = 23 -> base model - 193: ka = 2.24 with 193 intensity
;       set = 21 -> high T model (base model with high external temp.)
;       set = 22 -> smooth model (base model with smooth density prof.)
;       set = 24 -> hgres model (base model with high spatial resolution)
;       set = 25 -> hgres2d model-171 (2D base model with high spatial
;       set = 26 -> hgres2d model-193 (2D base model with high spatial
;       resolution and 193 intensity)
;       set = 3 -> long lambda model (longer wavelength, ka = 1.25)
set = 2

; Proceed with intensity calculation:

;!!!! .r ch_synthetic  -> changed parenthesis for brackets in line 1021

; define line-of-sight angles:
mua_d = [0.,30.,45.,60.]
; define ion:
if set eq 2 then ion = 'fe_9'
if set eq 23 then ion = 'fe_12'
if set eq 3 then ion = 'fe_9'
if set eq 25 then ion ='fe_9'
if set eq 26 then ion = 'fe_12'

; choose whether imaging (=1) or spectroscopic (=0) data:
imaging = 0
; wrapper for intensity calculation:
prl_slices, set=set, mua_d=mua_d, ion=ion, imaging=imaging

end

