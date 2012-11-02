
; GENERATE DATACUBES OF A CYLINDER WITH A FAST SAUSAGE WAVE

rotmat, dimx, dimy, dimz, ro, re, va, vae, co, ce, waka_ini

tubemodes, rho_int=ro, rho_ext=re, valfv_int=va, valv_ext=vae, cs_int=co, cs_ext=ce, waka_ini=waka_ini, waka_root, ka_root
dimx=204*3 & dimy=204*3

vel_modes, waka_root=waka_root, ka_root=ka_root, dimx=dimx, dimy=dimy, reg0, reg1, reg2, reg3, reg4, gridx, vr_md, vt_md, vz_md, pr_md, rr_md
dimt=30

;specify first the wavelength in velmod_wt.pro
velmod_t, waka_root=waka_root, ka_root=ka_root, gridx=gridx, dimt=dimt, reg3=reg3, vr_md=vr_md, aa, wk_rt, ka_rt, k_ind, dimz, tarr, gridz, vr_t, vz_t, rr_t, rtot_t, ptot_t, dispr, te_t

; if smooth profile case then run smprof_cubes.pro:
smprof_cubes, dimz, dimt, gridx, rtot_t, te_t, rtot_sm_t, te_sm_t
; specify first the wavelength in datacubes_wt.pro (add /smooth keyword for smooth profile)
; with no smoothing:
datacubes, rho_int=ro, rho_ext=re, valfv_int=va, valv_ext=vae, cs_int=co, cs_ext=ce, radius = aa, gridx=gridx, gridz=gridz, dimt=dimt, tarr=tarr, ka_rt=ka_rt, wk_rt=wk_rt, vr_t=vr_t, vz_t=vz_t, rtot_t=rtot_t, te_t=te_t, vr_cube, vz_cube, te_cube, rh_cube, /save_cubes)

; with smoothing:
;datacubes, rho_int=ro, rho_ext=re, valfv_int=va, valv_ext=vae, cs_int=co, cs_ext=ce, radius = aa, gridx=gridx, gridz=gridz, dimt=dimt, tarr=tarr, ka_rt=ka_rt, wk_rt=wk_rt, vr_t=vr_t, vz_t=vz_t, rtot_t=rtot_sm_t, te_t=te_sm_t, vr_cube, vz_cube, te_cube, rh_cube

divcubes,dimt=dimt,dimx=dimx,dimy=dimy,model=model,[/smooth]

----------------
; GENERATE INTENSITY CUBES FROM DATA

; A call to these 2 routines is first needed to load functions in memory (otherwise a strange error comes out)
.r params_r
.r lineongrid_goft_r
set=24
ix=0
it=0

params, set, it, ix, ro, re, va, vae, co, ce, aa, r0, gridx, gridy, gridz, dimx, dimy, dimz, tarr, kafix, ka_rt, wk_rt, te, rho, vr, vz, velx, vely, velz

lineongrid, rho, te, wayemi=2

retall

; Proceed with intensity calculation:

.r params
.r lineongrid_goft_tab
.r integrateemission
.r gridlos
.r prl_slices

prl_slices

end

