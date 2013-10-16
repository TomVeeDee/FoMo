
FUNCTION dispt1,X

  common vars1, A1, A2, A3, A4, A5, B1, B2, B3, B4, B5, n, ka
  common vars2, amp, rad, wk, moaa, meaa, t_k, th_l, z_j, r_i, aa1

  co = A1 & va = A2 & ct = A3 & ro = A4 & bo = A5
  ce = B1 & vae = B2 & cte = B3 & re = B4 & be = B5
  aa0 = 1. & sigi = -1. & sigk = 1.
  aa = rad
  merr = meaa*(r_i/aa)^2
  disp = X-amp*(vae^2+ce^2)*(wk^2-cte^2)/(wk^2*(wk^2-vae^2))*aa1*n/r_i*beselk(sqrt(sigk*merr),n,/double)/wk*aa/ka*sin(wk*ka*t_k/aa)*sin(n*X)*sin(ka*z_j/aa)/r_i-th_l

  return,disp
end
