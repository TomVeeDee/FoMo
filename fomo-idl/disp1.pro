
FUNCTION disp1,X

  common vars1, A1, A2, A3, A4, A5, B1, B2, B3, B4, B5, n, ka
  common vars2, amp, rad, wk, moaa, t_k, th_l, z_j, r_i
  
  co = A1 & va = A2 & ct = A3 & ro = A4 & bo = A5
  ce = B1 & vae = B2 & cte = B3 & re = B4 & be = B5
  aa0 = 1. & sigi = -1. & sigk = 1.
  aa = rad
  morr = moaa*(X/aa)^2
  disp = X+amp*((va^2+co^2)*(wk^2-ct^2)/(wk^2*(wk^2-va^2))*aa0*sqrt(sigi*moaa/aa^2)*dbeselj(sqrt(sigi*morr),n)*(aa/ka)^2)/wk*aa/ka*sin(wk*ka*t_k/aa)*cos(n*th_l)*sin(ka*z_j/aa)-r_i
  return,disp
end
