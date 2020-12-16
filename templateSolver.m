

syms tp dtp d2tp tm dtm d2tm  % to solve for
syms k1 k2 dx r T2m Tm Tp T2p % known

temp = tp == tm;                                % temperature continuity
flux = k1*(dtm + r*d2tm) == k2*(dtp + r*d2tp);  % heat flux continuity
t1 = Tm == tm - dx*dtm + 0.5*dx^2 * d2tm;       % Taylor expansion
t2 = T2m == tm - 2*dx*dtm + 2*dx^2 * d2tm;      %   "        "
t3 = Tp == tp + dx*dtp + 0.5*dx^2 * d2tp;       %   "        "
t4 = T2p == tp + 2*dx*dtp + 2*dx^2 * d2tp;      %   "        "

soln = solve([temp, flux, t1, t2, t3, t4], [tp dtp d2tp tm dtm d2tm]);

clear tp dtp d2tp tm dtm d2tm k1 k2 dx T2m Tm Tp T2p temp flux t1 t2 t3 t4
% flux = k1*dtm == k2*dtp;
