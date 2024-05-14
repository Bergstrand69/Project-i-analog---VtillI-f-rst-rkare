Rf = 100;
Rs = 50;
I1 = 0.010;
I2 = 0.004;
Vt = 0.025;
Bf = 200;
rpi1prim = 4*Bf*Vt/I1;
rpi2 = Bf*Vt/I2;
at_inf_0 = -1/Rf;   

cpi2 = 0.0000022;
cpi1 = 100*10^(-9);
AB_0 = -Bf*Bf*Rf/(Rf + Rs + rpi1prim);
p1 = -1/(rpi2*cpi2);

p2 = -(Rf + Rs + rpi1prim)/((Rf + Rs)*cpi1*rpi1prim);

disp(p1)

disp(p2)



phaseShift = pi/180 * 65.5;
cpi2_comp = sin(phaseShift)/tan(pi/2-phaseShift)*(AB_0)/(p2*rpi2);
disp('cpi2_comp');
disp(cpi2_comp);
p1_comp = -1/(rpi2*cpi2_comp);
disp('ny pol');
disp(p1_comp);

w0_comp = sqrt(p1_comp*p2*(1-AB_0))/(2*pi);
disp('f0_comp_calc:')
disp(w0_comp);



disp('AB(0)');
disp(20*log10(-AB_0));



s = zpk('s');
AB_s = AB_0*1/((1-s/(p1))*(1-s/(p2)));



AB_s_comp = AB_0*1/((1-s/(p1_comp))*(1-s/(p2)));

At_s = at_inf_0*AB_s/(1 - AB_s);

disp(bandwidth(At_s))

At_s_comp = at_inf_0*AB_s_comp/(1 - AB_s_comp);
disp('f0_comp:');
disp(bandwidth(At_s_comp)/(2*pi));

w0 = sqrt(p1*p2*(1-AB_0));

n_ph = -w0*w0/(sqrt(2)*w0 + p1 + p2);
disp('fanthom nolla:');
disp(n_ph);

L_ph = -Rf/n_ph;
c_ph = -1/(Rs*n_ph);

disp('L_ph');
disp(L_ph);
disp('c_ph');
disp(c_ph);

B_L_add = (1 - s/(-(Rf)/L_ph))/(1 - s/(-(Rf + Rs)/L_ph));

B_C_add = (1+ Rs*c_ph*s)/(1+ (Rf*Rs*c_ph*s)/(Rf + Rs));

ab_s_comp_L_ph = B_L_add*AB_0/((1-s/(p1))*(1-s/(p2)));
at_inf_s = -1/(Rf + L_ph*s);
At_s_comp_L_ph = at_inf_s*ab_s_comp_L_ph/(1 - ab_s_comp_L_ph);
disp('Systempoles');
disp(pole(At_s_comp));

ab_s_comp_C_ph = B_C_add * AB_0/((1-s/(p1))*(1-s/(p2)));
At_s_comp_C_ph = (-1/(Rf*(1+c_ph*Rs*s)))*ab_s_comp_C_ph/(1 - ab_s_comp_C_ph);

[gm, pm] = margin(-AB_s_comp);
disp('Phase Margin');
disp(pm);

figure(4)
hold on
rlocus(B_C_add/((1-s/(p1))*(1-s/(p2))));
hold off
disp(pole(At_s_comp_L_ph));

figure(1)

hold on
bode(AB_s_comp, 'y')
bode(AB_s, 'g')
bode(ab_s_comp_L_ph, 'b')
bode(ab_s_comp_C_ph, 'r--')
hold off
lgd = legend('Narrow banding', 'Okompenserad', 'Phantom kompenserad med spole','Phantom kompenserad med kondensator', 'Location','southwest');
fontsize(lgd,'decrease')
title('Sling förstärkningen AB');


figure(2)
hold on
bode(At_s_comp, 'y')
bode(At_s, 'g')
bode(At_s_comp_L_ph, 'b')
bode(At_s_comp_C_ph, 'r--')
title('Slutna överföringen');
lgd = legend('Narrow banding','Okompenserad', 'Phantom kompenserad med spole','Phantom kompenserad med kondensator', 'Location','southwest');
hold off


figure(3)
hold on
step(At_s_comp, 'y')
step(At_s, 'g')
step(At_s_comp_L_ph, 'b')
step(At_s_comp_C_ph, 'r--')
title('Steg svar för den slutna överföringen');
legend('Narrow banding','Okompenserad', 'Phantom kompenserad med spole','Phantom kompenserad med kondensator');
hold off




