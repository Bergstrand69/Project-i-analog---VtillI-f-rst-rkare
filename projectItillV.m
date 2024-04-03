Rf = 100;
Rs = 50;
I1 = 0.0096;
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

cpi2_comp = 1.5*(AB_0)/(p2*rpi2);
disp('cpi2_comp');
disp(cpi2_comp);
p1_comp = -1/(rpi2*cpi2_comp);

w0_comp = sqrt(p1_comp*p2*(1-AB_0))/(2*pi);
disp('f0_comp:')
disp(w0_comp);
disp('AB(0)')
disp(20*log10(-AB_0));

s = zpk('s');
AB_s = AB_0*1/((1-s/(p1))*(1-s/(p2)));

AB_s_comp = AB_0*1/((1-s/(p1_comp))*(1-s/(p2)));

At_s = at_inf_0*AB_s/(1 - AB_s);

At_s_comp = at_inf_0*AB_s_comp/(1 - AB_s_comp);

w0 = sqrt(p1*p2*(1-AB_0));

n_ph = -w0*w0/(sqrt(2)*w0 + p1 + p2);
disp('fanthom nolla:')
disp(n_ph)

L_ph = -Rf/n_ph;

B_L_add = (1 - s/(-(Rf)/L_ph))/(1 - s/(-(Rf + Rs)/L_ph));

ab_s_comp_L_ph = B_L_add*AB_0/((1-s/(p1))*(1-s/(p2)));
at_inf_s = -1/(Rf + L_ph*s);
At_s_comp_L_ph = at_inf_s*ab_s_comp_L_ph/(1 - ab_s_comp_L_ph);

figure(4)
hold on
rlocus(B_L_add/((1-s/(p1))*(1-s/(p2))));
hold off
disp(pole(At_s_comp_L_ph));

figure(1)
hold on
bode(AB_s_comp)
bode(AB_s)
bode(ab_s_comp_L_ph)
hold off

figure(2)
hold on
bode(At_s_comp)
bode(At_s)
bode(At_s_comp_L_ph)
hold off

figure(3)
hold on
step(At_s_comp)
step(At_s)
step(At_s_comp_L_ph)
hold off

