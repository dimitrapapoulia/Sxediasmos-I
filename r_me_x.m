clear;clc;


T = 633; %K

bco = exp(-22.25574)*exp(5629.059418/T);
bco2 = exp(-25.67799)*exp(7421.217224/T);
bh2o_bh2_05 = exp(-24.6281291)*exp(10103.43998/T);
k1 = exp(6.192362489)*exp(-13591.53236/T);
k2 = exp(21.8378944)*exp(-18390.66635/T);
k3 = exp(0.086177696)*exp(-12745/T);
K1 = exp(-6.04031207)*exp(11835.117/T);
K2 = exp(-4.672787)*exp(-4774.119/T);
K3 = exp(-1.36752507)*exp(7060.998/T);
k1bco = exp(-16.06338022)*exp(-7962.472937/T);
k1bco_K1 = exp(-10.02306815)*exp(-19797.58994/T);
bcobh2o_bh2_05 = exp(-46.88387181)*exp(15732.4994/T);
bc02bh2o_bh2_05 = exp(-50.3061191)*exp(17524.6572/T);
k2bco2 = exp(-3.840095601)*exp(-10969.44912/T);
k2bco2_K2 = exp(0.832691399)*exp(-6195.330122/T);
k3bco2 = exp(-25.59181123)*exp(-5323.782776/T);
k3bco2_K3 = exp(-24.22428723)*exp(-12384.78078/T);


%Initial values for vectors  -------------------------------------
ya0 = 0.22;
yb0 = 0.74;
ye0 = 1.e-14;
yc0 = 1.e-14;
yd0 = 1.e-14;

P=50.*10^5;
Pa0=P*ya0;
Pb0=P*yb0;
Pc0=P*yc0;
Pd0=P*yd0;
Pe0=P*ye0;

xa=[0.05 0.1 0.15 0.2 0.25 0.3 0.35];
xe=[0.1 0.2 0.3 0.4 0.5 0.6 0.7];

for i=1:7
Pa(i)=ya0*P*(1-0.2*xa(i)/(1-2*0.2*xa(i)*ya0) +...
    ya0*P*(1-0.8*xa(i)));
Pb(i)=ya0*P*((yb0/ya0)-3*0.2*xa(i))/(1-2*0.2*xa(i)*ya0) +...
    ye0*P*((yb0/ye0)-2*xe(i))/(1-2*xe(i)*ye0)+...
    ya0*P*((yb0/ya0)-0.8*xa(i));
Pc(i)=ya0*P*0.2*xa(i)/(1-2*0.2*xa(i)*ya0)+...
    ye0*P*xe(i)/(1-2*xe(i)*ye0);
Pd(i)=ya0*P*0.2*xa(i)/(1-2*0.2*xa(i)*ya0)+...
    ya0*P*0.8*xa(i);
Pe(i)=ye0*P*(1-xe(i))/(1-2*xe(i)*ye0)+...
    ya0*P*0.8*xa(i);

r1(i)= (k1bco*Pe(i)*Pb(i)^1.5-k1bco_K1*Pc(i)*Pb(i)^(-0.5))/(Pb(i)^0.5+bh2o_bh2_05*Pd(i)+bco*Pe(i)*Pb(i)^0.5+bcobh2o_bh2_05*Pe(i)*Pd(i)+bc02bh2o_bh2_05*Pa(i)*Pd(i));
r2(i)= (k2bco2*Pa(i)*Pb(i)-k2bco2_K2*Pe(i)*Pd(i))/(Pb(i)^0.5+bh2o_bh2_05*Pd(i)+bco*Pe(i)*Pb(i)^0.5+bcobh2o_bh2_05*Pe(i)*Pd(i)+bc02bh2o_bh2_05*Pa(i)*Pd(i));
r3(i)= (k3bco2*Pa(i)*Pb(i)^1.5-k3bco2_K3*Pc(i)*Pd(i)*Pb(i)^(-1.5))/(Pb(i)^0.5+bh2o_bh2_05*Pd(i)+bco*Pe(i)*Pb(i)^0.5+bcobh2o_bh2_05*Pe(i)*Pd(i)+bc02bh2o_bh2_05*Pa(i)*Pd(i));
ra(i)= -r2(i)-r3(i);
end 

format long
FA0 = 1191.; %mol/s
R = 8.314; %SI
CA0 = Pa0/(R*T) %mol/cum
CB0 = Pb0/(R*T)
CAR0 = (P*0.04)/(R*T)
Q = FA0/CA0 %cum/s
FB0 = Q*CB0
FAR0 = Q*CAR0

r_1 = -FA0*1./r1;
r_2 = -FA0*1./r2;
r_3 = -FA0*1./r3;
r_a = -FA0*1./ra;

%----------SIMPSON 3/8-----------------------------------------------------
p=2;
n=3*p+1;
s=xa(4)-xa(1);
for i=1:3:n-3
    I(i) = s*(r_a(i)+3*r_a(i+1)+3*r_a(i+2)+r_a(i+3))/8;
end

V=sum(I) %cum
SV=Q/V
%--------------------------------------------------------------------------

plot(xa,r_a)
xlabel('xa')
ylabel('-F,A0/ra')

%subplot(2,2,1)
%plot(xa, r_1)
%xlabel('xa')
%ylabel('-F,A0/r1')

%subplot(2,2,2)
%plot(xa, r_2)
%xlabel('xa')
%ylabel('-F,A0/r2')

%subplot(2,2,3)
%plot(xa,r_3)
%xlabel('xa')
%ylabel('-F,A0/r3')

%subplot(2,2,4)
%plot(xa,r_a)
%xlabel('xa')
%ylabel('-F,A0/ra')