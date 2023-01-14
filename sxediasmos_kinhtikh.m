 %Forward Euler 
%Solve the system of ODEs
clear;clc
%format long

%times statherwn-----------------------------------------------------------

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


t0=0;  % inital value of independent variable
tf=250;  % final value of independent variable

%Initial values for vectors  ----------------------------------------------
t(1)=t0;
ya0 = 0.143034;
yb0 = 0.66498;
ye0 = 0.16783;
yc0 = 0.;
yd0 = 0.;

P=50*10^5;
Pa0=P*ya0;
Pb0=P*yb0;
Pc0=P*yc0;
Pd0=P*yd0;
Pe0=P*ye0;

Pa(1)=Pa0;
Pb(1)=Pb0;
Pc(1)=Pc0;
Pd(1)=Pd0;
Pe(1)=Pe0;



%SYSTHMA RYTHMWN-----------------------------------------------------------

odefuna=@(t, Pa, Pb, Pc, Pd, Pe) -(k2bco2*Pa*Pb-k2bco2_K2*Pe*Pd)/(Pb^0.5+bh2o_bh2_05*Pd+bco*Pe*Pb^0.5+bcobh2o_bh2_05*Pe*Pd+bc02bh2o_bh2_05*Pa*Pd)...
                                 -(k3bco2*Pa*Pb^1.5-k3bco2_K3*Pc*Pd*Pb^(-1.5))/(Pb^0.5+bh2o_bh2_05*Pd+bco*Pe*Pb^0.5+bcobh2o_bh2_05*Pe*Pd+bc02bh2o_bh2_05*Pa*Pd);

odefunb=@(t, Pa, Pb, Pc, Pd, Pe) -2*(k1bco*Pe*Pb^1.5-k1bco_K1*Pc*Pb^(-0.5))/(Pb^0.5+bh2o_bh2_05*Pd+bco*Pe*Pb^0.5+bcobh2o_bh2_05*Pe*Pd+bc02bh2o_bh2_05*Pa*Pd)...
                                 -(k2bco2*Pa*Pb-k2bco2_K2*Pe*Pd)/(Pb^0.5+bh2o_bh2_05*Pd+bco*Pe*Pb^0.5+bcobh2o_bh2_05*Pe*Pd+bc02bh2o_bh2_05*Pa*Pd)...
                                 -3*(k3bco2*Pa*Pb^1.5-k3bco2_K3*Pc*Pd*Pb^(-1.5))/(Pb^0.5+bh2o_bh2_05*Pd+bco*Pe*Pb^0.5+bcobh2o_bh2_05*Pe*Pd+bc02bh2o_bh2_05*Pa*Pd);

odefunc=@(t, Pa, Pb, Pc, Pd, Pe) (k1bco*Pe*Pb^1.5-k1bco_K1*Pc*Pb^(-0.5))/(Pb^0.5+bh2o_bh2_05*Pd+bco*Pe*Pb^0.5+bcobh2o_bh2_05*Pe*Pd+bc02bh2o_bh2_05*Pa*Pd)...
                                 +(k3bco2*Pa*Pb^1.5-k3bco2_K3*Pc*Pd*Pb^(-1.5))/(Pb^0.5+bh2o_bh2_05*Pd+bco*Pe*Pb^0.5+bcobh2o_bh2_05*Pe*Pd+bc02bh2o_bh2_05*Pa*Pd);

odefund=@(t, Pa, Pb, Pc, Pd, Pe) (k2bco2*Pa*Pb-k2bco2_K2*Pe*Pd)/(Pb^0.5+bh2o_bh2_05*Pd+bco*Pe*Pb^0.5+bcobh2o_bh2_05*Pe*Pd+bc02bh2o_bh2_05*Pa*Pd)...
                                 +(k3bco2*Pa*Pb^1.5-k3bco2_K3*Pc*Pd*Pb^(-1.5))/(Pb^0.5+bh2o_bh2_05*Pd+bco*Pe*Pb^0.5+bcobh2o_bh2_05*Pe*Pd+bc02bh2o_bh2_05*Pa*Pd);

odefune=@(t, Pa, Pb, Pc, Pd, Pe) (k1bco*Pe*Pb^1.5-k1bco_K1*Pc*Pb^(-0.5))/(Pb^0.5+bh2o_bh2_05*Pd+bco*Pe*Pb^0.5+bcobh2o_bh2_05*Pe*Pd+bc02bh2o_bh2_05*Pa*Pd)...
                                 +(k2bco2*Pa*Pb-k2bco2_K2*Pe*Pd)/(Pb^0.5+bh2o_bh2_05*Pd+bco*Pe*Pb^0.5+bcobh2o_bh2_05*Pe*Pd+bc02bh2o_bh2_05*Pa*Pd);


%FORWARD EULER ------------------------------------------------------------

%Step size
dt=1;

i=1;
while true

  t(i+1) = t(i) + dt;
  Pa(i+1)=Pa(i)+odefuna(t(i),Pa(i),Pb(i),Pc(i),Pd(i),Pe(i))*dt;
  Pb(i+1)=Pb(i)+odefunb(t(i),Pa(i),Pb(i),Pc(i),Pd(i),Pe(i))*dt;
  Pc(i+1)=Pc(i)+odefunc(t(i),Pa(i),Pb(i),Pc(i),Pd(i),Pe(i))*dt;
  Pd(i+1)=Pd(i)+odefund(t(i),Pa(i),Pb(i),Pc(i),Pd(i),Pe(i))*dt;
  Pe(i+1)=Pe(i)+odefune(t(i),Pa(i),Pb(i),Pc(i),Pd(i),Pe(i))*dt;
 
  i=i+1;
  if(t(i)>=tf);break;end;
    
end

Pa=Pa';
Pb=Pb';
Pc=Pc';
Pd=Pd';
Pe=Pe';


for i=1:251
r1(i)= (k1bco*Pe(i)*Pb(i)^1.5-k1bco_K1*Pc(i)*Pb(i)^(-0.5))/(Pb(i)^0.5+bh2o_bh2_05*Pd(i)+bco*Pe(i)*Pb(i)^0.5+bcobh2o_bh2_05*Pe(i)*Pd(i)+bc02bh2o_bh2_05*Pa(i)*Pd(i));
r2(i)= (k2bco2*Pa(i)*Pb(i)-k2bco2_K2*Pe(i)*Pd(i))/(Pb(i)^0.5+bh2o_bh2_05*Pd(i)+bco*Pe(i)*Pb(i)^0.5+bcobh2o_bh2_05*Pe(i)*Pd(i)+bc02bh2o_bh2_05*Pa(i)*Pd(i));
r3(i)= (k3bco2*Pa(i)*Pb(i)^1.5-k3bco2_K3*Pc(i)*Pd(i)*Pb(i)^(-1.5))/(Pb(i)^0.5+bh2o_bh2_05*Pd(i)+bco*Pe(i)*Pb(i)^0.5+bcobh2o_bh2_05*Pe(i)*Pd(i)+bc02bh2o_bh2_05*Pa(i)*Pd(i));
ra(i)= r2(i) + r3(i);
rb(i)= 2*r1(i) +r2(i) +3*r3(i);
rc(i)= r1(i) + r3(i);
rd(i)= r2(i) + r3(i);
re(i)= r2(i) - r1(i);
end 

subplot(2,2,1)
plot(t, r1)
xlabel('t')
ylabel('r1')

subplot(2,2,2)
plot(t, r2)
xlabel('t')
ylabel('r2')

subplot(2,2,3)
plot(t, r3)
xlabel('t')
ylabel('r3')




