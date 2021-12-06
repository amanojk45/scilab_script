clear;
clc;
c = 0.9187;
k1 = 0.00017;
k=k1/(c*(1-c));
C=1.050*3600;
SOC1(1)=1;
y1(1)=c*C; //Charge
y2(1)=(1-c)*C;
Cav(1)=C/3600;
Rl=12;
a1=-5.331;
a2=0.3053;
a3=7.722*10^(-7);
a4=3.891*10^(-5);
a5=0.003777;
a6=11.09;
b0=0.1381;
b1=-0.6113;
b2=2.4469;
b3=-4.4802;
b4=3.8154;
b5=-1.2281;
c0=0.1092;
c1=-0.7727;
c2=3.6710;
c3=-8.2696;
c4=8.7117;
c5=-3.4389;
d0=170;
d1=-2115;
d2=16148;
d3=-40267;
d4=42161;
d5=-15922;
e0=-0.0005e7;
e1=0.0200e7;
e2=-0.1992e7;
e3=0.9682e7;
e4=-2.5396e7;
e5=3.6604e7;
e6=-2.7227e7;
e7=0.8154e7;

dt=1;
N=4000;
r = exp(-k * dt);
exec('C:\Users\bhagyashree.hajeri\Desktop\Scilab\Basic battery model\unit.sci', -1);
for t=1:0.1:N
    i(t)=2.21*(U(t)-U(t-60)+U(t-120)-U(t-180)+U(t-240)-U(t-300)+...
    U(t-360)-U(t-420)+U(t-480)-U(t-540)+U(t-600)-U(t-660)+U(t-720)-...
    U(t-780)+U(t-840)-U(t-900)+U(t-960)-U(t-1020)+U(t-1080)-U(t-1140)+..
    U(t-1200)-U(t-1260)+U(t-1320)-U(t-1380)+U(t-1440)-U(t-1500)+...
    U(t-1560)-U(t-1620)+U(t-1680)-U(t-1740)+U(t-1800)-U(t-1860)+...
    U(t-1920)-U(t-1980)+U(t-2040)-U(t-2100)+U(t-2160)-U(t-2220)+...
    U(t-2280)-U(t-2340)+U(t-2400)-U(t-2460)+U(t-2520)-U(t-2580)+...
    U(t-2640)-U(t-2700)+U(t-2760)-U(t-2820)+U(t-2880)-U(t-2940)+...
    U(t-3000)-U(t-3060)+U(t-3120)-U(t-3180)+U(t-3240)-U(t-3300)+...
    U(t-3360)-U(t-3420)+U(t-3480)-U(t-3540)+U(t-3600)-U(t-3660)+...
    U(t-3720)-U(t-3780)+U(t-3840)-U(t-3900)+U(t-3960)-U(t-4020));
end
for t=1:0.1:N    
    y(t) = y1(t) + y2(t);
    y1(t+1) = y1(t) * r + ((y(t) * k * c - i(t)) * (1 - r) - i(t) * c * (k * dt - 1 + r)) / k;
    y2(t+1) = y2(t) * r + y(t) * (1 - c) * (1 - r) - i(t) * (1 - c) * (k * dt - 1 + r) / k;

    del(t)=(y2(t)/(1-c)-y1(t)/c)/3600;
    cun(t)=(1-c)*del(t);
    Cav(t)=y1(t)/3600;
    Cunav(t)=y2(t)/3600;
    SOC1(t+1)=SOC1(t)-((i(t)*dt)+cun(t))/C;    
    SOC(t)=SOC1(t)*100;
//    if(SOC1(t+1)<0) then
//        x=cun(t);
//        break
//    end
    Voc(t)=a1*exp(-a2*SOC(t)) + a3*SOC(t)^3 + a4*SOC(t)^2 + a5*SOC(t) + a6;
    Ri(t)=b0+b1*SOC1(t)^1 + b2*SOC1(t)^2 + b3*SOC1(t)^3 + b4*SOC1(t)^4 + b5*SOC1(t)^5;
    Rts(t)=c0+c1*SOC1(t)^1 + c2*SOC1(t)^2 + c3*SOC1(t)^3 + c4*SOC1(t)^4 + c5*SOC1(t)^5;
    Rtl(t)=Rts(t);
   // Rtl(t)=Rts(t);
    Cts(t)=d0+d1*SOC1(t)^1 + d2*SOC1(t)^2 + d3*SOC1(t)^3 + d4*SOC1(t)^4 + d5*SOC1(t)^5;
    Ctl(t)=e0+e1*SOC1(t)^1 + e2*SOC1(t)^2 + e3*SOC1(t)^3 + e4*SOC1(t)^4 + e5*SOC1(t)^5+...
    e6*SOC1(t)^6+e7*SOC1(t)^7;
    VRi(t)=i(t)*Ri(t);
    if (i(t)>0) then
    Vts(t)=Rts(t)*i(t)*((1-r)/(Rts(t)*Cts(t)));
    Vtl(t)=Rtl(t)*i(t)*((1-r)/(Rtl(t)*Ctl(t)));
    end
    if (i(t)==0) then
    Vts(t)=Vts(t-1)*((1-r)/(Rts(t)*Cts(t)));
    Vtl(t)=Vtl(t-1)*((1-r)/(Rtl(t)*Ctl(t)));
    Vtsl(t)=Vts(t)+Vtl(t);
    //disp(t);
    //disp(Vts(t));
    //disp(Vtl(t));
    //disp(Vtsl(t));
    end
    Vbat(t)=Voc(t)-VRi(t)-Vts(t)-Vtl(t);
    if(Vbat(t)<10.5) then
        x=cun(t);
        break
    end
    I(t)=Vbat(t)/Rl;
end
N1=t;
disp(N1);
for t=N1+1:0.1:N
    cun(t)=x;
end
t=1:0.1:N1-1;
figure (0);
plot(t,Vbat(t)');
xgrid;
xlabel('Time (in Sec)');
ylabel('Terminal voltage (in Volt)');
title('Battery terminal Voltage'); 
// 
figure (2)
plot(t,SOC(t)');
xgrid;
xlabel('Time (in Sec)');
ylabel('State of charge (in %)');
title('Battery state of charge');
//t=380:400;
//figure (3)
//plot(t,Vtsl(t)');
//xgrid;
//xlabel('Time (in Sec)');
//ylabel('RC nw voltage (in Volt)');
//title('RC nw Voltage');
figure (4)
plot(SOC(t),Voc(t)');
xgrid;
xlabel('State of charge (in %)');
ylabel('Open circuit voltage (in Volt)');
title('Voc vs SOC');
//figure (6)
//plot(t,(Cav(t)*3600)');
//xgrid;
//xlabel('Time (in Sec)');
//ylabel('available charge (in As)');
//title('Battery available charge');
//figure (7)
//plot(t,I(t)');
//xgrid;
//xlabel('Time (in Sec)');
//ylabel('Output current (in Amp)');
//title('Battery Output current');
t=1:N;
figure (1)
plot(t,(cun(t)*3600)');
xgrid;
xlabel('Time (in Sec)');
ylabel('unavailable charge (in As)');
title('Battery unavailable charge');
figure (5)
plot(t,i(t)');
xgrid;
xlabel('Time (in Sec)');
ylabel('Current (in Amp)');
title('Battery input current');
