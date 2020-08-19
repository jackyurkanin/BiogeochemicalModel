% Main body function for biogeochemical model Cooley et al 2015
% using Georges Bank (GB)
% begin by initilizing all initial variables
TOTAL = ones(9,991); %r1 =t1 r2 =t2 r3 =s1 r4 = s2 r5 = DICs r6 = DICd r7 = TAs r8 = TAd
T1 = 23.29; T2 = 8.52;  % Temperature
S1 = 31.67; S2 = 32.95; % Salinity
ResT1 = 9.51; ResT2 = 10.00;
ResS1 = 32.47; ResS2 = 32.44;
H1 = 25.00; H2 = 45.00;
AT1 = 9.00; AT2 = 3.50;
AS1 = -0.75;
TauT1 = 0.09;
Tau2 = 0.15;
TauS1 = 0.27;
PhaseS1 = -0.45;
PhaseT1 = -0.97;
PhaseT2 = -1.54;
Kd = 10^-5;
Fratio = 0.25;
PICPOC = 0.04;
Remin = 0.80;
TAs = 2160.6;
DICs = 1890.1;
TAd = 2125.7;
DICd = 2019.0;
for t= 1:0.1:100
    display(t+2000)
    WindSpeedAve = 54.23*0.1-33.84*cos(2*pi*t+4.14);
    Sc = 668; % Schmidt number
    Kair = 0.251*WindSpeedAve^2*(Sc/660)^-0.5; % wanninkhof
    NPP = 1.04+0.58*cos(2*pi*t-2.11);
    NCP = NPP*Fratio;
    
    % Calculations
    
    dT1dt = (ResT1 - AT1*cos(2*pi*t+PhaseT1)-T1)/TauT1; %Eq. SI1
    dT2dt = (ResT2 - AT2*cos(2*pi*t+PhaseT2)-T2)/Tau2; %Eq. SI2
    dS1dt = (ResS1 - AS1*cos(2*pi*t+PhaseS1)-S1)/TauS1; %Eq. SI3
    dS2dt = (ResS2 - S2)/Tau2; %Eq. SI4
    % calculate density for surfacce box
    den11=999.842594+T1*(6.793952*10^-2)+T1^2*(-9.09529*10^-3)+...
        T1^3*(1.001685*10^-4);
    den21=den11+T1^4*(-1.120083*10^-6)+T1^5*(6.536332*10^-9);
    den31=den21+S1*(0.824493+T1*(-4.0899*10^-3)+T1^2*(7.6438*10^-5)+...
        T1^3*(-8.2467*10^-7)+T1^4*(5.3875*10^-9));
    density1=den31+S1^(3/2)*(-5.72466*10^-3+T1*1.0227*10^-4+T1^2*...
        (-1.6546*10^-6))+S1^2*(4.8314*10^-4);
    
    % calculate density for deeper box
    den12=999.842594+T2*(6.793952*10^-2)+T2^2*(-9.09529*10^-3)+...
        T2^3*(1.001685*10^-4);
    den22=den12+T2^4*(-1.120083*10^-6)+T2^5*(6.536332*10^-9);
    den32=den22+S2*(0.824493+T2*(-4.0899*10^-3)+T2^2*(7.6438*10^-5)+...
        T2^3*(-8.2467*10^-7)+T2^4*(5.3875*10^-9));
    density2=den32+S2^(3/2)*(-5.72466*10^-3+T2*1.0227*10^-4+T2^2*...
        (-1.6546*10^-6))+S2^2*(4.8314*10^-4);
    
    T1 = T1+dT1dt*0.1; T2 = dT2dt*0.1+T2; S1 = dS1dt*0.1+S1; S2 = dS2dt*0.1+S2;
    
    if density1 > density2 % check if the boxes mix
        T1 = (H1*T1+H2*T2)/(H1+H2);
        T2 = T1;
        S1 = (H1*S1+H2*S2)/(H1+H2);
        S2 =S1;
    end
    Press = 10.1325*(1 + H1/10); 
    [result] = CO2SYS(TAs,DICs,1,2,S1,T1,T1, Press,Press, 0,0,1,4,1);
    pCO2surf = result(19); pCO2atm = result(32);
    dDICdz = DICs-DICd; dTAdz = TAs-TAd;
    aCO2 = exp(-58.0931+90.5069*(100/(273+T1))+22.2940*log((T1+273)/100)+S1*...
        (0.027766-0.025888*((273+T1)/100)+0.0050578*((273+T1)/100)^2));
    Jair = Kair*aCO2*(pCO2atm-pCO2surf);% establish Jair
    dDIC1dt = (Jair+Kd*dDICdz-NCP-NCP*PICPOC)/H1; %Eq. SI7
    dDIC2dt = (Remin*NCP+NCP*PICPOC-Kd*dDICdz)/H2; %Eq. SI8
    dTA1dt = (Kd*dTAdz-2*NCP*PICPOC)/H1; %Eq. SI9
    dTA2dt = (-Kd*dTAdz+2*NCP*PICPOC)/H2; %Eq. SI10
    DICs = DICs+dDIC1dt*0.1; DICd = DICd+dDIC2dt*0.1; TAs = TAs+dTA1dt*0.1; TAd = TAd+dTA2dt*0.1;
    ValT = round(t*10)-9;
    TOTAL(1,ValT) = T1; TOTAL(2,ValT) = T2; TOTAL(3,ValT) = S1; TOTAL(4,ValT) = S2;
    TOTAL(5,ValT) = DICs; TOTAL(6,ValT) = DICd; TOTAL(7,ValT) = TAs; TOTAL(8,ValT) = TAd;
    TOTAL(9, ValT) = pCO2surf;
end
figure(1);
subplot(2,4,1); plot([2001:0.1:2100], TOTAL(1,:)); title('Surf Temp')
subplot(2,4,2); plot([2001:0.1:2100], TOTAL(2,:)); title('Deep Temp')
subplot(2,4,3); plot([2001:0.1:2100], TOTAL(3,:)); title('Surf Salinity')
subplot(2,4,4); plot([2001:0.1:2100], TOTAL(4,:)); title('Deep Salinity')
subplot(2,4,5); plot([2001:0.1:2100], TOTAL(5,:)); title('Surf DIC')
subplot(2,4,6); plot([2001:0.1:2100], TOTAL(6,:)); title('Deep DIC')
subplot(2,4,7); plot([2001:0.1:2100], TOTAL(7,:)); title('Surf TA')
subplot(2,4,8); plot([2001:0.1:2100], TOTAL(8,:)); title('Deep TA')
figure(2);
plot([2001:0.1:2100], TOTAL(9,:)); title('pCO2surf')