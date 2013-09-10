function MichaelisMenten
clear
clc
%--------------------------------------------------------------------------
T = 10; %time

%----------------------execution-------------------------------------------
% substrate = y(1), enzyme2 = y(2), intermediate = y(3), product = y(4)

 options = odeset('AbsTol',0.000001);
[tt,z] = ode15s(@lumped,[0 T],ones(24,1),options);%,option

%---------------------- graphical output ----------------------------------

% plot(t,y(:,1:4));
% hold on;
% plot(tt,z(:,1),'--',tt,s0-z(:,1)+p0,'--');
for i=1:size(z,2)
plot(tt,z(:,i));hold on;
end
% legend ('substrate','enzyme','intermediate','product','lumped substrate','lumped product');
xlabel('time'); ylabel('concentration');
grid;

%---------------------- functions ----------------------------------------
function dcdt = lumped(t,z)
dzdt=zeros(24,1);

% G6P=z(1);F6P=z(2);FBP=z(3);DHAP=z(4);G3P=z(5);BPG=z(6);PG3=z(7);PG2=z(8);PEP=z(9);Pyr=z(10);AcCoA=z(11);Citrate=z(12);Isocitrate=z(13);Ketoglutarate=z(14);SuccinylCoA=z(15);Succinate=z(16);
% Fumarate=z(17);Malate=z(18);Oxaloacetate=z(19);PGL=z(20);PGC=z(21);Ru5P=z(22);X5P=z(23);R5P=z(24);S7P=z(25);E4P=z(26);ATP=z(27);ADP=z(28);NAD=z(29);NADH=z(30);NADP=z(31);NADPH=z(32);DAHP=z(33);



G6P=z(1);
F6P=z(2);
FBP=z(3);DHAP=z(4);G3P=z(5);BPG=z(6);PG3=z(7);PG2=z(8);PEP=z(9);Pyr=z(10);
Oxaloacetate=z(11);PGL=z(12);Ru5P=z(13);X5P=z(14);R5P=z(15);S7P=z(16);E4P=z(17);ATP=z(18);ADP=z(19);NAD=z(20);NADH=z(21);NADP=z(22);NADPH=z(23);DAHP=z(24);

%dzdt(1) = -k(1)*k(2)*(e0+i0)*z/(k(2)+k(3)+k(1)*z);
Vf=11.09;Vr=0.31;KmA=60.25;KmP=60.57;%pgi for G6P
dzdtpgi= (Vf*G6P-Vr*F6P)/(1+G6P/KmA+F6P/KmP)
 
Vf=135.66;Vr=16.07;KmA=7.06;KmB=0.32;KmP=0.06;KmQ=1.29;KiA=0.84;KiB=0.38;KiP=3.11;KiQ=58.60;%pfk for F6P and ATP
 dzdtpfk=((Vf*F6P*ATP)/(KiA*KmB)-(Vr*FBP*ADP)/(KiQ*KmP))/(1+F6P/KiA+(KmA*ATP)/(KiA*KmB)+FBP*KmQ/KmP*KiQ+ADP/KiQ+(F6P*ATP)/(KiA*KmB)+(F6P*FBP*KmQ)/(KiA*KmP*KiQ)+(KmA*ATP*ADP)/(KiA*KmB*KiQ)+(FBP*ADP)/(KmP*KiQ)+(F6P*FBP*ATP)/(KiA*KmB*KiP)+(ATP*FBP*ADP)/(KiB*KmP*KiQ));

 Vf=69.91;Vr=1.00;KmA=1.11;KmP=10.87;KmQ=0.15;KiA=11.06;KiP=12.28;KiQ=11.95;%Ald F FBP
  dzdtald=((Vf*FBP)/(KiA*KmB)-(Vr*DHAP*G3P)/(KiQ*KmP))/(1+FBP/KiA+(KmQ*DHAP)/(KiQ*KmP)+G3P/KiQ+FBP/KiA*KmA+(FBP*DHAP*KmQ)/(KiA*KmP*KiQ)+(DHAP*G3P)/(KiQ*KmP)+(FBP*DHAP)/(KiA*KmA*KiP)+G3P/(KiQ*KiA));

 Vf=56.87;Vr=38.34;KmA=0.99;KmP=1.78;%Tim For DHAP
dzdttim=(Vf*DHAP-Vr*G3P)/(1+DHAP/KmA+G3P/KmP);
 
Vf=671.72;Vr=370.47;KmA=1.25;KmB=1.38;KmP=0.03;KmQ=0.22;KiA=2.89;KiB=2.43;KiP=7.60;KiQ=127.17;%Gapdh For GAP and NAD
 dzdtgapdh=((Vf*G3P*NAD)/(KiA*KmB)-(Vr*BPG*NADH)/(KiQ*KmP))/(1+G3P/KiA+(KmA*NAD)/(KiA*KmB)+BPG*KmQ/KmP*KiQ+NADH/KiQ+(G3P*NAD)/(KiA*KmB)+(G3P*BPG*KmQ)/(KiA*KmP*KiQ)+(KmA*NAD*NADH)/(KiA*KmB*KiQ)+(BPG*NADH)/(KmP*KiQ)+(G3P*BPG*NAD)/(KiA*KmB*KiP)+(NAD*BPG*NADH)/(KiB*KmP*KiQ));

 Vf=2225.00;Vr=589.80;KmA=1.40;KmB=2.00;KmP=0.3;KmQ=0.9;KiA=1.9;KiB=3.1;KiP=0.5;KiQ=4.00;%pgk For BPG and ADP
 dzdtpgk=((Vf*BPG*ADP)/(KiA*KmB)-(Vr*PG3*ATP)/(KiQ*KmP))/(1+BPG/KiA+(KmA*ADP)/(KiA*KmB)+PG3*KmQ/KmP*KiQ+ATP/KiQ+(BPG*ADP)/(KiA*KmB)+(BPG*PG3*KmQ)/(KiA*KmP*KiQ)+(KmA*ADP*ATP)/(KiA*KmB*KiQ)+(PG3*ATP)/(KmP*KiQ)+(F6P*PG3*ADP)/(KiA*KmB*KiP)+(ADP*PG3*ATP)/(KiB*KmP*KiQ));

 Vf=90.55;Vr=6.21;KmA=1.59;KmP=5.17;%pgm For3PG
dzdtpgm=(Vf*PG3-Vr*PG2)/(1+PG3/KmA+PG2/KmP);

 Vf=355.79;Vr=4.98;KmA=4.54;KmP=0.87;%eno For 2PG
dzdteno=(Vf*PG2-Vr*PEP)/(1+PG2/KmA+PEP/KmP);

Vf=109.60;Vr=0.80;KmA=2.80;KmB=1.60;KmP=33.60;KmQ=2.60;KiA=10.00;KiB=0.90;KiP=2629.40;KiQ=356.40;%G6pdh For G6P and NADH
 dzdtg6pdh=((Vf*G6P*NADP)/(KiA*KmB)-(Vr*PGL*NADPH)/(KiQ*KmP))/(1+G6P/KiA+(KmA*NADP)/(KiA*KmB)+PGL*KmQ/KmP*KiQ+NADPH/KiQ+(G6P*NADP)/(KiA*KmB)+(G6P*PGL*KmQ)/(KiA*KmP*KiQ)+(KmA*NADP*NADPH)/(KiA*KmB*KiQ)+(PGL*NADPH)/(KmP*KiQ)+(G6P*PGL*NADP)/(KiA*KmB*KiP)+(NADP*PGL*NADPH)/(KiB*KmP*KiQ));

 Vf=75.65;Vr=564.17;KmA=0.29;KmB=0.01;KmP=1.03;KmQ=9.24;KiA=5.16;KiB=0.18;KiP=1.39;KiQ=1.27;%Gnd for PGL and NADP
 dzdtgnd=((Vf*PGL*NADP)/(KiA*KmB)-(Vr*Ru5P*NADPH)/(KiQ*KmP))/(1+PGL/KiA+(KmA*NADP)/(KiA*KmB)+Ru5P*KmQ/KmP*KiQ+NADPH/KiQ+(PGL*NADP)/(KiA*KmB)+(PGL*Ru5P*KmQ)/(KiA*KmP*KiQ)+(KmA*NADP*NADPH)/(KiA*KmB*KiQ)+(Ru5P*NADPH)/(KmP*KiQ)+(PGL*Ru5P*NADP)/(KiA*KmB*KiP)+(NADP*Ru5P*NADPH)/(KiB*KmP*KiQ));

 Vf=21.05;Vr=10.62;KmA=0.41;KmP=1.03;%Rpe For Ru5P
dzdtrpe=(Vf*Ru5P-Vr*X5P)/(1+Ru5P/KmA+X5P/KmP);

Vf=15.95;Vr=4.31;KmA=0.09;KmP=0.09;%Rpi For Ru5P
dzdtrpi=(Vf*Ru5P-Vr*R5P)/(1+Ru5P/KmA+R5P/KmP);

 Vf=55.57;Vr=48.08;KmA=2.49;KmB=0.25;KmP=32.69;KmQ=1.33;KiA=3.18;KiB=0.58;KiP=80.45;KiQ=2.09;%TktAB1 For X5P and R5P
 dzdttktab1=((Vf*X5P*R5P)/(KiA*KmB)-(Vr*S7P*G3P)/(KiQ*KmP))/(1+X5P/KiA+(KmA*R5P)/(KiA*KmB)+G3P*KmQ/KmP*KiQ+S7P/KiQ+(X5P*R5P)/(KiA*KmB)+(X5P*G3P*KmQ)/(KiA*KmP*KiQ)+(KmA*R5P*S7P)/(KiA*KmB*KiQ)+(G3P*S7P)/(KmP*KiQ)+(X5P*R5P*G3P)/(KiA*KmB*KiP)+(G3P*R5P*S7P)/(KiB*KmP*KiQ));

  Vf=15.45;Vr=2.04;KmA=1.02;KmB=0.38;KmP=0.82;KmQ=0.97;KiA=0.14;KiB=0.4;KiP=11.55;KiQ=13.83;%TktAB2 For X5P and E4P
 dzdttktab2=((Vf*X5P*E4P)/(KiA*KmB)-(Vr*G3P*F6P)/(KiQ*KmP))/(1+X5P/KiA+(KmA*E4P)/(KiA*KmB)+F6P*KmQ/KmP*KiQ+G3P/KiQ+(X5P*E4P)/(KiA*KmB)+(X5P*F6P*KmQ)/(KiA*KmP*KiQ)+(KmA*E4P*G3P)/(KiA*KmB*KiQ)+(F6P*G3P)/(KmP*KiQ)+(X5P*E4P*F6P)/(KiA*KmB*KiP)+(F6P*E4P*G3P)/(KiB*KmP*KiQ));

 Vf=16.57;Vr=11.23;KmA=0.94;KmB=0.57;KmP=2.10;KmQ=10.35;KiA=0.72;KiB=1.22;KiP=9.18;KiQ=23.91;%Tal For S7P and GAP
 dzdttal=((Vf*G3P*S7P)/(KiA*KmB)-(Vr*E4P*F6P)/(KiQ*KmP))/(1+G3P/KiA+(KmA*S7P)/(KiA*KmB)+F6P*KmQ/KmP*KiQ+E4P/KiQ+(G3P*S7P)/(KiA*KmB)+(G3P*F6P*KmQ)/(KiA*KmP*KiQ)+(KmA*S7P*E4P)/(KiA*KmB*KiQ)+(F6P*E4P)/(KmP*KiQ)+(G3P*S7P*F6P)/(KiA*KmB*KiP)+(F6P*S7P*E4P)/(KiB*KmP*KiQ));

  Vf=10.17;Vr=0.59;KmA=0.29;KmB=0.39;KmP=1.48;KmQ=0.17;KiA=0.02;KiB=0.02;KiP=20.91;KiQ=4.59;%pk For PEP and ADP
 dzdtpk=((Vf*PEP*ADP)/(KiA*KmB)-(Vr*Pyr*ATP)/(KiQ*KmP))/(1+PEP/KiA+(KmA*ADP)/(KiA*KmB)+Pyr*KmQ/KmP*KiQ+ATP/KiQ+(PEP*ADP)/(KiA*KmB)+(PEP*Pyr*KmQ)/(KiA*KmP*KiQ)+(KmA*ADP*ATP)/(KiA*KmB*KiQ)+(Pyr*ATP)/(KmP*KiQ)+(PEP*Pyr*ADP)/(KiA*KmB*KiP)+(ADP*Pyr*ATP)/(KiB*KmP*KiQ));

 Vf=13.47;Vr=1.41;KmA=2.11;KmP=24.97;%Ppc For PEP
dzdtppc=(Vf*PEP-Vr*Oxaloacetate)/(1+PEP/KmA+Oxaloacetate/KmP);

 Vf=1.00;Vr=0.001;KmA=0.12;KmB=0.002;KmP=0.02;KiA=0.0001;KiB=0.002;KiP=83.35;%AroG For PEP and E4P
 dzdtarog=((Vf*PEP*E4P)/(KiA*KmB)-(Vr*DAHP)/(KiP*KmP))/(1+PEP/KiA+(KmA*E4P)/(KiA*KmB)+DAHP/KiP+(PEP*E4P)/(KiA*KmB)+PEP/(KiA*KiP)+(KmA*E4P*DAHP)/(KiA*KmB*KiP)+DAHP/(KmP*KiP)+(DAHP*E4P)/(KiB*KmP*KiP));

 
 
dzdtsersynth=1.75*PG3; % SerSynth For 3PG

dzdtsynth1=1.41*PEP;% Synth1 For PEP

dzdtsynth2=5.36*Pyr; % Synth2 For PYR

 dzdtpdh=18.80*Pyr;%Pdh For PYR

dzdtrpkk=1.03*R5P; % Rpkk For R5P

dzdtdahpout=0.69*DAHP;%DAHP out For DAHP

dzdtoaaout=4.27*Oxaloacetate;%OAA out For OAA
 
% Mass Balance
dcdtG6P=dzdtpgi+dzdtg6pdh
dzdtpgi
dzdtg6pdh
dcdtF6P=dzdtpfk-dzdtpgi-dzdttktab1-dzdttal;
dcdtFBP=dzdtald-dzdtpfk;
dcdtDHAP=dzdttim-dzdtald;
dcdtG3P=dzdtgapdh+dzdttal-dzdtald-dzdttim-dzdttktab1-dzdttktab2;
dcdtBPG=dzdtpgk-dzdtgapdh;
dcdtPG3=dzdtpgm-dzdtpgk;
dcdtPG2=dzdteno+dzdtsersynth-dzdtpgm;
dcdtPEP=dzdtppc+dzdtarog+dzdtpk+dzdtsynth1-dzdteno;
dcdtPYR=dzdtsynth2+dzdtpdh-dzdtpk;
dcdtOAA=dzdtoaaout-dzdtppc;
dcdtPGL=dzdtgnd-dzdtg6pdh;
dcdtRU5P=dzdtrpe+dzdtrpi-dzdtgnd;
dcdtX5P=dzdttktab1+dzdttktab2-dzdtrpe;
dcdtR5P=dzdtrpkk+dzdttktab1-dzdtrpi;
dcdtS7P=dzdttal-dzdttktab1;
dcdtE4P=dzdttktab2+dzdtarog-dzdttal;
dcdtDAHP=dzdtdahpout-dzdtarog;
dcdtATP=dzdtpfk-dzdtpgk-dzdtpk;
dzdtADP=-dzdtpfk+dzdtpgk+dzdtpk;
dcdtNAD=dzdtgapdh;
dcdtNADH=-dzdtgapdh;
dcdtNADP=dzdtg6pdh+dzdtgnd;
dcdtNADPH=-dzdtg6pdh-dzdtgnd;

dcdt=[dcdtG6P;
dcdtF6P;
dcdtFBP;
dcdtDHAP;
dcdtG3P;
dcdtBPG;
dcdtPG3;
dcdtPG2;
dcdtPEP;
dcdtPYR;
dcdtOAA;
dcdtPGL;
dcdtRU5P;
dcdtX5P;
dcdtR5P;
dcdtS7P;
dcdtE4P;
dcdtDAHP;
dcdtATP;
dzdtADP;
dcdtNAD;
dcdtNADH;
dcdtNADP;
dcdtNADPH];

