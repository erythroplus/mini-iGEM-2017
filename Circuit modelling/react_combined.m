% Version contaning all the differential equations
function dydt = react_combined(t,y)


dydt = zeros(size(y));
% Progesterone sensing block parameters
k1 = 0.3465;
c_cmv = 0.6; % 0.6 volt
km1 = 8.0; % 0.1 volt
n = 2.0;
k2 = 7.921;
k3 = 1.0;
k_3 = 10.0;
d1 = 0.01980;
d2 = 0.01980;
d3 = 0.01980;
k4 = 0.3465;
km2 = 0.1;
n2 = 2.0;
d5 = 0.1386;
k5 = 7.921;
d6 = 0.01980;

%Estrogene sensing block parameterse
e_k1 = 0.3465; 
e_c_cmv = 10.0; % 0.6 volt
e_km1 = 10.0; % 0.1 volt
e_n = 2.0;
e_k2 = 7.921;
e_k3 = 1.0;
e_k_3 = 10.0;
e_d1 = 0.01980;
e_d2 = 0.01980;
e_d3 = 0.01980;
e_k4 = 0.3465;
e_km2 = 0.1;
e_n2 = 2.0;
e_d5 = 0.1386;
e_k5 = 7.921;
e_d6 = 0.01980;
e_d7 = 0.0; % bit random, active deg 15 was

% Ratiometric block parameters
% s_k1 s_km1 s_n1
s_k1 = 0.3465; %0.3465
s_n1 = 2.0; % 2
s_km1 = 0.0005; % 0.1 %0.001 szupersÈges, de az elejÈn lÈvı szar az szar 
s_d1 = 0.3465; %mrna
s_k2 = 7.921;
s_d2 = 0.01980; % protein
s_k3 = 0.3465; % ezt megturbÛzzuk
s_km2 = 0.1;
s_n2 = 2.0;
s_d3 = 0.3465;
s_k4 = 7.921;
s_d4 = 0.01980;
s_k5 = 0.3465; % 0.3465
s_km3 = 0.1; % 0.1
s_n3 = 2.0; 
s_d5 = 0.3465;
s_k6 = 7.921; % 7.921
s_d6 = 0.3465; 

s_c = 0.01;

% Species list
mRNA_LEX = y(1);  
LEX_PR = y(2); 
LEX_PPR = y(3);
mRNA_TEV = y(4);
TEV = y(5);
pro = y(6);
mRNA_GAL = y(7);  
GAL_ER = y(8); 
GAL_EER = y(9);
mRNA_TALAV = y(10);
TALAV = y(11);
est = y(12);
mRNA_TALD = y(13);
TALD = y(14);
mRNA_TALB = y(15);
TALB = y(16);
mRNA_RDF = y(17);
RDF = y(18);

% Differential equations
dydt(1) = (k1*(LEX_PPR + c_cmv)^n)/(km1^n + (LEX_PPR + c_cmv)^n) - d1*mRNA_LEX - k2*mRNA_LEX; % dMRNA_LEX/dt
dydt(2) = k2*mRNA_LEX - k3*LEX_PR*pro + k_3*LEX_PPR - d2*LEX_PR; % d_LEX_PR/dt
dydt(3) = k3*LEX_PR*pro - k_3*LEX_PPR - d3*LEX_PPR; % d_LEX_PPR/dt
dydt(4) = (k4*LEX_PPR^n2)/(km2^n2 + LEX_PPR^n2) - d5*mRNA_TEV - k5*mRNA_TEV; % d_MRNA_TEV/dt
dydt(5) = k5*mRNA_TEV - d6*TEV; % dTEV/dt
dydt(6) = (2.*7.897.*(2.746 .* 10^4 - t).*exp(-(2.746 .* 10^4 - t).^2/(2697)^2))/(2697)^2 + (2.*39.04.*(3.252 .* 10^4 - t) .*exp(-(3.252.*10^4 - t).^2/6783^2))/6783^2; % PROGESTERONE DERIVATIVE
dydt(7) = (e_k1*(GAL_EER + e_c_cmv)^e_n)/(e_km1^e_n + (GAL_EER + e_c_cmv)^e_n) - e_d1*mRNA_GAL - e_k2*mRNA_GAL; % dMRNA_GAL/dt
dydt(8) = e_k2*mRNA_GAL - e_k3*GAL_ER*est + e_k_3*GAL_EER - e_d2*GAL_ER; % d_GAL_ER/dt
dydt(9) = e_k3*GAL_ER*est - e_k_3*GAL_EER - e_d3*GAL_EER; % d_GAL_EER/dt
dydt(10) = (e_k4*GAL_EER^e_n2)/(e_km2^e_n2 + GAL_EER^e_n2) - e_d5*mRNA_TALAV - e_k5*mRNA_TALAV; % d_MRNA_TALAVP16/dt
dydt(11) = e_k5*mRNA_TALAV - e_d6*TALAV - e_d7*TALAV*TEV; % dTALAVP16/dt
dydt(12) = (2.*0.7007.*(1.929.*10^4 - t) .*exp(-(1.929*10^4  - t).^2./4768^.2))/4768.^2 + (2 .*0.5464.*(1.929.*10^4 - t) .*exp(-(1.929.*10^4 - t).^2./8017.^2))./8017.^2;
%Talav active degradation added
dydt(13) = (s_k1*TALAV^s_n1)/(s_km1^s_n1 + TALAV^s_n1) - s_d1*mRNA_TALD - s_k2*mRNA_TALD; % dTALD_MRNA/dt
dydt(14) =  s_k2*mRNA_TALD - s_d2 * TALD;% dTALD/dt
dydt(15) = (s_k3*s_km2^s_n2)/(s_km2^n + TALD^n) - s_d3*mRNA_TALB - s_k4*mRNA_TALB; % dTALB_MRNA/dt
dydt(16) =  s_k4*mRNA_TALB - s_d4*TALB; % dTALB/dt
dydt(17) = (s_k5*TALB^s_n3)/(s_km3+TALB^s_n3)- s_d5*mRNA_RDF - s_k6*mRNA_RDF;% dRDF_MRNA/dt
dydt(18) = s_k6*mRNA_RDF - s_d6*RDF;% dRDF/dt
