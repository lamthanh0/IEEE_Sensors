close all; clear all;  clc;
tic
PdB            = -15:2:5;
xR             = [0.25 0.85 0.6]; 
ORR            = 2;
Eta            = 0.8;
AP             = 0.2;
Cth            = 0.5;
PL             = 3;
%bit_frame      = 10^6;
% C-PRS
%OP_CPRS_SIM1=CPRS_INID_SIM(PdB,xR,ORR,Eta,AP,Cth,PL,bit_frame);
OP_CPRS_EXACT1=CPRS_INID_EXACT(PdB,xR,ORR,Eta,AP,Cth,PL);
% S-PRS
%OP_SPRS_SIM1=SPRS_INID_SIM(PdB,xR,ORR,Eta,AP,Cth,PL,bit_frame);
OP_SPRS_EXACT1=SPRS_INID_EXACT(PdB,xR,ORR,Eta,AP,Cth,PL);
% C-FRS
%CFRS_INID_SIM(PdB,xR,ORR,Eta,AP,Cth,PL,bit_frame);
% ORS
%OP_ORS_SIM1=ORS_INID_SIM(PdB,xR,ORR,Eta,AP,Cth,PL,bit_frame);
OP_ORS_EXACT1=ORS_INID_EXACT(PdB,xR,ORR,Eta,AP,Cth,PL);
toc
figure(1)
h1 = semilogy(PdB,OP_CPRS_SIM1,'r--','LineWidth',2); hold on;
h2 = semilogy(PdB,OP_CPRS_EXACT1,'rs','LineWidth',2); hold on;
%h2 = semilogy(PdB,OP_PRS_AF_IND_EXACT_1,'r*','LineWidth',2); hold on;
h3 = semilogy(PdB,OP_SPRS_SIM1,'r--','LineWidth',2); hold on;
h4 = semilogy(PdB,OP_SPRS_EXACT1,'b*','LineWidth',2); hold on;
%h4 = semilogy(PdB,OP_PRS_AF_IND_EXACT_2,'ro','LineWidth',2); hold on;
h5 = semilogy(PdB,OP_ORS_SIM1,'r--','LineWidth',2); hold on;
h6 = semilogy(PdB,OP_ORS_EXACT1,'mo','LineWidth',2); hold on;
%h6 = semilogy(PdB,OP_PRS_AF_IND_EXACT_3,'rs','LineWidth',2); hold on;
grid on;
title(['OP versus transmit SNR with \rho=0.2, M=3 and \eta=0.8']);
xlabel('Transmit SNR \Psi (dB)');
ylabel('Outage Probability (OP)');
legend([h1,h2,h4,h6],'Monte Carlo Simulation','CPRS-Theo','SPRS-Theo','ORS-Theo');
