function OP_CPRS_SIM=CPRS_INID_SIM(PdB,xR,ORR,Eta,AP,Cth,PL,bit_frame)
%
OP_CPRS_SIM   = zeros(1,length(PdB));
%
for aa = 1 : length(PdB)
    fprintf('Running %d per %d \n',aa,length(PdB));
    PP      = 10^(PdB(aa)/10);
    MM      = length(xR);
    kap     = Eta*AP/(1-AP);
    LSR     = xR.^PL;
    LRD     = (1-xR).^PL;
    gth     = 2^(Cth/(1-AP))-1;
    %
    for bitnum   =  1 : bit_frame
        %
        SNR_SR = zeros(1,MM);
        SNR_RD = zeros(1,MM);
        SNR_RR = zeros(1,MM);
        for bb = 1 : MM
            hSR        = sqrt(1/2/LSR(bb))*(randn(1,1) + 1i*randn(1,1));
            SNR_SR(bb) = abs(hSR)^2;
            hRD        = sqrt(1/2/LRD(bb))*(randn(1,1) + 1i*randn(1,1));
            SNR_RD(bb) = abs(hRD)^2;
            hRR        = sqrt(1/2/ORR)*(randn(1,1) + 1i*randn(1,1));
            SNR_RR(bb) = abs(hRR)^2;
        end
        SNR_SR_max = max(SNR_SR);
        ID         = find(SNR_SR == SNR_SR_max);
        XX         = SNR_RD(ID);
        YY         = SNR_RR(ID);
        SNR1       = (1-kap.*YY)./kap./YY;
        SNR2       = kap*PP*SNR_SR_max*XX/(1-kap*YY);
        SNR_DF     = min(SNR1,SNR2);
        if (SNR_DF < gth)
            OP_CPRS_SIM(aa) = OP_CPRS_SIM(aa) + 1;           
        end                                                
    end
end
%
OP_CPRS_SIM= OP_CPRS_SIM/bit_frame;
OP_CPRS_SIM;
%semilogy(PdB,OP_CPRS,'bo'); grid on;hold on;
end





