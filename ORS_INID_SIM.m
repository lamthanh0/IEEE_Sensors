function OP_ORS_SIM=ORS_INID_SIM(PdB,xR,ORR,Eta,AP,Cth,PL,bit_frame)
%
OP_ORS_SIM    = zeros(1,length(PdB));
%
for aa = 1 : length(PdB)
    fprintf('Running %d per %d \n',aa,length(PdB));
    PP      = 10^(PdB(aa)/10);
    MM      = length(xR);
    kap     = Eta*AP/(1-AP);
    LSR     = xR.^PL;
    LRD     = (1-xR).^PL;
    gth     = 2^(Cth/(1-AP)) - 1;
    %
    for bitnum   =  1 : bit_frame
        %
        SNR_SRD_max = zeros(1,MM);
        for bb = 1 : MM
            hSR        = sqrt(1/2/LSR(bb))*(randn(1,1) + 1i*randn(1,1));
            SNR_SR     = abs(hSR)^2;
            hRD        = sqrt(1/2/LRD(bb))*(randn(1,1) + 1i*randn(1,1));
            SNR_RD     = abs(hRD)^2;
            hRR        = sqrt(1/2/ORR)*(randn(1,1) + 1i*randn(1,1));
            SNR_RR     = abs(hRR)^2;
            XX         = SNR_RD;
            YY         = SNR_RR;
            SNR1       = (1-kap.*YY)./kap./YY;
            SNR2       = kap*PP*SNR_SR*XX/(1-kap*YY);
            SNR_DF     = min(SNR1,SNR2);
            if (SNR_DF > SNR_SRD_max)
                SNR_SRD_max = SNR_DF;
            end            
        end       
        if (SNR_SRD_max < gth)
            OP_ORS_SIM(aa) = OP_ORS_SIM(aa) + 1;
        end
    end
end
%
OP_ORS_SIM = OP_ORS_SIM/bit_frame;
OP_ORS_SIM;
%semilogy(PdB,OP_ORS,'rs'); grid on;hold on;
end





