function OP_CPRS_EXACT=CPRS_INID_EXACT(PdB,xR,ORR,Eta,AP,Cth,PL)
OP_CPRS_EXACT   = zeros(1,length(PdB));
for aa = 1 : length(PdB)
    OP_CPRS_EXACT(aa) = ham(PdB(aa),xR,ORR,Eta,AP,Cth,PL);
end
%
OP_CPRS_EXACT;
%semilogy(PdB,OP_PRS,'b-'); grid on;hold on;
end
%
function OP = ham(PdB,xR,ORR,Eta,AP,Cth,PL)
%
PP      = 10^(PdB/10);
MM      = length(xR);
kap     = Eta*AP/(1-AP);
LSR     = xR.^PL;
LRD     = (1-xR).^PL;
gth     = 2^(Cth/(1-AP)) - 1;
%
OP      = 1;
for aa = 1 : MM
    Setaa = setdiff(1:MM,aa);
    for bb = 1 : MM - 1
        SET1 = nchoosek(Setaa,bb);
        [row col] = size(SET1);
        for cc = 1 : row
            SET2 = SET1(cc,:);
            LSR_sum = 0;
            for dd = 1 : col
                LSR_sum = LSR_sum + LSR(SET2(dd));
            end
            %               
            %ff1    = @(xx) ORR*exp(-ORR*xx).*exp(-LRD(aa)*(1-kap*xx*(1 + gth))/kap).*(exp(-LSR(aa)*(gth*(1-kap*xx))./(PP*(1-kap*(1+gth)*xx)))-LSR(aa)/(LSR(aa) + LSR_sum)*exp(-(LSR(aa) + LSR_sum)*(gth*(1-kap*xx))./(PP*(1-kap*(1+gth)*xx))));                                              
            %OP1    = integral(ff1,0,1/kap/(1+gth));   
            ff1    = @(xx)ORR.*exp(-ORR.*xx).*sqrt(4.*LSR(aa).*LRD(aa).*(gth-kap.*gth.*xx)./kap./PP).*besselk(1,sqrt(4.*LSR(aa).*LRD(aa).*(gth-kap.*gth.*xx)./kap./PP));
            OP1    = integral(ff1,0,1./kap./(1+gth));
            %ff2    = @(xx,yy,zz) ORR*exp(-ORR*xx)*LRD(aa).*exp(-LRD(aa)*yy).*(exp(-LSR(aa)*gth*(1-kap*xx)/kap/PP./yy) - LSR(aa)/(LSR(aa) + LSR_sum)*exp(-(LSR(aa) + LSR_sum)*gth*(1-kap*xx)/kap/PP./yy));
            %ymax   = @(xx)(1-kap*(1+gth)*xx)/kap;                        
            %OP2    = integral2(ff2,0,1/kap/(1+gth),0,ymax);  
            ff2    = @(xx) LSR(aa)./(LSR_sum+LSR(aa)).*ORR.*exp(-ORR.*xx).*...
                sqrt(4.*(LSR_sum+LSR(aa)).*LRD(aa).*(gth-kap.*gth.*xx)./kap./PP).*besselk(1,sqrt(4.*(LSR_sum+LSR(aa)).*LRD(aa).*(gth-kap.*gth.*xx)./kap./PP));
            OP2    = integral(ff2,0,1./kap./(1+gth));
            OP     = OP - (-1)^(bb+1)*(OP1 - OP2);
        end
    end
end
end






