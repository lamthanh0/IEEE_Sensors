function OP_ORS_EXACT=ORS_INID_EXACT(PdB,xR,ORR,Eta,AP,Cth,PL)
OP_ORS_EXACT   = zeros(1,length(PdB));
for aa = 1 : length(PdB)
    OP_ORS_EXACT(aa) = ham(PdB(aa),xR,ORR,Eta,AP,Cth,PL);
end
%
OP_ORS_EXACT;
%semilogy(PdB,OP_ORS_DF,'r-'); grid on;hold on;
end
%
function OP = ham(PdB,xR,ORR,Eta,AP,Cth,PL)
%
PP      = 10^(PdB/10);
MM      = length(xR);
kap     = Eta*AP/(1-AP);
LSR     = xR.^PL;
LRD     = (1-xR).^PL;
gth     = 2^(Cth/(1-AP))-1;
%
OP      = 1;
for aa = 1 : MM     
    %ff1 = @(xx) ORR*exp(-ORR*xx).*exp(-LRD(aa)*(1-kap*xx*(1 + gth))/kap).*exp(-LSR(aa)*(gth*(1-kap*xx))./(PP*(1-kap*(1+gth)*xx)));    
    %OP1 = integral(ff1,0,1/kap/(1+gth));    
    %ff2 = @(xx,yy)ORR*exp(-ORR*xx)*LRD(aa).*exp(-LRD(aa)*yy).*exp(-LSR(aa)*gth*(1-kap*xx)/kap/PP./yy);
    %ymax =@(xx)(1-kap*(1+gth)*xx)/kap; 
    %OP2 = integral2(ff2,0,1/kap/(1+gth),0,ymax);               
    %OP  = OP*(1 - OP1 - OP2);
    ff1  = @(xx)2.*ORR*exp(-ORR*xx).*sqrt(LSR(aa).*LRD(aa).*(gth-kap.*gth.*xx)./kap./PP).*besselk(1,2.*sqrt(LSR(aa).*LRD(aa).*(gth-kap.*gth.*xx)./kap./PP));
    OP1  = integral(ff1,0,1./kap./(1+gth));
    OP   = OP*(1-OP1);
end 
end






