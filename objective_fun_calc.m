%**************************************************************************
function [KGE,LOA] = objective_fun_calc(oQ,pQ,qv,minQV,maxQV)
e=eps;
%--------------------------------------------------------------------------
%kling_gupta efficiency
KGE                = 1 - sqrt(  (corr(pQ,oQ)-1).^2 ...
                              + (std(pQ)./(std(oQ)+e)-1).^2 ...
                              + (mean(pQ)./(mean(oQ)+e)-1).^2   ) ;
%--------------------------------------------------------------------------
LOA                = 1;

% dQ                 = max(abs([0;diff(oQ)]));
% UB                 = oQ + dQ; UB(UB>1)=1;
% LB                 = oQ - dQ; LB(LB<0)=0;
% 
% cond               = pQ>UB | pQ<LB;
% if nnz(cond)>0.1*length(oQ) || qv<minQV || qv>maxQV
%     LOA            = 0;
% end


if qv<minQV || qv>maxQV
    LOA            = 0;
end