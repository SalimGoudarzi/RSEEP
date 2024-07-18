%**************************************************************************
function [W,nb,idb] = LikelihoodWeightCal(KGE,LOA)

cond             = isnan(KGE) | LOA==0 | KGE<-100;

% target         = 0.99;
% cond           = KGE<target | isnan(KGE) | LOA==0;
% 
% Num            = floor(0.01*length(KGE));
% while nnz(1-cond)<Num && target>-10
%     target     = target-0.005;
%     cond       = KGE<target | isnan(KGE)| LOA==0;
% end

KGE(cond)        = min(KGE(:));
%--------------------------------------------------------------------------
%                   evaluate weights for all simulations
%--------------------------------------------------------------------------
% W               = zeros(length(KGE),1);
nb               = floor(0.01*length(KGE));
[~,id]           = sort(KGE,'descend');   
KGE(id(nb+1:end))= min(KGE(:));
idb              = id(1:nb);
% nb              = length(KGE)-nnz(cond==1); 
%calculate performacen metric weights
W                = ofs_weight(KGE);

%**************************************************************************
function [w,ofs] = ofs_weight(ofs)
worst            = min(ofs);
ofs(isnan(ofs))  = worst;
w                = (ofs-worst).^2;
w                = w./sum(w);