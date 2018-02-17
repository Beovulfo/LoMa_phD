function [k,sp] = compuspe(kx,ky,kz,var); 

nn = length(kx(:)); 
k2 = kx(:).^2 + ky(:).^2 + kz(:).^2; 
k = min(sqrt(k2(:))):max(2./3.*sqrt(k2(:))); 
sp = zeros(size(k)); 
ii = interp1(k,1:length(k),sqrt(kx(:).^2 + ky(:).^2 + kz(:).^2),'nearest'); 

kkk = var(:).*conj(var(:))./nn.^2; 

for ik = 1:length(k); 

    lista = ismember(ii,ik); 
    sp(ik) = sum(kkk(lista)); 

end 
