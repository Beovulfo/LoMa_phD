flag    = true;%splitting for chebfun
Q = struct('maxLength',2^18+1);
chebtech.techPref(Q);
cheboppref.setDefaults('display','iter')
cheboppref.setDefaults('errTol',5e-7)

chebfunpref.setDefaults('splitting',flag,'maxLength',40000)
chebfunpref.setDefaults('eps',1e-13)
%chebfunpref.setDefaults('splitMaxLength',10000)
%cheboppref.setDefaults('plotting','on')
cheboppref.setDefaults('vectorize',true)
