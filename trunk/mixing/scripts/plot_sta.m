U  = um/nacum;
V  = vm/nacum;
W  = wm/nacum;
uu = wk(:,1)/nacum;
vv = wk(:,2)/nacum;
ww = wk(:,3)/nacum;
uv = wk(:,7)/nacum;
%
urms = sqrt(uu-U.^2);
vrms = sqrt(vv-V.^2);
wrms = sqrt(ww-W.^2);
%
uv = uv - U.*V;
%
plot(y,urms,y,vrms,y,wrms)