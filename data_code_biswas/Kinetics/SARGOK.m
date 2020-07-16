%%Fits SAR data to a single saturating exponential%%
%%Solves for all aliquots at once, Georgina King, August 2015%%
%%georgina.king@uni-kolen.de%%

function out = SARGOK(beta0, D)
global a
I=beta0(1); D0=beta0(2);
Fgrowth = I.*(1-(1+(D./D0).*(a-1)).^(-1/(a-1)));
out=[Fgrowth];

