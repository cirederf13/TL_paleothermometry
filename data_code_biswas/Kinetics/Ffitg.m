%%Fits fading data to derive rho' after Huntley (2006)%%

function out1 = Ffitg(beta0,t)
I=beta0(1); m=beta0(2);
fit = m*log10(t)+I;
out1=fit;

