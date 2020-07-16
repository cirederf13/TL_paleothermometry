
    function f = single_peak(pars,T)
        k=8.617*10^-5;
        Im = pars(1);  
        Tm = pars(2);  
        E = pars(3);  
        b = pars(4);  
        f = Im*b^(b/(b-1))*exp(E/k./T.*(T-Tm)/Tm)....
        .*(1+(b-1)*2*k*Tm/E+(b-1)*(1-2*k*T/E).*(T.^2/Tm^2.*exp(E/k./T.*(T-Tm)/Tm))).^(-b/(b-1));
    end