function p = power_control(xuser, yuser, xAntenna, yAntenna, v, desv, theta)
    dref = 1;
    Lp_ref = 1;
    d = sqrt((xuser - xAntenna)^2 + (yuser - yAntenna)^2); %distance user-antenna
    Lp_lin = Lp_ref * (dref/d)^v; % pathloss in lineal
    Lp_dB = 10*log10(Lp_lin); % pathloss in dB
    shdw = normrnd(0,desv); % in dB
    interference = Lp_dB + shdw; % in dB
    p = 1 / ((10^(interference/10))^theta);
end

