function p = power_control(xuser, yuser, xAntenna, yAntenna, v, desv, theta)
    dref = 1; % Reference distance
    Lp_ref = 1; % Pathloss on reference distance
    d = sqrt((xuser - xAntenna)^2 + (yuser - yAntenna)^2); %distance user - antenna (expected its antenna)
    Lp_lin = Lp_ref * (dref/d)^v; % compute pathloss in lineal
    Lp_dB = 10*log10(Lp_lin); % pathloss in dB to be able to add shadowing
    shdw = normrnd(0,desv); % shadowing in dB is normal distributed
    interference = Lp_dB + shdw; % interference in dB, pathloss + shadowing
    p = 1 / ((10^(interference/10))^theta); % Compute the power control factor
end

