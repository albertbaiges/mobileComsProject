function gain = gain(xuser, yuser, v, desv, Lp_ref, dref)
    d = sqrt((xuser - 0)^2 + (yuser - 0)^2); %distance user-antenna
	Lp_lin = Lp_ref * (dref/d)^v; % pathloss in lineal
	Lp_dB = 10*log10(Lp_lin); % pathloss in dB
	shdw = normrnd(0,desv); % in dB
	interference = Lp_dB + shdw; % in dB
	gain = 10^(interference/10); % in lineal
end

