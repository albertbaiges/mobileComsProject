function gain = gain(xuser, yuser, v, desv, Lp_ref, dref)
    d = sqrt((xuser - 0)^2 + (yuser - 0)^2); % distance user - central Antenna
	Lp_lin = Lp_ref * (dref/d)^v; % pathloss in lineal
	Lp_dB = 10*log10(Lp_lin); % pathloss in dB to be able to add the shadowing
	shdw = normrnd(0,desv); % shadowing, normal distribution in dB
	interference = Lp_dB + shdw; % interference, pathloss + dhadowing in dB
	gain = 10^(interference/10); % go back to lineal and then return this value
end

