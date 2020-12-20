function  powerControl(v)
%close all;
hold on;

type = 1/3;
MAX = 1000;
draw_cells = 0; % to draw the cells
bandwidth = 100*10^3;
SIRGap = 4;

SIR_dB = zeros(1, MAX);
throughputs = zeros(1, MAX);

doPowerControl = true;


for theta = 0.05
    for loop = 1:MAX
        if(draw_cells == 1); figure(1); end
        radius = 3;
        %freq = 2*10^9; % Carrier Frequency
        %c = 3*10^8; % Speed of light
        t = 0:pi/3:2*pi;
        n = 1*10^0; % numero de puntos por sector
        ap = radius * sqrt(3)/2; % apotema
        desv = 8; % in dB
        %v = 3.8; % Pathloss exponent
        
        
        
        % Tier 1 Ring
        xt1 = ones(1,6) * (radius + radius/2) .* [0 1 1 0 -1 -1];
        yt1 = radius * ones(1,6) * sqrt(3)/2 .* [2 1 -1 -2 -1 1];
        %Tier 2 Ring
        xt2 = ones(1, 12) * (radius + radius/2) .* [0 1 2 2 2 1 0 -1 -2 -2 -2 -1];
        yt2 = radius * ones(1, 12) * sqrt(3)/2 .* [4 3 2 0 -2 -3 -4 -3 -2 0 2 3];
        
        % Central Cell
        x = radius*cos(t);
        y = radius*sin(t);
        if(draw_cells == 1)
            if(draw_cells == 1); plot(x, y, '-k'); end %hexagono central
            if(draw_cells == 1); plot(0, 0, "ok", 'MarkerSize', 3); end %cell central point
        end
        % Central Random User per sector
        for i=0:2
            sector = (2*pi)/3 * rand + i*(2*pi)/3;
            r = ap*sqrt(rand(n,1));
            xuser = 0 + r.*cos(sector);
            yuser = 0 + r.*sin(sector);
            if(draw_cells == 1); plot(xuser, yuser, "*b", 'MarkerSize', 3); end
            if i == 1
                xuser_ref = xuser;
                yuser_ref = yuser;
            end
        end
        
        % Pathloss of Central Cell
        d = sqrt((xuser_ref - 0)^2 + (yuser_ref - 0)^2); %distance user-antenna
        dref = 1;
        Lp_ref = 1;
        Lp_lin = Lp_ref * (dref/d)^v; % pathloss in lineal
        Lp_dB = 10*log10(Lp_lin); % pathloss in dB
        shdw = normrnd(0,desv); % in dB
        interference = Lp_dB + shdw; % in dB
        p = 1;
        if(doPowerControl)
            p = 1 / ((10^(interference/10))^theta);
        end
        tier_center_int = p*10^(interference/10); % in lineal
        %POWER = 1/(sqrt(tier_center_int));
        %tier_center_dB = 10*log10(tier_center_int); % in dB
        
        % Lineas discontinuas from the Center cell
        for j=1:2:6
            if(draw_cells == 1); plot([0 x(j)], [0 y(j)], '--k'); end
            if (j == 3)
                if(draw_cells == 1); plot([0 x(j)*5], [0 y(j)*5], '--c'); end
            elseif (j == 5)
                if(draw_cells == 1); plot([0 x(j)*5], [0 y(j)*5], '--c'); end
            end
        end
        
        % First Ring
        tier1r_Int = 0; % in lineal
        tier1r_Int_t13 = 0;
        for i=1:6
            if(draw_cells == 1)
                plot(x + xt1(i), y + yt1(i),'-k');
                plot(xt1(i), yt1(i), "ok", 'MarkerSize', 3);
            end
            for j=1:2:6
                if(draw_cells == 1); plot([xt1(i) xt1(i)+x(j)], [yt1(i) yt1(i)+y(j)], '--k'); end
            end
            % Tier 1 Random users per sector
            for k=0:2
                sector = (2*pi)/3 * rand + k*(2*pi)/3;
                r = ap*sqrt(rand(n,1));
                xuser = xt1(i) + r.*cos(sector);
                yuser = yt1(i) + r.*sin(sector);
                if (type == 1/3 && (i == 5 || i == 6) && k == 1)
                    d = sqrt((xuser - 0)^2 + (yuser - 0)^2); %distance user-antenna
                    Lp_lin = Lp_ref * (dref/d)^v; % pathloss in lineal
                    Lp_dB = 10*log10(Lp_lin); % pathloss in dB
                    shdw = normrnd(0,desv); % in dB
                    interference = Lp_dB + shdw; % in dB
                    p = 1;
                    if(doPowerControl)
                        p = power_control(xuser, yuser, xt1(i), yt1(i), v, desv, theta);
                    end
                    tier1r_Int = tier1r_Int + p*10^(interference/10); % in lineal
                    if(draw_cells == 1); plot(xuser, yuser, "*k", 'MarkerSize', 3); end
                    continue
                end
                if(draw_cells == 1); plot(xuser, yuser, "*r", 'MarkerSize', 3); end
            end        
        end
        
        
        tier2r_Int = 0; % in lineal
        for i = 1:12
            if(draw_cells == 1); plot(x + xt2(i), y+yt2(i), '-k'); end
            if(draw_cells == 1); plot(xt2(i), yt2(i), "ok", 'MarkerSize', 3); end
            for j=1:2:6
                if(draw_cells == 1); plot([xt2(i) xt2(i)+x(j)], [yt2(i) yt2(i)+y(j)], '--k'); end
            end
            % Tier 2 Random Users per sector
            for k=0:2
                sector = (2*pi)/3 * rand + k*(2*pi)/3;
                r = ap*sqrt(rand(n,1));
                xuser = xt2(i) + r.*cos(sector);
                yuser = yt2(i) + r.*sin(sector);
                if type == 1/3 && (i == 8 || i == 9 || i == 10 || i == 11 || i == 12) && k == 1
                    d = sqrt((xuser - 0)^2 + (yuser - 0)^2); %distance user-antenna
                    Lp_lin = Lp_ref * (dref/d)^v; % pathloss in lineal
                    Lp_dB = 10*log10(Lp_lin); % pathloss in dB
                    shdw = normrnd(0,desv); % in dB
                    interference = Lp_dB + shdw; % in dB
                    p = 1;
                    if(doPowerControl)
                        p = power_control(xuser, yuser, xt2(i), yt2(i), v, desv, theta);
                    end
                    tier2r_Int = tier2r_Int + p*10^(interference/10); % in lineal
                    %tier2r_Int_dB = 10*log10(tier2r_Int); % in dB
                    if(draw_cells == 1); plot(xuser, yuser, "*k", 'MarkerSize', 3); end
                    continue
                end
                if(draw_cells == 1); plot(xuser, yuser, "*g", 'MarkerSize', 3); end
            end
        end
        SIR = tier_center_int/(tier1r_Int+tier2r_Int); % in lineal + tier2r_Int
        SIR_dB(loop) = 10*log10(SIR); % in dB
        throughputs(loop) = throughput(bandwidth, SIR, SIRGap);
    end


    cdfplot(SIR_dB);
    hold on;
end
figure(1);
xlabel('SIR (dB)');
ylabel('Probability');
title('CDF');
figure(2);
xlabel('SIR (dB)');
ylabel("Throughput")
%legend('type 1/3 power control');
xline(-5);
yline(0.03);
end

