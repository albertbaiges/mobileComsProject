function plotSNR(type, v, doPowerControl, doThroughput)
global legendString throughputString
if isempty(legendString)
    legendString = "";
end
if isempty(throughputString)
    throughputString = "";
end
if(type ~= 1)
    frac = rat(type);
    fracName = frac(4:end);
else
    fracName = 1;
end
str = "type " + fracName + " with v = " + v;
if(doPowerControl == true)
    str = str + " with Power Control";
end
legendString(end+1) = str;
MAX = 1000;
draw_cells = 0; % to draw the cells
Lp_ref = 1;
dref = 1;
SIR_dB = zeros(1, MAX);
throughputs = zeros(1, MAX);
theta = 0.5;
bandwidth = 100*10^6;
SIRGap = 4; %dB
for loop = 1:MAX
    if(draw_cells == 1); figure(1); end
    radius = 3;
    %freq = 2*10^9; % Carrier Frequency
    %c = 3*10^8; % Speed of light
    t = 0:pi/3:2*pi;
    n = 1*10^0; % numero de puntos por  sector
    ap = radius * sqrt(3)/2; % apotema
    desv = 8; % in dB
    

    
    
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
    p = 1;
	if(doPowerControl)
        p = power_control(yuser_ref, yuser_ref, 0, 0, v, desv, theta);
	end
    desired_signal = p*gain(xuser_ref, yuser_ref, v, desv, Lp_ref, dref); % in lineal
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
            
            if type == 1 && (i == 5 || i == 6)
                tier1r_Int = tier1r_Int + gain(xuser, yuser, v, desv, Lp_ref, dref);
                if(draw_cells == 1); plot(xuser, yuser, "*k", 'MarkerSize', 3); end
                continue
            elseif (type == 1/3 && (i == 5 || i == 6) && k == 1)
                p = 1;
                if(doPowerControl)
                    p = power_control(xuser, yuser, xt1(i), yt1(i), v, desv, theta);
                end
                tier1r_Int = tier1r_Int + p*gain(xuser, yuser, v, desv, Lp_ref, dref); % in lineal
                if(draw_cells == 1); plot(xuser, yuser, "*k", 'MarkerSize', 3); end
                continue
            elseif (type == 1/9)
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
            
            if type == 1 && (i == 8 || i == 9 || i == 10 || i == 11 || i == 12)
                if(i == 8 && k == 0 && sector >= pi/3)
                    tier2r_Int = tier2r_Int + gain(xuser, yuser, v, desv, Lp_ref, dref); % in lineal
                    %tier2r_Int_dB = 10*log10(tier2r_Int); % in dB
                    if(draw_cells == 1); plot(xuser, yuser, "*k", 'MarkerSize', 3); end
                    continue
                elseif(i == 8 && k == 1) || (i == 9 || i == 10 || i == 11)
                    tier2r_Int = tier2r_Int + gain(xuser, yuser, v, desv, Lp_ref, dref); % in lineal
                    %tier2r_Int_dB = 10*log10(tier2r_Int); % in dB
                    if(draw_cells == 1); plot(xuser, yuser, "*k", 'MarkerSize', 3); end
                    continue
                elseif(i == 12) && (k == 1 || (k == 2 && sector < 5*pi/3 ))
                    tier2r_Int = tier2r_Int + gain(xuser, yuser, v, desv, Lp_ref, dref); % in lineal
                    %tier2r_Int_dB = 10*log10(tier2r_Int); % in dB
                    if(draw_cells == 1); plot(xuser, yuser, "*k", 'MarkerSize', 3); end
                    continue
                end
                
                if(draw_cells == 1); plot(xuser, yuser, "*g", 'MarkerSize', 3); end
                
            elseif type == 1/3 && (i == 8 || i == 9 || i == 10 || i == 11 || i == 12) && k == 1
                p = 1;
                if(doPowerControl)
                    p = power_control(xuser, yuser, xt2(i), yt2(i), v, desv, theta);
                end
                tier2r_Int = tier2r_Int + p*gain(xuser, yuser, v, desv, Lp_ref, dref); % in lineal
                %tier2r_Int_dB = 10*log10(tier2r_Int); % in dB
                if(draw_cells == 1); plot(xuser, yuser, "*k", 'MarkerSize', 3); end
                continue
                
            elseif type == 1/9 && (i == 8 || i == 10 || i == 12) && k == 1
                tier2r_Int = tier2r_Int + gain(xuser, yuser, v, desv, Lp_ref, dref); % in lineal
                %tier2r_Int_dB = 10*log10(tier2r_Int); % in dB
                if(draw_cells == 1); plot(xuser, yuser, "*k", 'MarkerSize', 3); end
                continue
            end
            if(draw_cells == 1); plot(xuser, yuser, "*g", 'MarkerSize', 3); end
        end
    end
    SIR = desired_signal/(tier1r_Int+tier2r_Int); % in lineal + tier2r_Int
    SIR_dB(loop) = 10*log10(SIR); % in dB
    %disp(SIR_dB);
    if(doThroughput)
        throughputs(loop) = throughput(bandwidth, SIR, SIRGap);
    end
end
figure(1);
cdfplot(SIR_dB);
hold on;
xlabel('SIR (dB)');
ylabel('Probability');
title('CDF');
legend(legendString(2:end));
if(doThroughput)
    figure(2);
    cdfplot(throughputs)
    xlabel("Throughput");
    ylabel("Probability");
    throughputString(end+1) = "type " + fracName;
    legend(throughputString(2:end));
    hold on;
end
end

