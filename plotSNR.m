function plotSNR(type, v, doPowerControl, doThroughput)
global legendString throughputString % Get the global variablees in which we acumulate the legends
if isempty(legendString) % Check if it has been cleanned or is "first run"
    legendString = ""; % Assign an empty string to later accummulate
end
if isempty(throughputString) % Check if it has been cleanned or its "first run"
    throughputString = ""; % Assign an empty string to later accummulate
end
if(type ~= 1) % Get the reuse factor in fraction to then use it in the legend
    frac = rat(type);
    fracName = frac(4:end);
else
    fracName = 1;
end
str = "type " + fracName + " with v = " + v; % Creating string that does into legend
if(doPowerControl == true)
    str = str + " with Power Control"; % Adding power control if it was specified
end
legendString(end+1) = str; % Append the new string to the legend
MAX = 1000; % Number of iterations in the montecarlo
draw_cells = false; % Enable this option to draw the cells and users. (Not recommended, slows the app)
Lp_ref = 1; % Pathloss on reference distance
dref = 1; % Reference distance
SIR_dB = zeros(1, MAX); % Number of samples we will compute of SIR
throughputs = zeros(1, MAX); % Number of samples we will compute of throughput
theta = 0.5; % Value of theta that will be used for the power control
bandwidth = 100*10^6; % Bandwith used then for the throughput
SIRGap = 4; % SIR Gap in dB
for loop = 1:MAX
    if(draw_cells == 1); figure(1); end % Draw central cell
    radius = 3; % Radius of the cells
    t = 0:pi/3:2*pi; % Angles of vertices hexagon
    n = 1; % Number of points by sector
    ap = radius * sqrt(3)/2; % Apotema
    desv = 8; % standard deviation of shadowing in dB
    

    
    
    % Tier 1 Ring of cells
    xt1 = ones(1,6) * (radius + radius/2) .* [0 1 1 0 -1 -1]; % Coefficients of x coordinates  cells
    yt1 = radius * ones(1,6) * sqrt(3)/2 .* [2 1 -1 -2 -1 1]; % Coefficients of y coordinates  cells
    % Tier 2 Ring of cells
    % Coefficients of x coordinates cells
    xt2 = ones(1, 12) * (radius + radius/2) .* [0 1 2 2 2 1 0 -1 -2 -2 -2 -1];
    %Coefficients of y coordinates cells
    yt2 = radius * ones(1, 12) * sqrt(3)/2 .* [4 3 2 0 -2 -3 -4 -3 -2 0 2 3];
    
    % Central Cell
    x = radius*cos(t); % X coordinates of axis central cell
    y = radius*sin(t); % Y coordinates of axis central cell
    if(draw_cells == 1) % Draw hexagon and sector lines
        if(draw_cells == 1); plot(x, y, '-k'); end %hexagono central
        if(draw_cells == 1); plot(0, 0, "ok", 'MarkerSize', 3); end %cell central point
    end
    % Central Random User per sector
    for i=0:2 % Sectors in cell, used to place users
        sector = (2*pi)/3 * rand + i*(2*pi)/3; % Compute random angle, and desfase based on current sector
        r = ap*sqrt(rand(n,1)); % Random distance from the center
        xuser = 0 + r.*cos(sector); % x coordinate of user computed from polars
        yuser = 0 + r.*sin(sector); % y coordinate of user computed from polars
        if(draw_cells == 1); plot(xuser, yuser, "*b", 'MarkerSize', 3); end % Draw the point if enabled
        if i == 1 % User used for computations
            xuser_ref = xuser; % Store the x coordinate of the left sector central cell user
            yuser_ref = yuser; % Store the y coordinate of the left sector central cell user
        end
    end
    
    % Pathloss of Central Cell
    p = 1; % In case there is not power control, control factor is 1 ("None")
	if(doPowerControl) % If there is power control, compute the value of p
        p = power_control(xuser_ref, yuser_ref, 0, 0, v, desv, theta);
    end
    % Compute the strengh of the signal by finding its gain to the central
    % antenna and multiplying by the power control factor
    desired_signal = p*gain(xuser_ref, yuser_ref, v, desv, Lp_ref, dref); % in lineal
    
    % Lineas discontinuas from the Center cell, only done when drawing is
    % enabled
    for j=1:2:6
        if(draw_cells == 1); plot([0 x(j)], [0 y(j)], '--k'); end
        if (j == 3)
            if(draw_cells == 1); plot([0 x(j)*5], [0 y(j)*5], '--c'); end
        elseif (j == 5)
            if(draw_cells == 1); plot([0 x(j)*5], [0 y(j)*5], '--c'); end
        end
    end
    
    % First Ring, cells "surrounding" the central cell
    tier1r_Int = 0; % Accumulated total interference from cells in first ring in lineal
    for i=1:6 % Number of cells, do computations for each of them
        if(draw_cells == 1) % Draw cells if enabled
            plot(x + xt1(i), y + yt1(i),'-k');
            plot(xt1(i), yt1(i), "ok", 'MarkerSize', 3);
        end
        for j=1:2:6 % Draw sectors if enabled
            if(draw_cells == 1); plot([xt1(i) xt1(i)+x(j)], [yt1(i) yt1(i)+y(j)], '--k'); end
        end
        % Tier 1 Random users per sector
        for k=0:2 % Number of sectors, done for each sector
            sector = (2*pi)/3 * rand + k*(2*pi)/3; % Compute random angle for location user
            r = ap*sqrt(rand(n,1)); % Compute random distance from the center
            xuser = xt1(i) + r.*cos(sector); % X coordinate of the user computed with polars
            yuser = yt1(i) + r.*sin(sector); % Y coordinate of the user computed with polars
            
            % Compute and acummulate the interferences, only ones that
            % interfere with reference user
            if type == 1 && (i == 5 || i == 6) % If universl bandwith, only cells 5 and 6 interfere
                % Compute the interference by finding the gain
                % and accumulating it
                % No support for power control added
                tier1r_Int = tier1r_Int + gain(xuser, yuser, v, desv, Lp_ref, dref);
                if(draw_cells == 1); plot(xuser, yuser, "*k", 'MarkerSize', 3); end
                continue
            elseif (type == 1/3 && (i == 5 || i == 6) && k == 1)
                % Compute the interference by finding the gain
                % multiplying it by the power control factor
                % and accummulating it
                p = 1; % Power control factor in case there is not
                if(doPowerControl) % If power control enabled, find power control factor
                    p = power_control(xuser, yuser, xt1(i), yt1(i), v, desv, theta);
                end
                % Compute the gain and multiply by power control to then
                % accumulate it
                tier1r_Int = tier1r_Int + p*gain(xuser, yuser, v, desv, Lp_ref, dref); % in lineal
                if(draw_cells == 1); plot(xuser, yuser, "*k", 'MarkerSize', 3); end
                continue
            elseif (type == 1/9) % Reuse factor 1/9
                continue % No users from this ring of cells creaters inteferences on 1/9
            end
            % Draw user if enabled
            if(draw_cells == 1); plot(xuser, yuser, "*r", 'MarkerSize', 3); end
        end        
    end
    
    
    % Secong Ring, cells "surrounding" the central cell
    tier2r_Int = 0; % Accumulated total interference from cells in second ring in lineal
    for i = 1:12 % Number of cells in second ring, computations done for each of them
        if(draw_cells == 1); plot(x + xt2(i), y+yt2(i), '-k'); end % Draw cells if enabled
        if(draw_cells == 1); plot(xt2(i), yt2(i), "ok", 'MarkerSize', 3); end % Draw antenna if enabled
        for j=1:2:6 % Number of sectors, line from center every two axis
            if(draw_cells == 1); plot([xt2(i) xt2(i)+x(j)], [yt2(i) yt2(i)+y(j)], '--k'); end
        end
        
        % Tier 2 Random Users per sector and computing of the interference
        % All the gains from these users to the central cell will account
        % as intereference
        for k=0:2 % Number of users per cell, one per sector
            sector = (2*pi)/3 * rand + k*(2*pi)/3; % Find random angle in which user is placed
            r = ap*sqrt(rand(n,1)); % Find random distance from the center in which user is placed
            xuser = xt2(i) + r.*cos(sector); % Find x coordinate of user using polars
            yuser = yt2(i) + r.*sin(sector); % Find y coordinate of user using polars
            
            % Find the interference with respect, reference user so only
            % certain users interact
            
            % Universal bandwith and cells that affect when universal
            % bandwith
            if type == 1 && (i == 8 || i == 9 || i == 10 || i == 11 || i == 12)
                if(i == 8 && k == 0 && sector >= pi/3) % Certain sectors
                    % Compute gain of the user and accummulate it
                    % No option for power control
                    tier2r_Int = tier2r_Int + gain(xuser, yuser, v, desv, Lp_ref, dref); % in lineal
                    %tier2r_Int_dB = 10*log10(tier2r_Int); % in dB
                    % Draw if enabled
                    if(draw_cells == 1); plot(xuser, yuser, "*k", 'MarkerSize', 3); end
                    continue
                elseif(i == 8 && k == 1) || (i == 9 || i == 10 || i == 11) % Certain sectors
                    % Compute gain of the user and accummulate it 
                    % No option for power control
                    tier2r_Int = tier2r_Int + gain(xuser, yuser, v, desv, Lp_ref, dref); % in lineal
                    %tier2r_Int_dB = 10*log10(tier2r_Int); % in dB
                    if(draw_cells == 1); plot(xuser, yuser, "*k", 'MarkerSize', 3); end
                    continue
                elseif(i == 12) && (k == 1 || (k == 2 && sector < 5*pi/3 )) % Certain sector and depending
                    % Compute gain of the user and accumulate it
                    % No option for power control
                    tier2r_Int = tier2r_Int + gain(xuser, yuser, v, desv, Lp_ref, dref); % in lineal
                    %tier2r_Int_dB = 10*log10(tier2r_Int); % in dB
                    if(draw_cells == 1); plot(xuser, yuser, "*k", 'MarkerSize', 3); end
                    continue
                end
                
                if(draw_cells == 1); plot(xuser, yuser, "*g", 'MarkerSize', 3); end
                
            % Reuse factor 1/3 and certain cells with only left sector
            elseif type == 1/3 && (i == 8 || i == 9 || i == 10 || i == 11 || i == 12) && k == 1
                % In case power control is not enabled p power control
                % factor is 1
                p = 1;
                if(doPowerControl) % If power control is enabled, compute power control factor
                    p = power_control(xuser, yuser, xt2(i), yt2(i), v, desv, theta);
                end
                % Compute interference by computing gain and multiplying by
                % factor then accummulating it on the total interferences
                tier2r_Int = tier2r_Int + p*gain(xuser, yuser, v, desv, Lp_ref, dref); % in lineal
                %tier2r_Int_dB = 10*log10(tier2r_Int); % in dB
                % Draw user if enabled
                if(draw_cells == 1); plot(xuser, yuser, "*k", 'MarkerSize', 3); end
                continue
            
            % Reuse factor 1/9 and certain cells with only left sector    
            elseif type == 1/9 && (i == 8 || i == 10 || i == 12) && k == 1
                % Compute interferences from this user and accummulate it
                % No option for power control
                tier2r_Int = tier2r_Int + gain(xuser, yuser, v, desv, Lp_ref, dref); % in lineal
                %tier2r_Int_dB = 10*log10(tier2r_Int); % in dB
                if(draw_cells == 1); plot(xuser, yuser, "*k", 'MarkerSize', 3); end
                continue
            end
            if(draw_cells == 1); plot(xuser, yuser, "*g", 'MarkerSize', 3); end
        end
    end
    % Computing the SIR in lineal
    SIR = desired_signal/(tier1r_Int+tier2r_Int); % Desired signal / sum of signals not desired
    % Computing the SIR in dB and store it to samples array
    SIR_dB(loop) = 10*log10(SIR); % in dB
    %disp(SIR_dB);
    if(doThroughput) % If enabled compute the throughput and store it to samples array
        throughputs(loop) = throughput(bandwidth, SIR, SIRGap);
    end
end
figure(1); % Figure 1 where CDF are drawn
cdfplot(SIR_dB); % Compute and plot the CDF
hold on; % Hold to enable comparations
% Tagging the labels and title
xlabel('SIR (dB)');
ylabel('Probability');
title('CDF');
legend(legendString(2:end)); % Appending legend using the global variable 
if(doThroughput) % If computing the throughput is enabled
    figure(2); % Draw in figure 2, separated from other plot
    cdfplot(throughputs); % Compute the CDF of the throughput
    % Taging labels and adding title
    xlabel("Throughput");
    ylabel("Probability");
    throughputString(end+1) = "type " + fracName; % Acumulating string in legend
    legend(throughputString(2:end)); % Ploting the legend
    hold on; % Hold plot to enable comparations
end
end

