
disp("This script offers a quick way to execute the functions");
type = 1; % Default reuse factor
disptype = 1; % Default reuse factor print on screen
v = 3.8; % Default pathloss exponent
doPowerControl = false; % Default, do not do power control
doThroughtput = false; % Default, do no compute throughput
while(true)
    disp("(1). Run with current configuration" + newline + ...
         "(2). Switch Reuse Factor(Currently: " + disptype + ")" + newline + ...
         "(3). Turn on/off Power Control (Currently: " + doPowerControl+")" + newline + ...
         "(4). Compute Troughtput (Currently: " + doThroughtput + ")" + newline + ...
         "(5). Change Pathloss Exponent (Currently: " + v + ")" + newline + ...
         "(6). Run this option everytime plot windows are closed" + newline + ...
         "(7). Exit");
      
    option = input("Select an option: ");
    switch(option)
        case 1
            plotSNR(type, v, doPowerControl, doThroughtput); % Execute SNR plot
        case 2
            disp("(1). Universal Bandiwth" + newline + ...
                 "(2). 1/3 of Reuse Factor" + newline + ...
                 "(3). 1/9 of Reuse Factor");
             rfactor = input("Select the reuse factor: ");
             switch rfactor
                 case 1
                     type = 1; % Universal bandwith
                     disptype = 1;
                 case 2
                     type = 1/3; % 1/3 of reuse factor
                     disptype = rat(type);
                     disptype = disptype(4:end);
                 case 3
                     type = 1/9; %1/9 of reuse factor
                     disptype = rat(type);
                     disptype = disptype(4:end);
             end
             doPowerControl = false; % Deactivate power control so cannot be launch with other factor
        case 3
            doPowerControl = ~doPowerControl; % Change state of doing reuse factor
            if(type ~= 1/3) % If enabled power control, force to use 1/3 of reuse factor
                disp("Option only available for reuse factor 1/3. Changing reuse factor...");
                type = 1/3;
                disptype = rat(type);
                disptype = disptype(4:end);
                pause(4);
            end
        case 4
            doThroughtput = ~doThroughtput; % Change state of computing throughput or not
        case 5
            disp("(1). 3" + newline + ...
                 "(2). 3.8" + newline + ...
                 "(3). 4.5");
             expo = input("Select the pathloss exponent: "); 
             switch expo % Select between the pathloss exponents (asked in exercises)
                 case 1
                     v = 3;
                 case 2
                    v = 3.8;
                 case 3
                    v = 4.5;
             end
        case 6 % Clear the global variables (used to store the legends)
            clear global;
        case 7 % Exit
            break;
    end
    clc;
end
