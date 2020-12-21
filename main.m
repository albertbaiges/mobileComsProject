
disp("This script offers a quick way to execute the functions");
type = 1;
disptype = 1;
v = 3.8;
doPowerControl = false;
doThroughtput = false;
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
            plotSNR(type, v, doPowerControl, doThroughtput);
        case 2
            disp("(1). Universal Bandiwth" + newline + ...
                 "(2). 1/3 of Reuse Factor" + newline + ...
                 "(3). 1/9 of Reuse Factor");
             rfactor = input("Select the reuse factor: ");
             switch rfactor
                 case 1
                     type = 1;
                     disptype = 1;
                 case 2
                     type = 1/3;
                     disptype = rat(type);
                     disptype = disptype(4:end);
                 case 3
                     type = 1/9;
                     disptype = rat(type);
                     disptype = disptype(4:end);
             end
        case 3
            doPowerControl = ~doPowerControl;
            if(type ~= 1/3)
                disp("Option only available for reuse factor 1/3. Changing reuse factor...");
                type = 1/3;
                disptype = rat(type);
                disptype = disptype(4:end);
                pause(4);
            end
        case 4
            doThroughtput = ~doThroughtput;
        case 5
            disp("(1). 3" + newline + ...
                 "(2). 3.8" + newline + ...
                 "(3). 4.5");
             expo = input("Select the pathloss exponent: "); 
             switch expo
                 case 1
                     v = 3;
                 case 2
                    v = 3.8;
                 case 3
                    v = 4.5;
             end
        case 6
            clear global;
        case 7
            break;
    end
    clc;
end
