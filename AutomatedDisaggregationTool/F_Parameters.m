
%% Parameters

% All the parameters and their values, i.e. specific values (obtained for
% household H1, household H2, household H3) and general values are grouped
% in the current function and recalled by the main code before the
% application functions for water use disaggregation.

function [Par] = F_Parameters(Par)

    % Habits 
    global tLmin tLmax
    global tDmin tDmax
        
    % Dishwasher parameters
    global tDWmin tDWmax
    global VDWmin VDWmax
    global dDWmin dDWmax
    global DDWmax
    global pDW12min pDW12max
    global pDW23min pDW23max
    global nDWmax
    
    % Kitchen sink parameters
    global DKSmin DKSmax

    % Washing machine parameters
    global tWMmin tWMmax
    global VWMmin VWMmax
    global dWMmin dWMmax
    global DWMmax
    global pWM12min pWM12max
    global pWM23min pWM23max
    global nWMmax
    
    % Shower parameters
    global tSmin tSmax
    global DSmin DSmax
    global VSmin VSmax
    global pSmax
    
    % Toilet parameters
    global VTHmin VTHmax
    global VTFmin VTFmax


    %% H1 specific values
    if Par == 1

       % Habits
       tLmin = 11;                % (hh)
       tLmax = 14;                % (hh)
       tDmin = 19;                % (hh)
       tDmax = 22;                % (hh)
        
       % Kitchen sink
       DKSmin = 4;                % (min)
       DKSmax = 10;               % (min)

       % Washing machine
       tWMmin = 18;               % (hh)
       tWMmax = 0;                % (hh)
       VWMmin = 8;                % (L)
       VWMmax = 15;               % (L)
       dWMmin = 2;                % (min)
       dWMmax = 4;                % (min)
       pWM12min = 15;             % (min)
       pWM12max = 95;             % (min)
       pWM23min = 10;             % (min)
       pWM23max = 20;             % (min)
       DWMmax = 130;              % (min)
       nWMmax = 1;                % (-)
        
       % Shower
       tSmin = 12;                % (hh)
       tSmax = 0;                 % (hh)
       VSmin = 20;                % (L)
       VSmax = 90;                % (L)
       DSmin = 8;                 % (min)
       DSmax = 16;                % (min)
       pSmax = 3;                 % (min)
        
       % Toilet
       VTHmin = 3;                % (L)
       VTHmax = 3;                % (L)
       VTFmin = 8;                % (L)
       VTFmax = 9;                % (L) 
    
       
    %% H2 specific values
    elseif Par == 2       

       % Habits
       tLmin = 11;                % (hh)
       tLmax = 14;                % (hh)
       tDmin = 19;                % (hh)
       tDmax = 22;                % (hh)
        
       % Dishwasher
       tDWmin = 0;                % (hh)
       tDWmax = 24;               % (hh)
       VDWmin = 2;                % (L)
       VDWmax = 5;                % (L)
       dDWmin = 2;                % (min)
       dDWmax = 5;                % (min)
       pDW12min = 20;             % (min)
       pDW12max = 110;            % (min)
       pDW23min = 10;             % (min)
       pDW23max = 15;             % (min)
       DDWmax = 150;              % (min)
       nDWmax = 2;                % (-)
       
       % Kitchen sink
       DKSmin = 4;                % (min)
       DKSmax = 10;               % (min)
       
       % Washing machine
       tWMmin = 9;                % (hh)
       tWMmax = 20;               % (hh)
       VWMmin = 10;               % (L)
       VWMmax = 20;               % (L)
       dWMmin = 2;                % (min)
       dWMmax = 5;                % (min)
       pWM12min = 10;             % (min)
       pWM12max = 120;            % (min)
       pWM23min = 10;             % (min)
       pWM23max = 20;             % (min)
       DWMmax = 140;              % (min)
       nWMmax = 2;                % (-)
       
       % Shower
       tSmin = 7;                 % (hh)
       tSmax = 20;                % (hh)
       VSmin = 25;                % (L)
       VSmax = 50;                % (L)
       DSmin = 3;                 % (min)
       DSmax = 6;                 % (min)
       pSmax = 2;                 % (min)
       
       % Toilet
       VTHmin = 0;                % (L)
       VTHmax = 0;                % (L) 
       VTFmin = 3;                % (L)
       VTFmax = 3;                % (L)
 
       
   %% H3 specific values
   elseif Par == 3
     
       % Habits
       tLmin = 11;                % (hh)
       tLmax = 14;                % (hh)
       tDmin = 19;                % (hh)
       tDmax = 22;                % (hh)

       % Dishwasher
       tDWmin = 0;                % (hh)
       tDWmax = 24;               % (hh)
       VDWmin = 2;                % (L)
       VDWmax = 4;                % (L)
       dDWmin = 2;                % (min)
       dDWmax = 4;                % (min)
       pDW12min = 20;             % (min)
       pDW12max = 80;             % (min)
       pDW23min = 15;             % (min)
       pDW23max = 85;             % (min)
       DDWmax = 150;              % (min)
       nDWmax = 2;                % (-)
       
              
       % Kitchen sink
       DKSmin = 2;                % (min)
       DKSmax = 15;               % (min)
       
       % Washing machine
       tWMmin = 0;                % (hh)
       tWMmax = 24;               % (hh)
       VWMmin = 5;                % (L)
       VWMmax = 18;               % (L)
       dWMmin = 2;                % (min)
       dWMmax = 4;                % (min)
       pWM12min = 15;             % (min)
       pWM12max = 95;             % (min)
       pWM23min = 10;             % (min)
       pWM23max = 20;             % (min)
       DWMmax = 120;              % (min)
       nWMmax = 2;                % (-)
              
       % Shower
       tSmin = 8;                 % (hh)
       tSmax = 21;                % (hh)
       VSmin = 20;                % (L)
       VSmax = 90;                % (L)
       DSmin = 5;                 % (min)
       DSmax = 20;                % (min)
       pSmax = 2;                 % (min)
        
       % Toilet
       VTHmin = 0;                % (L)
       VTHmax = 0;                % (L)
       VTFmin = 9;                % (L)
       VTFmax = 12;               % (L)
        
       
    %% General parameter values
    elseif Par == 0

       % Habits
       tLmin = 11;                % (hh)
       tLmax = 14;                % (hh)
       tDmin = 19;                % (hh)
       tDmax = 22;                % (hh)

       % Dishwasher
       tDWmin = 0;                % (hh)
       tDWmax = 24;               % (hh)
       VDWmin = 1;                % (L)
       VDWmax = 5;                % (L)
       dDWmin = 1;                % (min)
       dDWmax = 5;                % (min)
       pDW12min = 10;             % (min)
       pDW12max = 120;            % (min)
       pDW23min = 10;             % (min)
       pDW23max = 120;            % (min)
       DDWmax = 180;              % (min)
       nDWmax = 3;                % (-)
       
       % Kitchen sink
       DKSmin = 2;                % (min)
       DKSmax = 15;               % (min)
       
       % Washing machine
       tWMmin = 0;                % (hh)
       tWMmax = 24;               % (hh)
       VWMmin = 5;                % (L)
       VWMmax = 20;               % (L)
       dWMmin = 1;                % (min)
       dWMmax = 5;                % (min)
       pWM12min = 10;             % (min)
       pWM12max = 120;            % (min)
       pWM23min = 10;             % (min)
       pWM23max = 120;            % (min)
       DWMmax = 180;              % (min)
       nWMmax = 5;                % (-)

       % Shower
       tSmin = 0;                 % (hh)
       tSmax = 24;                % (hh)
       VSmin = 20;                % (L)
       VSmax = 150;               % (L)
       DSmin = 3;                 % (min)
       DSmax = 30;                % (min)
       pSmax = 3;                 % (min)
       
    else
    disp('Number is not valid!');
    
    end

end