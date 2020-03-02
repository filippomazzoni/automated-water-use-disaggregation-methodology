
%% Automated Household Water End-Use Disaggregation Through A Rule-Based Automated Methodology
% Filippo Mazzoni, Stefano Alvisi, Marco Franchini, Marco Ferraris, and Zoran Kapelan

% The methodology for automated water use disaggregation detects,
% disaggregates and classifies individual water uses one end-use category
% at a time. Specifically, the disaggregation of the aggregate water use
% time series collected at the domestic water inlet point is done by means
% of a set of functions. Such functions are applied in an order, starting
% with functions aimed at detecting end-uses that, given their very nature,
% are generally more regular in terms of water use, and ending with the
% most irregular ones.

% At first, electronic appliances use is investigated through function
% dishwasher (F_DW) and function washing machine (F_WM). Then, shower uses
% are searched through function shower (F_S). At last, toilet and tap uses
% (which are generally less-recognisable at one-minute resolution or often
% overlapped in time with other uses) are detected and classified by means
% of the function tap and toilet (F_TT). In conclusion, the water use time
% series for each end-use category is available.

% Operatively, the automated methodology for water use disaggregation was
% developed by using MATLAB® programming software and consists of a main
% code, where the Microsoft Excel® datasheet including the collected
% aggregate water use is loaded and functions for water use disaggregation
% are applied in turn.

clear all
close all
clc


%% Choice of household, period, and parameter values

disp('Household selection');
disp('1=H1, 2=H2, 3=H3, 4=H4');
HH = input('');

if HH == 1 || HH == 2 || HH == 3
    
   disp('Period selection')
   disp('1=Calibration, 2=Validation')
   Period = input('');   
   if Period ~= 1 && Period ~= 2
      disp('Not valid!')
   end
   
   disp('Parameter values selection')
   disp('1=H1 specific, 2=H2 specific, 3=H3 specific, 0=General')
   Par = input('');
   if Par ~= 0 && Par ~= 1 && Par ~= 2 && Par ~= 3
      disp('Not valid!')
   end
   
elseif HH == 4
   disp('Test household (H4) chosen.')
   disp('Disaggregation will be done with general parameter values and for the overall period.') 
   Par = 0;
else
   disp('Not valid!')
end


%% Parameter values recall
% All the required parameters values (specific or general) are recalled
% from the library through the function below.
[Par] = F_Parameters(Par);


%% Water use time series recall
% Aggregate water use time series and individual end-uses time series
% are recalled from the Microsoft Excel® datasheets.
if HH == 1 && Period == 1
    [Baseline, txt, ~] = xlsread('H1data','Baseline','A2:H44641');
elseif HH == 1 && Period == 2
    [Baseline, txt, ~] = xlsread('H1data','Baseline','A44642:H80641');
elseif HH == 2 && Period == 1
    [Baseline, txt, ~] = xlsread('H2data','Baseline','A2:H44641');
elseif HH == 2 && Period == 2
    [Baseline, txt, ~] = xlsread('H2data','Baseline','A44642:H80641');
elseif HH == 3 && Period == 1
    [Baseline, txt, ~] = xlsread('H3data','Baseline','A2:H40321');
elseif HH == 3 && Period == 2
    [Baseline, txt, ~] = xlsread('H3data','Baseline','A40322:H76321');
elseif HH == 4
    [Baseline, txt, ~] = xlsread('H4data','Baseline','A2:H60481');
end

% Length of the selected period (min)
n = size(Baseline,1);

% Dates are turned into numbers
for i = 1:size(txt,1)
    if length(txt{i,1})==10
       txt{i,1} = strcat(txt{i,1},' 00:00:00');
    end
end
for i = 1:n
    Datetime(i,:) = datevec(txt(i,1),'dd/mm/yyyy HH:MM:SS'); %#ok<*SAGROW>
end

% Errors in volumes (empty cells, negative numbers) are fixed
for i=1:size(Baseline,1)
    for j = 1:size(Baseline,2)
        if isnan(Baseline(i,j)) || Baseline(i,j)<0
           Baseline(i,j)=0;
        end
    end
end

% Aggregate water use matrix (A). The first column includes a progressive
% number, the second-fifth columns include date-time references, the sixth
% column includes the aggregate water use volume.
A = [[1:1:n]' Datetime(:,2:5) Baseline(:,1)]; %#ok<*NBRAK>

% Aggregate water use matrix only including time steps with water use (A0)
A0 = A(A(:,6)>0,:);
a0 = size(A0,1);

% Total observed water use in the selected period
TotalConsumption = sum(A0(:,6));


%% Individual end-uses time series initialization
% Initialization of a vector for each end-use category.
DW = [zeros(n,1)];  % Dishwasher uses time series
KS = [zeros(n,1)];  % Kitchen sink uses time series
WM = [zeros(n,1)];  % Washing machine uses time series
S =  [zeros(n,1)];  % Shower uses time series
BT = [zeros(n,1)];  % Bathroom taps uses time series
T =  [zeros(n,1)];  % Toilet uses time series


%% F_ToiletParameterValuesIdentification
% Function including the automated methodology for the obtainment of toilet
% parameter values in case of disaggregation by means of general parameter
% values.
if Par == 0
   [A,A0] = F_ToiletParameterValuesIdentification(A,A0);
end


%% F_DW
% Function for dishwasher use detection and classification. Household H1 is
% excluded since it does not include a dishwasher.
if HH ~= 1 
   [A,A0,DW,~] = F_DW(A,A0,DW);
end


%% F_WM
% Function for washing machine use detection and classification
[A,A0,WM,~] = F_WM(A,A0,WM);


%% F_S
% Function for shower use detection and classification
[A,A0,S,~] = F_S(A,A0,S,TotalConsumption);


%% F_TT
% Function for tap and toilet use detection and classification
[A,A0,KS,BT,T,~] = F_TT(A,A0,KS,BT,T);


%% Results
% Aggregate water use and disaggregated end-uses time series matrix
fprintf('Disaggregation complete.\n');
Automated = [Baseline(:,1) DW KS WM S BT T];
Results(1,:) = sum(Baseline,1);
Results(2,:) = sum(Automated,1);

if Results(1,2) == 0 % No dishwasher
   RVol = array2table(Results(:,[1 3:end]),'VariableNames',{'Aggregate','Kitchen_Sink','Washing_Machine','Shower','Bathroom_Taps','Toilet'},...
                'RowNames',{'Baseline data (L) ','Disaggregated data (L)'});
else
   RVol = array2table(Results,'VariableNames',{'Aggregate','Dishwasher','Kitchen_Sink','Washing_Machine','Shower','Bathroom_Taps','Toilet'},...
                'RowNames',{'Baseline data (L)','Disaggregated data (L)'});
end
disp(RVol);


%% Metrics calculation
% Water Contribution Accuracy (WCA) and Normalized Root-Mean-Square Error
% (NRMSE) are calculated in the following.

for i = 1:7
    
    % WCA (Water Contribution Accuracy)  
    WCA(i) = 100 - abs(Results(1,i)-Results(2,i))/Results(1,1)*100;
    
    % SE (Square Error)    
    for j = 1:n
        SE(j) = (Baseline(j,i)-Automated(j,i))^2;
    end
    
    % RMSE (Root-Mean-Square Error)
    RMSE(i) = sqrt(1/n*sum(SE));

    % NRMSE (Normalized Root-Mean-Square Error)
    NRMSE(i) = RMSE(i)/(max(Baseline(:,i))-min(Baseline(:,i)));  
        
end

% Metrics output
WCA = round(WCA,1);
NRMSE = round(NRMSE,3);

if Results(1,2) == 0 % No dishwasher
   MM = array2table([WCA(3:end); NRMSE(3:end)],'VariableNames',{'Kitchen_Sink','Washing_Machine','Shower','Bathroom_Taps','Toilet'},...
                'RowNames',{'WCA (%)','NRMSE (-)'});
else
   MM = array2table([WCA(2:end); NRMSE(2:end)],'VariableNames',{'Dishwasher','Kitchen_Sink','Washing_Machine','Shower','Bathroom_Taps','Toilet'},...
                'RowNames',{'WCA (%)','NRMSE (-)'});
end
disp(MM)