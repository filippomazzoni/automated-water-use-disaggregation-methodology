
%% F_ToiletParameterValuesIdentification
% This function includes the methodology for the automated identification
% of toilet parameter values (VTHmin, VTHmax, VTFmin, VTFmax) if
% disaggregation with general parameter values is required. Parameters
% values are fund by only by analysing the aggregate water use time series.
% Values are applied to global variables VTHmin, VTHmax, VTFmin, VTFmax.

function [A,A0] = F_ToiletParameterValuesIdentification(A,A0)

    global VTHmin VTHmax
    global VTFmin VTFmax

    %% Function beginning
       
    % All the water uses lasting between DTmin and DTmax are selected.
    % extracted from matrix A0. Those events can be related to single
    DTmin = 1;  % (min)
    DTmax = 2;  % (min)
    
    a0 = size(A0,1);
    cont = 1;
    cont_vector = 0; 
    
    while cont < a0
        Event = [];
        cont_event = cont;
        while cont_event < a0 && A0(cont_event,1)==A0(cont_event+1,1)-1
            cont_event = cont_event+1;
        end
        Event = A0(cont:cont_event,:);
        Duration = Event(end,1)-Event(1,1)+1;
        if Duration >= DTmin && Duration <= DTmax
           cont_vector = cont_vector+1;
           PossToilet_D(cont_vector) = Duration;
           PossToilet_V(cont_vector) = sum(Event(:,6));
        end
        cont = cont_event+1;
    end
        
    % The above water uses are divided into classes based on the volume of
    % water used. The lower and the upper thresholds are equalled to
    % parameters VTmin and VTmax. Matrix Freq includes the volume of each
    % class (first row) and the number of observed water uses (second row).
    VTmin = 3;   % (L)
    VTmax = 13;  % (L)
    
    for i = VTmin:VTmax
        Freq(1,i-2) = i;
        Freq(2,i-2) = length(find(PossToilet_V==i));
    end
    
    % Only water uses whose frequency of occurrence is at lest the 50% of
    % the mean frequency of occurrence are considered.
    tres = 0.50*mean(Freq(2,:));  
    
    % Volume classes associated with the highest frequency of occurrence
    % are given to parameters VTHmin,VTHmax,VTFmin,VTFmax.
    cont = 1;
    
    for i = 1:size(Freq,2)
          if Freq(2,i)>=tres
             PossVolume(cont) = Freq(1,i);
             cont = cont+1;
          end      
    end
    
    % First-attempt values (L)
    PossHF = PossVolume(PossVolume<=7);
    PossFF = PossVolume(PossVolume>7);
    
    % If there is a possible half flush without a possible full flush,
    % vector PossFF is substituted by vector PossHF
    if isempty(PossFF) && ~isempty(PossHF)     
       PossFF = PossHF;
       PossHF = [];           
    end
      
    
    %% Toilet full-flush values assignment
    
    pff = length(PossFF);
    
    % If pff only includes one element, such a value is given to both
    % VTFmin and VTFmax
    if pff ==1 
       VTFmin = PossFF;
       VTFmax = PossFF;
       
    % If PossFF includes two consequential elements, the lowest of them is
    % given to VTFmin and the highest of them is given to VTFmax.
    elseif pff==2 && PossFF(2)==PossFF(1)+1
       VTFmin = PossFF(1);
       VTFmax = PossFF(2);
       
    % If PossFF includes two non-consequential elements, the most frequent
    % of them is given to both VTFmin and VTFmax.   
    elseif pff==2 && PossFF(2)~=PossFF(1)+1
       pos1 = Freq(1,:)==PossFF(1);
       pos2 = Freq(1,:)==PossFF(2);
       Freq1 = Freq(2,pos1);
       Freq2 = Freq(2,pos2);  
       max_freq = max(Freq1,Freq2);
       max_pos = find(Freq(2,:)==max_freq);  
       VTFmin = PossFF(Freq(1,max_pos));
       VTFmax = PossFF(Freq(1,max_pos));
       
    % If PossFF includes three or more elements, the two-litre range
    % including the most frequent element is given to VTFmin and VTFmax. 
    elseif pff>2
        
       for i = 1:pff
           pos(i) = find(Freq(1,:)==PossFF(i));
           frequency(i) = Freq(2,pos(i));
       end
       [~, max_pos] = max(frequency);

       % If the most-frequent value of PossFF is the first of the vector,
       % that value is given to VTFmin. If the second element of vector
       % PossFF is related to a volume which is consequential to VTFmin,
       % that value is given to VTFmax. Otherwise, the value of the first
       % element is given to both VTFmin and VTFmax.
       if max_pos == 1       
          if PossFF(max_pos+1)==PossFF(max_pos)+1
             VTFmin = PossFF(max_pos);
             VTFmax = PossFF(max_pos+1);
          else
             VTFmin = PossFF(max_pos);
             VTFmax = PossFF(max_pos);
          end
         
        % If the most-frequent value of PossFF is in an intermediate
        % position, the previous and the following element are considered.
        % If both the previous and the following elements are equal to the
        % most-frequent value minus and plus one respectively, the one
        % whose frequency of occurrence is higher is considered.
       elseif max_pos>1 && max_pos<pff
          if PossFF(max_pos+1)==PossFF(max_pos)+1 && PossFF(max_pos)~=PossFF(max_pos-1)+1
             VTFmin = PossFF(max_pos);
             VTFmax = PossFF(max_pos+1);     
          elseif PossFF(max_pos+1)~=PossFF(max_pos)+1 && PossFF(max_pos)==PossFF(max_pos-1)+1
             VTFmin = PossFF(max_pos-1);
             VTFmax = PossFF(max_pos);    
          elseif PossFF(max_pos+1)==PossFF(max_pos)+1 && PossFF(max_pos)==PossFF(max_pos-1)+1
             freq_left = Freq(2,Freq(1,:)==PossFF(max_pos-1));
             freq_right = Freq(2,Freq(1,:)==PossFF(max_pos+1));
             higher_freq = max(freq_left,freq_right);
             other_value = Freq(1,Freq(2,:)==higher_freq);
             % If more than one value is found in matrix Freq, only the
             % values associated with a volume close to the most-frequent
             % is considered.
             if length(other_value)>1
                cont_to_delete = 0;
                for j = 1:length(other_value)
                    if other_value(j)<Freq(1,max_pos-1) || other_value(j)>Freq(1,max_pos+1)
                       cont_to_delete = cont_to_delete + 1;
                    end
                end
                if size(cont_to_delete,1)>0
                   other_value(cont_to_delete) = [];
                end
             end
             % If the previous and the following elements in vector PossFF
             % have the same frequency of occurrence, only one is
             % considered.
             if length(other_value)>1
                other_value(2) = [];
             end
             if other_value > PossFF(max_pos)
                VTFmin = PossFF(max_pos);
                VTFmax = other_value;
             else
                VTFmin = other_value;
                VTFmax = PossFF(max_pos);
             end   
          end
          
       % If the most-frequent value of PossFF is the last of the vector,
       % that value is given to VTFmax. If the previous element of vector
       % PossFF is related to a volume which is equal to VTFmin minus one,
       % that value is given to VTFmin. Otherwise, the value of the last
       % element is given to both VTFmin and VTFmax.    
       elseif max_pos == pff
          if PossFF(max_pos)==PossFF(max_pos-1)+1
             VTFmin = PossFF(max_pos-1);
             VTFmax = PossFF(max_pos);
          else
             VTFmin = PossFF(max_pos);
             VTFmax = PossFF(max_pos);
          end
       end
    end
    
    
    %% Toilet full-flush values assignment
    
    % This part of the function is considered only if vector PossHF is not
    % null.
    
    if isempty(PossHF)
       VTHmin = 0;
       VTHmax = 0;
       
    else
       
        phf = length(PossHF);
        
        % If phf only includes one element, such a value is given to both
        % VTHmin and VTHmax
        if length(PossHF) ==1
           VTHmin = PossHF;
           VTHmax = PossHF;

        % If PossHF includes two consequential elements, the lowest of them
        % is given to VTHmin and the highest of them is given to VTHmax.
        elseif phf==2 && PossHF(2)==PossHF(1)+1
           VTHmin = PossHF(1);
           VTHmax = PossHF(2);

        % If PossHF includes two non-consequential elements, the most
        % frequent of them is given to both VTHmin and VTHmax.   
        elseif phf==2 && PossHF(2)~=PossHF(1)+1
           pos1 = Freq(1,:)==PossHF(1);
           pos2 = Freq(1,:)==PossHF(2);
           Freq1 = Freq(2,pos1);
           Freq2 = Freq(2,pos2);
           max_freq = max(Freq1,Freq2);
           max_pos = find(Freq(2,:)==max_freq);
           VTHmin = PossHF(Freq(1,max_pos));
           VTHmax = PossHF(Freq(1,max_pos));

        % If PossHF includes three or more elements, the two-litre range
        % including the most frequent element is given to VTHmin and
        % VTHmax.
        elseif phf>2

           for i = 1:phf
               pos(i) = find(Freq(1,:)==PossHF(i));
               frequency(i) = Freq(2,pos(i));
           end
           [~, max_pos] = max(frequency);
           
           % If the most-frequent value of PossHF is the first of the
           % vector, that value is given to VTHmin. If the second element
           % of vector PossHF is related to a volume which is
           % consequential to VTHmin, tat value is given to VTHmax.
           % Otherwise, the value of the first element is given to both
           % VTHmin and VTHmax.
           if max_pos == 1
              if PossHF(max_pos+1)==PossHF(max_pos)+1
                 VTHmin = PossHF(max_pos);
                 VTHmax = PossHF(max_pos+1);
              else
                 VTHmin = PossHF(max_pos);
                 VTHmax = PossHF(max_pos);
              end

           % If the most-frequent value of PossHF is in an intermediate
           % position, the previous and the following element are
           % considered. If both the previous and the following elements
           % are equal to the most-frequent value minus and plus one 
           % respectively, the one whose frequency of occurrence is higher
           % is considered.
           elseif max_pos>1 && max_pos<phf
              if PossHF(max_pos+1)==PossHF(max_pos)+1 && PossHF(max_pos)~=PossHF(max_pos-1)+1
                 VTHmin = PossHF(max_pos);
                 VTHmax = PossHF(max_pos+1);
              elseif PossHF(max_pos+1)~=PossHF(max_pos)+1 && PossHF(max_pos)==PossHF(max_pos-1)+1
                 VTHmin = PossHF(max_pos-1);
                 VTHmax = PossHF(max_pos);
              elseif PossHF(max_pos+1)==PossHF(max_pos)+1 && PossHF(max_pos)==PossHF(max_pos-1)+1                                   
                 freq_left = Freq(2,find(Freq(1,:)==PossHF(max_pos-1)));
                 freq_right = Freq(2,find(Freq(1,:)==PossHF(max_pos+1)));
                 higher_freq = max(freq_left,freq_right);
                 other_value = Freq(1,Freq(2,:)==higher_freq);
                 % If more than one value is found in matrix Freq, only the
                 % values associated with a volume close to the
                 % most-frequent is considered.
                 if length(other_value)>1
                    cont_to_delete = 0;
                    for j = 1:length(other_value)
                        if other_value(j)<Freq(1,max_pos-1) || other_value(j)>Freq(1,max_pos+1)
                           cont_to_delete = cont_to_delete + 1;
                        end
                    end
                   if size(cont_to_delete,1)>0
                       other_value(cont_to_delete) = [];
                   end
                 end
                 % If the previous and the following elements in vector
                 % PossHF have the same frequency of occurrence, only one
                 % is considered.    
                 if length(other_value)>1
                    other_value(2) = [];
                 end
                 if other_value > PossHF(max_pos)
                    VTHmin = PossHF(max_pos);
                    VTHmax = other_value;
                 else
                    VTHmin = other_value;
                    VTHmax = PossHF(max_pos);
                 end
              end

           % If the most-frequent value of PossHF is the last of the
           % vector, that value is given to VTHmax. If the previous
           % element of vector PossHF is related to a volume which is
           % equal to VTHmax minus one, that value is given to VTHmin.
           % Otherwise, the value of the last element is given to both
           % VTHmin and VTHmax.    
           elseif max_pos == phf
               if PossHF(max_pos)==PossHF(max_pos-1)+1
                   VTHmin = PossHF(max_pos-1);
                   VTHmax = PossHF(max_pos);
               else
                   VTHmin = PossHF(max_pos);
                   VTHmax = PossHF(max_pos);
               end
           end
        end
        
    end

%% Function end
end