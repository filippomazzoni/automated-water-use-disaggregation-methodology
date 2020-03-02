
%% F_WM (function for washing machine water use detection and classification)

% This function is set to recognize and classify washing machine uses.
% These uses are searched within the specified period [tWMmin;tWMmax]
% according with the parameters chosen by the user. The function is set up
% to detect washing machine cycles whose withdrawals are in range [3;5].


function [A,A0,WM,Vol_WM] = F_WM(A,A0,WM)

global tWMmin tWMmax
global VWMmin VWMmax
global dWMmin dWMmax
global DWMmax
global pWM12min pWM12max
global pWM23min pWM23max
global nWMmax


%% Function beginning
a0 = size(A0,1);
cont = 1;

for i = 1:(a0-dWMmin+1)  
    
    % Moving window of duration dWMmin
    for k = 1:dWMmin
        MovingWindow(k,:)=A0(i+k-1,:);  
    end
       
    % If there are no NaN
    if ~isnan(MovingWindow(1:end,6))
        % If duration is <= dWMmax
        if (MovingWindow(end,1)-MovingWindow(1,1))<=dWMmax
            % If volume is between VWMmin and VWMmax
            if sum(MovingWindow(:,6))>=VWMmin && sum(MovingWindow(:,6))<=VWMmax
               % The use could be a washing machine withdrawal
               PossWM(cont:(cont+dWMmin-1),:)=MovingWindow;
               cont = cont+dWMmin;
            end
        end
    end
         
end

PossWM = unique(PossWM,'rows');
pWM = size(PossWM,1);

% Water uses are removed if the period is not between tWMmin and tWMmax
rowstodelete = [];
cont = 1;

if tWMmin <= tWMmax 
    for i = 1:pWM
        if PossWM(i,4)< tWMmin || PossWM(i,4)> tWMmax
           rowstodelete(cont) = i;
           cont = cont+1;
        end
    end  
elseif tWMmin > tWMmax
    
     for i = 1:pWM
        if PossWM(i,4)> tWMmax && PossWM(i,4)< tWMmin
           rowstodelete(cont) = i;
           cont = cont+1;
        end
     end
end
        
PossWM(rowstodelete,:)=[];
pWM = size(PossWM,1);


% If size(PossWM,1)=0, it means there are no washing machine uses in that
% day. Otherwise, withdrawals number and time of occurrence are studied.

if pWM >= 1

   pWM = size(PossWM,1);

   % All possible withdrawals are inserted in matrix PossFlush.
   PossFlush = []; 

   for i = 1:(pWM-1)
       SubEvent = [];
       
       if pWM >= 1
          if PossWM(end,1)-PossWM(1,1) == size(PossWM,1)-1
             SubEvent = PossWM;
          else
             cont = 1;
             SubEvent(cont,:) = PossWM(cont,:);
             while PossWM(cont+1,1)-PossWM(cont,1)==1
                   SubEvent(cont+1,:) = PossWM(cont+1,:);
                   cont = cont+1;
             end
          end

          PossWM(1:size(SubEvent,1),:) = [];
          pWM = size(PossWM,1); 

          % If a withdrawal lasts more than dWMmax, it is reduced
          while size(SubEvent,1) > dWMmax
                SubEvent(end,:) = [];
          end

          % If the volume is greater than VWMmax, it has to be reduced
          while sum(SubEvent(1:end,6)) > VWMmax
                SubEvent(end,:) = [];
          end

          % If the volume is lower than VWMmin, the use is removed
          if sum(SubEvent(1:end,6)) >= VWMmin
             PossFlush = [PossFlush; SubEvent];
          end
          
       end
   end
    
   % All possible washing machine withdrawals are in matrix PossFlush  
   pf = size(PossFlush,1);

   % If the number of minutes of withdrawal is lower than 3*dWMmin, there
   % cannot be a washing machine cycle. 
   if pf >= (3*dWMmin)

        % If distances between a number of possible withdrawal (i.e.
        % between 3 and 5) are in a range, a washing machine cycle is found.
        
        Cycles = []; 
        cont = 1;
        cont_WM = 0;
 
        F1I = [];        F1F = [];        F2I = [];        F2F = [];        
        F3I = [];        F3F = [];        F4I = [];        F4F = [];
        F5I = [];        F5F = [];
        
        
        for i = 1:(pf-5)
            for j = (i+1):(pf-4)
                 if PossFlush(i+1,1)-PossFlush(i,1) >= (pWM12min-1) && PossFlush(i+1,1)-PossFlush(i,1) <= (pWM12max+1)
                     if PossFlush(j+1,1)-PossFlush(j,1) >= (pWM23min-1) && PossFlush(j+1,1)-PossFlush(j,1) <= (pWM23max+1)
                        
                        % Case of three withdrawals     
                          if PossFlush(j,1)-PossFlush(i+1,1) <= dWMmax
                             F1F(cont,:) = PossFlush(i,:);   % First withdrawal end                           
                             F2I(cont,:) = PossFlush(i+1,:); % Second withdrawal beginning                          
                             F2F(cont,:) = PossFlush(j,:);   % Second withdrawal end
                             F3I(cont,:) = PossFlush(j+1,:); % Third withdrawal beginning
                             
                             % Case of four withdrawals
                             for k = (j+1):(pf-3)
                                 if PossFlush(k+1,1)-PossFlush(k,1) >= (pWM23min-1) && PossFlush(k+1,1)-PossFlush(k,1) <= (pWM23max+1)
                                    if PossFlush(k,1)-PossFlush(j+1,1) <= dWMmax
                                       F3F(cont,:) = PossFlush(k,:);   % Third withdrawal end
                                       F4I(cont,:) = PossFlush(k+1,:); % Fourth withdrawal beginning                                                 
                                       
                                       % Case of five withdrawals
                                       for m = (k+1):(pf-2)
                                           if PossFlush(m+1,1)-PossFlush(m,1) >= (pWM23min-1) && PossFlush(m+1,1)-PossFlush(m,1) <= (pWM23max+1)
                                              if PossFlush(m,1)-PossFlush(k+1,1) <= dWMmax
                                                 F4F(cont,:) = PossFlush(m,:);   % Fourth withdrawal end 
                                                 F5I(cont,:) = PossFlush(m+1,:); % Fifth withdrawal beginning 
                                              end
                                           end
                                       end
                                       
                                    end
                                 end                         
                             end
                                                                                                 
                             cont = cont+1;                      

                           end
                      end
                  end
             end
        end               

      % Withdrawal cycles are selected

      if size(F1F,1)>0
         dim = size(F1F,1);   
         
         % Due to implementation reason, all the matrices between F1F and
         % F3I (or F4I or F5I according to the case) should have the same
         % dimension of F1F, so null rows are added.

         F2I = [F2I; zeros(dim-size(F2I,1),6)];         F2F = [F2F; zeros(dim-size(F2F,1),6)];
         F3I = [F3I; zeros(dim-size(F3I,1),6)];         F3F = [F3F; zeros(dim-size(F3F,1),6)];
         F4I = [F4I; zeros(dim-size(F4I,1),6)];         F4F = [F4F; zeros(dim-size(F4F,1),6)];
         F5I = [F5I; zeros(dim-size(F5I,1),6)];
         
         % Initialization of F1I and F5F
         F1I = zeros(dim,6);
         F5F = zeros(dim,6);
         
         for s = 1:dim             
             
             % If in a row F3I>0 but F4I=0, a 3-withdrawal cycle is found.
             % If in a row F4I>0 but F5I=0, a 4-withdrawal cycle is found.
             % If in a row F5I>0, a 5-withdrawal cycle is found.
                         
             if F3I(s,1)>0 && F3F(s,1)>0
                 
                if F4I(s,1)>0 && F4F(s,1)>0    

                         % Possible 5-withdrawal cycle
                         % F1I e F5F have to be found
                         contF1 = find(PossFlush(:,1)==F1F(s,1));
                         contF5 = find(PossFlush(:,1)==F5I(s,1));
                         
                         % First (F1) and last (F5) withdrawal selection
                         if contF1 >= 2
                            while contF1>1 && PossFlush(contF1,1)-PossFlush(contF1-1,1)==1
                                  contF1 = contF1-1;
                            end     
                         end
                         if contF5 <= (pf-1)  
                            while contF5<pf && PossFlush(contF5+1,1)-PossFlush(contF5,1)==1
                                  contF5 = contF5+1;
                             end   
                         end
                                                  
                         % First and last minute of the cycle
                         F1I(s,:) = PossFlush(contF1,:);         
                         F5F(s,:) = PossFlush(contF5,:);   
                                       
                else
                   
                         % Possible 4-withdrawal cycle
                         % F1I e F4F have to be found
                         contF1 = find(PossFlush(:,1)==F1F(s,1));
                         contF4 = find(PossFlush(:,1)==F4I(s,1));

                         % First (F1) and last (F4) withdrawal selection
                         if contF1 >= 2
                            while contF1>1 && PossFlush(contF1,1)-PossFlush(contF1-1,1)==1
                                  contF1 = contF1-1;
                            end     
                         end
                         if contF4 <= (pf-1)  
                            while contF4<pf && PossFlush(contF4+1,1)-PossFlush(contF4,1)==1
                                  contF4 = contF4+1;
                             end   
                         end

                         % First and last minute of the cycle
                         F1I(s,:) = PossFlush(contF1,:);         
                         F4F(s,:) = PossFlush(contF4,:);
                
                end
                
             else                
                 
                 % Possible 3-withdrawal cycle
                 % F1I e F3F have to be found
                 contF1 = find(PossFlush(:,1)==F1F(s,1));
                 contF3 = find(PossFlush(:,1)==F3I(s,1));
                 
                 % First (F1) and last (F3) withdrawal selection
                 if contF1 >= 2
                    while contF1>1 && PossFlush(contF1,1)-PossFlush(contF1-1,1)==1
                          contF1 = contF1-1;
                    end     
                 end
                 if contF3 <= (pf-1)  
                    while contF3<pf && PossFlush(contF3+1,1)-PossFlush(contF3,1)==1
                          contF3 = contF3+1;
                     end   
                 end
                
                 % First and last minute of the cycle
                 F1I(s,:) = PossFlush(contF1,:);         
                 F3F(s,:) = PossFlush(contF3,:);
                 
             end
         end
                   
        % If the cycle lasts more than DWMmax, it cannot be a washing machine
        % cycle, so it has to be shorter.

        if dim >= 1
            rowstodelete = [];
            cont = 1;          
            for s = 1:dim  
                if F5I(s,1)==0 && F4I(s,1)==0 && (F3F(s,1)-F1I(s,1) > (DWMmax+5))
                   % Too long 3-withdrawal cycle
                   rowstodelete(cont) = s;
                   cont = cont+1;
                elseif F5I(s,1)==0 && F4I(s,1)>0 && (F4F(s,1)-F1I(s,1) > (DWMmax+5))
                   % Too long 4-withdrawal cycle
                   rowstodelete(cont) = s;
                   cont = cont+1;
                elseif F5I(s,1)>0 && (F5F(s,1)-F1I(s,1) > (DWMmax+5))
                   % Too long 5-withdrawal cycle
                   rowstodelete(cont) = s;
                   cont = cont+1;
                end 
            end
            
            % Too-long cycles removal from matrix F1I to F5F
            F1I(rowstodelete,:) = [];            F1F(rowstodelete,:) = [];
            F2I(rowstodelete,:) = [];            F2F(rowstodelete,:) = [];
            F3I(rowstodelete,:) = [];            F3F(rowstodelete,:) = [];
            F4I(rowstodelete,:) = [];            F4F(rowstodelete,:) = [];
            F5I(rowstodelete,:) = [];            F5F(rowstodelete,:) = [];
        
            % New number of possible washing machine cycles. 
            dim = size(F1I,1);    
        
        end
        
        % It is also assumed that time overlapped cycles cannot occur.
        
        if dim >= 2
           rowstodelete = [];
           cont = 1;

           for s = 1:(dim-1)

                % If the number of the last minute (F3F, F4F, or F5F value)
                % is greater than the number of the first minute (F1I) of
                % the following cycle, the second possible cycle is removed.
                                
                if F5I(s,1)==0 && F4I(s,1)==0 && F1I(s+1,1)<=F3F(s,1)
                   rowstodelete(cont) = s+1;
                   cont = cont+1;
                elseif F5I(s,1)==0 && F4I(s,1)>0 && F1I(s+1,1)<=F4F(s,1)
                   rowstodelete(cont) = s+1;
                   cont = cont+1;
                elseif F5I(s,1)>0 && F1I(s+1,1)<=F5F(s,1)
                   rowstodelete(cont) = s+1;
                   cont = cont+1;
                end
                
           end
            
           % Time-overlapped  cycles removal
           F1I(rowstodelete,:) = [];           F1F(rowstodelete,:) = [];
           F2I(rowstodelete,:) = [];           F2F(rowstodelete,:) = [];
           F3I(rowstodelete,:) = [];           F3F(rowstodelete,:) = [];
           F4I(rowstodelete,:) = [];           F4F(rowstodelete,:) = [];
           F5I(rowstodelete,:) = [];           F5F(rowstodelete,:) = [];

           % Number of possible washing machine cycles. 
           dim = size(F1I,1);
 
        end
        
    if dim >= 1        
        
        % Initialization of the matrix containing all the washing machine
        % events and the cell containing information about the volume. 
        Cycles = [];
        Vol_Cycles = 0;
        
        % Now, for every cycle, all the withdrawals are selected and the
        % volume of the whole cycle is calculated.
        
        for s = 1:dim
            
            F1 = [];            F2 = [];            F3 = [];
            F4 = [];            F5 = [];            
            
            if F5F(s,1)==0 && F4F(s,1)==0
               F1 = PossFlush(find(PossFlush(:,1)==F1I(s,1)):find(PossFlush(:,1)==F1F(s,1)),:);
               F2 = PossFlush(find(PossFlush(:,1)==F2I(s,1)):find(PossFlush(:,1)==F2F(s,1)),:);
               F3 = PossFlush(find(PossFlush(:,1)==F3I(s,1)):find(PossFlush(:,1)==F3F(s,1)),:);
               
               cont_WM = cont_WM+1;
               Cycles = [ Cycles; ...
                          F1 cont_WM*ones(size(F1,1),1); ...
                          F2 cont_WM*ones(size(F2,1),1); ...
                          F3 cont_WM*ones(size(F3,1),1)];

            elseif F5F(s,1)==0 && F4F(s,1)>0      
               F1 = PossFlush(find(PossFlush(:,1)==F1I(s,1)):find(PossFlush(:,1)==F1F(s,1)),:);
               F2 = PossFlush(find(PossFlush(:,1)==F2I(s,1)):find(PossFlush(:,1)==F2F(s,1)),:);
               F3 = PossFlush(find(PossFlush(:,1)==F3I(s,1)):find(PossFlush(:,1)==F3F(s,1)),:);
               F4 = PossFlush(find(PossFlush(:,1)==F4I(s,1)):find(PossFlush(:,1)==F4F(s,1)),:);
               
               cont_WM = cont_WM+1;
               Cycles = [ Cycles; ...
                          F1 cont_WM*ones(size(F1,1),1); ...
                          F2 cont_WM*ones(size(F2,1),1); ...
                          F3 cont_WM*ones(size(F3,1),1); ...
                          F4 cont_WM*ones(size(F4,1),1)];
                          
            elseif F5F(s,1)>0
               F1 = PossFlush(find(PossFlush(:,1)==F1I(s,1)):find(PossFlush(:,1)==F1F(s,1)),:);
               F2 = PossFlush(find(PossFlush(:,1)==F2I(s,1)):find(PossFlush(:,1)==F2F(s,1)),:);
               F3 = PossFlush(find(PossFlush(:,1)==F3I(s,1)):find(PossFlush(:,1)==F3F(s,1)),:);
               F4 = PossFlush(find(PossFlush(:,1)==F4I(s,1)):find(PossFlush(:,1)==F4F(s,1)),:);
               F5 = PossFlush(find(PossFlush(:,1)==F5I(s,1)):find(PossFlush(:,1)==F5F(s,1)),:);
               
               cont_WM = cont_WM+1;
               Cycles = [ Cycles; ...
                          F1 cont_WM*ones(size(F1,1),1); ...
                          F2 cont_WM*ones(size(F2,1),1); ...
                          F3 cont_WM*ones(size(F3,1),1); ...
                          F4 cont_WM*ones(size(F4,1),1); ...
                          F5 cont_WM*ones(size(F5,1),1);];
              
            end
        end
        
        % A maximum of nWMmax cycles is allowed per day
        days = unique(Cycles(:,2:3),'rows');
        to_delete = [];
        
        for i = 1:size(days,1)
            day_analyzed = [];
            for j = 1:size(Cycles,1)
                if Cycles(j,2:3)==days(i,1:2)
                   day_analyzed = [day_analyzed; Cycles(j,:)];
                end
            end
            if size(day_analyzed,1)>0 && day_analyzed(end,end)-day_analyzed(1,end)>(nWMmax-1)
               lim_sup = day_analyzed(1,end)+1;
               counter = find(day_analyzed(:,end)>lim_sup);
               to_delete = [to_delete; day_analyzed(counter,:)];
            end       
        end
            
        if size(to_delete,1)>0
            for i = 1:size(to_delete,1)
                Cycles(Cycles(:,1)==to_delete(i,1),:)=[];
            end
        end
        
        % Total volume associated with washing machine cycles.
        Vol_Cycles = sum(Cycles(:,6));

        % The selected uses are allocated in washing machine time series and
        % removed from the aggregate time series.
        for k = 1:size(Cycles,1)
            A(Cycles(k,1),6) = NaN;
            A0(find(A0(:,1)==Cycles(k,1)),6) = NaN;
            WM(Cycles(k,1)) = Cycles(k,6);
        end
        
     end
    end
   end
end

% The overall volume associated with the washing machine is sent out through Vol_WM.
if Vol_Cycles == 0
   Vol_WM = 0;
else
   Vol_WM = Vol_Cycles;
end

%% Function end
end