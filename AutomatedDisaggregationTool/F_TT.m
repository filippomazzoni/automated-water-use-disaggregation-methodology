
%% F_TT (function for tap and toilet water use detection and classification)

% This function is set to recognize and classify tap uses, toilet uses and
% combined (i.e. tap + toilet) uses. Besides, taps are divided between
% kitchen sink (KS) and bathroom taps (BT). Disaggregation is led by means
% of the parameter values chosen by the user.


function [A,A0,KS,BT,T,Vol_TT] = F_TT(A,A0,KS,BT,T)

global tLmin tLmax
global tDmin tDmax
global DKSmin DKSmax
global VTHmin VTHmax
global VTFmin VTFmax

%% Function beginning

% Matrix including all the residual water uses to be disaggregated
a0 = size(A0,1);
cont = 1;

for i = 1:a0
    if ~isnan(A0(i,6))
       ResidualEvents(cont,:) = A0(i,:);
       cont = cont+1;
    end
end

re = size(ResidualEvents,1);
Vol_TT = sum(ResidualEvents(:,6));


%% Water uses related with kitchen sink at meal time
cont = 1;
for i = 1:re
    % If the water use occurs at meal time it is likely to be related to
    % the kitchen sink
    if (A0(i,4)>=tLmin && A0(i,4)<=tLmax) || (A0(i,4)>=tDmin && A0(i,4)<=tDmax)
       PossibleKitchenSink(cont,:) = A0(i,:);
       cont = cont+1;
    end
end

pm = size(PossibleKitchenSink,1);

conttodelete = [];
for i = 1:pm
    if isnan(PossibleKitchenSink(i,6))
       conttodelete = [conttodelete; i];
    end    
end
PossibleKitchenSink(conttodelete,:) = [];
pm = size(PossibleKitchenSink,1);

% If the duration of a water use occurring at meal time is between DKSmin
% and DKSmax, the water use is likely to be related with a kitchen sink
% (e.g. meal preparation activity). Uses of matrix PossibleKitchenSink are
% analyzed in turn and, if the above conditions are respected, such uses
% are classified as kitchen sink uses.
Meal = [];
cont = 1;

while cont < pm
      Event = [];
      cont_event = cont;
      while cont_event < pm && PossibleKitchenSink(cont_event,1)==PossibleKitchenSink(cont_event+1,1)-1
            cont_event = cont_event+1;
      end
      Event = PossibleKitchenSink(cont:cont_event,:);
      Duration = Event(end,1)-Event(1,1)+1;
      if Duration >= DKSmin && Duration <= DKSmax
         Meal = [Meal; Event];
      end
      cont = cont_event+1;
end
m = size(Meal,1);

% The disaggregated above uses are allocated in vector KS and removed from
% the aggregate water use time series (i.e. matrices A and A0).
for i = 1:m
    KS(Meal(i,1))= Meal(i,6);
    A(Meal(i,1),6)= NaN;
    A0(find(A0(:,1)==Meal(i,1)),6) = NaN;
    ResidualEvents(find(ResidualEvents(:,1)==Meal(i,1)),6) = NaN;
end

% The above uses are also removed from matrix ResidualEvents.
rowstodelete = [];
cont = 1;
for i = 1:re
    if isnan(ResidualEvents(i,6))
        rowstodelete(cont) = i;
        cont = cont+1;
    end  
end

ResidualEvents(rowstodelete,:) = [];
re = size(ResidualEvents,1);
ResidualEvents;


%% Short (i.e. one-minute) residual uses selection and classification

% Selection
 cont = 1;
 for i = 1
     if ResidualEvents(i+1,1)>ResidualEvents(i,1)+1
        Event(cont,:)=ResidualEvents(i,:);
        cont = cont+1;
     end
 end
  for i = 2:(re-1)
     if ResidualEvents(i,1)>ResidualEvents(i-1,1)+1 && ResidualEvents(i+1,1)>ResidualEvents(i,1)+1
        Event(cont,:)=ResidualEvents(i,:);
        cont = cont+1;
     end
  end
 for i = re
    if ResidualEvents(i,1)>ResidualEvents(i-1,1)+1
       Event(cont,:)=ResidualEvents(i,:);
       cont = cont+1;
    end
 end
 
% Classification
e = size(Event,1);

for i = 1:e
    
    if Event(i,6) < VTHmin
        % One-minute uses including small volumes and occurring at meal time
        % are likely to be associated with the use of kitchen sink.
        % Otherwise, they are likely to be associated with the use of a
        % bathroom tap.
        if Event(i,6)<=3 && ((Event(i,4)>=tLmin && Event(i,4)<=tLmax) || (Event(i,4)>=tDmin && Event(i,4)<=tDmax))
            KS(Event(i,1)) = Event(i,6);
        else
            BT(Event(i,1)) = Event(i,6);
        end
        
    elseif Event(i,6) >= VTHmin && Event(i,6) <= VTHmax
        % One-minute uses including a volume between VTHmin and VTHmax are
        % likely to be associated with the use of toilet (half-flush)
        T(Event(i,1)) = Event(i,6);

    elseif Event(i,6) > VTHmax && Event(i,6) < VTFmin
        % One-minute uses including a volume between VTHmax and VTFmin are
        % likely to be associated with the use of a tap (kitchen sink or
        % bathroom tap) based on the time of the day.
        if ((Event(i,4)>=tLmin && Event(i,4)<=tLmax) || (Event(i,4)>=tDmin && Event(i,4)<=tDmax))
            KS(Event(i,1)) = Event(i,6);
        else
            BT(Event(i,1)) = Event(i,6);
        end
        
    % One-minute uses including a volume between VTFmin and VTFmax are
    % likely to be associated with the use of toilet (full-flush)    
    elseif Event(i,6) >= VTFmin && Event(i,6) <= VTFmax
        T(Event(i,1)) = Event(i,6);
     
    % One-minute uses including a volume greater than VTFmax are likely to
    % be associated with the use of toilet (full-flush). The rest is
    % associated with a tap based on the time of the day.
    elseif Event(i,6) > VTFmax
        T(Event(i,1)) = VTFmax;
        if ((Event(i,4)>=tLmin && Event(i,4)<=tLmax) || (Event(i,4)>=tDmin && Event(i,4)<=tDmax))
            KS(Event(i,1)) = Event(i,6)-VTFmax;
        else
            BT(Event(i,1)) = Event(i,6)-VTFmax;
        end
    end
    
    % The water use is removed from the aggregate water use time series
    % (i.e. matrices A and A0).
    Event(i,6) = NaN;
    A(Event(i,1),6) = NaN;
    A0(find(A0(:,1)==Event(i,1)),6) = NaN;
    ResidualEvents(find(ResidualEvents(:,1)==Event(i,1)),6) = NaN;
    
end

% Water uses previously allocated in the individual end-uses time series
% (i.e. KS, BT, T) are also removed from matrix ResidualEvents.
rowstodelete = [];
cont = 1;
for i = 1:re
    if isnan(ResidualEvents(i,6))
        rowstodelete(cont) = i;
        cont = cont+1;
    end
end
ResidualEvents(rowstodelete,:) = [];
re = size(ResidualEvents,1);
ResidualEvents;


%% Longer (i.e. multiple-minute) residual uses selection and classification

for i = 1:(re-1)

        if re>=1
             Event = [];
             
             % Longer residual uses selection
             if ResidualEvents(end,1)-ResidualEvents(1,1) == size(ResidualEvents,1)-1
                 Event = ResidualEvents;
             else
                 cont = 1;
                 Event(cont,:) = ResidualEvents(cont,:);
                 while ResidualEvents(cont+1,1)==ResidualEvents(cont,1)+1
                     cont = cont+1;
                     Event(cont,:) = ResidualEvents(cont,:);
                 end  
                 Event;
             end

             % Two-minute uses classification
             
             if size(Event,1) == 2
                 
                   % Two-minute uses including a volume lower than VTHmin
                   % are associated with a tap (kitchen sink or bathroom
                   % tap) based on the time of the day.
                   if sum(Event(:,6)) < VTHmin
                      if (Event(1,4)>=tLmin && Event(1,4)<=tLmax) || (Event(1,4)>=tDmin && Event(1,4)<=tDmax)
                         KS(Event(1,1)) = Event(1,6);
                         KS(Event(2,1)) = Event(2,6);         
                      else
                         BT(Event(1,1)) = Event(1,6);
                         BT(Event(2,1)) = Event(2,6);                    
                      end
                      
                   % Two-minute uses including a volume between VTHmin and
                   % VTHmax are likely to be associated with the use of
                   % toilet (half-flush)   
                   elseif sum(Event(:,6)) >= VTHmin && sum(Event(:,6)) <= VTHmax
                      [Vmax, Posmax] = max(Event(:,6));
                      if (Event(1,4)>=tLmin && Event(1,4)<=tLmax) || (Event(1,4)>=tDmin && Event(1,4)<=tDmax)
                          if Posmax == 1
                             Posmin = 2;
                          elseif Posmax == 2
                             Posmin = 1;
                          end
                          if Vmax >= VTHmin
                             T(Event(Posmax,1)) = VTHmin;
                             KS(Event(Posmax,1)) = Vmax-VTHmin;
                             KS(Event(Posmin,1)) = Event(Posmin,6);
                          elseif Vmax < VTHmin                             
                             T(Event(Posmax,1)) = Vmax;
                             T(Event(Posmin,1)) = VTHmin-Vmax;
                             KS(Event(Posmin,1)) = Event(Posmin,6)-T(Event(Posmin,1));
                          end       
                      else
                          if Posmax == 1
                             Posmin = 2;
                          elseif Posmax == 2
                             Posmin = 1;
                          end
                          if Vmax >= VTHmin
                             T(Event(Posmax,1)) = VTHmin;
                             BT(Event(Posmax,1)) = Vmax-VTHmin;
                             BT(Event(Posmin,1)) = Event(Posmin,6);
                          elseif Vmax < VTHmin                             
                             T(Event(Posmax,1)) = Vmax;
                             T(Event(Posmin,1)) = VTHmin-Vmax;
                             BT(Event(Posmin,1)) = Event(Posmin,6)-T(Event(Posmin,1));
                          end   
                      end 
                   
                   % One-minute uses including a volume between VTHmax and
                   % VTFmin are likely to be associated with the use of a
                   % tap (kitchen sink or bathroom tap) based on the time
                   % of the day.
                   elseif sum(Event(:,6)) > VTHmax && sum(Event(:,6)) < VTFmin
                      if ((Event(1,4)>=tLmin && Event(1,4)<=tLmax) || (Event(1,4)>=tDmin && Event(1,4)<=tDmax))
                         KS(Event(1,1)) = Event(1,6);
                         KS(Event(2,1)) = Event(2,6);
                      else  
                         BT(Event(1,1)) = Event(1,6);
                         BT(Event(2,1)) = Event(2,6);
                      end 

                   % Two-minute uses including a volume between VTFmin and
                   % VTFmax are likely to be associated with the use of
                   % toilet (half-flush)
                   elseif sum(Event(:,6)) >= VTFmin && sum(Event(:,6)) <= VTFmax
                      [Vmax, Posmax] = max(Event(:,6));
                      if (Event(1,4)>=tLmin && Event(1,4)<=tLmax) || (Event(1,4)>=tDmin && Event(1,4)<=tDmax)  
                         if Posmax == 1
                            Posmin = 2;
                         elseif Posmax == 2
                            Posmin = 1;
                         end
                         if Vmax >= VTFmin
                            T(Event(Posmax,1)) = VTFmin;
                            KS(Event(Posmax,1)) = Vmax-VTFmin;
                            KS(Event(Posmin,1)) = Event(Posmin,6);
                         elseif Vmax < VTFmin                             
                            T(Event(Posmax,1)) = Vmax;
                            T(Event(Posmin,1)) = VTFmin-Vmax;
                            KS(Event(Posmin,1)) = Event(Posmin,6)-T(Event(Posmin,1));
                         end        
                      else
                         [Vmax, Posmax] = max(Event(:,6));
                         if Posmax == 1
                            Posmin = 2;
                         elseif Posmax == 2
                            Posmin = 1;
                         end
                         if Vmax >= VTFmin
                            T(Event(Posmax,1)) = VTFmin;
                            BT(Event(Posmax,1)) = Vmax-VTFmin;
                            BT(Event(Posmin,1)) = Event(Posmin,6);
                         elseif Vmax < VTFmin                             
                            T(Event(Posmax,1)) = Vmax;
                            T(Event(Posmin,1)) = VTFmin-Vmax;
                            BT(Event(Posmin,1)) = Event(Posmin,6)-T(Event(Posmin,1));
                         end   
                      end 
                   
                   % One-minute uses including a volume greater than VTFmax
                   % are likely to be associated with the use of toilet
                   % (full-flush). The rest is associated with a tap based
                   % on the time of the day.
                   elseif sum(Event(:,6)) > VTFmax
                      [Vmax, Posmax] = max(Event(:,6));  
                      if Posmax == 1
                          Posmin = 2;
                      elseif Posmax == 2
                          Posmin = 1;
                      end
                      if Vmax > VTFmax
                         if (Event(1,4)>=tLmin && Event(1,4)<=tLmax) || (Event(1,4)>=tDmin && Event(1,4)<=tDmax)
                             T(Event(Posmax,1)) = VTFmax;
                             KS(Event(Posmax,1)) = Vmax-VTFmax;
                             KS(Event(Posmin,1)) = Event(Posmin,6);
                         else
                             T(Event(Posmax,1)) = VTFmax;
                             BT(Event(Posmax,1)) = Vmax-VTFmax;
                             BT(Event(Posmin,1)) = Event(Posmin,6);
                         end
                      elseif Vmax >= VTFmin && Vmax <= VTFmax
                          if (Event(1,4)>=tLmin && Event(1,4)<=tLmax) || (Event(1,4)>=tDmin && Event(1,4)<=tDmax)
                             T(Event(Posmax,1)) = Vmax;
                             KS(Event(Posmin,1)) = Event(Posmin,6);   
                          else                         
                             T(Event(Posmax,1)) = Vmax;
                             BT(Event(Posmin,1)) = Event(Posmin,6);   
                          end
                      elseif Vmax < VTFmin 
                          if (Event(1,4)>=tLmin && Event(1,4)<=tLmax) || (Event(1,4)>=tDmin && Event(1,4)<=tDmax)
                             T(Event(Posmax,1)) = Vmax;
                             T(Event(Posmin,1)) = VTFmin-Vmax;
                             KS(Event(Posmin,1)) = Event(Posmin,6)-T(Event(Posmin,1));
                          else                         
                             T(Event(Posmax,1)) = Vmax;
                             T(Event(Posmin,1)) = VTFmin-Vmax;
                             BT(Event(Posmin,1)) = Event(Posmin,6)-T(Event(Posmin,1)); 
                          end
                      end

                   end
                        
                   % The water use is removed from the aggregate water use
                   % time series
                   for j = 1:2
                       A(Event(j,1),6) = NaN;
                       A0(find(A0(:,1)==Event(j,1)),6) = NaN;
                       ResidualEvents(find(ResidualEvents(:,1)==Event(j,1)),6) = NaN;
                   end
                   Event(1:2,:)=[];
                   
                   
             % Three-minute (or longer) uses classification. The use of
             % kitchen sink is not considered for the above water uses,
             % since long kitchen sink uses have already been found in the
             % previous.
             elseif size(Event,1) >= 3
                 
             % A maximum number of allowable toilet flushes per
             % three-minute (or longer) water use (i.e. xTmax) is
             % considered.
             xTmax = 3;
             cont_Toilet = 0;
                
             % The three-minute (or longer) water use is progressively
             % disaggregated.
             E = Event;
             
             while size(E,1) >= 1
                 
                 % If the water use only includes one minute (since the
                 % previous ones have already been allocated), it is
                 % classified as follows.
                 if size(E,1) == 1 || (E(2,1)>E(1,1)+1)
                     if E(1,6) < VTHmin
                         BT(E(1,1)) = E(1,6); 
                     elseif E(1,6) >= VTHmin && E(1,6) <= VTHmax && cont_Toilet <= xTmax
                         T(E(1,1)) = E(1,6);
                         cont_Toilet = cont_Toilet+1;  
                     elseif E(1,6) >= VTHmin && E(1,6) <= VTHmax && cont_Toilet > xTmax
                         BT(E(1,1)) = E(1,6); 
                     elseif E(1,6) > VTHmax && E(1,6) < VTFmin
                         BT(E(1,1)) = E(1,6);  
                     elseif E(1,6) >= VTFmin && E(1,6) <= VTFmax && cont_Toilet <= xTmax
                         T(E(1,1)) = E(1,6);
                         cont_Toilet = cont_Toilet+1;    
                     elseif E(1,6) >= VTFmin && E(1,6) <= VTFmax && cont_Toilet > xTmax
                         BT(E(1,1)) = E(1,6);     
                     elseif E(1,6) > VTFmax && cont_Toilet <= xTmax
                         T(E(1,1)) = VTFmax;
                         BT(E(1,1)) = E(1,6)-VTFmax;
                         cont_Toilet = cont_Toilet+1;
                     elseif E(1,6) > VTFmax && cont_Toilet > xTmax
                         BT(E(1,1)) = E(1,6);
                     end
                     
                     % The water use is removed from the aggregate water
                     % use time series
                     A(E(1,1),6) = NaN;
                     A0(find(A0(:,1)==E(1,1)),6) = NaN;
                     ResidualEvents(find(ResidualEvents(:,1)==E(1,1)),6) = NaN;
                     E(1,:)=[];
                  
                 % If the water use still includes more than one minute,
                 % the first two minutes are considered and allocated as
                 % follows.
                 elseif E(2,1) == E(1,1)+1
                     if E(1,6)+E(2,6) < VTHmin
                        BT(E(1,1)) = E(1,6);
                        BT(E(2,1)) = E(2,6); 
                     elseif E(1,6)+E(2,6) >= VTHmin && E(1,6)+E(2,6) <= VTHmax && cont_Toilet <= xTmax
                        T(E(1,1)) = E(1,6);
                        T(E(2,1)) = E(2,6);
                        cont_Toilet = cont_Toilet+1;
                     elseif E(1,6)+E(2,6) >= VTHmin && E(1,6)+E(2,6) <= VTHmax && cont_Toilet > xTmax
                        BT(E(1,1)) = E(1,6);
                        BT(E(2,1)) = E(2,6);
                     elseif E(1,6)+E(2,6) > VTHmax && E(1,6)+E(2,6) < VTFmin
                        BT(E(1,1)) = E(1,6);
                        BT(E(2,1)) = E(2,6);
                     elseif E(1,6)+E(2,6) >= VTFmin && E(1,6)+E(2,6) <= VTFmax && cont_Toilet <= xTmax
                        T(E(1,1)) = E(1,6);
                        T(E(2,1)) = E(2,6);
                        cont_Toilet = cont_Toilet+1;
                     elseif E(1,6)+E(2,6) >= VTFmin && E(1,6)+E(2,6) <= VTFmax && cont_Toilet > xTmax
                        BT(E(1,1)) = E(1,6);
                        BT(E(2,1)) = E(2,6);
                     elseif E(1,6)+E(2,6) > VTFmax && cont_Toilet <= xTmax
                        [Vmax, Posmax] = max(E(1:2,6));
                        if Posmax == 1
                           Posmin = 2;
                        elseif Posmax == 2
                           Posmin = 1;
                        end
                        if Vmax > VTFmax
                           T(E(Posmax,1)) = VTFmax;
                           BT(E(Posmax,1)) = Vmax-VTFmax;
                           BT(E(Posmin,1)) = E(Posmin,6);
                           cont_Toilet = cont_Toilet+1; 
                        elseif Vmax >= VTFmin && Vmax <= VTFmax
                           T(E(Posmax,1)) = Vmax;
                           BT(E(Posmin,1)) = E(Posmin,6);
                           cont_Toilet = cont_Toilet+1; 
                        elseif Vmax < VTFmin
                           T(E(Posmax,1)) = Vmax;
                           T(E(Posmin,1)) = VTFmin-Vmax;
                           BT(E(Posmin,1)) = E(Posmin,6)-T(E(Posmin,1));
                           cont_Toilet = cont_Toilet+1;
                        end 
                     elseif E(1,6)+E(2,6) > VTFmax && cont_Toilet > xTmax
                        BT(E(1,1)) = E(1,6);
                        BT(E(2,1)) = E(2,6);
                     end
                     
                     % The water use is removed from the aggregate water
                     % use time series
                     for j = 1:2
                         A(E(j,1),6) = NaN;
                         A0(find(A0(:,1)==E(j,1)),6) = NaN;
                         ResidualEvents(find(ResidualEvents(:,1)==E(j,1)),6) = NaN;
                     end
                     E(1:2,:)=[];
                     
                 end
             end
             
             Event = [];
             
             end
             
             % Water uses previously allocated in the individual end-uses
             % time series are also removed from matrix ResidualEvents.
             rowstodelete = [];
             cont = 1;
             for i = 1:re
                 if isnan(ResidualEvents(i,6))
                     rowstodelete(cont) = i;
                     cont = cont+1;
                 end
             end
             ResidualEvents(rowstodelete,:) = [];
             re = size(ResidualEvents,1);
               
        end     
end

%% Function end

end