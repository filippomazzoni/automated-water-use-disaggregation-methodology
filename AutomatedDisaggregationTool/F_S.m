
%% F_S (function for shower water use detection and classification)

% This function is set to recognize and classify shower uses. These
% uses are searched within the specified period [tSmin;tSmax] according
% with the parameters chosen by the user.

function [A,A0,S,Vol_S] = F_S(A,A0,S,TotalConsumption)

global tSmin tSmax
global DSmin DSmax
global VSmin VSmax
global pSmax

%% Function beginning
a0 = size(A0,1);
cont = 1;

for i = 1:(a0-DSmin+1)  
    
    % Moving window of length DSmin
    for k = 1:(DSmin)
        MovingWindow(k,:)=A0(i+k-1,:);
        Duration = MovingWindow(end,1)-MovingWindow(1,1)+1;
    end
        
    % If there are no NaN
    if ~isnan(MovingWindow(1:end,6))    
        % If duration is <= DSmax
        if Duration <= DSmax
            % If volume is <= VSmax
            if sum(MovingWindow(:,6))<=VSmax
               % The water use could be a shower
               PossShower(cont:(cont+DSmin-1),:)=MovingWindow;
               cont = cont+DSmin;
            end
        end
    end        
end

PossShower = unique(PossShower,'rows');
ps = size(PossShower,1);

% Water uses are removed if the period is not between tSmin and tSmax.
rowstodelete = [];
cont = 1;

if tSmin <= tSmax
    for i = 1:ps
        if PossShower(i,4)< tSmin || PossShower(i,4)> tSmax
           rowstodelete(cont) = i;
           cont = cont+1;
        end
    end
elseif tSmin > tSmax
     for i = 1:ps
        if PossShower(i,4)> tSmax && PossShower(i,4)< tSmin
           rowstodelete(cont) = i;
           cont = cont+1;
        end
     end
end

PossShower(rowstodelete,:) = [];
ps = size(PossShower,1);

% The initial volume associated with shower is zero. Moreover, a maximum
% number of shower uses per day (nSmax) is considered here based on a
% selected volume treshold of daily consumption (100L).
Vol_S = 0;
avgS = 100;
avgConsumption = TotalConsumption/size(A,1)*1440; % (L/day)
nSmax = round(avgConsumption/avgS);

% At first, possible shower uses are selected
for i = 1:(ps)
    
    if ps>0
       Event = [];
        
        % Water use selection
        if PossShower(end,3) == PossShower(1,3)
           Event = PossShower;
        else
            cont = 1;
            Event(cont,:) = PossShower(cont,:);
            while PossShower(cont+1,3)==PossShower(cont,3)
                cont = cont+1;
                Event(cont,:) = PossShower(cont,:);
            end
        end
        
        % Water use classification
        Shower = [];
        SubEvent = [];
        Duration = [];
        cont_shower = 0;  
        step = [];
        step(1) = 1;  
        
        for j = 2:size(Event,1)
            step(j) = Event(j,1)-Event(j-1,1);
        end
        
        % If the water use includes interruptions of flow lasting more than
        % pSmax, the water use is split.
        k = find(step(:)>pSmax);
                
        if size(k,1) == 1
            
            % First part
            SubEvent = Event(1:(k(1)-1),:);
            Duration = SubEvent(end,1)-SubEvent(1,1)+1;
            
            % If duration is greater than DSmin
            if Duration >= DSmin
                % If volume is greater than VSmin
                if sum(SubEvent(:,6))>= VSmin
                   % If the number of already-occurred showers is lower
                   % than the limit
                   if cont_shower < nSmax
                       % Too long possible showers are reduced
                       if Duration > DSmax
                           cont = 1;
                           while Duration <= DSmax
                               cont = cont+1;
                               Duration = SubEvent(cont,1)-SubEvent(1,1)+1;
                           end
                           % If volume is greater than VSmin
                           if sum(SubEvent(1:cont,6))>= VSmin
                               % The water use is classified as shower use
                               cont_shower = cont_shower+1;
                               Shower = [Shower; SubEvent(1:cont,:)];
                           end
                       else
                           % The water use is classified as shower use
                           cont_shower = cont_shower+1;
                           Shower = [Shower; SubEvent];
                       end
                   end
                end
            end
            
            % Last part
            SubEvent = Event(k(1):end,:);
            Duration = SubEvent(end,1)-SubEvent(1,1)+1;
            
            if Duration >= DSmin
                if sum(SubEvent(:,6))> VSmin
                    if cont_shower < nSmax
                        if Duration > DSmax
                           cont = 1;
                           while Duration <= DSmax
                                cont = cont+1;
                                Duration = SubEvent(cont,1)-SubEvent(1,1)+1;
                           end
                           if sum(SubEvent(1:cont,6))>= VSmin
                              cont_shower = cont_shower+1;
                              Shower = [Shower; SubEvent(1:cont,:)];
                           end
                        else
                            cont_shower = cont_shower+1;
                            Shower = [Shower; SubEvent];
                        end
                    end
                end
            end
        
            
        elseif size(k,1)>=2
            
            % First part
            SubEvent = Event(1:(k(1)-1),:);
            Duration = SubEvent(end,1)-SubEvent(1,1)+1;
            
            if Duration >= DSmin
                if sum(SubEvent(:,6))> VSmin
                    if cont_shower < nSmax
                        if Duration > DSmax
                           cont = 1;
                           while Duration <= DSmax
                                cont = cont+1;
                                Duration = SubEvent(cont,1)-SubEvent(1,1)+1;
                           end
                           if sum(SubEvent(1:cont,6))>= VSmin
                              cont_shower = cont_shower+1;
                              Shower = [Shower; SubEvent(1:cont,:)];
                           end
                        else
                            cont_shower = cont_shower+1;
                            Shower = [Shower; SubEvent];
                        end
                    end
                end
            end
            
            % Intermediate parts
            for j = 2:(size(k,1))
                SubEvent = Event(k(j-1):(k(j)-1),:);
                Duration = SubEvent(end,1)-SubEvent(1,1)+1;
                
                if Duration >= DSmin
                    if sum(SubEvent(:,6))> VSmin
                        if cont_shower < nSmax
                            if Duration >DSmax
                               cont = 1;
                               while Duration <= DSmax
                                     cont = cont+1;
                                     Duration = SubEvent(cont,1)-SubEvent(1,1)+1;
                               end
                               if sum(SubEvent(1:cont,6))>= VSmin
                                  cont_shower = cont_shower+1;
                                  Shower = [Shower; SubEvent(1:cont,:)];
                               end
                            else
                                cont_shower = cont_shower+1;
                                Shower = [Shower; SubEvent];
                            end
                        end
                    end
                end
                
            end
            
            % Final part
            SubEvent = Event(k(end):end,:);
            Duration = SubEvent(end,1)-SubEvent(1,1)+1;
            
            if Duration >= DSmin
                if sum(SubEvent(:,6))> VSmin
                    if cont_shower < nSmax
                        if Duration > DSmax
                            cont = 1;
                            while Duration <= DSmax
                                cont = cont+1;
                                Duration = SubEvent(cont,1)-SubEvent(1,1)+1;
                            end
                            if sum(SubEvent(1:cont,6))>= VSmin
                               cont_shower = cont_shower+1;
                               Shower = [Shower; SubEvent(1:cont,:)];
                            end
                        else
                            cont_shower = cont_shower+1;
                            Shower = [Shower; SubEvent];
                        end
                    end
                end
            end
          
         
        % If the water use does not include flow interruptions lasting more
        % the pSmax, it is not split.
        else

           Duration = Event(end,1)-Event(1,1)+1;
           
           if Duration >= DSmin
               if sum(Event(:,6))> VSmin
                   if cont_shower < nSmax
                       if Duration > DSmax
                           cont = 1;
                           while Duration <= DSmax
                                 cont = cont+1;
                                 Duration = Event(cont,1)-Event(1,1)+1;
                           end
                           if sum(Event(1:cont,6))>= VSmin
                              cont_shower = cont_shower+1;
                              Shower = [Shower; Event(1:cont,:)];
                           end
                       else
                           cont_shower = cont_shower+1;
                           Shower = [Shower; Event];
                       end
                   end
               end
           end
           
        end
        
        % The detected shower use is allocated in the shower uses time
        % series (i.e. matrix S) and removed from the aggregate water use
        % time series (i.e. matrices A and A0)
        for s = 1:size(Shower,1)
            S(Shower(s,1))= Shower(s,6);
            A(Shower(s,1),6)= NaN;
            A0(A0(:,1)==Shower(s,1),6) = NaN;
        end
                
        % The detected shower use is also removed from matrix PossShower
        rowlisttodelete = [];
        for j = 1:size(Event,1)
            rowlisttodelete(j) = find(PossShower(:,1)==Event(j,1));
        end
        
        PossShower(rowlisttodelete,:) = [];
        ps = size(PossShower,1);
        
        % The aggregate volume of shower is updated
        if size(Shower)>0
            Vol_S = Vol_S + sum(Shower(:,6));
        end
        
    end
    
end

%% Function end

end