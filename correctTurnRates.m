function [junct] = correctTurnRates(junct,downstr)
%correctTurnRates: Apply corrections in order to keep all routes open
% prev = junct.turn;
%case Barcelona:

%remove zeros from the end
for i=1:length(junct.turn(:,1))
    %replace all zeros with the average of the simulation
    if sum(junct.turn(i,:)==0)<length(junct.turn(i,:)) %not all zeros
        if downstr(junct.or_index(i),2)==1
            %locate zeros
            %replace zeros with the average of non zeros
            junct.turn(i, junct.turn(i,:)==0) = mean(junct.turn(i,junct.turn(i,:)~=0));
        elseif downstr(junct.or_index(i),2)==2
            ind = find(junct.turn(i,:)==0); %indices of zeros
            sec_ind = setdiff(find(junct.origin==junct.origin(i)),i); %second index
            for j=1:length(ind)
                %check intersection
                if junct.turn(i,ind(j)) + junct.turn(sec_ind,ind(j)) == 0
                    junct.turn(i,ind(j)) = mean(junct.turn(i,junct.turn(i,:)~=0));
                    junct.turn(sec_ind,ind(j)) = 1 - junct.turn(i,ind(j));
                end
            end
        end
        
    end
end

%*** manual corrections for Barcelona network ***** 
junct.turn(11,:) = 0;
junct.turn(687,:) = 0.5; %687, 688 were both 0 in the same section (17027)
junct.turn(688,:) = 0.5; 
junct.turn(337,:) = 0.33; 
junct.turn(338,:) = 0.33; 
junct.turn(339,:) = 0.33; 


%after 10th interval, keep constant until 16 and then average until end
for i=1:length(junct.turn(:,1))
    junct.turn(i,12:15) = mean(junct.turn(i,9:11));
    junct.turn(i,16:end) = mean(junct.turn(i,1:10));
end

%every link that has only 1 and 0s and 1 downstream link -> put always 1
for i=1:length(junct.turn(:,1))
    check1 = junct.turn(i,:)==0;
    check2 = junct.turn(i,:)==1;
    if sum(check1+check2)==length(junct.turn(i,:))
        if downstr(junct.or_index(i),2)==1
            junct.turn(i,:) = 1;
        elseif downstr(junct.or_index(i),2)==2
            %for links with 2 downstream we must ensure that the sum of turn ratios is equal to 1 in
            %unsignalized
            if junct.splan(i) == 0 %not signalized
                ind = find(junct.origin == junct.origin(i));
                if sum(junct.turn(ind(1),:))==0
                    junct.turn(ind(2),:)=1;
                elseif sum(junct.turn(ind(2),:))==0
                    junct.turn(ind(1),:)=1;
                end
            end
        end
    end
end

return

