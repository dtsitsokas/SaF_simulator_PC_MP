function [upair,counter_1,counter_2] = correct_outflows(ind_or_pos,ind_dest_pos,upair,k,satFlow,Links2,ind_or,ind_dest,dest_index,counter_1,counter_2)
    
    %adjust calculated flows 
    %check outflows 
    for i=ind_or_pos
        % adjust outflows of origin links   
        sum_ = sum(upair(ind_or(i,:),k));
        if satFlow(i)-sum_<-1
            disp(strcat('correction outflow, link = ', num2str(Links2(i,1)),', k = ', num2str(k)))
            upair(ind_or(i,:),k) = satFlow(i)*upair(ind_or(i,:),k)/sum_;
            counter_1 = counter_1 + 1; 
        end
    end
    
    %check inflows 
    for i=ind_dest_pos
        %fix inflows
        sum_ = sum(upair(ind_dest(i,:),k));
        if satFlow(indata.junct2.dest_index(find(ind_dest(i,:),1)))-sum_<-1
             %disp('correction')
            disp(strcat('correction inflow, link = ', num2str(Links2(i,1)),', k = ', num2str(k)))
            upair(ind_dest(i,:),k) = satFlow(dest_index(find(ind_dest(i,:),1)))* upair(ind_dest(i,:),k)/sum_; %divide according to analogy of calculated flows - not sure if this is the priority - maybe v. queues should wait for everyone else (mostly for nodes with centroids)
            counter_2 = counter_2 + 1; 
        end
        
    end
    
end

