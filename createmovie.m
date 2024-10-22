function [v] = createmovie(x,step,start_step,capacity,LINKS,Nodes,DT,fname)
% x = current queues
% kmax = end of interval for video
% capacity = link storage capacity
% Links
% Nodes
% DT
% step
% fname
global X1 ind

link_occupancies=zeros(size(x,1),ceil(size(x,2)/step));
%Occupancy matrix over time for all links: link_occupancies
k = 1;
for i=step:step:size(x,2)
    %mean link occupancy over the step time window
    link_occupancies(:,k) = mean(x(:,i-step+1:i),2)./capacity;
    k = k + 1; 
end

for ti=1:length(link_occupancies(1,:))%size(Data1,2)
    X1 = link_occupancies(:,ti); %occupancies at t=ti
    fig = figure('units','normalized','outerposition',[0 0 1 1],'color','white');
    
    axis off
    hold on
    
    for i=1:size(LINKS,1)
        ind = X1(i);
        gplotdc_colormap([0 1;0 0],Nodes([find(Nodes(:,1)==LINKS(i,2)),find(Nodes(:,1)==LINKS(i,3))],[2 3]),ind);
        hold on;
    end
    colorbar
    
    h1 = title('link occupancy [0 to 1]','fontweight','bold','fontsize',20);
    pbaspect([1 1 1]);
    ylim([4.5805*10^6 4.5836*10^6])
    xlim([4.277*10^5 4.33*10^5])
    
    nSecs = (start_step+ti*step)*DT*3600;
    v = axis;
    text(430200,4581000,['Time: ',datestr(nSecs/86400, 'HH:MM:SS')],'fontsize',20,'EdgeColor','k','Margin',5,'LineWidth',1);
    M(ti) = getframe(fig);
    ti
    close all
end


[h, w, p] = size(M(1,1).cdata);  % use 1st frame to get dimensions
hf = figure;
set(hf, 'position', [90 140 w h]);
axis off
% movie(hf,M,1,8);
save Mfinal.mat M
%movie2avi(M(1,:), 'myMovie.avi','FPS',6);

%make movie file
v = VideoWriter(strcat(fname,'_movie.avi'));

%v.CompressionRatio = 3;
v.FrameRate = 5;

open(v)
writeVideo(v,M(1,:));
close(v);


end

