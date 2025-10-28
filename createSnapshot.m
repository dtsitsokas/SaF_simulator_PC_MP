function [] = createSnapshot(snap_t,DT,metric,metric_n,metric_name,LINKS,Nodes)
%snap_t: time of snapshot 
%metric: the requested value to create colormap 
%metric_n: the normalization value 
%LINKS
%Nodes

%normalize the metric to correspond to 0-1 units (for the color)
metric = metric(:,floor(snap_t/DT));
X1 = metric./metric_n;
figure('units','normalized','outerposition',[0 0 1 1],'color','white');
axis off
hold on

for i=1:size(LINKS,1)
    ind = X1(i);
    gplotdc_colormap([0 1;0 0],Nodes([find(Nodes(:,1)==LINKS(i,2)),find(Nodes(:,1)==LINKS(i,3))],[2 3]),ind);
    hold on;
end
c = colorbar('TickLabelInterpreter' ,'latex');
c.Label.String = metric_name; %'link occupancy $x/c$';
c.Label.Interpreter = 'latex';
%c.TickLabels = round((0:metric_n(1)/10:metric_n(1)));

%title(strcat(metric_name,'[0 to 1]'),'fontweight','bold','fontsize',20);
% pbaspect([1 1 1]);
ylim([4.5805*10^6 4.5836*10^6])
xlim([4.277*10^5 4.33*10^5])

nSecs = snap_t*3600;
text(430200,4581000,['Time: ',datestr(nSecs/86400, 'HH:MM:SS')],'fontsize',20,'EdgeColor','k','Margin',5,'LineWidth',1);


% Plot lines to indicate borders between central region and periphery 
x1_l1 = Nodes(Nodes(:,1)==19076,2); 
y1_l1 = Nodes(Nodes(:,1)==19076,3);
x2_l1 = Nodes(Nodes(:,1)==45906,2);
y2_l1 = Nodes(Nodes(:,1)==45906,3);
x1_l2 = Nodes(Nodes(:,1)==46370,2);
y1_l2 = Nodes(Nodes(:,1)==46370,3);
x1_l3 = Nodes(Nodes(:,1)==46947,2);
y1_l3 = Nodes(Nodes(:,1)==46947,3);

%line_1
plot([x1_l1; x2_l1],[y1_l1; y2_l1],'Linewidth',6,'Color','k')

%line_2 
plot([x2_l1; x1_l2],[y2_l1; y1_l2],'Linewidth',6,'Color','k')

%line_3
plot([x1_l2; x1_l3],[y1_l2; y1_l3],'Linewidth',6,'Color','k')

%line_4 
plot([x1_l3; x1_l1],[y1_l3; y1_l1],'Linewidth',6,'Color','k')
end

