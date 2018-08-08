function [hop_xy, fly_xy] = ClusterFlyHop(datafly, datahop, clustersize, plotflag)

% make cluster of varying size for a data set
xfly = datafly.y(:,1);
yfly = datafly.y(:,2);
xhop = datahop.y(:,1);
yhop = datahop.y(:,2);
if plotflag>0
    figure(15)
    plot(xhop, yhop,'o')
    hold on
    plot(xfly, yfly,'o')
    drawnow
end

% find scaling factors for data based on max of all data
xscale = max([datafly.y(:,1); datahop.y(:,1)]); % for position
yscale = max([datafly.y(:,2); datahop.y(:,2)]); % for velocity
% find minima point during hop
minhop = find(datahop.y(:,1) == min(datahop.y(:,1)));
% find maxima data point during fly
maxfly = find(datafly.y(:,1) == max(datafly.y(:,1)));

% point(s) to find cluster(s) around
xpt_fly = datafly.y(maxfly,1)/xscale;
ypt_fly = datafly.y(maxfly,2)/yscale;

xpt_hop = datahop.y(minhop,1)/xscale;
ypt_hop = datahop.y(minhop,2)/yscale;



for ii = 1:length(clustersize)
    
    k = clustersize(ii);
    hop_xy{ii} = knnsearch([xhop/xscale,yhop/yscale],[xpt_hop, ypt_hop], 'K',k);
    fly_xy{ii} = knnsearch([xfly/xscale,yfly/yscale],[xpt_fly, ypt_fly], 'K',k);
    if plotflag>0
        figure(15)
        plot(xhop(hop_xy{ii}), yhop(hop_xy{ii}),'x')
        plot(xfly(fly_xy{ii}), yfly(fly_xy{ii}),'x')
        xlabel('position')
        ylabel('velocity')
        title('clustering in phase space')
    end
end