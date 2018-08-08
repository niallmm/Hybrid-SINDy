% Run Hopping Model 
% run this file first, then 
% Ex_Hopping_validation.m
% and
% Ex_Hopping_PLOTS.m

close all
clear all

p.plottag =1;
addpath('hopper'); addpath('SINDy'); addpath('validation');
% kappa = k*y0/m*g  % nondimensional param
p.kappa = 10;

p.eps = 1e-6; % noise level
rng(1)
% variables are: time, t, vertical position, y, vertical velocity, yp

tend = 5;
dt = 0.033;
p.tspan = [0:dt:tend];

% start with mass slightly below it's equilibrium position
% velocity downward
p.yinitvec = [0.8 -0.1
    0.78 -0.1
    0.82 -0.1];

p.options = odeset('RelTol',1e-12,'AbsTol',1e-12, 'Events', @events_hopper);

[data]  = RUN_Hopper(p);
t_out = data.tout; 
y_out = data.yout+ rand(size(data.yout))*p.eps;
a_out =data.aout;
 

%% Run CCM on y_out
clear x y
x= y_out(1:end-1,1); % height
y= y_out(1:end-1,2); % velocity

dimMan = 2; % dimension of the Manifold
numOfPts = length(x); % number of points in the data
%k = dimMan+2;
k = 30; % number of points in a cluster

idx_xy = knnsearch([x/max(x) y/max(y)],[x/max(x) y/max(y)],'K',k);
idx_xy = idx_xy(:,2:end);


%% Run SINDy on clusters
% delay coordinates are xdelay and ydelay
% clusters are defined by idx_x(delay,cluster_num) and
% idx_y(delay,cluster_num)

% initialize a structure library
II_Xi =[]; ind_all = [];
xt = [x; x(end)];
yt = [y; y(end)];
plottag =0;
for ii = 1:floor(numOfPts)
    
    % stack all the position values in cluster ii into a vector
    X = [x(idx_xy(ii,:)') y(idx_xy(ii,:)')];
    % stack all the velocity values in cluster ii into a vector
    XP = [a_out(idx_xy(ii,:)',1) a_out(idx_xy(ii,:)',2)];
    
    
    % create function library
    polyorder = 2;  % search space up to 2nd order polynomials
    usesine = 0;    % no trig functions
    laurent = 0;    % no laurent terms
    dyin =0;        % just regular SINDy
    dyorder = 0;    % no derivatives in library
    Theta = poolDatady(X,2,polyorder,usesine, laurent, dyin, dyorder);
    m = size(Theta,2);
    for rr = 1:m
        normTheta(rr) = mean(Theta(:,rr));
        ThetaN(:,rr) = Theta(:,rr)/mean(Theta(:,rr));
    end
    Thetalib.Theta = Theta;
    Thetalib.normTheta = 1;
    Thetalib.dx = XP;
    Thetalib.polyorder = polyorder;
    Thetalib.usesine = usesine;
    
    lambdavals.numlambda = 20;
    lambdavals.lambdastart = -5;
    lambdavals.lambdaend = 2;
    
    Xistruct =multiD_Lambda(Thetalib,lambdavals);
    Xicomb = Xistruct.Xicomb;
    
    % add new models to library and find the indices for this set
    [II_Xitemp, II_ind] = build_Xi_Library(Xicomb, II_Xi);
    II_Xi = II_Xitemp;
    ind_all = [ind_all; II_ind];

    Cluster{ii}.centroidx =  mean(X(:,1));
    Cluster{ii}.centroidy = mean(X(:,2));
     % store results by cluster number
    Cluster{ii}.Xistruct = Xistruct;
    Cluster{ii}.II_ind = II_ind;
    Cluster{ii}.X = X;
    Cluster{ii}.XP = XP;
end

dateformatout = 'mmddyyyy';
save([datestr(now, dateformatout) 'Hopping_data_eps1e_6.mat'])

