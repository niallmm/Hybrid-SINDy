% code to run an SIR model with varying transmission rates
% Run this file then
% Ex_SIR_validation
% and 
% Ex_SIR_PLOT

clear all
close all
changeplot

plottag = 1; % 2 - maximum details of switch point detection during validation
addpath('SIR')
addpath('SINDy')
addpath('clustering')
rmpath('hopper')
% parameters
p.omega = 1/8; % exposed period of 8 days
p.gamma = 1/5; % infectious period of 5 days
p.beta_hat =9.336; % mixing/transmission rage for 5-11 year olds
% p.B = 1/(55*365); % birth and death rate
p.B = 1/365;
p.b = 0.8; % seasonality amplitude b = [0, 1] in paper
p.N = 1000; % total population
p.d = p.B;
years = 5;

eps = 0; % noise level

% initial conditions
initial_cond = [12 13 975];
kicksize = 2;
[numic,n] = size(initial_cond);

rng(1)

[t, S, I, R, beta, dS, dI, dR] = RunSIR_wk(p, initial_cond, years, kicksize, eps);
tsave = t; S_save = S; I_save = I; R_save = R; dS_save = dS; dI_save = dI;
dR_save = dR; beta_save = beta;

%% Run CCM on y_out
clear dx dy
dx = dS_save;
dy = dI_save;

clear x y 
x = S_save;
y = I_save;

dimMan = 2; % dimension of the Manifold
numOfPts = length(x); % number of points in the data

k = 30; % number of points in a cluster

idx_xy = knnsearch([x/max(x) y/max(y)],[x/max(x) y/max(y)],'K',k);
idx_xy = idx_xy(:,2:end);

%% Run SINDy on clusters
% delay coordinates are xdelay and ydelay
% clusters are defined by idx_x(delay,cluster_num) and
% idx_y(delay,cluster_num)
dt =1;
% initialize a structure library
II_Xi =[]; ind_all = [];
xt = [x; x(end)];
yt = [y; y(end)];

n = 2;
for ii = 1:floor(numOfPts)
    disp(['Cluster num:' num2str(ii) ' of' num2str(numOfPts)])
    
      
    % stack all the position values in cluster ii into a vector
    X = [x(idx_xy(ii,:)') y(idx_xy(ii,:)')];% z(idx_xy(ii,:)')];
    % stack all the velocity values in cluster ii into a vector
    XP = [dx(idx_xy(ii,:)') dy(idx_xy(ii,:)')];% dz(idx_xy(ii,:)')];

    % create function library
    polyorder = 3;  % search space up to 3rd order polynomials
    usesine = 0;    % no trig functions
    laurent = 0;    % no laurent terms
    dyin =0;        % just regular SINDy
    dyorder = 0;    % no derivatives in library
    [Theta, Thetastring] = poolDatady(X,n,polyorder,usesine, laurent, dyin, dyorder);
    m = size(Theta,2);
    
    Thetalib.Theta = Theta;
    Thetalib.normTheta = 1;
    Thetalib.dx = XP;
    Thetalib.polyorder = polyorder;
    Thetalib.usesine = usesine;
    
    lambdavals.numlambda = 20;
    lambdavals.lambdastart = -8;
    lambdavals.lambdaend = 1;
    
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
t_out = t;

dateformatout = 'mmddyyyy';
save([datestr(now, dateformatout) 'SIR_modelselection_eps0.mat'])