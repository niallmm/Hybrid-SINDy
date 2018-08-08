function [t, S, I, R, beta, dS, dI, dR] = RunSIR_wk(p, initial_cond,years, kicksize, eps)
% code to run an SEIR model with seasonal forcing


    
tsave =[];
S_save = [];
I_save = [];
R_save = [];
beta_save =[];
dS_save = [];
dI_save = [];
dR_save = [];

options = odeset('RelTol',1e-6,'AbsTol',1e-6,'Events', @events_diverge); 

[numic,n] = size(initial_cond);

for ic = 1:numic
    yinit = initial_cond(ic,:);
    
    tsave =[tsave; 0];
    S_save =[S_save; yinit(1)];
    I_save =[I_save; yinit(2)];
    R_save =[R_save; yinit(3)];
    beta_save =[beta_save p.beta_hat/(1+p.b)];
    dy = SIR_1(0,yinit,p);
    dS_save = [dS_save; dy(1)];
    dI_save = [dI_save; dy(2)];
    dR_save = [dR_save; dy(3)];
    size(dS_save)
    for kk = 1:years
        % winter break
        tspan = 0:35;
        [t,y] = ode45(@(t,y)SIR_1(t,y,p),tspan,yinit,options);
        dy = SIR_1(t,y,p);
        tsave = [tsave; t(2:end)+(kk-1)*365];
        S_save = [S_save; y(2:end,1)];
        I_save = [I_save; y(2:end,2)];
        R_save = [R_save; y(2:end,3)];
        beta_save = [beta_save p.beta_hat/(1+p.b)*ones(1, length(tspan)-1)];
        dy = SIR_1(t,y,p);
        dS_save = [dS_save; dy(1,2:end)'];
        dI_save = [dI_save; dy(2,2:end)'];
        dR_save = [dR_save; dy(3,2:end)'];
        % spring term
        tspan = 35:155;
        yinit = y(end,:) + randi(kicksize*[-1 1], 1,3);
        [t,y] = ode45(@(t,y)SIR1(t,y,p),tspan,yinit,options);
        tsave = [tsave; t(2:end)+(kk-1)*365];
        S_save = [S_save; y(2:end,1)];
        I_save = [I_save; y(2:end,2)];
        R_save = [R_save; y(2:end,3)];
        beta_save = [beta_save p.beta_hat*(1+p.b)*ones(1, length(tspan)-1)];
        dy = SIR1(t,y,p);
        dS_save = [dS_save; dy(1,2:end)'];
        dI_save = [dI_save; dy(2,2:end)'];
        dR_save = [dR_save; dy(3,2:end)'];
        % summer break
        tspan = 155:225;
        yinit = y(end,:) + randi(kicksize*[-1 1], 1,3);
        [t,y] = ode45(@(t,y)SIR_1(t,y,p),tspan,yinit,options);
        tsave = [tsave; t(2:end)+(kk-1)*365];
        S_save = [S_save; y(2:end,1)];
        I_save = [I_save; y(2:end,2)];
        R_save = [R_save; y(2:end,3)];
        beta_save = [beta_save p.beta_hat/(1+p.b)*ones(1, length(tspan)-1)];
        dy = SIR_1(t,y,p);
        dS_save = [dS_save; dy(1,2:end)'];
        dI_save = [dI_save; dy(2,2:end)'];
        dR_save = [dR_save; dy(3,2:end)'];
        % fall term
        tspan = 225:365;
        yinit = y(end,:) + randi(kicksize*[-1 1], 1,3);
        [t,y] = ode45(@(t,y)SIR1(t,y,p),tspan,yinit,options);
        tsave = [tsave; t(2:end)+(kk-1)*365];
        S_save = [S_save; y(2:end,1)];
        I_save = [I_save; y(2:end,2)];
        R_save = [R_save; y(2:end,3)];
        beta_save = [beta_save p.beta_hat*(1+p.b)*ones(1, length(tspan)-1)];
        dy = SIR1(t,y,p);
        dS_save = [dS_save; dy(1,2:end)'];
        dI_save = [dI_save; dy(2,2:end)'];
        dR_save = [dR_save; dy(3,2:end)'];
        yinit = y(end,:) + randi(kicksize*[-1 1], 1,3);
    end
end
dbeta_save = zeros(size(beta_save));
clear t S I R beta dS dI dR
t = tsave; S = S_save+eps*randn(size(S_save)); I = I_save+eps*randn(size(I_save));
R = R_save; beta = beta_save;
dS = dS_save; dI = dI_save; dR = dR_save;

