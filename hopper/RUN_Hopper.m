% Run Hopping Model
function[data] = RUN_Hopper(p);
% variables in data are: time, t, vertical position, y, vertical velocity, yp and
% acceleration

kappa = p.kappa;
tspan_in = p.tspan;
plottag = p.plottag;
yinitvec = p.yinitvec;
options = p.options;

dt = tspan_in(2)-tspan_in(1);
tend = tspan_in(end);

current_t = 0; % starting value of the current time
t_out = []; y_out = []; a_out = [];

for mm = 1:size(yinitvec,1)
    mm
    yinit  = yinitvec(mm,:)
    current_t = 0; % starting value of the current time
    t_out = [t_out; 0];
    y_out = [y_out; yinit];
    a_out = [a_out; hopperspring(0,yinit',kappa)'];
    t= 0;
    tspan = tspan_in;
    
    % m*xpp = -m*g - k(x-x0), x<=x0
    % m*xpp = -m*g,            x>x0
    while current_t < tend
        if length(tspan)>1 % check that we haven't reached the end within error
            if abs(yinit(1)-1)<1e-12 % if we are within error of transition
                if yinit(2)<0 % springing if velocity is negative
                   disp('hopper springing')
                    yinit(1) = 1-1e-12;
                    [t,y]=ode45(@(t,y) hopperspring(t,y,kappa),tspan,yinit,options);
                    a= hopperspring(t,y',kappa)';
                    figure(1)
                    plot(t,y(:,1))
                    hold on
                    drawnow
                else % flying
                    yinit(1) = 1+1e-12;
                    disp('hooper flying')
                    [t,y]=ode45(@(t,y) hopperflight(t,y),tspan,yinit,options);
                    a = hopperflight(t,y')';
                 
                    plot(t,y(:,1))
                    hold on
                    drawnow
                end
            elseif (yinit(1)<1) % if we are below 1 we are springing
               disp('hopper springing')
                [t,y]=ode45(@(t,y) hopperspring(t,y,kappa),tspan,yinit,options);
                a= hopperspring(t,y',kappa)';
                plot(t,y(:,1))
                    hold on
                    drawnow
            elseif (yinit(1) > 1)
                disp('hooper flying')
                [t,y]=ode45(@(t,y) hopperflight(t,y),tspan,yinit,options);
                a = hopperflight(t,y')';
                plot(t,y(:,1))
                    hold on
                    drawnow
            else
                disp('simulation ended early')
                current_t
                return
            end
                t_out = [t_out; t(2:end)];
                y_out = [y_out; y(2:end,:)];
                a_out = [a_out; a(2:end,:)];
                yinit = y(end,:);
                tspan = t(end):dt:tend+1;
                current_t = t_out(end);

         end
        
    end
   
    
end
 
% save data
data.yout = y_out;
data.aout = a_out;
data.tout = t_out;

