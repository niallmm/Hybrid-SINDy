function dy = SIR1(t,y,p)
if length(y) == 3
    S = y(1);
    I = y(2);
    R = y(3);
elseif length(y)>3
    S = y(:,1);
    I = y(:,2);
    R = y(:,3);
else
    error('SIR1 model has wrong input');
end


beta = p.beta_hat*(1+p.b);
dS = p.B*p.N - beta*I.*S/p.N -p.d*S;
dI = beta*I.*S/p.N  - (p.gamma+p.d)*I;
dR = p.gamma*I - p.d*R;

dy = [dS dI dR]';