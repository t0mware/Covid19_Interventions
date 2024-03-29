%ODE SEIRD model code with waning immunity

function [Classes] = ODE_SIRDV_immunity_model_2(para,ICs,beta)

%Run ODE using ODE45
opts = odeset('RelTol',1e-5);
%Only give outputs at defined timepoints
[t, pop] = ode45(@diff_SIRDV_immunity_model_2, [para.T0:1:para.T], [ICs.S1 0 ICs.I1 0 0 0 0 0 beta], opts,para,beta);

%Convert output to structure
Classes = struct('S1',pop(:,1),'S2',pop(:,2),'I1',pop(:,3),'I2',pop(:,4),'R',pop(:,5),'D',pop(:,6),'V',pop(:,7),'CV',pop(:,8),'l',pop(:,9),'t',t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Diff equations

function dPop = diff_SIRDV_immunity_model_2(t,pop,para,beta)

%Assign the population matrix into the classes
S1=pop(1);
S2=pop(2);
I1=pop(3);
I2=pop(4);
R=pop(5);
D=pop(6);
V=pop(7);
CV=pop(8);
l=pop(9);

% Vaccination rates
if t>para.T_start && t<para.T_stop
    %Check whether vaccine pool is used up
    if CV<para.v_pool
        if S1+S2+R>para.v_max 
            v1=para.v_max*S1/(S1+S2+R);
            v2=para.v_max*S2/(S1+S2+R);
            v3=para.v_max*R/(S1+S2+R);
        else 
            v1=S1;
            v2=S2;
            v3=R;
        end
    else
        v1=0;
        v2=0;
        v3=0;
    end
    
else
    v1=0;
    v2=0;
    v3=0;
end

%Define beta so that we lockdown at highest rate possible then increase
%beta when possible to give the lowest beta that keeps R<0.9
if t<52
    b=beta;
elseif t>=52 && S1+S2 > 0 && 0.9*para.N*(para.gamma+para.delta)/(S1+S2)<beta
    b=max(beta*(1-0.75*(t-52)/t),0.9*para.N*(para.gamma+para.delta)/(S1+S2));
else
    b=beta;
end

%Write down the ODE system
dS1 = -b*S1*(I1+I2)/para.N - v1;
dS2 = -b*S2*(I1+I2)/para.N + para.kappa*(R+V) - v2;
dI1 = b*S1*(I1+I2)/para.N - (para.gamma+para.delta)*I1;
dI2 = b*S2*(I1+I2)/para.N - (para.gamma + para.mu*para.delta)*I2;
dR = para.gamma*(I1+I2) - para.kappa*R - v3;
dD = para.delta*(I1+para.mu*I2);
dV = v1 + v2 + v3 - para.kappa*V; 
dCV = v1 + v2 + v3;
%Define differential equation to be able to store b
dl = b;
%Reshape the derivatives into a column vector
dPop = [dS1;dS2;dI1;dI2;dR;dD;dV;dCV;dl];

end

end

