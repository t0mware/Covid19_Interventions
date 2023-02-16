%ODE SIRDV model code

function [Classes] = ODE_SIRDV_model(para,ICs,beta)

%Run ODE using ODE45
opts = odeset('RelTol',1e-5);
%Only give outputs at defined timepoints
[t, pop] = ode45(@diff_SIRDVmodel, [para.T0:1:para.T], [ICs.S ICs.I ICs.R ICs.D 0 0 0 ], opts,para,beta);

%Convert output to structure
Classes = struct('S',pop(:,1),'I',pop(:,2),'R',pop(:,3),'D',pop(:,4),'V',pop(:,5),'UV',pop(:,6),'CV',pop(:,7),'t',t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Diff equations

function dPop = diff_SIRDVmodel(t,pop,para,beta)

%Assign the population matrix into the classes
S=pop(1);
I=pop(2);
R=pop(3);
D=pop(4);
V=pop(5);       %Vaccinated successfully
UV=pop(6);      %Vaccinated unsuccessfully
CV=pop(7);      %Count total vaccinations given

%Vaccines
if t>para.T_start && t<para.T_stop
    %Check whether vaccine pool is used up
    if CV<para.v_pool
        if S>para.v_max 
            v=para.v_max;
        else 
            v=S;
        end
    else
        v=0;
    end
    
else
    v=0;
end

% Define beta
b = beta(t);

%Write down the ODE system
dS = -b*S*I/para.N - v;
dI = b*(S+UV)*I/para.N - para.gamma*I;
dR = (1-para.p_d)*para.gamma*I;
dD = para.p_d*para.gamma*I;

%Successfully vaccinated
dV = v*para.v_eff;          
%Unsuccessful vaccines given (will not be vaccinated again but could be infected)
dUV = (1-para.v_eff)*v - b*UV*I/para.N; 
%Cumulative count of vaccinated
dCV = v; 

%Reshape the derivatives into a column vector
dPop = [dS;dI;dR;dD;dV;dUV;dCV];

end

end
