%ODE SIRD model code

function [Classes] = ODE_SIRD_model(para,ICs,beta)


%Run ODE using ODE45
opts = odeset('RelTol',1e-5);
%Only give outputs at defined timepoints
[t, pop] = ode45(@diff_SIRD_model, [para.T0:1:para.T], [ICs.S ICs.I ICs.R ICs.D],opts,para,beta);

%Convert output to structure
Classes = struct('S',pop(:,1),'I',pop(:,2),'R',pop(:,3),'D',pop(:,4),'t',t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Diff equations

function dPop = diff_SIRD_model(t,pop,para,beta)

%Assign the population matrix into the classes
S=pop(1);
I=pop(2);
R=pop(3);
D=pop(4);

%Define beta
b=beta(t,I);

%Write down the ODE system
dS = -b*S*I/para.N;
dI = b*S*I/para.N - para.gamma*I;
dR = (1-para.p_d)*para.gamma*I;
dD = para.p_d*para.gamma*I;

%Reshape the derivatives into a column vector
dPop = [dS;dI;dR;dD];
end

end

