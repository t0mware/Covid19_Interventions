%% ========================================================================
% MIA Project - Thomas Ware

% Modeling of the Covid-19 pandemic in the UK with different intervention 
% stratagies and policy reccomendations.

%% ========================================================================
% Preamble

% Cleaning and Settings
clc
clear all
close all

%Figure Settings
set(0,'DefaultFigureColor',[1 1 1])
set(0, 'defaultaxesfontsize', 15)
set(0, 'defaultlinelinewidth', 1.5)

%Latex Font
set(0,'DefaultTextInterpreter', 'latex')
set(0,'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex')

%Define colour map
CMap = [0 0.4470 0.7410; 0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;...
    0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.3010 0.7450 0.9330;0.6350 0.0780 0.1840];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter Estimating

% We start by estimating the parameters for our models using methods and
% data from the UK government. Look at non-lockdown periods for base
% estimates.
% We will use likelihood functions and a MCMC approach

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simple SIR model with overshoot

% Start with a simple SIR model to see overshoot and  herd immunity. Define
% maxtime, the parameters for the ODE system and initial conditions
maxtime = 200;
para = struct('beta',1/3.5,'gamma',1/7,'N',67330000,'T0',0,'T',maxtime);
ICs = struct('S',para.N-1,'I',1,'R',0);
[Classes] = ODE_SIR_model(para,ICs);

herd_immunity = para.N*para.gamma/para.beta;

% Dynamics with overshoot
figure(1)
clf
subplot(1,2,1)
plot(Classes.t,Classes.I)
xlabel('Time (Days)')
ylabel('Infections')

subplot(1,2,2)
p=plot(Classes.t,Classes.R,Classes.t,Classes.R(end)*ones(maxtime+1,1),'--',...
    Classes.t,para.N*(1-para.gamma*ones(maxtime+1,1)/para.beta),'--')
ylabel('Recovered')
xlabel('Time (Days)')
legend([p(2),p(3)],{'Total recovered','Herd immunity threshold'},'Location','southeast')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lockdown estimations

% We start by modelling varying severities of lockdowns by changing the
% value for beta. Define a new parameter lt to be the lockdown tightness
% to give beta_t = (1-lt)*beta.

maxtime = 2*365;
para = struct('gamma',1/7,'p_d',0.0075,'N',67330000,'T0',0,'T',maxtime);
ICs = struct('S',para.N-1,'I',1,'R',0,'D',0);

lt = [0.1 0.2 0.3 0.4];

for i=1:length(lt) 
beta = @(t,I) (t<=52)/3.5 + (t>52)*(1-lt(i))/3.5;
[Classes] = ODE_SIRD_model(para,ICs,beta);
Mat{i}(:,1) = Classes.I;
Mat{i}(:,2) = Classes.R+Classes.D;
end

% Infections
figure(2)
clf
hold on
for i=1:length(lt)
plot(Classes.t,Mat{i}(:,1),'DisplayName',['lt = ' num2str(lt(i))])
xlabel('Time (Days)')
ylabel('Infections')
legend('show')
title('Effect of changing lockdown tightness')
end
xlim([0,500])

% Recoveries
figure(3)
clf
hold on
for i=1:length(lt)
plot(Classes.t,Mat{i}(:,2),'DisplayName',['lt = ' num2str(lt(i))])
xlabel('Time (Days)')
ylabel('Recovered')
title('')
legend('show','Location','southeast')
end
xlim([0 500])
y1 = yline(para.N-herd_immunity,'--','Herd immunity threshold','HandleVisibility','off');
y1.LabelVerticalAlignment = 'bottom';
y1.LineWidth = 1.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start day

% Now show lockdown estimates with same duration taken to be 5 months (~168
% days) starting from different days
Start_day = [52,80,113,140];

for i=1:length(Start_day)
    
beta = @(t,I) (t<=Start_day(i)|t>(Start_day(i)+140))/3.5 + (1-(t<=Start_day(i)|t>(Start_day(i)+140)))*0.6/3.5;
            [Classes] = ODE_SIRD_model(para,ICs,beta);
            Mat2{i}(:,1) = Classes.I;
            Mat2{i}(:,2) = Classes.D+Classes.R;
end

%Infections
figure(4)
clf
hold on
for i=1:length(Start_day)
plot(Classes.t,Mat2{i}(:,1),'DisplayName',['Start day = ' num2str(Start_day(i))])
xlabel('Time (Days)')
ylabel('Infections')
legend('show')
title('Effect of changing lockdown start date')
end
xlim([0 400])

% Recoveries
figure(5)
clf
hold on
for i=1:length(Start_day)
plot(Classes.t,Mat2{i}(:,2),'DisplayName',['Start Day = ' num2str(Start_day(i))])
xlabel('Time (Days)')
ylabel('Recovered')
title('')
legend('show','Location','southeast')
end
y1 = yline(para.N-herd_immunity,'--','Herd immunity threshold','HandleVisibility','off');
y1.LabelVerticalAlignment = 'bottom';
y1.LineWidth = 1.5;
xlim([0 400])
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now start our cost analysis for six different intervention stratagies.
% These six interventions will be:
% 1 - No intervention
% 2 - Lockdown that stops after 5 months
% 3 - Lockdown until end of pandemic
% 4 - Periodic lockdown
% 5 - ICU capacity lockdown
% 6 - Partial vaccination
% 7 - Full vacciantion

% Set up time horizon
NumYears = 3;
maxtime = 365*NumYears;
t_Yr=[1:365:maxtime];

% Set up discounting and disability weighting
high_dw = 0.15;
mid_dw = 0.01;
r = 0.03;
disc = (1/(1+r)).^([1:NumYears]-1);

% Introduce variation to parameters to simulate unknowns

% Pre-compute costs of the vaccine and lockdown
% Vacc_cost = gamrnd(5,3/5,1000,1);
% Lock_cost =

% Set up number of iterations and timestep
iterations = 1;
% Gives start point (index) of each year
t_Yr=[1:365:NumYears*365+1];

% Set up the number of strategies
NumStrats = 7;

% Preallocate matrices for storing DALYs and Costs
CostMat = zeros(iterations,NumStrats);
CostMat_disc = zeros(iterations,NumStrats);
DALYMat = zeros(iterations,NumStrats);
DALYMat_disc = zeros(iterations,NumStrats);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now run for 3 years
para = struct('gamma',1/7,'p_d',0.0075,'N',67330000,'T0',0,'T',maxtime,'T_start',312,'T_stop',maxtime,'v_pool',50000000,'v_max',133135,'v_eff',0.85);
ICs = struct('S',para.N-1,'I',1,'R',0,'D',0);

for r = 1:iterations
    for Strat = 1:NumStrats
        if Strat==1 %Do nothing
            % Define beta as a function of t and possibly I
            beta = @(t,I) 1/3.5;
            % Run model for 3 years using ICs as the initial conditions
            [Classes] = ODE_SIRD_model(para,ICs,beta);
            % Save outputs to matricies
            SMat{Strat}(r,:) = Classes.S;
            IMat{Strat}(r,:) = Classes.I;
            RMat{Strat}(r,:) = Classes.R;
            DMat{Strat}(r,:) = Classes.D;
            
        elseif Strat==2 %Lockdown for 5 months
            % Lockdown from 116 to 256    
            beta = @(t,I) (t<=116|t>256)/3.5 + (1-(t<=116|t>256))*0.6/3.5;
            [Classes] = ODE_SIRD_model(para,ICs,beta);
            SMat{Strat}(r,:) = Classes.S;
            IMat{Strat}(r,:) = Classes.I;
            RMat{Strat}(r,:) = Classes.R;
            DMat{Strat}(r,:) = Classes.D;
  
        elseif Strat==3 %Lockdown until end of pandemic for herd immunity
            % After 52 days, change beta to 0.7*beta
            beta = @(t,I) ((t<=52)/3.5 + (t>52)*0.7/3.5);
            % Run for 3 years
            [Classes] = ODE_SIRD_model(para,ICs,beta);
            % Save outputs to matrices
            SMat{Strat}(r,:) = Classes.S;
            IMat{Strat}(r,:) = Classes.I;
            RMat{Strat}(r,:) = Classes.R;
            DMat{Strat}(r,:) = Classes.D;

        elseif Strat==4 %Weekly periodic Lockdown
            %  Lockdown three weeks out of every four. Define beta by:
            beta = @(t,I) (0.6+0.4*((t<=52)|((t>52)&&(mod(t-52,28)>=21)&&(mod(t-52,28)<=27))))/3.5;
            % Run for 3 years
            [Classes] = ODE_SIRD_model(para,ICs,beta);
            % Save outputs to matrices
            SMat{Strat}(r,:) = Classes.S;
            IMat{Strat}(r,:) = Classes.I;
            RMat{Strat}(r,:) = Classes.R;
            DMat{Strat}(r,:) = Classes.D;       
            
        elseif Strat==5 % Monthly periodic lockdown
            % Lockdown three months out of every four.
            beta = @(t,I) (0.6+0.4*((t<=52)|((t>52)&&(mod(t-52,112)>=84)&&(mod(t-52,112)<=111))))/3.5;
            % Run for 3 years
            [Classes] = ODE_SIRD_model(para,ICs,beta);
            % Save outputs to matrices
            SMat{Strat}(r,:) = Classes.S;
            IMat{Strat}(r,:) = Classes.I;
            RMat{Strat}(r,:) = Classes.R;
            DMat{Strat}(r,:) = Classes.D;  
                
        elseif Strat==6 %Full Vaccination
            % First vaccine was given in the UK on the 8th December 2020
            % which accounts to day 312 in our model. Full vaccination is
            % giving people 2 doses to get a high level of protection.
            % Define beta
            beta = @(t) (t<=52|t>464)/3.5 + (1-(t<=52|t>464))*0.6/3.5;
            % Run for 3 years
            [Classes] = ODE_SIRDV_model(para,ICs,beta)
            % Save outputs to matrices
            SMat{Strat}(r,:) = Classes.S;
            IMat{Strat}(r,:) = Classes.I;
            RMat{Strat}(r,:) = Classes.R;
            DMat{Strat}(r,:) = Classes.D; 
            VMat{Strat}(r,:) = Classes.V;
            CVMat{Strat}(r,:) = Classes.CV;
            
        elseif Strat==7 %Partial Vaccination
            % Partial vaccination equates to giving one dose but with a 
            % lower level of protection.
            para.v_eff=0.65;
            beta = @(t) (t<=52|t>502)/3.5 + (1-(t<=52|t>502))*0.6/3.5;
            % Run for 3 years
            [Classes] = ODE_SIRDV_model(para,ICs,beta)
            % Save outputs to matrices
            SMat{Strat}(r,:) = Classes.S;
            IMat{Strat}(r,:) = Classes.I;
            RMat{Strat}(r,:) = Classes.R;
            DMat{Strat}(r,:) = Classes.D; 
            VMat{Strat}(r,:) = Classes.V;
            CVMat{Strat}(r,:) = Classes.CV;
              
        end
    end
end

% Infections
figure(6)
clf
hold on
for Strat = 1:7
plot(Classes.t,IMat{Strat}(1,:))
end
legend('No intervention','5 month lockdown','Full duration Lockdown','Week based lockdown','Month based lockdown','Full vaccination','Partial vaccination')
xlabel('Time (Days)')
ylabel('Infections')
xlim([0 500])
ylim([0 11000000])

% Recoveries and vaccinations
figure(7)
clf
hold on
for Strat = 1:5
plot(Classes.t,RMat{Strat}(1,:))
end
plot(Classes.t,RMat{6}(1,:)+VMat{6}(1,:))
plot(Classes.t,RMat{7}(1,:)+VMat{7}(1,:))
xlim([0 500])
legend('No intervention','5 month lockdown','Full duration Lockdown','Week based lockdown','Month based lockdown','Full vaccination','Partial vaccination')
xlabel('Time (Days)')
ylabel('Removed')
xlim([0 600])

% Deaths
figure(8)
clf
hold on
for Strat = 1:7
plot(Classes.t,DMat{Strat}(1,:))
end
legend('No intervention','5 month lockdown','Full duration Lockdown','Week based lockdown','Month based lockdown','Full vaccination','Partial vaccination')
xlabel('Time (Days)')
ylabel('Deaths')
xlim([0 600])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate durations and deaths of each method

for Strat=1:NumStrats
    Duration(Strat) = find(IMat{Strat}(1,:)<1,1,'first')
    Deaths(Strat) = DMat{Strat}(1,Duration(Strat))
end

% Plot as bar charts
figure(9)
clf
subplot(2,1,1)
b=bar(Deaths/1000);
title('Expected deaths (thousands)')
b.FaceColor='Flat';
b.CData=CMap;
xticklabels({''})

subplot(2,1,2)
b=bar(Duration);
title('Expected duration (days)')
b.FaceColor='Flat';
b.CData=CMap;
xticklabels({''})
