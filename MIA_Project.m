%% ========================================================================
% MIA Project - Thomas Ware

% Modeling of the Covid-19 pandemic in the Uk with different intervention 
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

%% ========================================================================
% Simple SIR model with overshoot
% Start with a simple SIR model to see overshoot and herd immunity
% Define maxime, parameters for the model and initial conditions
maxtime = 365;
para = struct('beta',2.6/12,'gamma',1/12,'N',67026292,'T0',0,'T',maxtime);
ICs = struct('S',para.N-1,'I',1,'R',0);
[Classes] = ODE_SIR_model(para,ICs);

% Calculate herd immunity threshold
herd_immunity = para.N*para.gamma/para.beta;

% Dynamics of infection with overshoot
figure(1)
clf
subplot(1,2,1)
plot(Classes.t,Classes.I)
xlabel('Time (Days)')
ylabel('Infections')

subplot(1,2,2)
p=plot(Classes.t,Classes.R,Classes.t,Classes.R(end)*ones(maxtime+1,1),'--',...
    Classes.t,para.N*(1-para.gamma*ones(maxtime+1,1)/para.beta),'--')
ylabel('Removed')
xlabel('Time (Days)')
legend([p(2),p(3)],{'Total recovered','Herd immunity threshold'},'Location','southeast')

%% ========================================================================
% Lockdown estimations
% We start by modelling varying severities of lockdowns by changing the
% value for beta. Define a new parameter q to be the quarantine severity.
% beta0 = R0*(gamma + delta)
beta0 = 2.6*1.009/12;
maxtime = 2*365;
para = struct('gamma',1/12,'delta',0.009/12,'N',67026292,'T0',1,'T',maxtime+1);
ICs = struct('S',para.N-1,'I',1,'R',0,'D',0);

q = [0 0.25 0.5 0.75];

% For each severity lockdown for full duration after 52 days
for i=1:length(q) 
beta = @(t) (t<=52)*beta0 + (t>52)*beta0*(1-q(i)*(t-52)/t);
[Classes] = ODE_SIRD_model(para,ICs,beta);
IMat1(i,:) = Classes.I;
RMat1(i,:) = Classes.R+Classes.D;
end

% Plot infections and recovered for each value of q
figure(2)
clf
subplot(1,2,1)
hold on
for i=1:length(q)
    plot(Classes.t,IMat1(i,:),'DisplayName',['q = ' num2str(q(i))])
xlabel('Time (Days)')
ylabel('Infections')
legend('show','Location','northeast')
xlim([0 450])
end

subplot(1,2,2)
hold on
for i=1:length(q)
plot(Classes.t,RMat1(i,:),'DisplayName',['q = ' num2str(q(i))])
xlabel('Time (Days)')
ylabel('Removed')
end
y1 = yline(para.N-herd_immunity,'--','Herd immunity threshold','HandleVisibility','off');
y1.LabelVerticalAlignment = 'bottom';
y1.LineWidth = 1.5;
xlim([0 450])

%% ========================================================================
% Lockdown start date
% Now show lockdown estimates with same duration taken to be 5 months (~140
% days) starting from different days
Start_day = [42,52,62,72];

for i=1:length(Start_day)
    
beta = @(t,I) (t<=Start_day(i)|t>(Start_day(i)+140))*3/7 + (1-(t<=Start_day(i)|t>(Start_day(i)+140)))*0.4*3/7;
            [Classes] = ODE_SIRD_model(para,ICs,beta);
            IMat2{i}(:,1) = Classes.I;
            RMat2{i}(:,2) = Classes.D+Classes.R;
end

% Plot infections and recoveries for each start date
figure(4)
clf
subplot(1,2,1)
hold on
for i=1:length(Start_day)
plot(Classes.t,IMat2{i}(:,1),'DisplayName',['Start day = ' num2str(Start_day(i))])
xlabel('Time (Days)')
ylabel('Infections')
legend('show')
title('Effect of changing lockdown start date')
end
xlim([0 365])

subplot(1,2,2)
hold on
for i=1:length(Start_day)
plot(Classes.t,RMat2{i}(:,2),'DisplayName',['Start Day = ' num2str(Start_day(i))])
xlabel('Time (Days)')
ylabel('Recovered')
title('')
legend('show','Location','southeast')
end
y1 = yline(para.N-herd_immunity,'--','Herd immunity threshold','HandleVisibility','off');
y1.LabelVerticalAlignment = 'bottom';
y1.LineWidth = 1.5;
xlim([0 365])
            
%% ========================================================================
% Alternative strategies
% Now start our simulations and cost analysis for seven different
% intervention stratagies. This code was run with different lockdown
% methods and the same lockdowns plus vaccination to show that vaccination
% is necessary.
% After this our seven interventions will be:
% 1 - No intervention
% 2 - No lockdowns with vaccination
% 3 - Lockdown that stops after 6 months
% 4 - Lockdown until vaccine becomes availible
% 5 - Lockdown until end of the pandemic
% 6 - Monthly periodic lockdown
% 7 - Lockdown based on R number

%Set up time horizon
NumYears = 3;
maxtime = 365*NumYears+1;
t_Yr=[1:365:maxtime];

%Set up discounting and disability weighting
dw = 0.133;
r = 0.03;
disc = (1/(1+r)).^([1:NumYears]-1);

%Set up number of iterations and timestep
iterations = 1000;

%Introduce variation to parameters to simulate unknowns
lambda = betarnd(5,15,1,iterations); 
gamma = (lambda*14+(1-lambda)*9).^-1; %~1/12 days^-1
mu = 1-betarnd(15,3,1,iterations); %~0.1
kappa = normrnd(224,20,1,iterations).^(-1); %~1/(8*30) days^-1


%Set up the number of strategies
NumStrats = 7;

%Pre-allocate costs
%Cost of lockdown is around £125 bn per year
lock_day = normrnd(125,10,1,iterations).*1e9/365;
%Vaccination in the UK cost around £11.1 bn and provided 133.5m doses 
% (267m vaccinations total)
vacc = normrnd(11.1,1,1,iterations).*1e3/133.5;

%Preallocate matrices for storing DALYs and Costs
CostMat = zeros(iterations,NumStrats);
CostMat_disc = zeros(iterations,NumStrats);
DALYMat = zeros(iterations,NumStrats);
DALYMat_disc = zeros(iterations,NumStrats);

% Now run for 3 years
for r = 1:iterations
    %Set parameters, initial conditions and non-lockdown contact rate    
    para = struct('gamma',gamma(r),'delta',0.009*gamma(r),'kappa',1/(8*30),'mu',mu(r),'N',67026292,'T0',1,'T',maxtime,'T_start',312,'T_stop',maxtime,'v_pool',75624410,'v_max',446800);
    %Fix R0 to be 2.6
    beta0 = 2.6*1.009*gamma(r);
    %Assume one initial infection
    ICs = struct('S1',para.N-1,'I1',1);
    
    for Strat = 1:NumStrats
        if Strat==1 %Do nothing
            % Define beta as a function of t
            beta = @(t) beta0;
            % Remove vaccination
            para.v_pool = 0;
            % Run model for 3 years using ICs as the initial conditions
            [Classes] = ODE_SIRDV_immunity_model(para,ICs,beta);
            % Save outputs to matricies
            SMat{Strat}(r,:) = Classes.S1 + Classes.S2;
            I1Mat{Strat}(r,:) = Classes.I1;
            I2Mat{Strat}(r,:) = Classes.I2;
            RMat{Strat}(r,:) = Classes.R;
            DMat{Strat}(r,:) = Classes.D;
            VMat{Strat}(r,:) = Classes.V;
            CVMat{Strat}(r,:) = Classes.CV;
            LRMat{Strat}(r,:) = lockdown_rate(beta,beta0,maxtime);
            
        elseif Strat==2 %No lockdowns with vaccination
            beta = @(t) beta0;
            %Add vaccination
            para.v_pool=151248820;
            [Classes] = ODE_SIRDV_immunity_model(para,ICs,beta);
            % Save outputs to matricies
            SMat{Strat}(r,:) = Classes.S1 + Classes.S2;
            I1Mat{Strat}(r,:) = Classes.I1;
            I2Mat{Strat}(r,:) = Classes.I2;
            RMat{Strat}(r,:) = Classes.R;
            DMat{Strat}(r,:) = Classes.D;
            VMat{Strat}(r,:) = Classes.V;
            CVMat{Strat}(r,:) = Classes.CV;
            LRMat{Strat}(r,:) = lockdown_rate(beta,beta0,maxtime);
            
        elseif Strat==3 %Lockdown for 6 months
            % Lockdown from 52 to 220  then after beta increases linearly to original value in two months  
            beta = @(t) (t<52|t>276)*beta0 + (t>=52&&t<220)*beta0*(1-0.75*(t-52)/t) + ...
                    (t>=220&&t<=276)*(beta0*(1-0.75*(168/220))+ beta0*0.75*(168/220)/56*(t-220)); 
            para.v_pool=151248820;
            [Classes] = ODE_SIRDV_immunity_model(para,ICs,beta);
            SMat{Strat}(r,:) = Classes.S1 + Classes.S2;
            I1Mat{Strat}(r,:) = Classes.I1;
            I2Mat{Strat}(r,:) = Classes.I2;
            RMat{Strat}(r,:) = Classes.R;
            DMat{Strat}(r,:) = Classes.D;
            VMat{Strat}(r,:) = Classes.V;
            CVMat{Strat}(r,:) = Classes.CV;
            LRMat{Strat}(r,:) = lockdown_rate(beta,beta0,maxtime);
              
       elseif Strat==4 % Lockdown until vaccine becomes availible
           % Similarly lockdown unitl vaccine then increase beta linearly to original value in two months 
            beta = @(t) (t<52|t>368)*beta0 + (t>=52&&t<312)*beta0*(1-0.75*(t-52)/t) + ...
                    (t>=312&&t<=368)*(beta0*(1-0.75*(260/312))+ beta0*0.75*(260/312)/56*(t-312));
            para.v_pool=151248820;
            % Run for 3 years
            [Classes] = ODE_SIRDV_immunity_model(para,ICs,beta);
            % Save outputs to matrices
            SMat{Strat}(r,:) = Classes.S1 + Classes.S2;
            I1Mat{Strat}(r,:) = Classes.I1;
            I2Mat{Strat}(r,:) = Classes.I2;
            RMat{Strat}(r,:) = Classes.R;
            DMat{Strat}(r,:) = Classes.D; 
            VMat{Strat}(r,:) = Classes.V;
            CVMat{Strat}(r,:) = Classes.CV;
            LRMat{Strat}(r,:) = lockdown_rate(beta,beta0,maxtime);
            
        elseif Strat==5 %Lockdown until end of pandemic
            %Lockdown constantly
            beta = @(t) (t<52)*beta0 + (t>=52)*beta0*(1-0.75*(t-52)/t);
            para.v_pool=151248820;
            % Run for 3 years
            [Classes] = ODE_SIRDV_immunity_model(para,ICs,beta);
            % Save outputs to matrices
            SMat{Strat}(r,:) = Classes.S1 + Classes.S2;
            I1Mat{Strat}(r,:) = Classes.I1;
            I2Mat{Strat}(r,:) = Classes.I2;
            RMat{Strat}(r,:) = Classes.R;
            DMat{Strat}(r,:) = Classes.D;
            VMat{Strat}(r,:) = Classes.V;
            CVMat{Strat}(r,:) = Classes.CV;
            LRMat{Strat}(r,:) = lockdown_rate(beta,beta0,maxtime);

        elseif Strat==6 %Monthly lockdown
            % Lockdown three months out of every four.
            beta = @(t) monthly_periodic(beta0,t)
            para.v_pool=151248820;
            % Run for 3 years
            [Classes] = ODE_SIRDV_immunity_model(para,ICs,beta);
            % Save outputs to matrices
            SMat{Strat}(r,:) = Classes.S1 + Classes.S2;
            I1Mat{Strat}(r,:) = Classes.I1;
            I2Mat{Strat}(r,:) = Classes.I2;
            RMat{Strat}(r,:) = Classes.R;
            DMat{Strat}(r,:) = Classes.D;  
            VMat{Strat}(r,:) = Classes.V;
            CVMat{Strat}(r,:) = Classes.CV;
            LRMat{Strat}(r,:) = lockdown_rate(beta,beta0,maxtime);
            
        elseif Strat==7 %Lockdown to keep R<0.9.
            para.v_pool=151248820;
            % Run for 3 years
            [Classes] = ODE_SIRDV_immunity_model_2(para,ICs,beta0);
            % Save outputs to matrices
            SMat{Strat}(r,:) = Classes.S1 + Classes.S2;
            I1Mat{Strat}(r,:) = Classes.I1;
            I2Mat{Strat}(r,:) = Classes.I2;
            RMat{Strat}(r,:) = Classes.R;
            DMat{Strat}(r,:) = Classes.D; 
            VMat{Strat}(r,:) = Classes.V;
            CVMat{Strat}(r,:) = Classes.CV;
            % Calculate estimate for beta
            beta = [Classes.l(2:end)-Classes.l(1:end-1);beta0];
            LRMat{Strat}(r,:) = lockdown_rate(beta,beta0,maxtime);

            elseif Strat==7 %Lockdown to keep R=0.9.
            para.v_pool=151248820;
            % Run for 3 years
            [Classes] = ODE_SIRDV_immunity_model_2(para,ICs,beta0);
            % Save outputs to matrices
            SMat{Strat}(r,:) = Classes.S1 + Classes.S2;
            I1Mat{Strat}(r,:) = Classes.I1;
            I2Mat{Strat}(r,:) = Classes.I2;
            RMat{Strat}(r,:) = Classes.R;
            DMat{Strat}(r,:) = Classes.D; 
            VMat{Strat}(r,:) = Classes.V;
            CVMat{Strat}(r,:) = Classes.CV;
            % Calculate estimate for beta
            beta = [Classes.l(2:end)-Classes.l(1:end-1);beta0];
            LRMat{Strat}(r,:) = lockdown_rate(beta,beta0,maxtime);
        end
    end
end
%%
%Plot on infection dynamics for all strategies
figure(5)
clf
hold on
for Strat = 1:NumStrats
plot(Classes.t,I1Mat{Strat}(1,:)+I2Mat{Strat}(1,:))
end
legend('No intervention','Vaccination','6 month lockdown','Lockdown until vaccine','Full duration','Month based','R number')
xlabel('Time (Days)')
ylabel('Infections')
xlim([0 1000])
x1 = xline(312,'--k','Vaccination start','HandleVisibility','off')
x1.LabelVerticalAlignment = 'top';
x1.LineWidth = 1.5;

%Plot of removed with includes recovered and vaccinated
figure(6)
clf
hold on
for Strat = 1:NumStrats
plot(Classes.t,RMat{Strat}(1,:)+VMat{Strat}(1,:))
end
legend('No intervention','Vaccination','6 month lockdown','Lockdown until vaccine','Full duration','Month based','R number','location','southeast')
xlabel('Time (Days)')
ylabel('Removed')
xlim([0 700])
x1 = xline(312,'--k','Vaccination start','HandleVisibility','off')
x1.LabelVerticalAlignment = 'bottom';
x1.LineWidth = 1.5;

%Plot for deaths from each
figure(7)
clf
hold on
for Strat = 1:NumStrats
plot(Classes.t,DMat{Strat}(1,:))
end
legend('No intervention','Vaccination','6 month lockdown','Lockdown until vaccine','Full duration','Month based','R number','location','southeast')
xlabel('Time (Days)')
ylabel('Deaths')
xlim([0 700])
x1 = xline(312,'--k','Vaccination start','HandleVisibility','off')
x1.LabelVerticalAlignment = 'bottom';
x1.LineWidth = 1.5;

% Calculate durations and deaths of each method
for r=1:iterations
    for Strat=1:NumStrats
        dur = [maxtime find((I1Mat{Strat}(r,:)+I2Mat{Strat}(r,:))<1,1,'first')];
        Duration(Strat,r) = min(dur);
        Deaths(Strat,r) = DMat{Strat}(r,Duration(Strat));
        [y,ix] = max(I1Mat{Strat}(r,:)+I2Mat{Strat}(r,:));
        Peak_time(Strat,r) = ix;
        Peak_inf(Strat,r) = y;
    end
end

%Plot as bar charts
figure(8)
clf
subplot(1,2,1)
b=bar(mean(Deaths,2)/1000);
title('Expected deaths (thousands)')
b.FaceColor='Flat';
b.CData=CMap;
xticklabels({''})

subplot(1,2,2)
b=bar(mean(Duration,2));
title('Expected duration (days)')
b.FaceColor='Flat';
b.CData=CMap;
xticklabels({''})

%% ========================================================================
% Cost analysis
% Compute costs and DALYs for each method
for r = 1:iterations
    for Strat = 1:NumStrats
        for i=1:length(t_Yr)-1
            %YLL
            %Calculate deaths per year and then multiply by 11.1607
            YLL(i) = 11.1607*(DMat{Strat}(r,t_Yr(i+1))-DMat{Strat}(r,t_Yr(i)));
            %YLD
            %Person years infected annualy is the integral over time of the
            %number of people infected  
            high_person_years_annual(i) = trapz(Classes.t([t_Yr(i):t_Yr(i+1)]), I1Mat{Strat}(r,[t_Yr(i):t_Yr(i+1)]))/365;
            low_person_years_annual(i) = trapz(Classes.t([t_Yr(i):t_Yr(i+1)]), I2Mat{Strat}(r,[t_Yr(i):t_Yr(i+1)]))/365;
            YLD(i) = dw*high_person_years_annual(i) + 0.2*dw*low_person_years_annual(i);
            %Costs
            %Vaccine is cost of vaccination per person plus a £600m cost in
            %first year for finding a vaccine
            Cost_vac(i) = vacc(r)*(CVMat{Strat}(r,t_Yr(i+1)) - CVMat{Strat}(r,t_Yr(i))) + isequal(i,1)*600e6;
            %Cost of lockdown is the cost of each day in lockdown
            %multiplied by the lockdown rate
            Cost_lock(i) = lock_day(r)*sum(LRMat{Strat}(r,[t_Yr(i):min(t_Yr(i+1),Duration(Strat))]));
            
        end
        % Calculate DALYs which is disability weighting x person years
        DALYs = YLL + YLD;
        % Save to matricies
        DALYMat(r,Strat) = sum(DALYs,2);
        % Compute discounted DALYs
        DALYMat_disc(r,Strat) = sum(DALYs.*disc,2);
        
        % CostMat is the cost over the 3 years without discounting
        CostMat(r,Strat) = sum(Cost_vac+Cost_lock);
        % CostMat_disc is costs with discounting
        CostMat_disc(r,Strat) = sum((Cost_vac+Cost_lock).*disc,2);
    end
end

%% ========================================================================
% Calculate mean DALYs and costs for each method
for Strat = 1:NumStrats
    Mean_DALYs(Strat) = mean(DALYMat(:,Strat));
    Mean_Cost(Strat) = mean(CostMat(:,Strat));
end

% Now create a plot of undiscounted DALYs and costs
% Start by calculating delta costs and DALYs
DCostMat = CostMat - repmat(CostMat(:,1),1,NumStrats);
DDALYMat = repmat(DALYMat(:,1),1,NumStrats) - DALYMat;

%Means using discounted
MeanDCost = mean(DCostMat,1);
MeanDDALY = mean(DDALYMat,1);

%Once calculated we reorder. No methods are dominated.
figure(9)
clf
hold on
%Plot outcomes as a scatter (with transparency)
for Strat = [1 2 3 6 7 4 5]
    scatter(DDALYMat(:,Strat)./1e3,DCostMat(:,Strat)./1e9,[],CMap(Strat,:),'filled','markerfacealpha',0.5)
%Mark on means
    hs(Strat)=scatter(MeanDDALY(Strat)/1e3,MeanDCost(Strat)/1e9,[],CMap(Strat,:),'filled');
end
%Draw on ICERs
plot(MeanDDALY([1 2 3 6 7 4 5])./1e3,MeanDCost([1 2 3 6 7 4 5])./1e9,'-k')
legend(hs,'No intervention','Vaccination','6 month lockdown','Lockdown until vaccine','Full duration','Month based','R number','location','southeast')
xlabel('Undiscounted DALYs averted (thousands)')
ylabel('Undiscounted additional costs ($\pounds$ B)')

%Compute Delta Costs/ Delta DALYs and compute NMBs with discounting
%Discounted
DCostMat_disc = CostMat_disc - repmat(CostMat_disc(:,1),1,NumStrats);
DDALYMat_disc = repmat(DALYMat_disc(:,1),1,NumStrats) - DALYMat_disc;

%Means using discounted
MeanDCost_disc = mean(DCostMat_disc);
MeanDDALY_disc = mean(DDALYMat_disc);

%ICERs
ICER(1)=0;
ICER(2)=MeanDCost_disc(2)/MeanDDALY_disc(2);
ICER(3)=(MeanDCost_disc(3)-MeanDCost_disc(2))/(MeanDDALY_disc(3)-MeanDDALY_disc(2));
ICER(4)=(MeanDCost_disc(6)-MeanDCost_disc(3))/(MeanDDALY_disc(6)-MeanDDALY_disc(3));
ICER(5)=(MeanDCost_disc(7)-MeanDCost_disc(6))/(MeanDDALY_disc(7)-MeanDDALY_disc(6));
ICER(6)=(MeanDCost_disc(4)-MeanDCost_disc(7))/(MeanDDALY_disc(4)-MeanDDALY_disc(7));
ICER(7)=(MeanDCost_disc(5)-MeanDCost_disc(4))/(MeanDDALY_disc(5)-MeanDDALY_disc(4));

ICERtable1 = array2table([MeanDCost_disc' MeanDDALY_disc' ICER']);
ICERtable1.Properties.VariableNames={'Delta Costs','Delta DALYs','ICER'};

%% ========================================================================
%Plot CE plane with discounting
figure(10)
clf
hold on
%Plot outcomes as a scatter (with transparency)
for Strat = 1:NumStrats
    scatter(DDALYMat_disc(:,Strat)./1e3,DCostMat_disc(:,Strat)./1e9,[],CMap(Strat,:),'filled','markerfacealpha',0.5)
%Mark on means
    hs(Strat)=scatter(MeanDDALY_disc(Strat)/1e3,MeanDCost_disc(Strat)/1e9,[],CMap(Strat,:),'filled');
end
%Draw on ICERs
plot(MeanDDALY_disc([1 2 3 6 7 4 5])./1e3,MeanDCost_disc([1 2 3 6 7 4 5])./1e9,'-k')

legend(hs,'No intervention','Vaccination','6 month lockdown','Lockdown until vaccine','Full duration','Month based','R number','location','southeast')
xlabel('DALYs averted (thousands)')
ylabel('Additional costs ($\pounds$ B)')

%NMB for $1000 increments
%N.B. using original order of strategies

w=0:1:1000; %WTP thresholds
for j=1:length(w)
    for Strat=1:NumStrats
        NMB(:,(Strat-1)*length(w)+j) = w(j)*DDALYMat_disc(:,Strat) - DCostMat_disc(:,Strat)/1e3;
    end
end

%Check if each strategy is optimal (lowest NBM for each WTP) per replicate
Optimal = zeros(iterations,length(w)*NumStrats);
for j=1:length(w)
    for r= 1:iterations
        [im,iy] = max(NMB(r,j:length(w):Strat*length(w)));
        Optimal(r,j+ (iy-1)*length(w)) = 1;
    end
end

%Compute the probability that each strategy is CE
for j=1:length(w)
    ProbCE(:,j) = (sum(Optimal(:,j:length(w):length(w)*NumStrats))/iterations)';
end

%Plot CEACs
figure(11)
clf
set(gca, 'ColorOrder', CMap, 'NextPlot', 'replacechildren');
%Plot each CEAC
p=plot(w,ProbCE);
ylabel('Probability cost-effective')
xlabel('Willingness to pay per DALY averted ($\pounds$ thousands)')
xlim([0 1000])
ylim([0 1])

legend(p,'No intervention','Vaccination','6 month lockdown','Lockdown until vaccine','Full duration','Month based','R number','location','NorthEast')
