function beta = monthly_periodic(beta0,t)
if t<52
    beta = beta0;
elseif t>=52 && (mod(t-52,112)>=84) && (mod(t-52,112)<=111)
    %Lockdown at given rate
    beta = beta0*(1-0.75*84/136)+ (beta0-beta0*(1-0.75*84/136))/(111-84)*(t-find(mod([1:1:t]-52,112)==84,1,'last'));    
else
    %Increase intact rate linearly in month that follows
    beta = beta0*(1-0.75*(t - find(mod([1:1:t]-52,112)==111,1,'last'))/(t+52-find(mod([1:1:t]-52,112)==111,1,'last')));
end
