function lr = lockdown_rate(beta,beta0,maxtime)
     l = zeros(1,maxtime);
    for i=1:maxtime
        l(i) = 1-beta(i)/beta0;
    end
    lr =l;
end