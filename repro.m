function [repo_pop,temp_popn] = repro(popn,N)
    d = (0.10)*N;
    repo_pop = popn(d+1:N,:);
    temp_popn = popn(1:d,:);
end

    