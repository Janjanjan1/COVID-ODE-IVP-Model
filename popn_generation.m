function popn = popn_generation(N)
    % Making the new Population 

    popn = zeros(N,4);

    for i = 1:N
        for j = 2:4
            popn(i,j) = rand();
        end
    end


