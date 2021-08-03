function repo_pop = secs(repo_pop)
    format long 

    temp_pop = zeros(size(repo_pop,1),4);
    temp_prob = zeros(); 
    som = sum(1:size(repo_pop,1));
    for i  = 1:size(repo_pop,1)
        temp_prob(i) = ((size(repo_pop,1)+1)-i)/som;
    end
   
    z = 1;
    
   while z ~= size(repo_pop,1)
        temp_rand = rand();
        temp_cp = [0,cumsum(temp_prob)];
        temp_ind = find(temp_rand>temp_cp, 1, 'last');
        p1 = repo_pop(temp_ind,:);

        temp_rand = rand();
        temp_cp = [0,cumsum(temp_prob)];
        temp_ind = find(temp_rand>temp_cp, 1, 'last');
        p2 = repo_pop(temp_ind,:);
        
        % Crossing Over
        temp_rand = rand();
        of1 = zeros(1,4);
        of2 = zeros(1,4);
        if temp_rand <= 0.9
            for i = 2:4
                temp_rand = rand();
                of1(i) = abs(p1(1,i) * (1 - temp_rand) - p2(1,i) * (temp_rand));
                of2(i) = abs(p2(1,i) * (1 - temp_rand) - p1(1,i) * (temp_rand));
            end
            of1 = p1;
            of2 = p2;
        end
            temp_rand = rand();
            if temp_rand <= 0.9
                temp_rand = randi([2,4],1);
                of1(1,temp_rand) = rand();     
                temp_pop(z,:) = of1;
                z = z + 1;
            else
                temp_pop(z,:) = of1;
                z = z+1 ;
            end
   
            temp_rand = rand();
            if temp_rand <= 0.9
                temp_rand = randi([2,4],1);
                of2(1,temp_rand) = rand();          
                temp_pop(z,:) = of2;
            else
                temp_pop(z,:) = of2;
            end
   end
    clear repo_pop;
    repo_pop = temp_pop;
    
end