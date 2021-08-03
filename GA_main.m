format long
rng shuffle
N = input('Enter Popn Size');
tot = sus2(1)+exposed(1)+infected(1)+removed(1);
% gov = gov/100;
% Popn Generation
popn = popn_generation(N);
min_fitness = zeros();
status = 0;
gen = 0;
fit_pop = zeros(N,4);

while status == 0
   gen = gen + 1;
    %Slicing Popn Matrix for B1,B2,gama
        tic;
        parfor p= 1:N
  
            B1 = popn(p,2);
            B2 = popn(p,3);
            gam= popn(p,4);


            % running the runge-kutta algo to make S,E,I,R based on those 
            [~,I,R]= model_6th(B1,B2,gam,sus2,infected,exposed,removed,tot);

            % running the fitness algorithm:

            % writing fitness value 
            fit_pop(p,1) = fit_check(I,R,removed,infected,tot);
            
        end
        toc
        popn(:,1) = fit_pop(:,1);
  
    
  
    if gen>1000
        status = 1;
        break;
    end
    
    popn = sortrows(popn);

    if popn(1,1) == 0
        B1 = popn(1,2);
        B2 = popn(1,3);
        gam= popn(1,4);
        [~,I,R] = model_6th(B1,B2,gam,sus2,infected,exposed,removed,tot);
        popn(1,1) = fit_check(I,R,removed,infected,tot);
    end
    
    % Recording the Average Fitness Of Population
    min_fitness(1,gen) = popn(1,    1);

    semilogy(min_fitness(1,:),'r'); drawnow

    % Sorting the popn:
    
   
    
    % Selecting the btm 50% and halving the matrix:
    [repo_pop,temp_popn] = repro(popn,N);
    clear popn;
   
    
    %Parents, Crossing Over and Replacing:
    repo_pop = secs(repo_pop);
   

    % Adding Reproduced Popn In To the Bigger Population:
  
    
    popn = [temp_popn;repo_pop];

    clear temp_popn;
    clear repo_pop;
  
    
    
    
end

B1 = popn(1,2);
B2 = popn(1,3);
gam= popn(1,4);
[T,I,R] = model_6th(B1,B2,gam,sus2,infected,exposed,removed,tot);
I = I.*tot;
R = R.*tot;



