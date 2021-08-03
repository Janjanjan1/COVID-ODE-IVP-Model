function fitness = fit_check(I,R,removed,infected,tot)
   
    I = I .* tot;
    R = R .* tot;

    mean_I = mean((infected - I).^ 2);
    mean_R = mean((removed - R).^2);
    f1 = 1.5 *(mean_I + 10*((infected(1)-I(1))^2 + ((infected(end) - I(end))^2)));
    f2 = 0.5 *(mean_R + 10*(((removed(1)-R(1))^2)+((removed(end)-R(end))^2)));
    fitness = f1+f2;
   
 
       
   
    
 
end
