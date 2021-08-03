    % RK-4 Method 
    % These Three have to be generated thru the genetic algorit
    function [T,I,R]= model1(B1,B2,gam,sus2,infected,exposed,removed)  
        format long
        tot = 29136810;
        t0 = 1;
        h = 1;
        tEnd = size(infected,1);
        N = (tEnd-t0)/h ;
        T = [t0:h:tEnd]';
        S = zeros(N+1,1);
        E = zeros(N+1,1);
        I = zeros(N+1,1);
        R = zeros(N+1,1);
        c = zeros(N+1,1);
        k = zeros(4,1);
        l = zeros(4,1);
        m = zeros(4,1);
        n = zeros(4,1);
        S(1) = sus2(1)/tot;
        E(1) = exposed(1)/tot;
        I(1) = infected(1)/tot;

        R(1) = removed(1)/tot;

        for i = 1:N
             k1 = h* fun_s(T,B1,B2,S(i),E(i),I(i),tot);
             l1 = h* fun_e(T,B1,B2,S(i),E(i),I(i),tot);
             m1 = h* fun_i(T,B1,B2,S(i),E(i),I(i),tot,gam);
   

             k2 = h* fun_s(T+(h/2),B1,B2,S(i) + (k1/2),E(i) + (l1/2),I(i)+(m1/2),tot);
             l2 = h* fun_e(T+(h/2),B1,B2,S(i) + (k1/2),E(i) + (l1/2),I(i)+(m1/2),tot);
             m2 = h* fun_i(T+(h/2),B1,B2,S(i) + (k1/2),E(i) + (l1/2),I(i)+(m1/2),tot,gam);
        

             k3 = h* fun_s(T+(h/2),B1,B2,S(i) + (k2/2),E(i) + (l2/2),I(i)+(m2/2),tot);
             l3 = h* fun_e(T+(h/2),B1,B2,S(i) + (k2/2),E(i) + (l2/2),I(i)+(m2/2),tot);
             m3 = h* fun_i(T+(h/2),B1,B2,S(i) + (k2/2),E(i) + (l2/2),I(i)+(m2/2),tot,gam);
           

             k4 = h* fun_s(T+(h),B1,B2,S(i) + (k3),E(i) + (l3),I(i)+(m3),tot);
             l4 = h* fun_e(T+(h),B1,B2,S(i) + (k3),E(i) + (l3),I(i)+(m3),tot);
             m4 = h* fun_i(T+(h),B1,B2,S(i) + (k3),E(i) + (l3),I(i)+(m3),tot,gam);
            

             S(i+1) = S(i) + (1/6) * (k1 + 2*k2 + 2*k3 + k4);
             E(i+1) = E(i) + (1/6) * (l1 + 2*l2 + 2*l3 +l4);
             I(i+1) = I(i) + (1/6) * (m1 + 2*m2 + 2*m3 +m4);
             R(i+1) = 1 - (S(i+1) + I(i+1) + E(i+1));

        end
    
        
    end
    


