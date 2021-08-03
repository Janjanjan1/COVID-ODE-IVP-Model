function [Time, I , R] = model_6th_STR(B1,B2,gam,sus2,infected,exposed,removed,tot,gov)
    format long
        t0 = 1;
        h = 1;
        tEnd = size(infected,1);
        N = (tEnd-t0)/h ;
        Time = (t0:h:tEnd);
        S = zeros(N+1,1);
        E = zeros(N+1,1);
        I = zeros(N+1,1);
        R = zeros(N+1,1);
        S(1) = sus2(1)/tot;
        E(1) = exposed(1)/tot;
        I(1) = infected(1)/tot;
        R(1) = removed(1)/tot;
        
        for i = 1:N
            T = Time(i);
            infection_func = iF(gov(i),B1,B2,S(i),E(i),I(i));
            k1 = h* fun_s_STR(infection_func);
            l1 = h* fun_e_STR(E(i),tot,infection_func);
            m1 = h* fun_i_STR(E(i),I(i),tot,gam);
            n1 = h* fun_r_STR(I(i),gam);
            
            k2 = h* fun_s_STR(infection_func);
            l2 = h* fun_e_STR(E(i)+(1/4)*l1*h,tot,infection_func);
            m2 = h* fun_i_STR(E(i)+(1/4)*l1*h,I(i)+(1/4)*m1*h,tot,gam);
            n2 = h* fun_r_STR(I(i)+(1/4)*m1*h,gam);
            
            k3 = h* fun_s_STR(infection_func);
            l3 = h* fun_e_STR(E(i)+(1/8)*l1*h+(1/8)*l2*h,tot,infection_func);
            m3 = h* fun_i_STR(E(i)+(1/8)*l1*h+(1/8)*l2*h,I(i)+(1/8)*m1*h+(1/8)*m2*h,tot,gam);
            n3 = h* fun_r_STR(I(i)+(1/8)*m1*h+(1/8)*m2*h,gam);
            
            k4 = h* fun_s_STR(infection_func);
            l4 = h* fun_e_STR(E(i)-(1/2)*l2*h+l3*h,tot,infection_func);
            m4 = h* fun_i_STR(E(i)-(1/2)*l2*h+l3*h,I(i)-(1/2)*m2*h+m3*h,tot,gam);
            n4 = h* fun_r_STR(I(i)-(1/2)*m2*h+m3*h,gam);
            
            k5 = h* fun_s_STR(infection_func);
            l5 = h* fun_e_STR(E(i)+(3/16)*l1*h+(9/8)*l4*h,tot,infection_func);
            m5 = h* fun_i_STR(E(i)+(3/16)*l1*h+(9/8)*l4*h,I(i)+(3/16)*m1*h+(9/8)*m4*h,tot,gam);
            n5 = h* fun_r_STR(I(i)+(3/16)*m1*h+(9/8)*m4*h,gam);
            
            k6 = h* fun_s_STR(infection_func);
            l6 = h* fun_e_STR(E(i)-(3/7)*l1*h+(2/7)*l2*h+(12/7)*l3*h-(12/7)*l4*h+(8/7)*l5*h,tot,infection_func);
            m6 = h* fun_i_STR(E(i)-(3/7)*l1*h+(2/7)*l2*h+(12/7)*l3*h-(12/7)*l4*h+(8/7)*l5*h,I(i)-(3/7)*m1*h+(2/7)*m2*h+(12/7)*m3*h-(12/7)*m4*h+(8/7)*m5*h,tot,gam);
            n6 = h* fun_r_STR(I(i)-(3/7)*m1*h+(2/7)*m2*h+(12/7)*m3*h-(12/7)*m4*h+(8/7)*m5*h,gam);
            
             S(i+1) = S(i) + (1/90) * (7*k1 + 32*k2 + 12*k4 + 32*k5 + 7*k6);
             E(i+1) = E(i) + (1/90) * (7*l1 + 32*l2 + 12*l4 + 32*l5 + 7*l6);
             I(i+1) = I(i) + (1/90) * (7*m1 + 32*m2 + 12*m4 + 32*m5 + 7*m6);
             R(i+1) = R(i) + (1/90) * (7*n1 + 32*n2 + 12*n4 + 32*n5 + 7*n6);
        end
        
            