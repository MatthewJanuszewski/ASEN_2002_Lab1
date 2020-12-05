%% Set Up

clc
clear

% assign variables 

L = 27.25; 
W_not = 2001; 
p  = 50; 

% Aluminum 6061

rho = 0.098; %[lb/in^3]
bendingYS = 35*1000; %[psi] %????
% shear_YS = 27; %[ksi]
cost = 8.03; %[$/lb]
t_min = 0.25; %[in]

% equations

% W_w= -rho*A; 

%% Calculate Area

% circle

 r_Vec = linspace(1/8,6,1000);
 A_good_1 = [];

    for i = 1:1000
        
          r = r_Vec(i);

          A = pi*r^2; 
          I = (pi*r^4)/4;
          y = r; % max distance from centroid
          
          W_w = -rho*A; 
          M = moments(L,W_not,W_w);
          n_1 = safetyfactor(I,y,M,bendingYS);
          
          if n_1 >= 1.5 && n_1 <= 1.53
              
              A_good_1 = [A_good_1,A];
              
          end
              

    end
    
    r_good = sqrt(A_good_1/pi);
    minA_1 = min(A_good_1)
    
% shape 2

l_Vec = linspace(1/4,12,500);
w_Vec = linspace(1/4,12,500);

 A_good_2 = [];
 l_good_2 = [];
 w_good_2 = [];

    for i = 1:500
        
        for j = 1:500
            
            l = l_Vec(i);
            w = w_Vec(j);

            A = l*w;
            I = (w*l^3)/12;
            y = l/2;
            
            W_w = -rho*A; 
            M = moments(L,W_not,W_w);
            n_2 = safetyfactor(I,y,M,bendingYS);
            
           if n_2 >= 1.5 && n_2 <= 1.53
            
            A_good_2 = [A_good_2,A];
            l_good_2 = [l_good_2,l];
            w_good_2 = [w_good_2,w];
            
           end 
          
        end
        
    end
    
  [minA_2,I2] = min(A_good_2)
  min_l_2 = l_good_2(I2)
  min_w_2 = w_good_2(I2)
    
% shape 3

A_good_3 = [];
l_good_3 = [];
l2_good_3 = [];

l_3_Vec = linspace(0.5,12,500);
l2_3_Vec = linspace(0,11.5,500);
        
    for i = 1:500
        
        l = l_3_Vec(i);
         
       for j = 1:500
           
         l2 = l2_3_Vec(j); 
           
           if l2 <= l - 0.5
               
              A = l^2 - (l2^2); 
              I = ((l^4)/12) - ((l2^4)/12);
              y = l/2;

              W_w = -rho*A; 
              M = moments(L,W_not,W_w);  

              n_3 = safetyfactor(I,y,M,bendingYS);
             
           else
               
             n_3 = nan; %skip where shape inverts 
             
           end
           
           
           if n_3 >= 1.5 && n_3 <= 1.53
               
             A_good_3 = [A_good_3,A];
             l_good_3 = [l_good_3,l];
             l2_good_3 = [l2_good_3,l2];
             
           end

       end
       
    end
    
 [minA_3,I3] = min(A_good_3,[],'omitnan')
 min_l_3 = l_good_3(I3)
 min_l2_3 = l2_good_3(I3)
    
% shape 4 

% area of top +bottom piece will always be the same 

l_4_Vec = linspace(1/4,11.5,100);
w_4_Vec = linspace(1/4,12,100);
t_4_Vec = linspace(1/4,12,100);
h_4_Vec = linspace(1/4,6,100);

A_good_4 = [];
l_good_4 = [];
w_good_4 = [];
h_good_4 = [];
t_good_4 = [];

    for i = 1:100
        
        l = l_4_Vec(i);
        
        for j = 1:100
            
            w = w_4_Vec(j);
            
            for k = 1:100
                
                t = t_4_Vec(k);
                
                for z = 1:100
                    
                    h = h_4_Vec(z);
                    
                    if (l + 2*h) <= 12
                        
                        A = (l*t)+2*(h*w);
                        I_prime = ((t*l^3)/12) + (2*(w*h^3)/12); 

                        y_prime = l/2 + h/2; %top + bottom piece centroid y distance from actual centroid

                        A_t_b = h*w;

                        I = I_prime + 2*(A_t_b*(y_prime^2));
                        y = (l/2) + h;

                        W_w = -rho*A; 
                        M = moments(L,W_not,W_w);

                        n_4 = safetyfactor(I,y,M,bendingYS);

                    else

                         n_4 = nan; %skip where shape inverts 
             
                    end
                    
                    if n_4 >= 1.5 && n_4 <= 1.53
               
                        A_good_4 = [A_good_4,A];
                        l_good_4 = [l_good_4,l];
                        w_good_4 = [w_good_4,w];
                        h_good_4 = [h_good_4,h];
                        t_good_4 = [t_good_4,t];
             
                    end
           
                end
            end
        end
    end
    
 [minA_4,I4] = min(A_good_4,[],'omitnan')
 min_l_4 = l_good_4(I4)
 min_w_4 = w_good_4(I4)
 min_h_4 = h_good_4(I4)
 min_t_4 = t_good_4(I4)
    

% shape 5

l_5_Vec = linspace(1/4,11.75,100);
w_5_Vec = linspace(1/4,12,100);
t_5_Vec = linspace(1/4,12,100);
h_5_Vec = linspace(1/4,12,100);

A_good_5 = [];
l_good_5 = [];
w_good_5 = [];
h_good_5 = [];
t_good_5 = [];

    for i = 1:100
        
        l = l_5_Vec(i);
        
        for j = 1:100
            
            w = w_5_Vec(j);
            
            for k = 1:100
                
                t = t_5_Vec(k);
                
                for z = 1:100
                    
                    h = h_5_Vec(z); 
                    
                    if (l + h) <= 12

                        A = (w*h)+(l*t);

                        I_prime = ((t*l^3)/12) + ((w*h^3)/12); 

                        A_top = h*w;
                        A_middle = l*t;

                        y_c = (((l+(h/2))*(A_top)) + ((l/2)*(A_middle)))/(A_top + A_middle); %distance from bottom

                        y_prime_top = abs((l+(h/2)) - y_c);
                        y_prime_middle = abs(y_c - (l/2));

                        I = I_prime + (A_top*(y_prime_top^2)) + (A_middle*(y_prime_middle^2)); 

                            if y_prime_middle > y_prime_top

                                y =  y_prime_middle;

                            else

                                y = y_prime_top;
                            end

                        W_w = -rho*A; 
                        M = moments(L,W_not,W_w);
                        n_5 = safetyfactor(I,y,M,bendingYS);

                    else

                        n_5 = nan; %skip where shape inverts
                        
                    end
                        
                    if n_5 >= 1.5 && n_5 <= 1.53
               
                        A_good_5 = [A_good_5,A];
                        l_good_5 = [l_good_5,l];
                        w_good_5 = [w_good_5,w];
                        h_good_5 = [h_good_5,h];
                        t_good_5 = [t_good_5,t];
             
                    end
                end
            end 
        end
    end
    
    
 [minA_5,I5] = min(A_good_5,[],'omitnan')
 min_l_5 = l_good_5(I5)
 min_w_5 = w_good_5(I5)
 min_h_5 = h_good_5(I5)
 min_t_5 = t_good_5(I5)


%% Moment Function

function [Moment] = moments(L,W_not,W_w)

x = 0; % max bending moment occurs at origin

% Calculate moment using distributed load equation
    
M_dist = (0.5*L*W_not)*(1/3*L); %distributed load at x = 0;

%M_dist = (W_not/6)*(L-x)^2*(1-(x/L)); 


M_weight = W_w*((L-x)^2)/2; %calculate moment caused by weight 

Moment = M_dist + M_weight; 

end

%% FS Function

function n_bending = safetyfactor(I,y,M,bendingYS)

sigma_bending = -M*y/I; 

n_bending = abs(bendingYS/sigma_bending); 

end
