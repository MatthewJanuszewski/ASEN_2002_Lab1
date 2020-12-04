%% Set Up

clc
clear

% assign variables 

L = 27.25; 
W_not = 2001; 
p  = 50; 

% Aluminum 6061

rho = 0.098; %[lb/in^3]
bendingYS = 35; %[ksi] %????
% shear_YS = 27; %[ksi]
cost = 8.03; %[$/lb]
t_min = 0.25; %[in]

% equations

% W_w= -rho*A; 

%% Calculate Area

% circle

 r_Vec = linspace(1/8,6,100);

    for i = 1:100
        
          r = r_Vec(i);

          A = pi*r^2; 
          I = (pi*r^4)/4;
          y = r; % max distance from centroid
          
          W_w = -rho*A; 
          M = moments(L,W_not,W_w);
          n_1_Vec(i) = safetyfactor(I,y,M,bendingYS);

    end

    [min_safety,I_1] = min(n_1_Vec);
    r_min = r_Vec(I_1); 
    
% shape 2

l_Vec = linspace(1/4,12,100);
w_Vec = linspace(1/4,12,100);

    for i = 1:100
        
        for j = 1:100
            
            l = l_Vec(i);
            w = w_Vec(j);

            A = l*w;
            I = (w*l^3)/12;
            y = l/2;
            
            W_w = -rho*A; 
            M = moments(L,W_not,W_w);
            n_2_Vec = safetyfactor(I,y,M,bendingYS); % matrix???
          
        end
        
    end
    
    
% shape 3
        
    for l = minl:maxl
        
       for l2 = minl2:maxl2

             A = l^2 - (l - 2*t)^2; 
             I = ((l^4)/12) - ((l2^4)/12);
             y = l/2;
             
             W_w = -rho*A; 
             M = moments(L,W_not,W_w);
             n_3 = safetyfactor(I,y,M,bendingYS);

       end
       
    end
    
% shape 4 

% assume width + thickness of top/bottom pieces are the same 

    for l = minl:maxl
        
        for w = minw:maxw
            
            for t = mint:maxt
                
                for h = minh:maxh

                    A = (l*t)+2*(h*w);
                    I_prime = ((t*l^3)/12) + (2*(w*h^3)/12); 

                    y_prime = l/2 + h/2; %top + bottom piece centroid y distance from actual centroid

                    A_t_b = h*w;

                    I = I_prime + (A_t_b*(y_prime^2));
                    y = (l/2)+(h/2);

                    W_w = -rho*A; 
                    M = moments(L,W_not,W_w);
                    n_4 = safetyfactor(I,y,M,bendingYS);
           
                end
            end
        end
    end

% shape 5


    for l = minl:maxl
        
        for w = minw:maxw
            
            for t = mint:maxt
                
                for h = minh:maxh

                    A = (w*h)+(l*t);

                    I_prime = ((t*l^3)/12) + ((w*h^3)/12); 

                    A_top = h*w;
                    A_middle = l*t;

                    y_c = (((l+(h/2))*(A_top)) + ((l/2)*(A_middle)))/(A_top + A_middle); %distance from bottom

                    y_prime_top = (l+(h/2)) - y_c;
                    y_prime_middle = y_c - (l/2);

                    I = I_prime + (A_top*(y_prime_top^2)) + (A_middle*(y_prime_middle^2)); 

                        if y_prime_middle > y_prime_top

                            y =  y_prime_middle;

                        else

                            y = y_prime_top;
                        end

                    W_w = -rho*A; 
                    M = moments(L,W_not,W_w);
                    n_5 = safetyfactor(I,y,M,bendingYS);

                end
            end 
        end
    end

%% Moment Function

function [Moment] = moments(L,W_not,W_w)

x = 0; % max bending moment occurs at origin

% Calculate moment using distributed load equation
    
M_dist = (W_not/6)*(L-x)^2*(1-(x/L)); %calculate moment caused by distributed load

M_weight = W_w*((L-x)^2)/2; %calculate moment caused by weight 

Moment = M_dist + M_weight; 

end

%% FS Function

function n_bending = safetyfactor(I,y,M,bendingYS)

sigma_bending = -M*y/I; 

n_bending = sigma_bending/bendingYS; 

end
