
time =  linspace(0,30,100)  
Cp = Cp10(time)

ktrans = 0.25
ve     = 0.50
kep    = ktrans / ve

Ct = Tofts([ktrans, kep],time,Cp)

% ######################################             DEFINE FUNCTIONS TO BE USED ABOVE

function [Cp_out] = Cp10(t)
    %Input:  time t in minutes
    %THIS FUNCTION CALCULATES AN AIF WITH A SIMULATED INJECTION TIME OF 10
    %SECONDS
    %Injection of 10 seconds
    A= 30.0 ; %mM/min
    B= 1.0  ;
    C= 4.0  ; %min^-1
    D= 0.65 ; %mM
    E= 5.0  ; %min
    F= 0.04 ;  %min-1
    
    Cp_out=A.*(t.^B).*exp(-t.*C)+ D.*(1-exp(-t.*E)).* exp(-t.*F);%inject = 10sec
    
    end

function [c_toi] = Tofts(pars,time,Cp)
        % Simulate Concentration in a tissue of interest usign the Tofs model
        %
        % [c_toi] =myToftsCt2(ktrans,kep,t,Cp)
        %
        % ktrans=pars(1);
        % kep=pars(2);
        % t=X(:,1);
        % Cp=X(:,2);
        
        % Authors:
        % Joey DeGrandchamp                 Julio Cardenas-Rodriguez
        % University of Arizona             University of Arizona
        % jdegrandchamp@email.arizona.edu   cardenaj@email.arizona.edu
        %
        %                       www.cardenaslab.org
        
        ktrans=pars(1);
        kep=pars(2);
        
        n_points=length(time);
        % expo=zeros(1,n_points);
        % crpexp=expo;
        % c_toi=crpexp;
        
        c_toi=zeros(n_points,1);
        
        for k = 2:n_points
            int_t = time(k);
            
            for j = 1:k
                dummy_time = time(j);
                expo(j) =exp(-((kep.*(int_t-dummy_time)))); %#ok<*AGROW>
                crpexp(j) = Cp(j)*expo(j);
            end
            
            t2 = time(1:k);
            %%
            crpexp_integral = trapz(t2,crpexp);
            c_toi(k) = ktrans*crpexp_integral; 
        end
        
        end
        