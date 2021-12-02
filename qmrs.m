function [sol,numIts,touch] = qmrs(A,r,target,epsilon,maxiter,K,eigv0,shift,touch)                                                                                    
%% QMRS   Symmetric Quasi-Minimum-Residual (Solves A*t = -r)
%
% [t,k,hist] = QMRS(A,r,Etolerance,maxiter,K,eigv0,shift,hist)
%
% A             n by m matrix
% r             n by 1 vector
% normA         Estimated norm of A (not A^TA)
% epsilon       User target tolerance
% maxiter       maximum number of iterations
% K             preconditioner for A
% eigv0         initial eigenvalue estimate
% shift         shift needed for eigenresidual estimate
% hist          eigenresidual convergence history (appended to input)
%
% Preconditioner options for QMRS include:
%       K = {P1,P2} where P1 and P2 are function handles that approximate
%       inv(A)*x = P1(P2(x));
%
%       K = P1 where P1 is a function handle that approximates inv(A)*x
%       = P1(x);
%
%       K = P1 where P1 is a matrix that approximates A\x = P1\x
%
%       K = 1 for no preconditioning
LTolerance_factor = 1.8^(-touch);
ETolerance_factor = 1.8^(-touch);
LTolerance = norm(r)*eps;
ETolerance = norm(r)*0.1;

g = r;

%Preconditioning
if isempty(K)
    d = r;
elseif isa(K,'cell')
    d = K{1}(K{2}(r));
elseif isa(K,'function_handle')
    d = K(r);
elseif ismatrix(K)
    d = K\r;
else
    d = r;
end

rho_prev = g'*d;

Theta_prev = 0.0;

tau_init = norm(r);
tau_prev = tau_init;

Beta = 0.0; Beta_prev = 0.0;
Delta = 0.0; Delta_prev = 0.0;
Psi = 0.0; Psi_prev = 0.0;

%delta and sol are vectors of size n initially 0
delta = zeros(size(r));
sol = zeros(size(r));

eval_prev = eigv0; 

eres_updated = 0.0;
eres_prev = 0.0;

Gamma = 0.0; Gamma_prev = 0.0;
Phi = 0.0; Phi_prev = 0.0;

numIts = 0;

while numIts < maxiter
    w = A(d);
    sigma_prev = d'*w;
    
    if (sigma_prev == 0.0)
        %disp('sigma == 0')
        break;
    end
    
    alpha_prev = rho_prev/sigma_prev;
    if (abs(alpha_prev) < eps || abs(alpha_prev) > 1/eps)
        %disp('alpha_prev < eps || > 1/eps')
        break;
    end
    
    g = g - alpha_prev*w;
    
    Theta = g'*g;
    Theta = sqrt(Theta);
    Theta = Theta/tau_prev;
    
    c = 1.0/sqrt(1+Theta^2);
    
    tau = tau_prev*Theta*c;
    
    gamma = c^2*Theta_prev^2;
    
    eta = alpha_prev*c^2;
    
    delta = gamma*delta + eta*d;
    sol = delta + sol;
    
    numIts = numIts + 1;
    
    if rho_prev == 0
        %disp('rho_prev == 0')
        break;
    end
    
    if numIts > 1 && tau < LTolerance
            %disp('Met LTolerance')
            break;
    end
    
    Delta = gamma*Delta_prev + eta*rho_prev;
    Beta = Beta_prev - Delta;
    Phi = gamma^2*Phi_prev + eta^2*sigma_prev;
    Psi = gamma*Psi_prev + gamma*Phi_prev;
    Gamma = Gamma_prev + 2*Psi + Phi;
    
    dot_sol = sol'*sol;
    eval_updated = shift + (eigv0 - shift + 2*Beta + Gamma)/(1 + dot_sol);
    eres2_updated = (tau^2)/(1+dot_sol) + (((eigv0 - shift + Beta)^2)/(1+dot_sol) ...
        - (eval_updated - shift)^2);
    eres_prev = eres_updated;
    if eres2_updated < 0
        eres_updated = sqrt( (tau^2)/(1 + dot_sol) );
    else
        eres_updated = sqrt(eres2_updated);
    end
    
    if (ETolerance > 0.0 || ETolerance_factor > 0)
        if (numIts > 20 && (tau_prev <= eres_updated || eres_prev <= tau))
            %disp('tau_prev <= eresupdated || eresprev <= tau')
            break;
        end
        
        
        %Not really sure what this is for....
        %looking for smallest
        if eval_updated > eval_prev && strcmpi(target,'S')
            %disp('eval_updated > eval_prev');
            break;

        %looking for largest
        elseif eval_updated < eval_prev && strcmpi(target,'L')
                %disp(eval_updated < eval_prev');
                break;
        
        %looking for interior
        elseif abs(eigv0-eval_updated) > tau_init+eres_updated && target ~= 0 && target ~= inf
                %disp('abs(eigv0-eval_updated) > tau_init+eres_updated');
                break;   
        end
        
        
        if numIts > 1 && eres_updated < ETolerance
            %disp('ETolerance met')
            break;
        end
        
        tol = min(tau/LTolerance_factor, eres_updated/ETolerance_factor);
        if tol < epsilon
            %disp('Convergence test met')
            touch = touch + 1;
            break;
        end
        
        eval_prev = eval_updated;
        
    end
    if numIts < maxiter
        
        %Preconditioning
        if isempty(K)
            w = g;
        elseif isa(K,'cell')
            w = K{1}(K{2}(g));
        elseif isa(K,'function_handle')
            w = K(g);
        elseif ismatrix(K)
            w = K\g;
        else
            w = g;
        end
        
        
        rho = g'*w;
        beta = rho/rho_prev;
        d = w + beta*d;
        
        rho_prev = rho;
        tau_prev = tau;
        Theta_prev = Theta;
        Delta_prev = Delta;
        Beta_prev = Beta;
        Phi_prev = Phi;
        Psi_prev = Psi;
        Gamma_prev = Gamma;
    end
    if numIts == 1
            ETolerance = tau*0.1;
    end
end
if numIts == 0
    sol = d;
end

end