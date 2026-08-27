% Evaluate all fourteen zone solvers on a fixed set of random parameter vectors
% and write the results to CSV, so the Python translation can be checked against
% the MATLAB original it came from.
%
% Run from this folder with the MATLAB sources on the path:
%   matlab -batch "addpath('<Matlab dir>'); dump_matlab_reference"

rng(20260827);

zones = {'zone111','zone211','zone212','zone221','zone222','zone311','zone312', ...
         'zone321','zone322','zone411','zone412','zone421','zone422','zone4222'};

nSets = 12;
% column vector on purpose. The zone files loop `for i = 1:size(beta_array)`,
% and size() on a ROW vector returns [1 N], so 1:[1 N] collapses to 1:1 and only
% the first entry is filled. A column vector gives [N 1] and the loop is correct.
betas = [1.05; 2; 5; 12; 30; 60; 120; 250];

fid = fopen('matlab_reference.csv','w');
fprintf(fid, 'zone,set,beta,k,M\n');

for s = 1:nSets
    % random but physically sensible section
    L     = 1000 + 3000*rand;
    b     = 100  + 200*rand;
    h     = 150  + 350*rand;
    alpha = 0.70 + 0.25*rand;
    E     = 20000 + 30000*rand;
    epcr  = 1.0e-4 + 1.0e-4*rand;
    beta_1 = 2 + 10*rand;
    beta_2 = beta_1 + 10 + 60*rand;
    beta_3 = beta_2 + 50 + 300*rand;
    mu_1 = 0.2 + 1.2*rand;  mu_2 = 0.1 + 0.6*rand;  mu_3 = 0.02 + 0.3*rand;
    eta_1 = (mu_1 - 1)/(beta_1 - 1);
    eta_2 = (mu_2 - mu_1)/(beta_2 - beta_1);
    eta_3 = (mu_3 - mu_2)/(beta_3 - beta_2);
    xi    = 1.0 + 0.5*rand;
    omega = 6 + 20*rand;
    mu_c  = 0.6 + 0.4*rand;
    lambda_cu = omega + 5 + 30*rand;
    eta_c = omega*(mu_c - 1)/(lambda_cu - omega);
    n     = 3 + 6*rand;
    kappa = 8 + 30*rand;
    mu_s  = 1.0 + 0.5*rand;
    chi_su = kappa + 20 + 200*rand;
    eta_s = kappa*(mu_s - 1)/(chi_su - kappa);
    rho_c = 0.002 + 0.008*rand;
    rho_t = 0.004 + 0.020*rand;
    iota  = 0.5*rand*10;
    psi   = 2 + 6*rand;
    rho_x = 0.001 + 0.006*rand;

    % M_cr exactly as the driver builds it, from zone111 at beta = 1
    [kcr, ~] = zone111(1, L, b, h, alpha, E, epcr, beta_1, beta_2, beta_3, ...
                       eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, ...
                       rho_c, rho_t, iota, psi, rho_x);
    M_cr = (epcr * E * b * h^2) / (12 * (1 - kcr));

    % echo the parameter set so Python uses identical inputs
    if s == 1
        pf = fopen('matlab_params.csv','w');
        fprintf(pf, ['set,L,b,h,alpha,E,epcr,beta_1,beta_2,beta_3,eta_1,eta_2,eta_3,', ...
                     'xi,omega,eta_c,n,kappa,eta_s,rho_c,rho_t,iota,psi,rho_x,M_cr\n']);
    else
        pf = fopen('matlab_params.csv','a');
    end
    fprintf(pf, '%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n', ...
            s, L, b, h, alpha, E, epcr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, ...
            xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t, iota, psi, rho_x, M_cr);
    fclose(pf);

    for z = 1:numel(zones)
        zn = zones{z};
        f  = str2func(zn);
        try
            if strcmp(zn, 'zone111')
                [kk, MM] = f(betas, L, b, h, alpha, E, epcr, beta_1, beta_2, beta_3, ...
                             eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, ...
                             rho_c, rho_t, iota, psi, rho_x);
            else
                [kk, MM] = f(M_cr, betas, L, b, h, alpha, E, epcr, beta_1, beta_2, beta_3, ...
                             eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, ...
                             rho_c, rho_t, iota, psi, rho_x);
            end
            for i = 1:numel(betas)
                fprintf(fid, '%s,%d,%.17g,%.17g,%.17g\n', zn, s, betas(i), kk(i), MM(i));
            end
        catch ME
            fprintf('  %s set %d FAILED: %s\n', zn, s, ME.message);
        end
    end
end

fclose(fid);
fprintf('wrote matlab_reference.csv and matlab_params.csv\n');
