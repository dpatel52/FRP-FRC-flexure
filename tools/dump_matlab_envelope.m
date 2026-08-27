% Build the full fourteen zone envelope in MATLAB for three cases, exactly the way
% the driver does, and write it to CSV so the Python state machine can be checked
% end to end rather than zone by zone.
%
%   matlab -batch "addpath('<Matlab dir>'); dump_matlab_envelope"

cases = struct('name', {}, 'p', {});

% 1. Shabani et al. 2025, beam 4#4-F0.5-45, no bonded skin
p1.b = 250; p1.h = 400; p1.cover = 45.5; p1.L = 3836;
p1.E = 31528.6; p1.epcr = 1.27658e-4;
p1.mu_1 = 0.146; p1.mu_2 = 0.146; p1.mu_3 = 0.146;
p1.beta_1 = 2; p1.beta_2 = 50; p1.beta_3 = 300;
p1.xi = 1.001; p1.omega = 11.18; p1.mu_c = 1; p1.ecu = 0.005;
p1.Es = 61200; p1.kappa_strain = 0.0174; p1.mu_s = 1.0; p1.chi_su_strain = 0.0175;
p1.rho_t = 507/(250*400); p1.rho_c = 0;
p1.iota = 0; p1.psi = 0; p1.rho_x = 0;
cases(1).name = 'shabani_noskin'; cases(1).p = p1;

% 2. same section, with a bonded skin active from the start
p2 = p1; p2.psi = 6.6; p2.rho_x = 0.0011; p2.iota = 0;
cases(2).name = 'skin_iota0'; cases(2).p = p2;

% 3. same section, skin activating after cracking
p3 = p1; p3.psi = 6.6; p3.rho_x = 0.0011; p3.iota = 10.0;
cases(3).name = 'skin_iota10'; cases(3).p = p3;

fid = fopen('matlab_envelope.csv','w');
fprintf(fid,'case,i,beta,k,M\n');
pf = fopen('matlab_envelope_params.csv','w');
fprintf(pf,['case,b,h,alpha,L,E,epcr,beta_1,beta_2,beta_3,eta_1,eta_2,eta_3,xi,', ...
            'omega,eta_c,n,kappa,eta_s,rho_c,rho_t,iota,psi,rho_x,M_cr\n']);

for c = 1:numel(cases)
    p = cases(c).p;  nm = cases(c).name;
    alpha = (p.h - p.cover)/p.h;
    eta_1 = (p.mu_1 - 1)/(p.beta_1 - 1);
    eta_2 = (p.mu_2 - p.mu_1)/(p.beta_2 - p.beta_1);
    eta_3 = (p.mu_3 - p.mu_2)/(p.beta_3 - p.beta_2);
    lambda_cu = p.ecu/p.epcr;
    eta_c = p.omega*(p.mu_c - 1)/(lambda_cu - p.omega);
    n = p.Es/p.E;
    kappa = p.kappa_strain/p.epcr;
    chi_su = p.chi_su_strain/p.epcr;
    eta_s = kappa*(p.mu_s - 1)/(chi_su - kappa);

    [kcr, ~] = zone111(1, p.L, p.b, p.h, alpha, p.E, p.epcr, p.beta_1, p.beta_2, p.beta_3, ...
                       eta_1, eta_2, eta_3, p.xi, p.omega, eta_c, n, kappa, eta_s, ...
                       p.rho_c, p.rho_t, p.iota, p.psi, p.rho_x);
    M_cr = (p.epcr*p.E*p.b*p.h^2)/(12*(1 - kcr));

    fprintf(pf,'%s,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n', ...
        nm, p.b, p.h, alpha, p.L, p.E, p.epcr, p.beta_1, p.beta_2, p.beta_3, ...
        eta_1, eta_2, eta_3, p.xi, p.omega, eta_c, n, kappa, eta_s, ...
        p.rho_c, p.rho_t, p.iota, p.psi, p.rho_x, M_cr);

    beta_z1 = linspace(0, 1, 200)';
    beta_z2 = linspace(1, p.beta_1, 500)';
    beta_z3 = linspace(p.beta_1, p.beta_2, 500)';
    beta_z4 = linspace(p.beta_2, p.beta_3, 500)';

    a = {p.L, p.b, p.h, alpha, p.E, p.epcr, p.beta_1, p.beta_2, p.beta_3, ...
         eta_1, eta_2, eta_3, p.xi, p.omega, eta_c, n, kappa, eta_s, ...
         p.rho_c, p.rho_t, p.iota, p.psi, p.rho_x};

    % exactly as the driver wires it, each zone group over its own beta segment
    [k111,M111]   = zone111(beta_z1, a{:});
    [k211,M211]   = zone211(M_cr, beta_z2, a{:});
    [k212,M212]   = zone212(M_cr, beta_z2, a{:});
    [k221,M221]   = zone221(M_cr, beta_z2, a{:});
    [k222,M222]   = zone222(M_cr, beta_z2, a{:});
    [k311,M311]   = zone311(M_cr, beta_z3, a{:});
    [k312,M312]   = zone312(M_cr, beta_z3, a{:});
    [k321,M321]   = zone321(M_cr, beta_z3, a{:});
    [k322,M322]   = zone322(M_cr, beta_z3, a{:});
    [k411,M411]   = zone411(M_cr, beta_z4, a{:});
    [k412,M412]   = zone412(M_cr, beta_z4, a{:});
    [k421,M421]   = zone421(M_cr, beta_z4, a{:});
    [k422,M422]   = zone422(M_cr, beta_z4, a{:});
    [k4222,M4222] = zone4222(M_cr, beta_z4, a{:});

    Envelope = Envelope_Final(kappa, p.omega, p.epcr, beta_z1, beta_z2, beta_z3, beta_z4, ...
        k111,M111, k211,M211, k212,M212, k221,M221, k222,M222, ...
        k311,M311, k312,M312, k321,M321, k322,M322, ...
        k411,M411, k412,M412, k421,M421, k422,M422, k4222,M4222, ...
        p.beta_1, p.beta_2, p.beta_3, alpha);

    beta_all = [beta_z1; beta_z2; beta_z3; beta_z4];
    for i = 1:size(Envelope,1)
        fprintf(fid,'%s,%d,%.17g,%.17g,%.17g\n', nm, i, beta_all(i), Envelope(i,1), Envelope(i,2));
    end
    fprintf('%s done, %d rows\n', nm, size(Envelope,1));
end

fclose(fid); fclose(pf);
fprintf('wrote matlab_envelope.csv and matlab_envelope_params.csv\n');
