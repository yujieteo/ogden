clear all;

% Imports N = 1, eta = 1, lambda_max = 4.5, s-lambda loading curve.
dat1_imp    =   detectImportOptions('dat/2d-s-lambda-eta-1-N-1-loading.csv');
dat1_dbl    =   readmatrix('dat/2d-s-lambda-eta-1-N-1-loading.csv', dat1_imp);
load_x_lambda_Nr1     =   dat1_dbl(:,1);
load_y_lambda_Nr1     =   dat1_dbl(:,2);

% Imports smax-N curve, lambda_max = 3, decay curve.
dat2_imp    =   detectImportOptions('dat/2e-smax-N-lambdamax-3.csv');
dat2_dbl    =   readmatrix('dat/2e-smax-N-lambdamax-3.csv', dat2_imp);
decay_x     =   dat2_dbl(:,1);
decay_y     =   dat2_dbl(:,2)./dat2_dbl(1,1);   % Normalised decay.

% Imports N = 1, lambda_max = 4.5, s-lambda UN-loading curve.
dat3_imp    =   detectImportOptions('dat/2d-s-lambda-eta-1-N-1-unloading.csv');
dat3_dbl    =   readmatrix('dat/2d-s-lambda-eta-1-N-1-unloading.csv', dat3_imp);
load_x_lambda_Nu1   =   dat3_dbl(:,1);
load_y_sigma_Nu1    =   dat3_dbl(:,2);
lambdamax           =   4.5;

% Imports N = 2, lambda_max = 4.5, s-lambda UN-loading curve.
dat4_imp    =   detectImportOptions('dat/3d-s-lambda-eta-1-N-2-unloading.csv');
dat4_dbl    =   readmatrix('dat/3d-s-lambda-eta-1-N-2-unloading.csv', dat4_imp);
load_x_lambda_Nu2   =   dat4_dbl(:,1);
load_y_sigma_Nu2    =   dat4_dbl(:,2);
lambdamax           =   4.5;
N                   =   2;

% Imports N = 2, lambda_max = 4.5, s-lambda loading curve.
dat5_imp    =   detectImportOptions('dat/3d-s-lambda-eta-1-N-2-loading.csv');
dat5_dbl    =   readmatrix('dat/3d-s-lambda-eta-1-N-2-loading.csv', dat5_imp);
load_x_lambda_Nr2   =   dat5_dbl(:,1);
load_y_sigma_Nr2    =   dat5_dbl(:,2);
lambdamax           =   4.5;
N                   =   2;

% Fitting of Ogden model.

ft_N1r    =   fittype(  'mu*(x^(alpha-1) - x^(-alpha-1))',...
                'independent', 'x', ...
                'coefficients', {'mu','alpha'});

fo_N1r        =   fitoptions(ft_N1r);
fo_N1r.Lower  =   [3.97, 2.158];
fo_N1r.Upper  =   [3.97, 2.158];

f_N1r         = fit(load_x_lambda_Nr1,load_y_lambda_Nr1,ft_N1r,fo_N1r);
co_N1r        = coeffvalues(f_N1r);

mu              = co_N1r(1,1);
alpha           = co_N1r(1,2);

% Fitting of shakedown st   ress decay

ft_decay        =   fittype( 'kappa*(N)^(-gamma)',...
                    'independent', 'N', ...
                    'coefficients', {'kappa','gamma'});

fo_decay        =   fitoptions(ft_decay);
fo_decay.Lower  =   [0.85, 0.16];
fo_decay.Upper  =   [0.85, 0.16];

f_decay         = fit(decay_x,decay_y,ft_decay,fo_decay);
co_decay        = coeffvalues(f_decay);

kappa           = co_decay(1,1);
gma             = co_decay(1,2);

stress_anon = @(r, c0, beta, nu1, nu3, alpha, mu, lambda0, lambda) ...
    (1 - (1./r) .* erf((((mu./alpha).*((lambda0.^alpha) + 1 + (lambda0.^-alpha) - 3) + (mu./2).*(nu1.*(lambda0.^2 - 1) + nu3.*((lambda0.^(-2) - 1)))) - ((mu./alpha).*((lambda.^alpha) + 1 + (lambda.^-alpha) - 3)  + (mu./2).*(nu1.*(lambda.^2 - 1) + nu3.*((lambda.^(-2) - 1)))))./(c0 + (beta .* ((mu./alpha).*((lambda0.^alpha) + 1 + (lambda0.^-alpha) - 3)  + (mu./2).*(nu1.*(lambda0.^2 - 1) + nu3.*((lambda0.^(-2) - 1)))))))) .* mu .* (lambda.^(alpha - 1) ...
    - lambda .^(-alpha - 1)) + (1 - (1 - (1./r) .* erf((((mu./alpha).*((lambda0.^alpha) + 1 + (lambda0.^-alpha) - 3)  + (mu./2).*(nu1.*(lambda0.^2 - 1) + nu3.*((lambda0.^(-2) - 1)))) - ((mu./alpha).*((lambda.^alpha) + 1 + (lambda.^-alpha) - 3)  + (mu./2).*(nu1.*(lambda.^2 - 1) + nu3.*((lambda.^(-2) - 1)))))./(c0 + (beta .* ((mu./alpha).*((lambda0.^alpha) + 1 + (lambda0.^-alpha) - 3)  + (mu./2).*(nu1.*(lambda0.^2 - 1) + nu3.*((lambda0.^(-2) - 1))))))))) .* mu .* ...
    (nu1 .* lambda - nu3 .* (lambda .^ (-3)));

stress_N_u_anon = @(cu, r, beta, nu1, nu3, alpha, mu, lambda0, kappa, gma, N, lambda) ...
    kappa .* (N .^ (- gma)) .* ...
    (1 - (1./r) .* erf((((mu./alpha).*((lambda0.^alpha) + 1 + (lambda0.^-alpha) - 3) + (mu./2).*(nu1.*(lambda0.^2 - 1) + nu3.*((lambda0.^(-2) - 1)))) - ((mu./alpha).*((lambda.^alpha) + 1 + (lambda.^-alpha) - 3)  + (mu./2).*(nu1.*(lambda.^2 - 1) + nu3.*((lambda.^(-2) - 1)))))./(cu + (beta .* ((mu./alpha).*((lambda0.^alpha) + 1 + (lambda0.^-alpha) - 3)  + (mu./2).*(nu1.*(lambda0.^2 - 1) + nu3.*((lambda0.^(-2) - 1)))))))) .* mu .* (lambda.^(alpha - 1) ...
    - lambda .^(-alpha - 1)) + (1 - (1 - (1./r) .* erf((((mu./alpha).*((lambda0.^alpha) + 1 + (lambda0.^-alpha) - 3)  + (mu./2).*(nu1.*(lambda0.^2 - 1) + nu3.*((lambda0.^(-2) - 1)))) - ((mu./alpha).*((lambda.^alpha) + 1 + (lambda.^-alpha) - 3)  + (mu./2).*(nu1.*(lambda.^2 - 1) + nu3.*((lambda.^(-2) - 1)))))./(cu + (beta .* ((mu./alpha).*((lambda0.^alpha) + 1 + (lambda0.^-alpha) - 3)  + (mu./2).*(nu1.*(lambda0.^2 - 1) + nu3.*((lambda0.^(-2) - 1))))))))) .* mu .* ...
    (nu1 .* lambda - nu3 .* (lambda .^ (-3)));

stress_N_r_anon = @(cr, r, beta, nu1, nu3, alpha, mu, lambda0, kappa, gma, N, lambda) ...
    kappa .* (N .^ (- gma)) .* ...
    (1 - (1./r) .* erf((((mu./alpha).*((lambda0.^alpha) + 1 + (lambda0.^-alpha) - 3) + (mu./2).*(nu1.*(lambda0.^2 - 1) + nu3.*((lambda0.^(-2) - 1)))) - ((mu./alpha).*((lambda.^alpha) + 1 + (lambda.^-alpha) - 3)  + (mu./2).*(nu1.*(lambda.^2 - 1) + nu3.*((lambda.^(-2) - 1)))))./(cr + (beta .* ((mu./alpha).*((lambda0.^alpha) + 1 + (lambda0.^-alpha) - 3)  + (mu./2).*(nu1.*(lambda0.^2 - 1) + nu3.*((lambda0.^(-2) - 1)))))))) .* mu .* (lambda.^(alpha - 1) ...
    - lambda .^(-alpha - 1)) + (1 - (1 - (1./r) .* erf((((mu./alpha).*((lambda0.^alpha) + 1 + (lambda0.^-alpha) - 3)  + (mu./2).*(nu1.*(lambda0.^2 - 1) + nu3.*((lambda0.^(-2) - 1)))) - ((mu./alpha).*((lambda.^alpha) + 1 + (lambda.^-alpha) - 3)  + (mu./2).*(nu1.*(lambda.^2 - 1) + nu3.*((lambda.^(-2) - 1)))))./(cr + (beta .* ((mu./alpha).*((lambda0.^alpha) + 1 + (lambda0.^-alpha) - 3)  + (mu./2).*(nu1.*(lambda0.^2 - 1) + nu3.*((lambda0.^(-2) - 1))))))))) .* mu .* ...
    (nu1 .* lambda - nu3 .* (lambda .^ (-3)));


ft_N1u       =   fittype( stress_anon,...
                'coefficients', {'r','c0','beta','nu1','nu3'}, ...
                'independent', 'lambda', ...
                'problem', {'alpha', 'mu', 'lambda0'} ...
                );

fo_N1u        =   fitoptions(ft_N1u);
fo_N1u.Lower  =   [1, 5, 0, 0.1, 0.1];
fo_N1u.Upper  =   [Inf, Inf, Inf, Inf, Inf];

f_N1u         =    fit(load_x_lambda_Nu1, ...
                    load_y_sigma_Nu1, ...
                    ft_N1u, ...
                    fo_N1u, ...
                    'problem', {alpha, mu, lambdamax});

co_N1u      = coeffvalues(f_N1u);

r           = co_N1u(1,1);
c0          = co_N1u(1,2);
beta        = co_N1u(1,3);
nu1         = co_N1u(1,4);
nu3         = co_N1u(1,5);

ft_Nu2    =   fittype( stress_N_u_anon,...
                'coefficients', {'cu'}, ...
                'independent', 'lambda', ...
                'problem', {'r','beta','nu1','nu3','alpha', 'mu', 'lambda0', 'kappa', 'gma', 'N'} ...
                );

fo_Nu2          =   fitoptions(ft_Nu2);
fo_Nu2.Lower    =   25;
fo_Nu2.Upper    =   25;

f_Nu2           =   fit(load_x_lambda_Nu2, ...
                    load_y_sigma_Nu2, ...
                    ft_Nu2, ...
                    fo_Nu2, ...
                    'problem', {r, beta, nu1, nu3, alpha, mu, lambdamax, kappa, gma, N});

co_Nu2          =   coeffvalues(f_Nu2);
cu              =   co_Nu2;

ft_Nr2          =   fittype( stress_N_r_anon,...
                'coefficients', {'cr'}, ...
                'independent', 'lambda', ...
                'problem', {'r','beta','nu1','nu3','alpha', 'mu', 'lambda0', 'kappa', 'gma', 'N'} ...
                );

fo_Nr2          =   fitoptions(ft_Nr2);
fo_Nr2.Lower    =   45;
fo_Nr2.Upper    =   45;

f_Nr2           =   fit(load_x_lambda_Nr2, ...
                    load_y_sigma_Nr2, ...
                    ft_Nr2, ...
                    fo_Nr2, ...
                    'problem', {r, beta, nu1, nu3, alpha, mu, lambdamax, kappa, gma, N});

co_Nr2          =   coeffvalues(f_Nr2);
cr              =   co_Nr2;

res = get(0,'screensize');
fig = figure;
set(fig, 'position', res);

subplot(1,2,1);

hold on;

xlim([1 4.5]);
ylim([-0.5 23]);
N1u_expt    =   plot(load_x_lambda_Nu1,load_y_sigma_Nu1,'Ro');
N1u_fit     =   plot(f_N1u,'B-');
N1r_expt    =   plot(load_x_lambda_Nr1,load_y_lambda_Nr1,'Ro');
N1r_fit     =   plot(f_N1r,'B-');

plot1_legend = legend([N1u_expt N1u_fit N1r_expt, N1r_fit], ...
    'Unloading (experimental)', ...
    'Unloading (curve fit)', ...
    'Loading (experimental)', ...
    'Loading (curve fit)', ...
    'Location','southeast','interpreter','latex');

title1 = title('Reloading-unloading hystersis ($\lambda_\mathrm{max} = 4.5$, $N = 1$)','interpreter','latex');
xlabel1 = xlabel('Stretch $\lambda$','interpreter','latex');
ylabel1 = ylabel('Nominal stress $\sigma_{\mathrm{nom}}$ (kPa)','interpreter','latex');

ax = gca;
ax.FontSize = 12; 

set(plot1_legend,'FontSize',16);
set(title1,'FontSize',16);
set(xlabel1,'FontSize',16);
set(ylabel1,'FontSize',16);

xlim([1 4.75]);
ylim([-1 25]);
grid on

hold off;

subplot(1,2,2);

hold on;

xlim([1 4.5]);
ylim([-0.5 17]);

N2u_expt    =   plot(load_x_lambda_Nu2,load_y_sigma_Nu2, 'Ro');
N2u_fit     =   plot(f_Nu2, 'B-');
N2r_expt    =   plot(load_x_lambda_Nr2,load_y_sigma_Nr2, 'Ro');
N2r_fit     =   plot(f_Nr2, 'B-');

grid on

plot2_legend = legend([N2u_expt N2u_fit N2r_expt N2r_fit], ...
    'Unloading (experimental)', ...
    'Unloading (curve fit)', ...
    'Loading (experimental)', ...
    'Loading (curve fit)', ...
    'Location','southeast','interpreter','latex');

title2  = title('Reloading-unloading hystersis ($\lambda_\mathrm{max} = 4.5$, $N = 2$)','interpreter','latex');
xlabel2 = xlabel('Stretch $\lambda$','interpreter','latex');
ylabel2 = ylabel('Nominal stress $\sigma_{\mathrm{nom}}$ (kPa)','interpreter','latex');

ax = gca;
ax.FontSize = 12; 

set(plot2_legend,'FontSize',16);
set(title2,'FontSize',16);
set(xlabel2,'FontSize',16);
set(ylabel2,'FontSize',16);


xlim([1 4.75]);
ylim([-1 18]);

hold off;
