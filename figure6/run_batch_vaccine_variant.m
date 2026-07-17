clear; clc;

%% =========================
%  批量模拟设置
% =========================

% vaccine_types = {'nonRBD', 'RBD', 'inactivated'};
% vaccine_doses = [1e10, 1e12, 1e14, 1e16];
% alpha_list = 0:4;

vaccine_types = {'inactivated'};
vaccine_doses = [1e10];
alpha_list = 0:4;
% 保存目录
result_dir = 'batch_results_vaccine_variant';
if ~exist(result_dir, 'dir')
    mkdir(result_dir);
end

% ODE 选项
options_default = odeset('RelTol', 1e-4, ...
                         'AbsTol', 1e-8);

%% =========================
%  主循环
% =========================

for ivac = 1:length(vaccine_types)

    vaccine_type = vaccine_types{ivac};

    for idose = 1:length(vaccine_doses)

        dose = vaccine_doses(idose);

        fprintf('\n====================================================\n');
        fprintf('Primary vaccination: type = %s, dose = %.1e\n', vaccine_type, dose);
        fprintf('====================================================\n');

        %% ---------------------------------------------------------
        %  1. 初次疫苗接种模拟
        % ---------------------------------------------------------

        [p_vac, y0_vac, tspan_vac] = build_primary_vaccine_model(vaccine_type, dose);

        non_neg_indices = 1:length(y0_vac);
        options = odeset(options_default, 'NonNegative', non_neg_indices);

        tic;
        [t_vac, y_vac] = ode15s(@(t, y) sys_ode_vaccine_dispatch(t, y, p_vac), ...
                                tspan_vac, y0_vac, options);
        time_primary = toc;

        fprintf('Primary vaccination finished. Time = %.2f min\n', time_primary / 60);

        % 保存初次接种结果
        primary_file = fullfile(result_dir, ...
            sprintf('primary_%s_dose_%1.0e.mat', vaccine_type, dose));

        metadata_primary = struct();
        metadata_primary.vaccine_type = vaccine_type;
        metadata_primary.dose = dose;
        metadata_primary.stage = 'primary_vaccination';
        metadata_primary.runtime_seconds = time_primary;

        save(primary_file, 't_vac', 'y_vac', 'p_vac', 'metadata_primary', '-v7.3');

        %% ---------------------------------------------------------
        %  2. 对每一个 alpha 进行二次感染模拟
        % ---------------------------------------------------------

    plot_alpha_0
    end
end

fprintf('\nAll simulations finished.\n');