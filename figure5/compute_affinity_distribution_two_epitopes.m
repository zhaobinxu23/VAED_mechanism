function dist = compute_affinity_distribution_two_epitopes(t, y, t_grid, p)
% ========================================================================
% compute_affinity_distribution_two_epitopes
% ========================================================================
% Purpose:
%   Compute affinity-resolved antibody distributions for two epitope groups:
%
%   1) Neutralizing antibodies
%      - RBD-specific repertoire
%      - first epitope antibody block
%
%   2) Non-neutralizing antibodies
%      - non-RBD-specific repertoire
%      - second epitope antibody block
%
% This function is designed for the two-epitope model used in figure3a.m
% and sys_ode_two_epitopes_new.m.
%
% ------------------------------------------------------------------------
% State vector structure:
%
%   y = [T; I_1...I_N; V; Tc; X_1...X_3604; cumulative_death]
%
% According to sys_ode_two_epitopes_new.m:
%
%   X = y(N+4 : end-3)
%
% Therefore:
%
%   X(k) corresponds to y(:, p.N + 3 + k)
%
% because:
%   y(N+4) = X(1)
%
%
% ------------------------------------------------------------------------
% Antibody indices inside X:
%
% First epitope / RBD / neutralizing repertoire:
%   X(101:200)    = soluble IgM, neutralizing
%   X(301:400)    = soluble IgG, neutralizing
%
% Second epitope / non-RBD / non-neutralizing repertoire:
%   X(1901:2000)  = soluble IgM, non-neutralizing
%   X(2101:2200)  = soluble IgG, non-neutralizing
%
% ------------------------------------------------------------------------
% Inputs:
%   t       : ODE output time vector
%   y       : ODE output state matrix
%   t_grid  : uniform time grid for interpolation
%   p       : parameter structure, must contain p.N
%
% Output:
%   dist    : structure containing affinity-resolved distributions
% ========================================================================

%% ------------------------------------------------------------------------
% Basic checks
% -------------------------------------------------------------------------
if nargin < 4
    error('compute_affinity_distribution_two_epitopes requires inputs: t, y, t_grid, p');
end

if ~isfield(p, 'N')
    error('Parameter structure p must contain field p.N');
end

% Make sure t_grid is a row vector for consistent output
if iscolumn(t_grid)
    t_grid = t_grid';
end

nT = length(t_grid);

% Helper: convert X index to y-column index
% Since X(1) = y(:, p.N + 4), X(k) = y(:, p.N + 3 + k)
idxX = @(k) p.N + 3 + k;

%% ------------------------------------------------------------------------
% Preallocate 100-clone antibody distributions
% Rows: 100 antibody affinity clones
% Columns: interpolated time points
% -------------------------------------------------------------------------
Neu_IgM_100    = zeros(100, nT);
Neu_IgG_100    = zeros(100, nT);
NonNeu_IgM_100 = zeros(100, nT);
NonNeu_IgG_100 = zeros(100, nT);

%% ------------------------------------------------------------------------
% Extract and interpolate soluble antibodies
% -------------------------------------------------------------------------
for i = 1:100
    % ---------------------------------------------------------------------
    % First epitope: RBD-specific / neutralizing antibodies
    % ---------------------------------------------------------------------
    idx_Neu_IgM = idxX(100 + i);   % X(101:200)
    idx_Neu_IgG = idxX(300 + i);   % X(301:400)

    % ---------------------------------------------------------------------
    % Second epitope: non-RBD-specific / non-neutralizing antibodies
    % ---------------------------------------------------------------------
    idx_NonNeu_IgM = idxX(1900 + i);  % X(1901:2000)
    idx_NonNeu_IgG = idxX(2100 + i);  % X(2101:2200)

    % Interpolate onto t_grid
    Neu_IgM_100(i,:) = interp1(t, y(:, idx_Neu_IgM), ...
                               t_grid, 'linear', 'extrap');

    Neu_IgG_100(i,:) = interp1(t, y(:, idx_Neu_IgG), ...
                               t_grid, 'linear', 'extrap');

    NonNeu_IgM_100(i,:) = interp1(t, y(:, idx_NonNeu_IgM), ...
                                  t_grid, 'linear', 'extrap');

    NonNeu_IgG_100(i,:) = interp1(t, y(:, idx_NonNeu_IgG), ...
                                  t_grid, 'linear', 'extrap');
end

% Avoid tiny negative interpolation artifacts
Neu_IgM_100    = max(Neu_IgM_100, 0);
Neu_IgG_100    = max(Neu_IgG_100, 0);
NonNeu_IgM_100 = max(NonNeu_IgM_100, 0);
NonNeu_IgG_100 = max(NonNeu_IgG_100, 0);

%% ------------------------------------------------------------------------
% Aggregate 100 antibody clones into 19 K_D affinity groups
% -------------------------------------------------------------------------
% This follows the same grouping rule as your original
% compute_affinity_distribution.m:
%
%   kd_index = fix((i - 1)/10) - mod(i - 1, 10) + 10;
%
% This maps 100 antibody classes into 19 affinity groups.
% The corresponding K_D values are labeled as:
%
%   -31, -30, ..., -13
% -------------------------------------------------------------------------

nKD = 19;
kd_values = -31:1:-13;

Neu_IgM_kd    = zeros(nKD, nT);
Neu_IgG_kd    = zeros(nKD, nT);
NonNeu_IgM_kd = zeros(nKD, nT);
NonNeu_IgG_kd = zeros(nKD, nT);

for i = 1:100
    kd_index = fix((i - 1)/10) - mod(i - 1, 10) + 10;

    if kd_index < 1 || kd_index > nKD
        error('Computed kd_index = %d is outside valid range 1:%d', kd_index, nKD);
    end

    Neu_IgM_kd(kd_index,:)    = Neu_IgM_kd(kd_index,:)    + Neu_IgM_100(i,:);
    Neu_IgG_kd(kd_index,:)    = Neu_IgG_kd(kd_index,:)    + Neu_IgG_100(i,:);
    NonNeu_IgM_kd(kd_index,:) = NonNeu_IgM_kd(kd_index,:) + NonNeu_IgM_100(i,:);
    NonNeu_IgG_kd(kd_index,:) = NonNeu_IgG_kd(kd_index,:) + NonNeu_IgG_100(i,:);
end

%% ------------------------------------------------------------------------
% Combined antibody populations
% -------------------------------------------------------------------------
% For primary infection with ancestral virus, the key comparison is:
%
%   RBD-specific neutralizing antibodies
%       vs
%   non-RBD-specific non-neutralizing antibodies
%
% Therefore we combine IgM + IgG for each epitope group.
% -------------------------------------------------------------------------
Neu_total_kd    = Neu_IgM_kd    + Neu_IgG_kd;
NonNeu_total_kd = NonNeu_IgM_kd + NonNeu_IgG_kd;

Neu_total    = sum(Neu_total_kd, 1);
NonNeu_total = sum(NonNeu_total_kd, 1);

Neu_IgM_total    = sum(Neu_IgM_kd, 1);
Neu_IgG_total    = sum(Neu_IgG_kd, 1);
NonNeu_IgM_total = sum(NonNeu_IgM_kd, 1);
NonNeu_IgG_total = sum(NonNeu_IgG_kd, 1);

%% ------------------------------------------------------------------------
% Fold expansion relative to baseline
% -------------------------------------------------------------------------
% Use a small baseline floor to avoid division by zero.
% If initial abundance is extremely small, this prevents numerical problems.
% -------------------------------------------------------------------------
baseline_floor = 1;

Neu_fold = Neu_total ./ max(Neu_total(1), baseline_floor);
NonNeu_fold = NonNeu_total ./ max(NonNeu_total(1), baseline_floor);

Neu_IgM_fold = Neu_IgM_total ./ max(Neu_IgM_total(1), baseline_floor);
Neu_IgG_fold = Neu_IgG_total ./ max(Neu_IgG_total(1), baseline_floor);

NonNeu_IgM_fold = NonNeu_IgM_total ./ max(NonNeu_IgM_total(1), baseline_floor);
NonNeu_IgG_fold = NonNeu_IgG_total ./ max(NonNeu_IgG_total(1), baseline_floor);

%% ------------------------------------------------------------------------
% Optional summary metrics
% -------------------------------------------------------------------------
[Neu_peak, Neu_peak_idx] = max(Neu_total);
[NonNeu_peak, NonNeu_peak_idx] = max(NonNeu_total);

Neu_peak_time = t_grid(Neu_peak_idx);
NonNeu_peak_time = t_grid(NonNeu_peak_idx);

peak_ratio_NonNeu_to_Neu = NonNeu_peak / max(Neu_peak, eps);

%% ------------------------------------------------------------------------
% Store outputs
% -------------------------------------------------------------------------
dist = struct();

% Time and affinity axes
dist.t_grid = t_grid;
dist.kd_values = kd_values;

% 100-clone raw antibody distributions
dist.Neu_IgM_100 = Neu_IgM_100;
dist.Neu_IgG_100 = Neu_IgG_100;
dist.NonNeu_IgM_100 = NonNeu_IgM_100;
dist.NonNeu_IgG_100 = NonNeu_IgG_100;

% K_D aggregated antibody distributions: 19 x nT
dist.Neu_IgM_kd = Neu_IgM_kd;
dist.Neu_IgG_kd = Neu_IgG_kd;
dist.NonNeu_IgM_kd = NonNeu_IgM_kd;
dist.NonNeu_IgG_kd = NonNeu_IgG_kd;

% Combined IgM + IgG distributions
dist.Neu_total_kd = Neu_total_kd;
dist.NonNeu_total_kd = NonNeu_total_kd;

% Total abundance over all affinity groups
dist.Neu_total = Neu_total;
dist.NonNeu_total = NonNeu_total;

dist.Neu_IgM_total = Neu_IgM_total;
dist.Neu_IgG_total = Neu_IgG_total;
dist.NonNeu_IgM_total = NonNeu_IgM_total;
dist.NonNeu_IgG_total = NonNeu_IgG_total;

% Fold expansion
dist.Neu_fold = Neu_fold;
dist.NonNeu_fold = NonNeu_fold;

dist.Neu_IgM_fold = Neu_IgM_fold;
dist.Neu_IgG_fold = Neu_IgG_fold;
dist.NonNeu_IgM_fold = NonNeu_IgM_fold;
dist.NonNeu_IgG_fold = NonNeu_IgG_fold;

% Peak information
dist.Neu_peak = Neu_peak;
dist.NonNeu_peak = NonNeu_peak;
dist.Neu_peak_time = Neu_peak_time;
dist.NonNeu_peak_time = NonNeu_peak_time;
dist.peak_ratio_NonNeu_to_Neu = peak_ratio_NonNeu_to_Neu;

end