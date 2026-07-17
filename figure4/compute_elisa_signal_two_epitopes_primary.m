function elisa = compute_elisa_signal_two_epitopes_primary(simData, elisaTimeGrid, labelText)
% compute_elisa_signal_two_epitopes_primary
%
% Compute simulated ELISA signals for neutralizing and non-neutralizing
% antibodies during primary infection using the original ELISA kinetic model.
%
% Neutralizing antibodies:
%   RBD-specific repertoire
%   IgM: 1901:2000
%   IgG: 2101:2200
%
% Non-neutralizing antibodies:
%   non-RBD-specific repertoire
%   IgM: 101:200
%   IgG: 301:400
%
% ELISA settings:
%   antigen amount       = 1e15
%   antibody dilution    = 100-fold
%   assay duration       = 2
%   ELISA readout        = 1e15 - antigen remaining at t = 2

if nargin < 3
    labelText = 'Primary infection';
end

t = simData.t(:);
y = simData.y;

elisaTimeGrid = elisaTimeGrid(:);

% -------------------------------------------------------------------------
% Extract antibody repertoires using the correct primary-infection indices
% -------------------------------------------------------------------------
ab = extract_two_epitope_antibody_matrices(t, y);

% Interpolate each antibody matrix to ELISA time grid
neuIgM_t    = interpolate_antibody_matrix(t, ab.neu.IgM, elisaTimeGrid);
neuIgG_t    = interpolate_antibody_matrix(t, ab.neu.IgG, elisaTimeGrid);
nonNeuIgM_t = interpolate_antibody_matrix(t, ab.nonNeu.IgM, elisaTimeGrid);
nonNeuIgG_t = interpolate_antibody_matrix(t, ab.nonNeu.IgG, elisaTimeGrid);

% -------------------------------------------------------------------------
% Original ELISA kinetic parameters
% -------------------------------------------------------------------------
[para, para_new, para_new_1] = get_elisa_parameters();

% -------------------------------------------------------------------------
% Compute ELISA signal at each sampled infection time
% -------------------------------------------------------------------------
nT = numel(elisaTimeGrid);

neuIgMSignal    = zeros(nT, 1);
neuIgGSignal    = zeros(nT, 1);
nonNeuIgMSignal = zeros(nT, 1);
nonNeuIgGSignal = zeros(nT, 1);

for k = 1:nT
    neuIgMSignal(k) = run_single_elisa_assay( ...
        squeeze(neuIgM_t(k, :, :)), ...
        zeros(10, 10), ...
        para, para_new, para_new_1);

    neuIgGSignal(k) = run_single_elisa_assay( ...
        zeros(10, 10), ...
        squeeze(neuIgG_t(k, :, :)), ...
        para, para_new, para_new_1);

    nonNeuIgMSignal(k) = run_single_elisa_assay( ...
        squeeze(nonNeuIgM_t(k, :, :)), ...
        zeros(10, 10), ...
        para, para_new, para_new_1);

    nonNeuIgGSignal(k) = run_single_elisa_assay( ...
        zeros(10, 10), ...
        squeeze(nonNeuIgG_t(k, :, :)), ...
        para, para_new, para_new_1);
end

% -------------------------------------------------------------------------
% Package output
% -------------------------------------------------------------------------
elisa.label = labelText;
elisa.time = elisaTimeGrid;

elisa.neu.IgM = neuIgMSignal;
elisa.neu.IgG = neuIgGSignal;
elisa.neu.total = neuIgMSignal + neuIgGSignal;

elisa.nonNeu.IgM = nonNeuIgMSignal;
elisa.nonNeu.IgG = nonNeuIgGSignal;
elisa.nonNeu.total = nonNeuIgMSignal + nonNeuIgGSignal;

end

%% ========================================================================
% Run one ELISA assay
% ========================================================================

function elisaValue = run_single_elisa_assay( ...
    IgM_matrix, IgG_matrix, para, para_new, para_new_1)

IgM_matrix = max(IgM_matrix, 0);
IgG_matrix = max(IgG_matrix, 0);

IgM_matrix = reshape_to_10_by_10(IgM_matrix);
IgG_matrix = reshape_to_10_by_10(IgG_matrix);

x0 = zeros(1402, 1);

% -------------------------------------------------------------------------
% Antibody dilution before ELISA assay
% -------------------------------------------------------------------------
antibodyDilution = 100;

% In the original ELISA model:
%   101:200    IgM antibody species
%   301:400    IgG antibody species
%   1402       free antigen
%
% Antibodies are diluted 100-fold before entering the ELISA model.
x0(101:200) = IgM_matrix(:) / antibodyDilution;
x0(301:400) = IgG_matrix(:) / antibodyDilution;

% ELISA antigen amount
x0(1402) = 1e15;

odeOptions = odeset( ...
    'RelTol', 1e-6, ...
    'AbsTol', 1e-9, ...
    'NonNegative', 1:1402);

try
    [tElisa, zz] = ode15s( ...
        @(tt, yy) pathway_model_many_antibody_immune_include_plasma_elisa( ...
            tt, yy, para, para_new, para_new_1), ...
        [0 2], ...
        x0, ...
        odeOptions);

    virusAtEnd = interp1(tElisa, zz(:, 1402), 2);
    elisaValue = 1e15 - virusAtEnd;

    if ~isfinite(elisaValue) || elisaValue < 0
        elisaValue = 0;
    end

catch
    elisaValue = NaN;
end

end

%% ========================================================================
% Extract antibody matrices from primary infection simulation
% ========================================================================

function ab = extract_two_epitope_antibody_matrices(t, y)

nT = numel(t);

% Full model state:
%   y = [T; I_1...I_N; V; Tc; X_1...X_3604; cumulative_death]
%
% For p.N = 400:
%   X starts at y column 404.
xOffset = 403;

idxNeuIgM    = xOffset + (1901:2000);
idxNeuIgG    = xOffset + (2101:2200);
idxNonNeuIgM = xOffset + (101:200);
idxNonNeuIgG = xOffset + (301:400);

requiredMaxIndex = max([idxNeuIgM, idxNeuIgG, idxNonNeuIgM, idxNonNeuIgG]);

if size(y, 2) < requiredMaxIndex
    error(['Simulation output y has only %d columns, but full-state ELISA ', ...
           'extraction requires at least %d columns. If y contains only X, ', ...
           'use the no-offset index version instead.'], ...
           size(y, 2), requiredMaxIndex);
end

ab.neu.IgM = reshape_time_series(y(:, idxNeuIgM), nT);
ab.neu.IgG = reshape_time_series(y(:, idxNeuIgG), nT);

ab.nonNeu.IgM = reshape_time_series(y(:, idxNonNeuIgM), nT);
ab.nonNeu.IgG = reshape_time_series(y(:, idxNonNeuIgG), nT);

end

%% ========================================================================
% ELISA parameters
% ========================================================================

function [para, para_new, para_new_1] = get_elisa_parameters()

para = zeros(20, 1);

para(1) = 1e-20;
para(2) = 0.5;
para(3) = 2.2e17 + 1e13;
para(4) = 0.5e13;
para(5) = 0.01;
para(6) = 0.005;
para(7) = 1;
para(8) = 5e-7;
para(9) = 0.1;
para(10) = 0.05;
para(11) = 4.4e9;
para(12) = 1e9;
para(13) = 0.1;
para(14) = 0.05;
para(15) = 2e-7;
para(16) = 0.02;
para(17) = 0.1;
para(18) = 0.5;
para(19) = 0;
para(20) = 1e5;

para_new = zeros(10, 1);
para_new(1) = 1e-22;
para_new(2) = 1e-21;
para_new(3) = 1e-20;
para_new(4) = 1e-19;
para_new(5) = 1e-18;
para_new(6) = 1e-17;
para_new(7) = 1e-16;
para_new(8) = 1e-15;
para_new(9) = 1e-14;
para_new(10) = 1e-13;

para_new_1 = zeros(10, 1);
para_new_1(1) = 1e0;
para_new_1(2) = 1e1;
para_new_1(3) = 1e2;
para_new_1(4) = 1e3;
para_new_1(5) = 1e4;
para_new_1(6) = 1e5;
para_new_1(7) = 1e6;
para_new_1(8) = 1e7;
para_new_1(9) = 1e8;
para_new_1(10) = 1e9;

end

%% ========================================================================
% Reshape time x 100 into time x 10 x 10
% ========================================================================

function out = reshape_time_series(block, nT)

block = max(block, 0);

if size(block, 2) ~= 100
    error('Each antibody block must contain 100 species.');
end

out = zeros(nT, 10, 10);

for k = 1:nT
    out(k, :, :) = reshape(block(k, :), 10, 10);
end

end

%% ========================================================================
% Interpolate time x 10 x 10 antibody matrices
% ========================================================================

function matOut = interpolate_antibody_matrix(t, matIn, tGrid)

nGrid = numel(tGrid);
matOut = zeros(nGrid, 10, 10);

flatIn = reshape(matIn, numel(t), 100);
flatOut = interp1(t, flatIn, tGrid, 'linear', 'extrap');
flatOut = max(flatOut, 0);

for k = 1:nGrid
    matOut(k, :, :) = reshape(flatOut(k, :), 10, 10);
end

end

%% ========================================================================
% Force 10 x 10
% ========================================================================

function mat = reshape_to_10_by_10(mat)

if isequal(size(mat), [10 10])
    return;
end

if numel(mat) ~= 100
    error('ELISA antibody input must contain 100 values.');
end

mat = reshape(mat, 10, 10);

end