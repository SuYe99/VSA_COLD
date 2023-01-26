function main_COLD_brdf_map_tc(tile_name, n_row, CI, ...
    Tmax_cg, conse, nlines_read, varargin)
% Matlab code for Continous Change Detection (for standalone version)
% Example: main_COLD_brdf_map_tc_clean('h10v04', 1, 0.75, 1, 14, 10);

%% get parameters from inputs
p = inputParser;
% default values.
addParameter(p, 'date_start', 2013001); % Date range of the data
addParameter(p, 'date_end', 2021365); 
addParameter(p, 'col1', 1); % Dimension of the map columns
addParameter(p, 'col2', 2400); 
addParameter(p, 'save_code', 1); % save the record or not
addParameter(p, 'n_buf', 2); % 5x5 buffer window 
addParameter(p, 'prc_cld', 100); % exclude all cloud/snow buffer
addParameter(p, 'num_toler', 1); % tolerate 1 obs. within the consecutive anomaly observations
addParameter(p, 'vza_lmt_all', [0, 60]); % VZA range
addParameter(p, 'vzaints', [0 20; 20 40; 40 60]); % VZA intervals 

% request user's input
parse(p, varargin{:});
date_start = p.Results.date_start;
date_end = p.Results.date_end;
col1 = p.Results.col1;
col2 = p.Results.col2;
save_code = p.Results.save_code;
n_buf = p.Results.n_buf;
prc_cld = p.Results.prc_cld;
num_toler = p.Results.num_toler;
vza_lmt_all = p.Results.vza_lmt_all;
vzaints = p.Results.vzaints;

cmd_l = '/home/til19015/Black_Marble/COLD-master/';
addpath(cmd_l);  

n_irows = ceil(n_row./nlines_read); 
% %% Tuning variables for continuous change detection
% % probability for detecting surface change
try inputs = textread('COLD_Parameters.txt'); %#ok<DTXTRD>
    % need at least three varibles
    % change threshold
    T_cg = inputs(1); % 0.99
    % number of consecutive obs
    conse = inputs(2); % 6
    % maximum number of coefficients
    num_c = inputs(3); % 8
catch me %#ok<NASGU>
    % change threshold
    T_cg = CI;
    % maximum number of coefficients
    num_c = 4;
end
%% Constants:
% get image parameters automatically
% folder name of all COLD results
tcom_name = 'dly';
n_rst = sprintf('TSFitMap_%s_%s_con3', tile_name, tcom_name);
dir_l_new = '/shared/zhulab/Tian/VIIRS_NTL/COLDResult_l_new3';
dir_scn_BRDF_new = fullfile('BlackMarble_BRDFcorrected', tile_name);
dir_BRDF_new = fullfile(dir_l_new, dir_scn_BRDF_new);

% make TSFitMap folder for storing coefficients
if ~isfolder(fullfile(dir_BRDF_new, n_rst))
    mkdir(fullfile(dir_BRDF_new, n_rst));
end

% not exist or corrupt
fprintf('Processing row%d of the %s: CI = %f, conse = %d\n', n_row, ...
    tile_name, CI, conse);

% run COLD for the target row
TrendSeasonalFit_map_COLDContinousModel_row(date_start, date_end, ...
    col1, col2, save_code, n_buf, prc_cld, num_toler, vzaints, vza_lmt_all, ...
    tile_name, n_row, T_cg, Tmax_cg, conse, nlines_read, n_irows, ...
    dir_BRDF_new, num_c, n_rst)

% profile viewer
% exit
end % end of function 
