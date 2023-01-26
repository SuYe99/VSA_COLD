function [check_file, line_torg, line_torg_raw, line_tmq, line_tvza, line_tsaa, line_tm2,  line_QFDNB, line_tm_buf2] = ...
    createMergeRowdata(tile_name, n_irows, date_start, date_end, sdate, varargin)
% Example: createStackRowdata('h32v12', 1800, 1949)

%% get parameters from inputs
p = inputParser;
% default values.
addParameter(p, 'nlines_read', 10); % default 10 lines
addParameter(p, 'num_buf', 5); % default 5 (1-5 pixel buffers)
addParameter(p, 'replace_code', 0); 
addParameter(p, 'save_code', 0); 

% request user's input
parse(p, varargin{:});
nlines_read = p.Results.nlines_read;
num_buf = p.Results.num_buf;
replace_code = p.Results.replace_code;
save_code = p.Results.save_code;

% Set the upper directories for torage
% dir_l_stack = '/shared/cn452/Tian';
% dir_l_BRDF = '/scratch/til19015';
dir_l_stack = '/shared/zhulab/Tian';
dir_l_BRDF = dir_l_stack;

% dir_l_stack = fullfile(dir_l_stack, 'VIIRS_NTL_2');
% dir_l_BRDF = fullfile(dir_l_BRDF, 'VIIRS_NTL_BT');
% dir_BRDF = fullfile(dir_l_stack, 'BlackMarble_BRDFcorrected_BT', tile_name);
% dir_BRDF = fullfile(dir_l_BRDF, 'BlackMarble_BRDFcorrected_BT', tile_name);
dir_l_stack = fullfile(dir_l_stack, 'VIIRS_NTL');
dir_l_BRDF = fullfile(dir_l_BRDF, 'VIIRS_NTL');
dir_stack = fullfile(dir_l_stack, 'Stacked', tile_name);
dir_BRDF = fullfile(dir_l_BRDF, 'BlackMarble_BRDFcorrected', tile_name);

% Get all of the dates
num_t = length(sdate);

% Get the default position limits
[row1, row2, col1, col2] = deal(1, 2400, 1, 2400);
ncols = col2-col1+1;

% The start and end saving line array irows for the map
irows = row1+nlines_read.*(0:ceil((row2-row1)/nlines_read)-1);

% For the target row array
n_row1 = irows(n_irows);
n_row2 = n_row1+nlines_read-1;
n_rst_stack = fullfile(dir_stack, 'stackdata', sprintf('%04d%04d', n_row1, n_row2));
n_rst_BRDF = fullfile(dir_BRDF, 'stackdata', sprintf('%04d%04d', n_row1, n_row2));

% make the storage folder for the lines
if ~isfolder(n_rst_BRDF)
    mkdir(n_rst_BRDF);
end

filename_new = sprintf('map_tsdata_buffer_%04d%04d_stacked_%d_%d.mat', ...
    n_row1, n_row2, date_start, date_end);

% Skip if file existed and saving is needed
if replace_code == 0 && save_code == 1
    check_file = 0;
%     var_list = {'line_torg', 'line_torg_raw', 'line_tsf', 'line_tmq', ...
%         'line_tvza', 'line_tsaa', 'line_tm2', 'line_QFDNB', 'line_m12', ...
%         'line_m15', 'line_m16', 'line_tmask', 'line_tm_buf2', 'sdate'};
    var_list = {'line_torg', 'line_torg_raw', 'line_tsf', 'line_tmq', ...
        'line_tvza', 'line_tsaa', 'line_tm2', 'line_QFDNB', 'line_tmask', ...
        'line_tm_buf2', 'sdate'};
    var_null = {'nulldata'};
    
    try
        variableInfo = who('-file', fullfile(n_rst_BRDF, filename_new));
        check_file = sum(ismember(var_list, variableInfo)) == length(var_list);
        check_file = checkfile | isequal(variableInfo, var_null);
    catch
        
    end
    
    if check_file == 1
        [line_torg, line_tmq, line_tvza, line_tsaa, line_tm2,  line_QFDNB, line_tm_buf2] = ...
            deal(0, 0, 0, 0, 0, 0, 0);
        return
    end
else
    check_file = 0;
end

% try
%     load(fullfile(n_rst_BRDF, filename_new));
%     check_file = 1;
%     return
% catch   
% end

% if isfile(fullfile(n_rst_BRDF, filename_new))
%     check_file = 1;
%     return
% end

% Preparing the default values and casts for the stack variables
line_torg_new = zeros(nlines_read, ncols, num_t, 'uint16'); % BRDF-corrected DNB
line_torg_raw_new = zeros(nlines_read, ncols, num_t, 'uint16'); % Raw DNB
line_tsf_new = zeros(nlines_read, ncols, num_t, 'uint8'); % BRDF-corrected snow flag
line_tmq_new = zeros(nlines_read, ncols, num_t, 'uint8'); % BRDF-corrected mandetory quality flag
line_tvza_new = zeros(nlines_read, ncols, num_t, 'int16'); % Raw VZA
line_tsaa_new = zeros(nlines_read, ncols, num_t, 'int16'); % Raw VAA
line_tm2_new = zeros(nlines_read, ncols, num_t, 'uint8'); % Cloud/snow/poor-quality masked valid mask for BRDF-corrected DNB
line_QFDNB_new = zeros(nlines_read, ncols, num_t, 'uint16'); % BRDF-corrected DNB quality flag
% line_m12_new = zeros(nlines_read, ncols, num_t, 'uint16'); % Thermal band m12
% line_m15_new = zeros(nlines_read, ncols, num_t, 'uint16'); % Thermal band m15
% line_m16_new = zeros(nlines_read, ncols, num_t, 'uint16'); % Thermal band m16
line_tmask_new = zeros(nlines_read, ncols, num_t, 'uint16'); % Raw VCM
line_tm_buf2_new = zeros(nlines_read, ncols, num_t, num_buf, 'uint8'); % Cloud/snow buffer (1-5 dialog pixel)

% For each of the date
for n_t = 1:num_t
    num_date = datetime(sdate(n_t), 'ConvertFrom', 'datenum', ...
        'Format', 'yyyyDDD');
    filename_row = sprintf('map_tsdata_buffer_%04d%04d_%04d%03d.mat', ...
        n_row1, n_row2, year(num_date), day(num_date, 'DOY'));
    
    % Reread row data if no Black Marble data
    if ~isfile(fullfile(n_rst_stack, filename_row))
%         createRowdata(tile_name, num_date)
        continue
    end
    
    % Continue if full nulldata
    try
        variableInfo = who('-file', fullfile(n_rst_stack, filename_row));
    catch
        num_date
        continue
%         createRowdata(tile_name, num_date, 'replace_code', 1);
%         variableInfo = who('-file', fullfile(n_rst_stack, filename_row));
    end
    
    var_error = {'errordata'};
    var_null = {'nulldata'};
    if isequal(variableInfo, var_null)
        continue
    elseif isequal(variableInfo, var_error)
        num_date
%         createRowdata(tile_name, num_date)
        continue
    end
    
    try
%         load(fullfile(n_rst_stack, filename_row), 'line_torg', 'line_torg_raw', ...
%             'line_tsf', 'line_tmq', 'line_tvza', 'line_tsaa', 'line_tm2', ...
%             'line_QFDNB', 'line_m12', 'line_m15', 'line_m16', 'line_tmask', ...
%             'line_tm_buf2');
        load(fullfile(n_rst_stack, filename_row), 'line_torg', 'line_torg_raw', ...
            'line_tsf', 'line_tmq', 'line_tvza', 'line_tsaa', 'line_tm2', ...
            'line_QFDNB', 'line_tmask', 'line_tm_buf2');
    catch
        continue
        %%%%%%%%%%%%createRowdata(tile_name, num_date, 'replace_code', 1);
        
%         load(fullfile(n_rst_stack, filename_row), 'line_torg', 'line_torg_raw', ...
%             'line_tsf', 'line_tmq', 'line_tvza', 'line_tsaa', 'line_tm2', ...
%             'line_QFDNB', 'line_m12', 'line_m15', 'line_m16', 'line_tmask', ...
%             'line_tm_buf2');
        load(fullfile(n_rst_stack, filename_row), 'line_torg', 'line_torg_raw', ...
            'line_tsf', 'line_tmq', 'line_tvza', 'line_tsaa', 'line_tm2', ...
            'line_QFDNB', 'line_tmask', 'line_tm_buf2');
    end
    
    line_torg_new(:, :, n_t) = line_torg;
    line_torg_raw_new(:, :, n_t) = line_torg_raw;
    line_tsf_new(:, :, n_t) = line_tsf;
    line_tmq_new(:, :, n_t) = line_tmq;
    line_QFDNB_new(:, :, n_t) = line_QFDNB;
%     line_m12_new(:, :, n_t) = line_m12;
%     line_m15_new(:, :, n_t) = line_m15;
%     line_m16_new(:, :, n_t) = line_m16;
    line_tsaa_new(:, :, n_t) = line_tsaa;
    line_tm2_new(:, :, n_t) = line_tm2;
    line_tmask_new(:, :, n_t) = line_tmask;
    line_tvza_new(:, :, n_t) = line_tvza;
    line_tm_buf2_new(:, :, n_t, :) = line_tm_buf2;
    
%     clear line_torg line_torg_raw ine_tsf line_tmq line_tvza line_tsaa ...
%         line_tm2 line_QFDNB line_m12 line_m15 line_m16 ine_tmask line_tm_buf2
    clear line_torg line_torg_raw ine_tsf line_tmq line_tvza line_tsaa ...
        line_tm2 line_QFDNB line_tmask line_tm_buf2
end

line_torg = line_torg_new;
line_torg_raw = line_torg_raw_new;
line_tsf = line_tsf_new;
line_tmq = line_tmq_new;
line_QFDNB = line_QFDNB_new;
% line_m12 = line_m12_new;
% line_m15 = line_m15_new;
% line_m16 = line_m16_new;
line_tsaa = line_tsaa_new;
line_tmask = line_tmask_new;
line_tvza = line_tvza_new;
line_tm2 = line_tm2_new;
line_tm_buf2 = line_tm_buf2_new;

% save(fullfile(n_rst_BRDF, [filename_new, '.tmp']), 'sdate', 'line_torg', ...
%     'line_torg_raw', 'line_tsf', 'line_tmq', 'line_tvza', 'line_tsaa', ...
%     'line_tm2', 'line_QFDNB', 'line_m12', 'line_m15', 'line_m16', ...
%     'line_tmask', 'line_tm_buf2');

% Skip saving if needed
if save_code == 0
    return
end
save(fullfile(n_rst_BRDF, [filename_new, '.tmp']), 'sdate', 'line_torg', ...
    'line_torg_raw', 'line_tsf', 'line_tmq', 'line_tvza', 'line_tsaa', ...
    'line_tm2', 'line_QFDNB', 'line_tmask', 'line_tm_buf2');
movefile(fullfile(n_rst_BRDF, [filename_new, '.tmp']), fullfile(n_rst_BRDF, ...
    filename_new))
end