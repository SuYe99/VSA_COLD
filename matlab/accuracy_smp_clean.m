function [f1, omi, com] = accuracy_smp_clean(varargin)
% Example: accuracy_smp_clean(1);
%
% Objective:
% Counting the abrupt and transition changes and calculating the accuracy
% assessment indices (omission/ commission/ f1 score). Recording the
% unmatched mapped breaks' and references' sample ids for each type of the
% changes.
%
% Default setting:
% The change probability threshold applied: 0.75;
% Number of consecutive oberservations: 14;
% Selected model fitting method: 5x5 cold buffer with 100% of buffer
%     deleted valid obs. rate applied.;
% Judging methods to match the mapped breaks with the references: merged
%     1-year moving window breaks judged by reference 1-year moving windows. 
% Temporal composited period for the time series: Daily.
%
% Inputs:
% yr_s: start year for the accuracy assessment process;
% yr_e: end year for the accuracy assessment process;
% tlr_t: tolerance time for reference changes' start and end time.
% vld_typ: considered types of changes (1/2=abrupt recovered;
%     3/4=abrupt unrecovered; 5/6=transition).
% NTL_lim: threshold for dark pixel removal process (nW*cm-2*sr-1).


%% get parameters from inputs
p = inputParser;
% default values.
addParameter(p, 'yr_s', 2013); 
addParameter(p, 'yr_e', 2020); 
addParameter(p, 'tlr_t', 365.25*0.5); 
addParameter(p, 'vld_typ', [1, 2, 3, 4, 5, 6]); 
addParameter(p, 'NTL_lim', 1); 

% request user's input
parse(p, varargin{:});

yr_s = p.Results.yr_s;
yr_e = p.Results.yr_e;
tlr_t = p.Results.tlr_t;
vld_typ = p.Results.vld_typ;
NTL_lim = p.Results.NTL_lim;

%% Preparing for the accuracy measurement
% Get the directories
dir_l = '/*/VIIRS_NTL/';
n_rst = 'VZACOLDResult';
dir_save = fullfile(dir_l, 'AccuracyAssessment', 'result');
if ~isfolder(dir_save)
    mkdir(dir_save)
end
filename_save = 'CalibrationSample.mat';

% Read the mannual interpreted table
load(fullfile(dir_l, 'AccuracyAssessment', 'sample_interpretation_clean.mat'), 'ref_smp');

%% Loop checking for the samples

% Initialize reference change segments and mapped breaks
num_refchg = 0;
mtch_refchg = 0;
num_break = 0;
mtch_break = 0;
skip_ids = [];
for i = 1:length(ref_smp)
    mapbreak = [];

    %% Read reference segments
    % Obtain the timing and types of the reference change segments
    tmp_refbreak_s = ref_smp(i).date_start;
    tmp_refbreak_e = ref_smp(i).date_end;
    tmp_type = ref_smp(i).chg_type;
    tmp_conf = ref_smp(i).chg_conf;

    % Exclude references in the first and last years, and the types
    % that we do not consider so far (7=gradual; 8=cyclical)
    tmp_ids_vld = (tmp_refbreak_s >= datenum(yr_s, 1, 1) & ...
        tmp_refbreak_e <= datenum(yr_e, 12, 31)) &... % within the range of years
        ismember(tmp_type, vld_typ); % within the types we will focus on
    tmp_refbreak_s = tmp_refbreak_s(tmp_ids_vld);
    tmp_refbreak_e = tmp_refbreak_e(tmp_ids_vld);
    tmp_conf = tmp_conf(tmp_ids_vld);

    % Skip samples with low confidence(=1) transition changes
    if min(tmp_conf) < 2
        continue
    end

    %% Read mapped changes.
    file_name_rec = sprintf('record_change_smp%03d', i);
    % Get the IDs of failed samples
    if ~isfile(fullfile(dir_l, n_rst, file_name_rec))
        skip_ids = [skip_ids i]; 
        continue
    end
    load(fullfile(dir_l, n_rst, file_name_rec), 'rec_cg');

    % start and end time
    t_start = [rec_cg.t_start];
    t_end = [rec_cg.t_end];
    num_fit = size(rec_cg, 2);
    for n = 1:num_fit
        % Skip if no valid break mapped
        if rec_cg(n).change_prob < 1 || rec_cg(n).t_break == 0
            continue
        end
        
        %% Consistant dark pixel removal check
        % Get the break identified VZA interval.
        n_vza = mod(rec_cg(n).change_prob, 10);

        % Model overall values before breaks
        start_value = 0.1*(rec_cg(n).coefs(1, n_vza)+rec_cg(n).coefs(2, n_vza)*t_start(n));
        end_value = 0.1*(rec_cg(n).coefs(1, n_vza)+rec_cg(n).coefs(2, n_vza)*t_end(n));
        rad_bfr = max([start_value, end_value]);
        
        % Model overall values after breaks
        rad_aft = nan;
        if n < num_fit
            start_value = 0.1*(rec_cg(n+1).coefs(1, n_vza)+rec_cg(n+1).coefs(2, n_vza)*t_start(n+1));
            end_value = 0.1*(rec_cg(n+1).coefs(1, n_vza)+rec_cg(n+1).coefs(2, n_vza)*t_end(n+1));
            rad_aft = max([start_value, end_value]);
        end

        % Predicted change magnitude value
        chg_mag = 0.1*rec_cg(n).magnitude;

        % Remove if 'low confidence' break points detected
        if rad_bfr < NTL_lim && rad_aft < NTL_lim && chg_mag < NTL_lim
            continue
        end

        %% Collect mapped breaks of the stratified models
        mapbreak = [mapbreak, rec_cg(n).t_break];
    end

    % Total number of mapped breaks and reference changes
    num_break = num_break+length(mapbreak); % Mapped breaks    
    % Skip if no change identified in the reference
    if isempty(tmp_refbreak_s)
        continue
    end
    num_refchg = num_refchg+length(tmp_refbreak_s); % Reference changes
    
    %% Summerize number of matched reference time segments and mapped breaks
    check_mtchref = zeros(length(tmp_refbreak_s), 1); 
    for n = 1:length(mapbreak)
        check_match = true;
        n_ref = find(tmp_refbreak_s <= mapbreak(n) & mapbreak(n) <= tmp_refbreak_e);

        if isempty(n_ref)
            dif_ref = [abs(tmp_refbreak_s-mapbreak(n)); abs(tmp_refbreak_e-mapbreak(n))];
            tmp_dif = min(dif_ref(:));
            [~, n_ref] = find(dif_ref == tmp_dif);
            if length(n_ref)> 1
                n_ref = n_ref(1);
            end

            check_match = tmp_dif <= tlr_t;
        end

        if check_match == false
            continue
        end

        check_mtchref(n_ref) = 1;
        mtch_break = mtch_break+1;
    end

    mtch_refchg = mtch_refchg+sum(check_mtchref == 1);
end

%% Calculate accuracy estimator rates
% Omissions
unmtch_refchg = num_refchg-mtch_refchg;
omi = unmtch_refchg/num_refchg;
% Commissions
unmtch_break = num_break-mtch_break;
com = unmtch_break/num_break;
% F1 scores
f1 = 200*(1-com)*(1-omi)/((1-com)+(1-omi));

% Save the accuracy assessment result
save(fullfile(dir_save, filename_save), 'omi', 'com', 'f1');
end
