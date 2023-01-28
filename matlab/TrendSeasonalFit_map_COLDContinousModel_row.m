function TrendSeasonalFit_map_COLDContinousModel_row(date_start, date_end, ...
    col1, col2, save_code, n_buf, prc_cld, num_toler, vzaints, vza_lmt_all, ...
    tile_name, n_row, T_cg, Tmax_cg, conse, nlines_read, n_irows, ...
    dir_BRDF_new, num_c, n_rst)

warning ('off', 'all');

num_vza = size(vzaints, 1);
if num_vza == 1
    num_mod = 1;
else
    num_mod = num_vza+1;
end

min_num_c = num_c;
max_num_c = min_num_c;

% number of clear observation / number of coefficients
n_times = 12; % initial value for 16-day Landsat imagery
% number of days per year
num_yrs = 365.25;

% minimum year for model intialization
mini_yrs = 1;
% threshold (degree) of mean included angle
nsign = 45;

% consecutive number
def_conse = conse;

% Tmasking of noise
Tmax_cg = norminv(1-(1-Tmax_cg)/2);
% adjust threshold based on chi-squared distribution
def_pT_cg = T_cg;
def_T_cg = norminv(1-(1-def_pT_cg)/2);

% check if row already exist
filename_save = sprintf('%s_record_change%04d.mat', tile_name, n_row);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Get ready for Xs & Ys

% The year and doy for the start and end time
yr_s = floor(date_start/1000);
doy_s = mod(date_start, 1000);
yr_e = floor(date_end/1000);
doy_e = mod(date_end, 1000);
date_start = datenum(yr_s, 1, doy_s);
date_end = datenum(yr_e, 1, doy_e);
sdate = (date_start: date_end)';

[~, line_torg, ~, line_tmq, line_tvza, line_tsaa, line_tm2,  line_QFDNB, line_tm_buf2] = ...
    createMergeRowdata(tile_name, n_irows, date_start, date_end, sdate);
    
% Trial data:    
% load('/scratch/zhz18039/til19015/share_c/trail_data.mat', 'line_torg', 'line_tmq', 'line_tvza', 'line_tsaa', 'line_tm2',  'line_QFDNB', 'line_tm_buf2')

% Types of data normalizartion used
line_t = line_torg;
line_tm_buf = line_tm_buf2(:, :, :, n_buf);

irow_ids = n_row-(n_irows-1)*nlines_read;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Get ready for Xs & Ys
% initialize NUM of Functional Curves
num_fc = 0;
% initialize the struct data of RECording of ChanGe (rec_cg)
rec_cg = struct('t_start', [], 't_end', [], 't_break', [], 'coefs', [], ...
    'rmse', [], 'pos', [], 'change_prob', [], 'num_obs', [], 'category', ...
    [], 'magnitude', [], 'adj_rmse', [], 'rec_rmse', [], 'rec_normdif', []);
prc_bufferexclude = zeros(col2-col1+1, 1);

% Loop COLD for the columns
for icol_ids = col1:col2
    % tic

    qa = 0;
    fit_cft = nan(max_num_c, num_mod);
    rmse = nan(1, num_mod);
    
    % get default conse & T_cg
    conse = def_conse;
    T_cg = def_T_cg;
    
    % Recording the all data records
    x_all = sdate;
    y_all = line_t(irow_ids, icol_ids, :);
    y_all = y_all(:);
    MQF_all = line_tmq(irow_ids, icol_ids, :);
    QFDNB_all = line_QFDNB(irow_ids, icol_ids, :);
    vza_all = line_tvza(irow_ids, icol_ids, :);
    vza_all = vza_all(:);
    saa_all = line_tsaa(irow_ids, icol_ids, :);
    saa_all = saa_all(:);
    
    % Selecting the cloud & cloud buffer mask based on the buffer
    % deleted valid data rate
    mask_all1 = line_tm2(irow_ids, icol_ids, :);
    mask_all1 = mask_all1(:);
    
    mask_all2 =  line_tm_buf(irow_ids, icol_ids, :);
    mask_all2 = mask_all2(:);
    
    num_vld = sum(line_t(irow_ids, icol_ids, :) < 65535 ...
        & line_t(irow_ids, icol_ids, :) > 0);
    num_vld_buf = sum(line_tm_buf(irow_ids, icol_ids, :));
    num_bufdlt = num_vld-num_vld_buf;
    prc_bufferexclude(icol_ids) = num_bufdlt/num_vld;
    
    if prc_bufferexclude(icol_ids) > prc_cld
        mask_all =  mask_all1;
    else
        mask_all =  mask_all1 & mask_all2;
    end
    
    idrange_vzaint = zeros(length(vza_all), num_vza);
    % clear pixels
    idclr = zeros(length(vza_all), 1);
    
    for nvza = 1:num_vza
        vza_lmt = [vzaints(nvza, 1), vzaints(nvza, 2)];
        
        % Obtaining the valid obs. ids for the modeling
        [~, idrange_vzaint(:, nvza), idrange_all, idclr, ~, ~] = ...
            datarange_brdf(x_all, y_all, MQF_all, QFDNB_all, date_start, ...
            date_end, mask_all, vza_all, vza_lmt, vza_lmt_all, saa_all, ...
            line_tm2, irow_ids, icol_ids);
    end   
    
    %% Case 1: All time series data do not have enough clear observations for change detection
    if sum(idclr) < n_times*max_num_c % not enough for all data
        % TBD, but create a file to record
        %             save(fullfile(dir_BRDF_new, n_rst, filename_save), 'filename_save');
        continue;
    end   
    
    %% Case 2: Stratified zenith angle data to detect changes
    
    % clear and within physical range pixels based on all time series data
    idgood = idrange_all & idclr; % Modified
    
    % Xs & Ys for computation after removing noises such as cloud
    [clrx, clry, clrvza, ~] = find_ids(sdate, line_t, line_tvza, line_tsaa, ...
        irow_ids, icol_ids, idgood, 1);
    
    % continue if not enough clear pixels, all data time series
    if length(clrx) < n_times*min_num_c
        continue;
    end
    
    % find the vza intervals of each clear obs
    clrvzaints = zeros(size(clrx)); % 1, 2, 3 means different angle range's data
    for nvza = 1:num_vza
        index_tmp = clrvza >= vzaints(nvza, 1)*100 & clrvza < vzaints(nvza, 2)*100;
        clrvzaints = clrvzaints+nvza*double(index_tmp);
    end
    
    % !!! Start here !!!
    % !!!
    % !!!
    % !!!
    
    % one variable, all data list, which will benefit code organizing
    clr_data = [clrx, clry, clrvzaints];
    rec_normdif  = zeros(size(clr_data, 1), num_vza);
    rec_rmse  = zeros(size(clr_data, 1), num_vza, 2); % [mod_rmse, tmp_rmse]
    
    % adjust RMSE for each stratifed data
    adj_rmse = zeros(num_mod, 1); % +1 means all data
    for nvza = 1:num_vza
        clry_vza = clr_data(clr_data(:, end) == nvza, 2);
        if length(clry_vza) < 2
            continue
        end
        adj_rmse(nvza) = median(abs(clry_vza(2:end)-clry_vza(1:end-1)), 1);
    end
    
    % adjust RMSE for all data for Tmask
    adj_rmse_all = median(abs(clr_data(2:end, 2)-clr_data(1:end-1, 2)), 1);
    adj_rmse(end) = adj_rmse_all;
    clear clry_vza;
    
    % the first observation for TSFit for each zenith angle data
    i_start = ones(num_mod, 1); % the all data is 1
    for nvza = 1:num_vza
        try
            i_start(nvza) = find(clr_data(:, end) == nvza, 1);
        catch
        end
    end
    
    % record the start of the model initialization (0=>initial;1=>done)
    BL_train = zeros(num_mod, 1);
    i_count = zeros(num_mod, 1);
    % identified and move on for the next curve
    num_fc = num_fc+1; % NUM of Fitted Curves (num_fc)
    % record the num_fc at the beginning of each pixel
    rec_fc = num_fc;
    % initialize i_dense
    i_dense = ones(num_mod, 1);
    % start with the miminum requirement of clear obs for all data
    i = n_times*min_num_c; % first i  goes further
    
    % record the difference of conse observations, all records same format
    rec_dif = nan(4, conse, num_mod);
    % dimension 1: v_dif_mag % record the magnitude of change
    % dimension 2: v_dif     % value of difference for conse obs
    % dimension 3: vec_mag   % chagne vector magnitude
    % dimension 4: index of observation in orginal datalist
    
    % while loop - process till the last clear observation - conse
    while i <= size(clr_data, 1)-conse
        
        % find the ID of vza interval of the i's observation
        a_id = num_mod; % indicates all data
        if num_mod == 1
            n_id = a_id;
        else
            vza_id = clr_data(i, end); % indicates zenith data's demension in "clr_data"
            n_id = [a_id vza_id];
        end
        
        changeconfirmed = 0; % label change confirmed, and there is no change at default
        
        % "all data" and "zenith angle data" will be checked at the same time
        for n = n_id % n indicates the ID of stratified zenith data, and all data has higher priority
            
            % if change has been confirmed by all data, do not repeat zenith angle data
            if changeconfirmed > 0
                continue;
            end
            
            % set input dataset all the same between "all data" and
            % "zenith angle data", which can share same process code below
            if n > num_vza
                clr_data_n = clr_data;
                clr_data_n(:, end) = n; % all in the 4th stratifed data
            else
                clr_data_n = clr_data; % do not change the orginal zenith angle ID
            end
            
            % remaining observations: if remaining observations less than conse, next process
            if sum(clr_data_n(i+1:end, end) == n) < conse
                continue;
            end
            
            % basic requrirements: 1) enough observations; 2) enough time
            if ~basicRequire(clr_data_n, i_start, i, n, min_num_c, mini_yrs, n_times, num_yrs)
                continue; % if it fails, we will move to next observation
            end
            
            % initializing model
            if BL_train(n) == 0
                [clrx_n, clry_n] = findStratifiedData(clr_data_n, n, i_start(n), i);  % select the stratified data x, accordingly
                % check max time difference
                if max(clrx_n(2:end)-clrx_n(1:end-1)) > num_yrs % max time difference
                    i_start(n)  = moveStartIndex(clr_data_n, i_start, n); % start from next clear obs ONLY from the stratified data list
                    i_dense(n) = i_start(n); % i that is dense enough
                    continue;
                end

                % initial model fit
                [fit_cft(:, n), rmse(n), rec_v_dif, ~] = autoRobustFit2(clrx_n, clry_n, min_num_c);
                
                % record the ith RMSEs and normalized differences for check
                rec_rmse(i, n, 1)  = rmse(n); % [mod_rmse, tmp_rmse]
                [isstable, rec_normdif(i, n)] = isStableModel(rec_v_dif, fit_cft(:, n), adj_rmse(n), rmse(n), clrx_n, T_cg);
                
                % stable model requirement
                if isstable
                    BL_train(n) = 1;
                else
                    fit_cft(:, n) = nan(max_num_c, 1); % withdraw model
                    rmse(n) = nan; % withdraw model's rmse
                    i_start(n)  = moveStartIndex(clr_data_n, i_start, n); % start from next clear obs
                    continue;
                end
                
                % backward change detection if model ready
                i_count(n) = 0;
                
                % find the previous break point
                if num_fc == rec_fc
                    i_break = 1; % first curve
                else
                    i_break = find(clr_data_n(:, 1) >= rec_cg(num_fc-1).t_break, 1); % after the first curve
                end
                
                % backward change detection
                changeconfirmed_back = 0;
                if i_start(n) > i_break % if observations are available between i_break and i_start(n)
                    % loop observations from i_start(n) -1 back to i_break to confirm its change
                    % find out the remaining conse observations
                    ids_check = find(clr_data_n(:, end) == n);
                    ids_check = ids_check(i_break <= ids_check & ids_check <= i_start(n)-1); % back to real ID in all datalist, this code can select the real IDs that we need to check
                    if ~isempty(ids_check) % sometimes no data
                        ids_check = ids_check(end:-1:1); % backfarward e.g., 8, 7, 6, 5, 4, 3, 2, 1
                        for itmp = 1:length(ids_check) % loop the observation
                            i_ini = ids_check(itmp); % real ID in all datalist
                            ids_conse = ids_check(itmp:length(ids_check)); % observations that will be compared
                            ids_conse = ids_conse(1:min(conse, length(ids_conse))); % not exceed conse
                            
                            % compare observation and predication
                            rec_dif_back = diffConseObs(n, clr_data_n, fit_cft, rmse, adj_rmse, ids_conse, num_vza, conse);
                           
                            % confirm change
                            rec_rmse(i, n, 1)  = rmse(n); % [mod_rmse, tmp_rmse]
                            [changeconfirmed_back, v_dif_mag, rec_normdif(i, n)] = ...
                                confirmChange(n, rec_dif_back, num_toler, T_cg, nsign, Tmax_cg);
                            
                            if changeconfirmed_back > 0
                                break; % stop detecting change backward
                            elseif changeconfirmed_back == -1
                                clr_data(i_ini, :) = []; % remove noise
                                rec_normdif(i_ini, :) = [];
                                rec_rmse(i_ini, :, :) = []; 
                                
                                rec_dif(4, :, :) = rec_dif(4, :, :)-1; % remove index in the record of difference in the orginal datalist
                                i = i-1; % stay & check again after noise removal
                            end
                            % update new_i_start if i_ini is not a confirmed break
                            i_start(n) = i_ini; % may other stratified data varied
                            %% end of confirming change
                        end
                    end
                end
                
                
                % fit the fisrt curve if enough observations remaining
                % only fit first curve if 1) have more than
                % conse obs 2) previous obs is less than a year
                
                if changeconfirmed_back > 0 && (num_fc == rec_fc || clr_data_n(i_start(n)-1, 1)-clr_data_n(i_dense(n), 1) > num_yrs)
                    % find out the remaining conse observations
                    ids_check = find(clr_data_n(:, end) == n);
                    if sum(i_dense(n) <= ids_check & ids_check <= i_start(n)) >= conse
                        
                        fit_cft_disturb = nan(max_num_c, num_mod);
                        rmse_disturb = nan(1, num_mod);
                        qa = 10;
                        
                        [fit_cft_disturb(:, end), rmse_disturb(end)] = ...
                            autoRobustFit2(clr_data_n(i_dense(n):i_start(n)-1, 1), clr_data_n(i_dense(n):i_start(n)-1, 2), min_num_c);
                        
                        % check the data at each stratified
                        for nvza =1:num_vza
                            [clrx_n, clry_n] = findStratifiedData(clr_data_n, nvza, i_dense(n), i_start(n)-1);
                            if length(clrx_n) < conse % if less than conse obs, no model fitting
                                continue;
                            end
                            
                            [fit_cft_disturb(:, nvza), rmse_disturb(nvza)] = autoRobustFit2(clrx_n, clry_n, min_num_c);
                        end
                        % if model available, just record it or them
                        if sum(~isnan(fit_cft_disturb(:))) > 0
                            % record time of curve end
                            rec_cg(num_fc).t_end = clr_data_n(i_start(n)-1, 1);
                            % record postion of the pixel
                            rec_cg(num_fc).pos = icol_ids;
                            % record fitted coefficients
                            rec_cg(num_fc).coefs = fit_cft_disturb;
                            % record rmse of the pixel
                            rec_cg(num_fc).rmse = rmse_disturb;
                            % record break time
                            rec_cg(num_fc).t_break = clr_data(i_start(n), 1);
                            % record change probability
                            rec_cg(num_fc).change_prob = changeconfirmed_back;
                            rec_cg(num_fc).t_start = clr_data(i_dense(n), 1);
                            % record fit category
                            rec_cg(num_fc).category = qa+min_num_c;
                            % record number of observations
                            rec_cg(num_fc).num_obs = i_start(n)-i_dense;
                            % record change magnitude
                            rec_cg(num_fc).magnitude = -median(v_dif_mag, 'omitnan');
                            
                            % record the RMSEs and normalized differences
                            rec_cg(num_fc).adj_rmse = adj_rmse;
                            rec_cg(num_fc).rec_rmse = rec_rmse;
                            rec_cg(num_fc).rec_normdif = rec_normdif;
                            
                            % identified and move on for the next functional curve
                            num_fc = num_fc+1;
                            
                            i_start(:) = i_start(n); % unique i_start because of new change is confirmed
                            i_dense(:) = i; % update to record i_dense
                        end
                        clear fit_cft_disturb rmse_disturb;
                    end
                end
            end % end of inialization
            
            % continous change detection
            [clrx_n, clry_n] = findStratifiedData(clr_data_n, n, i_start(n), i);
            update_num_c = max_num_c;
            % initial model fit when there are not many obs
            if  i_count(n) == 0 || length(clrx_n) <= max_num_c*n_times
                % update i_count at each interation
                i_count(n) = sum(clr_data_n(i_start(n):i, end) == n);
                % update i_count at each interation               
                [fit_cft(:, n), rmse(n), ~] = autoRobustFit2(clrx_n, clry_n, max_num_c); % update the time series model                
                
                tmpcg_rmse = 0;
                rec_rmse(i, n, 2)  = tmpcg_rmse;
                
                %% start of updating record
                % record time of curve start
                rec_cg(num_fc).t_start = clr_data(min(i_start), 1);
                % record time of curve end
                rec_cg(num_fc).t_end = clr_data(i, 1);
                % record break time
                rec_cg(num_fc).t_break = 0; % no break at the moment
                % record postion of the pixel
                rec_cg(num_fc).pos = icol_ids;
                % record fitted coefficients
                rec_cg(num_fc).coefs = fit_cft;
                % record rmse of the pixel
                rec_cg(num_fc).rmse = rmse;
                % record change probability
                rec_cg(num_fc).change_prob = 0;
                % record number of observations
                rec_cg(num_fc).num_obs = i-min(i_start)+1;
                % record fit category
                rec_cg(num_fc).category = qa+0; % no change at default
                % record change magnitude
                rec_cg(num_fc).magnitude = 0;
                
                % record the RMSEs and normalized differences
                rec_cg(num_fc).adj_rmse = adj_rmse;
                rec_cg(num_fc).rec_rmse = rec_rmse;
                rec_cg(num_fc).rec_normdif = rec_normdif;

                %% end of updating record
                
                % find out the remaining conse observations
                ids_conse = find(clr_data_n(:, end) == n);
                ids_conse = ids_conse(find(ids_conse > i, conse)); % the ID after i at the current dataset
                
                if length(ids_conse) == conse % enough conse data
                    % compare observation and predication
                    rec_dif_n = diffConseObs(n, clr_data_n, fit_cft, rmse, adj_rmse, ids_conse, num_vza, conse);
                    rec_dif(:, :, n) = rec_dif_n(:, :, n); % update to record
                end
            else
                % update i_count at each interation
                i_count(n) = sum(clr_data_n(i_start(n):i, end) == n);
                % update the model fitting
                [fit_cft(:, n), rmse(n)] = autoRobustFit2(clrx_n, clry_n, update_num_c);
                
                % record fitted coefficients
                rec_cg(num_fc).coefs = fit_cft;
                % record rmse of the pixel
                rec_cg(num_fc).rmse = rmse;
                % record number of observations
                rec_cg(num_fc).num_obs =  i-min(i_start)+1;
                % record fit category
                rec_cg(num_fc).category = qa+0;
                % record time of curve end
                rec_cg(num_fc).t_end = clr_data_n(i, 1);
                
                % record the RMSEs and normalized differences
                rec_cg(num_fc).adj_rmse = adj_rmse;
                rec_cg(num_fc).rec_rmse = rec_rmse;
                rec_cg(num_fc).rec_normdif = rec_normdif;
                
                % move the ith col to i-1th col at zenith angle data
                i_end_n = rec_dif(4, end, n); % record the index of observation
                % find out next conse observation in the zenith angle data
                rec_dif(:, 1:end-1, n) = rec_dif(:, 2:end, n);
                rec_dif(:, end, n) = nan;
                if ~isnan(i_end_n) % no comparation at last
                    i_end_n = find(clr_data_n(:, 1) > clr_data_n(i_end_n, 1) & clr_data_n(:, end) == n, 1); % move the last observation in this zenith angle datalist to the next one
                end
                if ~isempty(i_end_n) && ~isnan(i_end_n) % to the end
                    % all data new comparision
                    % absolute difference
                    rec_dif(1, end, n) = clr_data_n(i_end_n, 2)-autoTSPred_OLS(clr_data_n(i_end_n, 1), fit_cft(:, clr_data_n(i_end_n, 3))); % i_conse may be different from n
                    % minimum rmse
                    mini_rmse = max([adj_rmse(clr_data_n(i_end_n, 3)), rmse(clr_data_n(i_end_n, 3)), tmpcg_rmse]); %  must be the same as i + conse
                    % z-scores
                    rec_dif(2, end, n) = rec_dif(1, end, n)/mini_rmse;
                    rec_dif(3, end, n) = norm(rec_dif(2, end, n)); %^2;
                    rec_dif(4, end, n) = i_end_n;
                end
            end
            
            % confirm change
            rec_rmse(i, n, 1)  = rmse(n); % [mod_rmse, tmp_rmse]
            [changeconfirmed, v_dif_mag, rec_normdif(i, n)] = ...
                confirmChange(n, rec_dif, num_toler, T_cg, nsign, Tmax_cg);
            
            if changeconfirmed > 0
                % record break time
                rec_cg(num_fc).t_break = clr_data_n(i+1, 1);
                % record change probability
                rec_cg(num_fc).change_prob = 10+changeconfirmed;
                %1X: change prob,
                % X1: confirmed by continous model,
                % X2: confirmed by the 1st zenith model
                % X3: confirmed by the 2st zenith model
                % X4: confirmed by the 3st zenith model
                % record change magnitude
                rec_cg(num_fc).magnitude = median(v_dif_mag, 'omitnan');
                
                % record the RMSEs and normalized differences
                rec_cg(num_fc).adj_rmse = adj_rmse;
                rec_cg(num_fc).rec_rmse = rec_rmse;
                rec_cg(num_fc).rec_normdif = rec_normdif;
                
                % identified and move on for the next functional curve
                num_fc = num_fc+1;
                % start from i+1 for the next functional curve
                i_start(:) = i+1; % all i_start new
                
                % start training again
                BL_train(1:num_mod) = 0;
                i_count(1:num_mod) = 0;
                % initialize the variables
                rec_dif = nan(4, conse, num_mod); % record:
                fit_cft = nan(max_num_c, num_mod);
                rmse = nan(1, num_mod);
                i_dense(:) = i; % update to record i_dense
            elseif changeconfirmed == -1
                % remove noise
                clr_data(i+1, :) = [];
                rec_normdif(i+1, :) = [];
                rec_rmse(i+1, :, :) = [];
                
                % remove index for diff variable
                rec_dif(4, :, :) = rec_dif(4, :, :)-1;% index remove to index -1
                i = i-1; % stay & check again after noise removal
            end
        end
        % move forward to the i+1th clear observation
        i = i+1;
    end % end of while iterative
    
    
    % Two ways for processing the end of the time series
    a_id = num_mod;
    if num_mod == 1
        n_id = a_id;
    else
        n_id = [a_id 1:num_vza];
    end
    
    
    for i_n = n_id % priority all data, small angle, midel angle, and large angle
        if BL_train(i_n) == 1 && i_n == num_mod %  prob based on all data
            % 1) if no break find at the end of the time series
            % define probability of change based on conse
            rec_dif_n_last = rec_dif(:, :, i_n);
            rec_dif_n_last(:, isnan(rec_dif_n_last(4, :))) = [];
            if isempty(rec_dif_n_last)
                continue;
            end
            id_last = 0;
            for i_conse = size(rec_dif_n_last, 2):-1:1
                % sign of change vector
                max_angl = mean(angl(rec_dif_n_last(2, i_conse:size(rec_dif_n_last, 2))));
                if rec_dif_n_last(3, i_conse) <= T_cg || max_angl >= nsign
                    % the last stable id
                    id_last = i_conse;
                    break;
                end
            end
            
            % update change probability
            rec_cg(num_fc).change_prob = (size(rec_dif_n_last, 2)-id_last)/size(rec_dif_n_last, 2);
            % update end time of the curve
            rec_cg(num_fc).t_end = clr_data(end-size(rec_dif_n_last, 2)+id_last, 1);
            
            % record the RMSEs and normalized differences
            rec_cg(num_fc).adj_rmse = adj_rmse;
            rec_cg(num_fc).rec_rmse = rec_rmse;
            rec_cg(num_fc).rec_normdif = rec_normdif;
            
            if size(rec_dif, 2) > id_last % > 1
                % update time of the probable change
                rec_cg(num_fc).t_break = clr_data(end-size(rec_dif_n_last, 2)+id_last+1, 1);
                % update magnitude of change
                rec_cg(num_fc).magnitude = median(rec_dif_n_last(1, id_last+1:size(rec_dif_n_last, 2))); % based on all data
            end
        elseif BL_train(i_n) == 0
            % 2) if break find close to the end of the time series
            % Use [conse,min_num_c*n_times+conse) to fit curve
            % end of fit qa = 20
            qa = 20;
            if num_fc == rec_fc
                % first curve
                i_start = 1;
            else
                i_start = find(clr_data(:, 1) >= rec_cg(num_fc-1).t_break);
                i_start = i_start(1);
            end
            
            if i_n > num_vza
                clr_data_n = clr_data;
                clr_data_n(:, end) = i_n; % all in the 4th stratifed data
            else
                clr_data_n = clr_data; % do not change the orginal zenith angle ID
            end
            [clrx_n, clry_n] = findStratifiedData(clr_data_n, i_n, i_start, size(clr_data, 1));
            clear clr_data_n;
            if length(clrx_n) >= conse
                [fit_cft(:, i_n), rmse(i_n)] = autoRobustFit2(clrx_n, clry_n, min_num_c);
                
                % record time of curve start
                rec_cg(num_fc).t_start = clrx(i_start);
                % record time of curve end
                rec_cg(num_fc).t_end=clrx(end);
                % record break time
                rec_cg(num_fc).t_break = 0;
                % record postion of the pixel
                rec_cg(num_fc).pos = icol_ids;
                % record fitted coefficients
                rec_cg(num_fc).coefs = fit_cft;
                % record rmse of the pixel
                rec_cg(num_fc).rmse = rmse;
                % record change probability
                rec_cg(num_fc).change_prob = 0;
                % record number of observations
                rec_cg(num_fc).num_obs = length(clrx(i_start:end));
                % record fit category
                rec_cg(num_fc).category = qa+min_num_c;
                % record change magnitude
                rec_cg(num_fc).magnitude = 0;
                
                % record the RMSEs and normalized differences
                rec_cg(num_fc).adj_rmse = adj_rmse;
                rec_cg(num_fc).rec_rmse = rec_rmse;
                rec_cg(num_fc).rec_normdif = rec_normdif;
            end
        end
    end
    
%     fprintf('Finishing in %0.2f s\n', toc);
end%end of for icol_ids loop

if save_code == 1
    save(fullfile(dir_BRDF_new, n_rst, [filename_save, '.temp']), 'rec_cg', 'prc_bufferexclude');
    movefile(fullfile(dir_BRDF_new, n_rst, [filename_save, '.temp']), fullfile(dir_BRDF_new, n_rst, filename_save));
    clear clr_vza
end

end % end of function

% function to caculate included angle between ajacent pair of change vectors
function y = angl(v_dif)
v_dif(isnan(v_dif)) = []; % remove nan data
[row, ~] = size(v_dif);
y = zeros(row-1, 1);

if row > 1
    for i = 1:row-1
        a = v_dif(i, :);
        b = v_dif(i+1, :);
        % y measures the opposite of cos(angle)
        y(i) = acos(a*b'/(norm(a)*norm(b)));
    end
else
    y = 0;
end

% convert angle from radiance to degree
y = y*180/pi;
end

% function to find the next i_start for the target VZA interval
function i_start_next = moveStartIndex(clr_data, i_start, n)
i_start_next = find(clr_data(:, 1) > clr_data(i_start(n), 1) & clr_data(:, end) == n, 1);
end

function [clrx_n, clry_n]= findStratifiedData(clr_data, n, i_start, i_end)
% This function is to find out the stratified data from the orginal data
% list, according to indexes of start and end
%
% clr_data: the 3-dimensional clear dataset
% n: the index of the stratified level, such as 1, 2, 3
% i_start: the starting index in all data list
% i_end: the ending index in all data list
%

% update selecting the stratified data x and y, accordingly
clrx_n = clr_data(clr_data(:, 3) == n, 1);
clry_n = clr_data(clr_data(:, 3) == n, 2);

if ~exist('i_start', 'var') && ~exist('i_end', 'var')
    return; % all data at stratified n
else
    ids_n = find(clr_data(i_start, 1) <= clrx_n & clrx_n <= clr_data(i_end, 1)); % IDs of the current period within the stratified data
    clrx_n = clrx_n(ids_n); % pick up the final data
    clry_n = clry_n(ids_n); % pick up the final data
end
end

function isenough = basicRequire(clr_data, istart, i, n, min_num_c, mini_yrs, n_times, num_yrs)
% examine enough observations
if n == length(istart)
    i_span = i-istart(n)+1;
    time_span = (clr_data(i, 1)-clr_data(istart(n), 1))./num_yrs;
else
    clrx_n = findStratifiedData(clr_data, n, istart(n), i);
    if isempty(clrx_n) % no enough data at stratried dataset
        isenough = 0;
        return;
    end
    i_span = length(clrx_n);
    time_span = (clrx_n(end)-clrx_n(1) )./num_yrs;
end
% 1) not enough observations;   % 2) not enough time
isenough = i_span >= n_times*min_num_c && time_span >= mini_yrs;
end

function [isstable, v_dif] = isStableModel(rec_v_dif, fit_cft, adj_rmse, rmse, clrx, T_cg)
% mini rmse
mini_rmse = max(adj_rmse, rmse);
% compare the first clear obs
v_start = rec_v_dif(1)/mini_rmse;
% compare the last clear observation
v_end = rec_v_dif(end)/mini_rmse;
% anormalized slope values
v_slope = fit_cft(2)*(clrx(end)-clrx(1))/mini_rmse;
% differece in model intialization
v_dif = abs(v_slope)+max(abs(v_start), abs(v_end));
v_dif = norm(v_dif); %^2;
isstable = v_dif <= T_cg; % true is stable model, false is not
end

%% function of comparing observations and predications according to conse
function rec_dif_back = diffConseObs(n, clr_data_n, fit_cft, rmse, adj_rmse, ids_conse, num_mod, conse)
% initialize the difference variables, that are same as normal change detection below
rec_dif_back = nan(4, conse, num_mod); % record:
% dimension 1: v_dif_mag % record the magnitude of change
% dimension 2: v_dif     % value of difference for conse obs
% dimension 3: vec_mag   % chagne vector magnitude
% dimension 4: index of observation in orginal datalist
for i = 1:length(ids_conse)
    id_now = ids_conse(i);
    v_dif_mag_tmp = clr_data_n(id_now, 2)-autoTSPred_OLS(clr_data_n(id_now, 1), fit_cft(:, n));
    v_dif_tmp= v_dif_mag_tmp/max(adj_rmse(n), rmse(n)); % z-scores
    vec_mag_tmp = norm(v_dif_tmp); %^2;
    % record zenith angle comparisions
    rec_dif_back(1, i, n) = v_dif_mag_tmp;
    rec_dif_back(2, i, n) = v_dif_tmp;
    rec_dif_back(3, i, n) = vec_mag_tmp;
    rec_dif_back(4, i, n) = id_now;
end
end

%% function of change confirmation
function [changeconfirmed_back, v_dif_mag, v_dif] = confirmChange(n, rec_dif_back, num_toler, T_cg, nsign, Tmax_cg)
% confirm change
changeconfirmed_back = 0;
v_dif_mag = [];
vec_mag_sort = rec_dif_back(3, :, n); % zenith angle data first
vec_mag_sort = sort(vec_mag_sort);
v_dif = rec_dif_back(3, 1, n);
if v_dif > T_cg && vec_mag_sort(min(1 + num_toler, length(vec_mag_sort))) > T_cg ...
        && mean(angl(rec_dif_back(2, :, n))) < nsign % change confirm
    changeconfirmed_back = n; % change identified by zenith angle
    v_dif_mag = rec_dif_back(1, :, n); % record change magnitude
elseif v_dif > Tmax_cg
    changeconfirmed_back = -1; % false change identified by all data
end
end
