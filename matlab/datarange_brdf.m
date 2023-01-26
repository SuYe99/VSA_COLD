function [idrange_phy, idrange, idrange_all, idclr, sn_pct, idsn] ...
    = datarange_brdf(x_all, y_all, MQF_all, QFDNB_all, date_start, date_end, ...
    mask_all, vza_all, vza_lmt, vza_lmt_all, vaa_all, line_tm, irow_ids, smp_ids)

idrange_phy = y_all > 0 & y_all < 65535;
idrange_DNB = idrange_phy & QFDNB_all(:) == 0;
idrange_date = abs(x_all) >= date_start & abs(x_all) <= date_end;
idrange_cld = mask_all == 1 & MQF_all(:) < 2; 
idrange_vza = vza_all >= vza_lmt(1)*100 & vza_all < vza_lmt(2)*100;
idrange_vza_all = vza_all >= vza_lmt_all(1)*100 & vza_all < vza_lmt_all(2)*100;
idrange = idrange_DNB & idrange_vza & idrange_date;% & idrange_vaa;
idrange_all = idrange_DNB & idrange_vza_all & idrange_date;

% # of clear observatons
idclr = idrange_cld;
% snow pixels
idsn = squeeze(line_tm(irow_ids, smp_ids, :)) == 3;
% percent of snow observations
sn_pct = sum(sum(idsn))/(sum(idclr)+sum(idsn)+0.01);

% x_all = x_all(idrange_DNB & idrange_vza);
% y_all = y_all(idrange_DNB & idrange_vza);
% vza_all = vza_all(idrange_DNB & idrange_vza);
% vaa_all = vaa_all(idrange_DNB & idrange_vza);
% 
% % [x_all, uniq_id, ~] = unique(x_all);
% % mean of repeated values
% % get the mean values
% tmp_y = y_all(uniq_id);
% tmp_vza = vza_all(uniq_id);
% tmp_vaa= vaa_all(uniq_id);
% 
% y_all = tmp_y;
% vza_all = tmp_vza;
% vaa_all = tmp_vaa;
end