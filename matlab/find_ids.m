function [clrx, clry, clrsza, clrsaa] = find_ids(sdate, line_t, line_tsza, ...
    line_tsaa, irow, icol, idgood, nbands)
clrx = sdate(idgood);
clry = line_t(irow,icol,idgood);
clry = double(clry(:));
clrsza = line_tsza(irow,icol,idgood);
clrsza = double(clrsza(:));
clrsaa = line_tsaa(irow,icol,idgood);
clrsaa = double(clrsaa(:));

% find repeated ids
[clrx,uniq_id,~] = unique(clrx);
% mean of repeated values
tmp_y = zeros(length(clrx),nbands);
tmp_sza = zeros(length(clrx),nbands);
tmp_saa = zeros(length(clrx),nbands);
% get the mean values
for i = 1:nbands
    tmp_y(:,i) = clry(uniq_id,i);
    tmp_sza(:,i) = clrsza(uniq_id,i);
    tmp_saa(:,i) = clrsaa(uniq_id,i);
end
clry = tmp_y;
clrsza = tmp_sza;
clrsaa = tmp_saa;

% x_all = clrx;
% y_all = clry;
% sza_all = clrsza;
% saa_all = clrsaa;

end