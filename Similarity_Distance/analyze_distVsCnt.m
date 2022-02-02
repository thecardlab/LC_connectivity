close all
clear all

% import distances are in nm, should convert to um for plotting
dist_rescale = 1e-3; 

% assumes that input files are sorted by the same LC ordering

infile_dist = './Input/Imaged_LC-LC_RawEuclidDist_nm.csv';
infile_vlpn = './Input/ImagedLCs_VLPNcount_2LC_syn(40)_frac(0.1).csv';

out_dir = './Output/';


epsfig = hgexport('factorystyle');
epsfig.Format = 'eps';
set(0, 'DefaultFigurePosition', [0 0 400 200]);
set(0, 'DefaultFigureColor', 'white');
set(0, 'DefaultFigurePaperPositionMode','auto');


% import distance matrix
opts = detectImportOptions(infile_dist);
opts.SelectedVariableNames = [2:11];

LC_names_distFile = opts.SelectedVariableNames;

LC_dist = readmatrix(infile_dist, opts) .* dist_rescale;

% import vlpn count matrix
opts = detectImportOptions(infile_vlpn);
opts.SelectedVariableNames = [2:11];

LC_names_vlpnFile = opts.SelectedVariableNames;

vlpn_cnt = readmatrix(infile_vlpn, opts);

% extract lower diagonal, or the unique entires for VLPN that integrates
% from 2 LC cell types
dist_tril = tril(LC_dist, -1);
ind = find(dist_tril);

dist_list = LC_dist(ind);
vlpn_list = vlpn_cnt(ind);

[~, ind_order] = sort(dist_list);

sorted_dist_list = dist_list(ind_order);
sorted_vlpn_list = vlpn_list(ind_order);


sorted_cumsum_vlpn_list = cumsum(sorted_vlpn_list);
cdf_vlpn = sorted_cumsum_vlpn_list ./ max(sorted_cumsum_vlpn_list);

tmp = ones( size(sorted_vlpn_list) );
tmp = cumsum( tmp );
null_cdf_vlpn = tmp ./ max(tmp);



figure; plot(sorted_dist_list, sorted_vlpn_list);
xlabel('distance'); ylabel('vlpn count');
hgexport(gcf, [out_dir, 'dist_vs_vlpn_count.eps']  ,epsfig,'Format','eps')
close

figure; plot(sorted_dist_list, cdf_vlpn, sorted_dist_list, null_cdf_vlpn, 'k');
xlabel('distance'); ylabel('cumulative probability');
hgexport(gcf, [out_dir, 'dist_vs_vlpn_cdf.eps']  ,epsfig,'Format','eps')
close


% export cdf
tmp_table = table( sorted_dist_list, cdf_vlpn, null_cdf_vlpn);
writetable(tmp_table, [out_dir, 'ImagedLCs_dist_vs_vlpnCnt_cdf.csv'],'WriteRowNames',true);


%%

% export data for Kolmogorov-Smirnov test in graphpad prism
% graphpad prism requires raw data counts, not frequency distributions


% replicate entries in sorted_dist_list based on the count in
% sorted_vlpn_list

raw_dist_data = NaN( sum(sorted_vlpn_list) , 1);

ridx = 1;
for k = 1:length(sorted_dist_list)
    val = sorted_dist_list(k);
    num_rep = sorted_vlpn_list(k);
    
    if num_rep > 0
        tmp = repmat(val, num_rep, 1);
        raw_dist_data( ridx:(ridx+num_rep-1) ) = tmp;
        ridx = ridx + num_rep;
    end
end


null_dist_data = sorted_dist_list; % equal probability for all dist, so one entry per dist

% export 'raw' data
tmp_table = table( raw_dist_data);
writetable(tmp_table, [out_dir, 'ImagedLCs_dist_data_raw.csv'],'WriteRowNames',true);

% export 'null' data
tmp_table = table( sorted_dist_list);
writetable(tmp_table, [out_dir, 'ImagedLCs_dist_data_null.csv'],'WriteRowNames',true);
