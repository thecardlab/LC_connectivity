close all
clear all


infile_ca = './Input/lc_round1_compiled_2019.08.06.mat';
infile_vlpn = './Input/ImagedLCs_VLPNcount_2LC_syn(40)_frac(0.1).csv';

out_dir = './Output/';

epsfig = hgexport('factorystyle');
epsfig.Format = 'eps';
set(0, 'DefaultFigurePosition', [0 0 400 200]);
set(0, 'DefaultFigureColor', 'white');
set(0, 'DefaultFigurePaperPositionMode','auto');





load(infile_ca);

% import vlpn count matrix
opts = detectImportOptions(infile_vlpn);
opts.SelectedVariableNames = [2:11];

LC_names_vlpnFile = opts.SelectedVariableNames;

vlpn_cnt = readmatrix(infile_vlpn, opts);


%%
% convert ca data into similarity matrix

% average data for each LC type
cell_types = unique(optic_glomerulus);

ca_avg = [];
for k = 1:length(cell_types)
    % find all entries for current cell type
    lc_match = strcmp(optic_glomerulus, cell_types(k));
    idx = find(lc_match);
    
    ca_avg(k, :) = mean( ca_resp(idx, :), 1);
end


% manually define LC order for calcium data so that it's in the same order
% as connectome vlpn count matrix
LC_order = [4,10,8,7,9,5,3,1,6,2];

LC_names_caFile = cell_types(LC_order);

% create ca resp similarity matrix, order by the new LC_order
sim = [];
for m = 1:length(LC_order)
    for n = 1:length(LC_order)
        r1 = LC_order(m);
        r2 = LC_order(n);
        sim(m,n) = sim_corr( ca_avg(r1, :), ca_avg(r2, :) );
    end
end

% export similarity matrix
tmp_table = array2table(sim, 'VariableNames', LC_names_caFile, 'RowNames', LC_names_caFile);
writetable(tmp_table, [out_dir, 'ImagedLCs_sim_matrix.csv'],'WriteRowNames',true);


%%
% extract unique values for entries that correspond to integrating from
% 2 LC cell types, then sort them

% extract lower diagonal, or the unique entires for VLPN that integrates
% from 2 LC cell types
sim_tril = tril(sim, -1);
ind = find(sim_tril);


sim_list = sim(ind);
vlpn_list = vlpn_cnt(ind);

[~, ind_order] = sort(sim_list);

sorted_sim_list = sim_list(ind_order);
sorted_vlpn_list = vlpn_list(ind_order);



sorted_cumsum_vlpn_list = cumsum(sorted_vlpn_list);
cdf_vlpn = sorted_cumsum_vlpn_list ./ max(sorted_cumsum_vlpn_list);

tmp = ones( size(sorted_vlpn_list) );
tmp = cumsum( tmp );
null_cdf_vlpn = tmp ./ max(tmp);



figure; plot(sorted_sim_list, sorted_vlpn_list);
xlabel('similarity'); ylabel('vlpn count');
hgexport(gcf, [out_dir, 'sim_vs_vlpn_count.eps']  ,epsfig,'Format','eps')
close

figure; plot(sorted_sim_list, cdf_vlpn, sorted_sim_list, null_cdf_vlpn, 'k');
xlabel('similarity'); ylabel('cumulative probability');
hgexport(gcf, [out_dir, 'sim_vs_vlpn_cdf.eps']  ,epsfig,'Format','eps')
close

% export cdf
tmp_table = table( sorted_sim_list, cdf_vlpn, null_cdf_vlpn);
writetable(tmp_table, [out_dir, 'ImagedLCs_sim_vs_vlpnCnt_cdf.csv'],'WriteRowNames',true);

%%

% export data for Kolmogorov-Smirnov test in graphpad prism
% graphpad prism requires raw data counts, not frequency distributions


% replicate entries in sorted_dist_list based on the count in
% sorted_vlpn_list

raw_sim_data = NaN( sum(sorted_vlpn_list) , 1);

ridx = 1;
for k = 1:length(sorted_sim_list)
    val = sorted_sim_list(k);
    num_rep = sorted_vlpn_list(k);
    
    if num_rep > 0
        tmp = repmat(val, num_rep, 1);
        raw_sim_data( ridx:(ridx+num_rep-1) ) = tmp;
        ridx = ridx + num_rep;
    end
end


null_sim_data = sorted_sim_list; % equal probability for all dist, so one entry per similarity

% export 'raw' data
tmp_table = table( raw_sim_data);
writetable(tmp_table, [out_dir, 'ImagedLCs_sim_data_raw.csv'],'WriteRowNames',true);

% export 'null' data
tmp_table = table( sorted_sim_list);
writetable(tmp_table, [out_dir, 'ImagedLCs_sim_data_null.csv'],'WriteRowNames',true);
