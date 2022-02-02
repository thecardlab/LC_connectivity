clear all; close all;


% clustered by synapses/VLPN 
LC_types = {'LC17', 'LPLC2', 'LC4', 'LC9', 'LPLC4', 'LC16', 'LC6', 'LC26', 'LC24', 'LC25', 'LC20', 'LC22', 'LPLC1', 'LC18', 'LC15', 'LC11', 'LC21', 'LC13', 'LC12'}; % new order based coInt_spi matrix, clustergram on half and half separately


file_dir = './neuPrint/';
filename_add = '_excludedOL';
out_dir = './processed/';


% user defined thresholds
threshold_min_num_syn = 40; % each VLP neuron must make at least # synapses across all LC cell types
threshold_min_LC_frac = 0.1; % proportion of VLPN<>LC synapses each LC cell type must make to be considered as connected


% define color map
% https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
custom_cmap = cmocean('ice',100);


epsfig = hgexport('factorystyle');
epsfig.Format = 'eps';
set(0, 'DefaultFigurePosition', [0 0 800 800]);
set(0, 'DefaultFigureColor', 'white');
set(0, 'DefaultFigurePaperPositionMode','auto');

%%
% clean up the table to remove LCs as possible "output_type"

% collate LC/LPLC body IDs
all_LC_bodyId = [];
num_LC_neurons = [];
for i = 1:length(LC_types)
    data = readtable([file_dir, LC_types{i},filename_add,'.csv']);

    % collate all LC/LPLC body IDs
    num_LC_neurons(i) = length( unique(data.LC_bodyId) );
    all_LC_bodyId = cat(1, all_LC_bodyId, unique(data.LC_bodyId));
end
writetable(table(LC_types', num_LC_neurons'),[out_dir, 'num_LC_neurons.csv'],'WriteRowNames',true)


% remove LC to LC connections by excluding any output_bodyID == LC body ID    
for i = 1:length(LC_types)
    
    data = readtable([file_dir, LC_types{i},filename_add,'.csv']);
    
    N_out_bodyId = unique(data.output_bodyId);
    
    % remove LC to LC connections, appear to be quite common in the PVLP also    
    ind_del = [];
    for k = 1:length(all_LC_bodyId)
        ind = find( N_out_bodyId == all_LC_bodyId(k) );
        ind_del = cat(1, ind_del, ind);
    end
    
    N_out_bodyId( ind_del ) = []; % delete entries with the any LC as output target
    
    
    for j = 1:length(N_out_bodyId)
        N_ind = find(data.output_bodyId==N_out_bodyId(j));
        N_out_name{j,1} = data.output_type(N_ind(1));
        num_LC_synapsing(j,1) = length(N_ind);
        num_synapses(j,1) = sum(data.c_weight(N_ind));
        avg_synapses(j,1) = num_synapses(j,1) / num_LC_synapsing(j,1);
    end
    
    temp_table = table(N_out_bodyId, N_out_name, num_LC_synapsing, num_synapses, avg_synapses);
    LC_tables{i} = sortrows(temp_table,'num_LC_synapsing','descend');
    
    clear data N* num* avg* temp*
    
    % Save the tables
    writetable(LC_tables{i},[out_dir, 'Collated_', LC_types{i},'.csv'],'WriteRowNames',true)
end
    


%% Make single matrix of LC types (row) and unique VLP neurons (col)

% Put together all 'unique' bodyId's for each LC and make one collated
% unique list
bodyId_list = [];
bodyId_names = [];
for i = 1:length(LC_tables)
    bodyId_list = [bodyId_list; LC_tables{i}.N_out_bodyId];
    bodyId_names = [bodyId_names; LC_tables{i}.N_out_name];
end

[bodyId_list_unique, IA, IC] = unique(bodyId_list);
bodyId_names_unique = bodyId_names(IA);


for i = 1:length(LC_tables)
    for j = 1:length(bodyId_list_unique)
        ind = find(LC_tables{i}.N_out_bodyId == bodyId_list_unique(j));
        if isempty(ind)
            mat_synapses(i,j) = 0;
            mat_numLC(i,j) = 0;
        else
            mat_synapses(i,j) = LC_tables{i}.num_synapses(ind);
            mat_numLC(i,j) = LC_tables{i}.num_LC_synapsing(ind);
        end
    end
end

bodyId_list_unique_STR = cellstr(num2str(bodyId_list_unique,'%010d'));
bodyId_list_unique_STR_valid = cellfun(@(c)['bodyId_', c ],bodyId_list_unique_STR,'UniformOutput',false);


%% filter data based on # of synapses

% enforce each VLP neuron must integrate from at least # of LC synapses
min_num_synapses = threshold_min_num_syn;

tot_synapses = sum(mat_synapses, 1);

tmp_ind = find( tot_synapses > min_num_synapses);
cnt_numSyn = mat_synapses(:, tmp_ind); % remove VLP neurons that dont integrate min # of synapses, each row is an LC cell type, col is a post-synaptic neuron

cnt_numLC = mat_numLC(:, tmp_ind); % each row is an LC cell type, col is a post-synaptic neuron
cnt_totSyn = sum(cnt_numSyn, 1);
cnt_fracSyn = cnt_numSyn ./ repmat(cnt_totSyn, [length(LC_types) 1]);

filtered_bodyId = bodyId_list_unique_STR_valid(tmp_ind);

% save output after min num syn thresholding
tmp_cnt_numSyn_table = array2table( cnt_numSyn );
tmp_cnt_numSyn_table.Properties.VariableNames = filtered_bodyId;
tmp_cnt_numSyn_table.Properties.RowNames = LC_types;

tmp_cnt_numLC_table = array2table( cnt_numLC );
tmp_cnt_numLC_table.Properties.VariableNames = filtered_bodyId;
tmp_cnt_numLC_table.Properties.RowNames = LC_types;


writetable(tmp_cnt_numSyn_table,[out_dir, 'LC2VLPN_synCount_minSynThreshold(', num2str(threshold_min_num_syn),').csv'],'WriteRowNames',true)
writetable(tmp_cnt_numLC_table,[out_dir, 'LC2VLPN_numInputLCneurons_minSynThreshold(', num2str(threshold_min_num_syn),').csv'],'WriteRowNames',true)


%%
% threshold data to define whether VLP neuron is integrating from an LC 
% cell type, note that each LC cell type can contain between 30-200 neurons
% this threshold effectively removes weakly connected LC cell types

% each LC cell type must contribute a minimum fraction of synapses
% to a VLP neuron (threshold defined by 'min_LC_frac')

% e.g. if a VLP neuron has 40 total of synapses across all LC types
% then each LC type must make min_LC_frac * 40 at minimum to be considered
% as connected
% if min_LC_frac = 0.1, then each LC type must make at least 4 synapses


min_LC_frac = threshold_min_LC_frac;



% pathway analysis


% cnt_numSyn: number of synapses between (many) pre-synaptic LC neurons and one post-synaptic neuron; each row is an LC cell type, each col is a post-synaptic neuron
% cnt_numLC: number of pre-syp neurons
% cnt_totSyn: total number of synapses from pre-syn to a given post-syn neuron
% cnt_fracSyn: cnt_numSyn / cnt_totSyn, for each pre-syn cell type, what proportion of connection does it account for
cnt_synDen = cnt_numSyn ./ cnt_numLC;


% filter data by requiring individual cnt_fracSyn >= X% of total connections
fcnt_numSyn = cnt_numSyn;
fcnt_numSyn( cnt_fracSyn < min_LC_frac ) = 0; % remove connections below threshold

fcnt_totSyn = sum( fcnt_numSyn, 1 );
fcnt_synDen = fcnt_numSyn ./ cnt_numLC; fcnt_synDen( isnan(fcnt_synDen) ) = 0;

fcnt_numPreTypes = sum( fcnt_numSyn > 0, 1 ); % how many pre-syn LC types are connected to a given post-syn neuron (each column)

fsuffix = ['_syn(', num2str(min_num_synapses), ')_frac(', num2str(min_LC_frac) ')'];

% LC to VLPN histogram for all connectivity data
tmp_cnt = hist(fcnt_numPreTypes, [1:length(LC_types)]);
temp_table = table([1:length(LC_types)]', tmp_cnt');
writetable(temp_table,[out_dir, 'conn_hist_all', fsuffix, '.csv'],'WriteRowNames',true);

figure; 
hist(fcnt_numPreTypes, [1:length(LC_types)]);
xlabel('Integrates from # of LC types'); ylabel('VLP neuron count');
hgexport(gcf, [out_dir,'conn_hist_all', fsuffix, '.eps'] ,epsfig,'Format','eps')
close

%%
% can generate full statistics for export
% post-synaptic neuron bodyIDs are stored in filtered_bodyId
connTypes = unique( fcnt_numPreTypes );
for k = 1:length(connTypes)
    
    idx = find( fcnt_numPreTypes == connTypes(k) );
    
    ex_numSyn = fcnt_numSyn(:, idx );
    ex_synDen = fcnt_synDen(:, idx );
    ex_numLC = cnt_numLC(:, idx );
    
    bodyId_post = filtered_bodyId( idx );
    
    % export to csv
    tbl_numSyn = array2table( transpose(ex_numSyn), 'VariableNames', LC_types, 'RowNames', bodyId_post);
    writetable(tbl_numSyn,[out_dir, 'VLP_(' ,num2str(connTypes(k)),')LC_numSyn', fsuffix, '.csv'],'WriteRowNames',true);

    tbl_synDen = array2table( transpose(ex_synDen), 'VariableNames', LC_types, 'RowNames', bodyId_post);
    writetable(tbl_synDen,[out_dir, 'VLP_(' ,num2str(connTypes(k)),')LC_synDen', fsuffix, '.csv'],'WriteRowNames',true);
   
    tbl_numLC = array2table( transpose(ex_numLC), 'VariableNames', LC_types, 'RowNames', bodyId_post);
    writetable(tbl_numLC,[out_dir, 'VLP_(' ,num2str(connTypes(k)),')LC_numInputCells', fsuffix, '.csv'],'WriteRowNames',true);
end




%%
% generate special connectivity matrix just for the two LC types to 1 VLPN on off-diagonals
% and 1 LC type to 1VLPN along diagonal
% use simple binary count for each connection, split ratio for 2 LC inputs
cumConn = zeros( length(LC_types) );

% off diagonals now
idx = find( fcnt_numPreTypes == 2 );
for k = 1:length(idx)
    conn = fcnt_synDen( :, idx(k) );
    LC_idx = find(conn);
    idx_a = LC_idx(1); idx_b = LC_idx(2); 

    % idx_a is always the stronger connection, idx_b is the weaker
    cumConn( idx_a, idx_b ) = cumConn( idx_a, idx_b ) + 1;
    cumConn( idx_b, idx_a ) = cumConn( idx_b, idx_a ) + 1;
end
tmp = min(cumConn); min_val = mean(tmp);
tmp = max(cumConn); max_val = mean(tmp); % smooth out outliers

% do the diagonals 
idx = find( fcnt_numPreTypes == 1 );
for k = 1:length(idx)
    conn = fcnt_synDen( :, idx(k) );
    LC_idx = find( conn );
    cumConn( LC_idx, LC_idx ) = cumConn( LC_idx, LC_idx ) + 1; % each conn has weight of 1
end


% clustergram(cumConn,'RowLabels',LC_types)
% clustergram(cumConn(1:10,1:10),'RowLabels',LC_types(1:10))
% clustergram(cumConn(11:20,11:20),'RowLabels',LC_types(11:20))


figure;
imagesc( cumConn )
colormap(custom_cmap)
xticks([1:length(LC_types)]); yticks([1:length(LC_types)]);
set(gca, 'XTickLabel', LC_types, 'YTickLabel', LC_types);
axis equal; colorbar;

caxis([min_val, max_val]);
hgexport(gcf, [out_dir, 'conn_2LC_connCount.eps']  ,epsfig,'Format','eps')
close

writetable( array2table(cumConn, 'VariableNames', LC_types,'RowNames', LC_types),[out_dir, 'conn_2LC_', fsuffix, '.csv'],'WriteRowNames',true);


%%
% functional imaging data is only for the first 10 LC cell types
% count how many single LC to VLPN and two LC to VLPN are accounted for
% cumConn_imaged = cumConn(1:10, 1:10);
% 
% cnt_ImgLC_single = sum( diag(cumConn_imaged) );
% cnt_ImgLC_double = sum( sum( tril(cumConn_imaged, -1) ) );
% cnt_ommittedLC_double = sum( sum( cumConn(11:end, 1:10) ) );

%%

% average synaptic input density, split into primary vs secondary inputs
% dominant or primary input is defined by total # of synapses, not synaptic
% density, and this is a decent metric because strong connections could be
% due to either high synaptic density or medium synaptic density but a lot
% of connections

% need to track how many synapses from LC type, and also how many LC
% neurons are inputs

% also which LC input is determined by the row index, column just specifies
% the unique pairing of LC inputs

cumNumLCsyn = zeros( length(LC_types) );
cumNumLC = zeros( length(LC_types) );

% off diagonals now
idx = find( fcnt_numPreTypes == 2 );
for k = 1:length(idx)
    numSyn = fcnt_numSyn( :, idx(k) );
    numLC = cnt_numLC(:, idx(k) );
    
    r_idx = find(numSyn); % LC cell type index
    c_idx = flip(r_idx); % LC-LC pairing index
    
    for m = 1:length(r_idx)
        cumNumLCsyn(r_idx(m), c_idx(m)) = numSyn(r_idx(m)) + cumNumLCsyn(r_idx(m), c_idx(m));
        cumNumLC(r_idx(m), c_idx(m)) = numLC(r_idx(m)) + cumNumLC(r_idx(m), c_idx(m));    
    end
end
cumAvgSynDen = cumNumLCsyn ./ cumNumLC;
cumAvgSynDen( isnan(cumAvgSynDen) ) = 0;

tmp = min(cumAvgSynDen); min_val = mean(tmp);
tmp = max(cumAvgSynDen); max_val = mean(tmp); % smooth out outliers

tmp = min(cumNumLCsyn); min_val_cnt = mean(tmp);
tmp = max(cumNumLCsyn); max_val_cnt = mean(tmp); % smooth out outliers


% % do the diagonals 
idx = find( fcnt_numPreTypes == 1 );
for k = 1:length(idx)
    numSyn = fcnt_numSyn( :, idx(k) );
    numLC = cnt_numLC(:, idx(k) );
    
    r_idx = find( numSyn ); % LC cell type index
    c_idx = flip(r_idx); % LC-LC pairing index

    for m = 1:length(r_idx)
        cumNumLCsyn(r_idx(m), c_idx(m)) = numSyn(r_idx(m)) + cumNumLCsyn(r_idx(m), c_idx(m));
        cumNumLC(r_idx(m), c_idx(m)) = numLC(r_idx(m)) + cumNumLC(r_idx(m), c_idx(m));    
    end
end

% calculate syn density
cumAvgSynDen = cumNumLCsyn ./ cumNumLC;
cumAvgSynDen( isnan(cumAvgSynDen) ) = 0;

figure;
imagesc( cumAvgSynDen )
colormap(custom_cmap)
xticks([1:length(LC_types)]); yticks([1:length(LC_types)]);
set(gca, 'XTickLabel', LC_types, 'YTickLabel', LC_types);

caxis([min_val, max_val]);
hgexport(gcf, [out_dir, 'conn_2LC_avgSynDen.eps']  ,epsfig,'Format','eps')
close

figure;
imagesc( cumNumLCsyn )
colormap(custom_cmap)
xticks([1:length(LC_types)]); yticks([1:length(LC_types)]);
set(gca, 'XTickLabel', LC_types, 'YTickLabel', LC_types);

caxis([min_val_cnt, max_val_cnt]);
hgexport(gcf, [out_dir, 'conn_2LC_totSynCnt.eps']  ,epsfig,'Format','eps')
close

% total # of synapses / # of VLPN
cumSynPerVLPN = cumNumLCsyn ./ cumConn;
cumSynPerVLPN( isnan(cumSynPerVLPN) ) = 0;

figure;
imagesc( cumSynPerVLPN );
colormap(custom_cmap)
xticks([1:length(LC_types)]); yticks([1:length(LC_types)]);
set(gca, 'XTickLabel', LC_types, 'YTickLabel', LC_types);

tmp = min(cumSynPerVLPN); min_val = mean(tmp);
tmp = max(cumSynPerVLPN); max_val = mean(tmp); % smooth out outliers
caxis([min_val, max_val]);
axis equal; colorbar;
hgexport(gcf, [out_dir, 'conn_2LC_synPerVLPN.eps']  ,epsfig,'Format','eps')

%clustergram(cumSynPerVLPN,'RowLabels',LC_types)

