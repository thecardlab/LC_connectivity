clear all; close all;



% clustered by synapses/VLPN 
LC_types = {'LC17', 'LPLC2', 'LC4', 'LC9', 'LPLC4', 'LC16', 'LC6', 'LC26', 'LC24', 'LC25', 'LC20', 'LC22', 'LPLC1', 'LC18', 'LC15', 'LC11', 'LC21', 'LC13', 'LC12'}; 


LC2plot = [1,2,3, 10, 13, 14, 15, 16, 17, 19]; % select subset that were imaged


file_dir = './neuPrint/';
filename_add = '_tBars_excludedOL';



epsfig = hgexport('factorystyle');
epsfig.Format = 'eps';
set(0, 'DefaultFigurePosition', [0 0 400 400]);
set(0, 'DefaultFigureColor', 'white');
set(0, 'DefaultFigurePaperPositionMode','auto');

out_dir = './processed/';

% load color map
cmap_file = 'LC_colors.csv';
color_table = readtable(cmap_file);


%%


% collate LC/LPLC body IDs
all_LC_bodyId = [];
num_LC_neurons = [];
for i = 1:length(LC_types)
    tmp = readtable([file_dir,LC_types{i},filename_add,'.csv']);

    num_LC_neurons(i) = length( unique(tmp.LC_bodyId) );
    all_LC_bodyId = cat(1, all_LC_bodyId, unique(tmp.LC_bodyId));
end
writetable(table(LC_types', num_LC_neurons'),[out_dir, 'num_LC_neurons.csv'],'WriteRowNames',true)


% remove LC to LC connections by excluding any output_bodyID == LC body ID    
for i = 1:length(LC_types)
    data = readtable([file_dir,LC_types{i},filename_add,'.csv']);
    
    N_out_bodyId = unique(data.output_bodyId); % neurons downstream of LC
    
    ind_del = [];
    for k = 1:length(all_LC_bodyId)
        ind = find( N_out_bodyId == all_LC_bodyId(k) ); % if no match (ie not LC), will return null
        ind_del = cat(1, ind_del, ind);
    end
    
    N_out_bodyId( ind_del ) = []; % delete entries with the any LC as output target
    
    
    % reconstruct data table, by using non LC output entries
    data_cleaned = table;
    
    for j = 1:length(N_out_bodyId)
        N_ind = find(data.output_bodyId==N_out_bodyId(j));
        
        % copy over entries to the new data_cleaned table
        tmp = data(N_ind, :);
        data_cleaned = [data_cleaned; tmp];
    end
    
    
  %  temp_table = table(N_out_bodyId, N_out_name, num_LC_synapsing, num_synapses, avg_synapses);
 %   LC_tables{i} = sortrows(temp_table,'num_LC_synapsing','descend');
    LC_tables{i} = data_cleaned;
    
    clear data N* num* avg* temp*
    
    % Save out the tables
    writetable(LC_tables{i},[out_dir, 'Collated_', LC_types{i},'.csv'],'WriteRowNames',true)
end


%%
% to visualize LC output T-bars, need to generate a list of unique T-bars
% one LC T-bar can be presynaptic to many post-synaptic neurons downstream of LC

% T-bars can be grouped either by individual LC cell or the post-synaptic
% target

% visualize by LC type, color each LC
tBars_arrays = {};
figure;

dot_scale = 5;
for i = 1:length(LC2plot)
    k = LC2plot(i);
    
    tBars = [LC_tables{k}.X, LC_tables{k}.Y, LC_tables{k}.Z];
    uniq_tBars = unique(tBars, 'rows');
    tBars_arrays{k} = uniq_tBars;
    
    tmp = strcmpi( LC_types{k}, color_table.LC_type );
    idx = find(tmp);
    
    rgb_mat = repmat([color_table.R(idx), color_table.G(idx), ...
        color_table.B(idx)]./255, size(uniq_tBars,1), 1);
    scatter3( uniq_tBars(:,1), uniq_tBars(:,2), uniq_tBars(:,3), dot_scale, rgb_mat)
    hold on;
end
view(0,180);
axis on
hgexport(gcf, [out_dir, 'LC_glom_withAxis.eps'] ,epsfig,'Format','eps')

axis off
hgexport(gcf, [out_dir, 'LC_glom.eps'] ,epsfig,'Format','eps')


%%

% calculate center positions, just use kmeans as shortcut, by assigning k=1

LC_centroids = [];
dot_scale = 200;
figure;
for m = 1:length(LC2plot)
    i = LC2plot(m);
    
    [~, centroid] = kmeans( tBars_arrays{i}, 1 ); % set k = 1, looking at a single LC cell type at a time
    LC_centroids = cat(1, LC_centroids, centroid);
    
    % plot the centroids
    tmp = strcmpi( LC_types{i}, color_table.LC_type );
    idx = find(tmp);
    
    rgb_val = [color_table.R(idx), color_table.G(idx), ...
        color_table.B(idx)]./255;
    scatter3( centroid(1), centroid(2), centroid(3), dot_scale, rgb_val, 'filled')
    hold on;
end
hold off;
view(0,180);
axis off
hgexport(gcf, [out_dir, 'LC_glom_centroids.eps'] ,epsfig,'Format','eps')

LC_dist = pdist2( LC_centroids, LC_centroids, 'euclidean' );

% according to Kazunori, 1 pixel = 8 nm in hemibrain dataset
LC_dist = LC_dist .* 8; % nm

% export the raw LC-LC distance table, should be in nanometers
tbl_LCrawDist = array2table(LC_dist, 'VariableNames', LC_types(LC2plot), 'RowNames', LC_types(LC2plot));
writetable(tbl_LCrawDist,[out_dir, 'LC-LC_RawEuclidDist_nm.csv'],'WriteRowNames',true);

