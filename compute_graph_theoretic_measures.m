function [graph_measures] = compute_graph_theoretic_measures(erptfd_assoc_matrix_name,outfilename,comps_defs_input,preproc_options,num_rand_null_networks),

% function [graph_measures] = compute_graph_theoretic_measures(erptfd_assoc_matrix,comps_defs_input,preproc_options,num_rand_null_networks),
%
%  this function calculates several network-based statistics via Sporn and Rubinov's Brain Connectivity Toolbox (https://sites.google.com/site/bctnet/)
%
%  network measures based on graph theory are calculated from connectivity matrices obtained from extract_averages_subsampling_all2all_connectivity
%  computed measures are stored in the structure graph_measures
%
%  regions-of-interest can be extracted from the time-frequency representation, which are then used to create association matrices
%  association matrices can be thresholded, and/or binarized, and then used to compute network statistics using BCT functions
%
%  OPTIONS:
%
%  preproc_options.threshold                 specifies whether or not to threshold matrix (0 = no, 1 = yes)
%                                                            1 = threhsolding performed
%                                                            0 = no threhsolding performed
%  preproc_options.thresh_type               specifies type of thresholding
%                                                            subsample_without_replacement: random sample unique sweeps (i.e., no duplicates)
%                                                            user_defined: uses sweeps from previous static_sets (define in static_sets)


%  function structure
%    1. creat symmetric grand average association matrix by collapsing across condition averages (catcodes)
%      - for each requested component
%      A. Decide to threshold... Yes [preproc_options.threshold == 1] No [preproc_options.threshold = 0]
%          A1. proportion of strongest weights (0 <= p <= 1) [preproc_options.thresh_type = 'proportion']
%              can be a scalar (applied to all components), or a vector of length(comps_defs_input) to
%              specify different proportions for each component. Values must be specified in [preproc_options.thresh_val]
%          A2. keeping edges with a value larger than or equal to one standard deviation above the median connectivity value
%              this value is derived from the grand average distribution (i.e., condition independent), and is used to threshold
%              the (MX Cohen, 2011, 2012, 2014)
%
%  EXAMPLE ROI:
%
%  the code will compute a grand average (collapsing catcodes) association matrix, which is used for:
%    1. generating null-model via a randomizded network (preserves weight/degree distributions, and strength if fully connected and weighted)
%       these null-models are used to compute null statistics to compute normalized network measures (e.g., small-worldness)
%    2. calculating condition-independent thresholds via 1 SD > median, if requested (MX Cohen, 2011, 2012, 2014)
%
%  list of computed graph theory measures ([X] denotes catcode, [ROI] denotes component name)
%
%    measures of integration
%      graph_measures.catcode_[X].[ROI]_char_path_length_[X]:  characteristic path length (scalar)
%      graph_measures.catcode_[X].[ROI]_Eglob_[X]:             global efficiency (scalar)
%      graph_measures.catcode_[X].[ROI]_nEglob_[X]:            normalized global efficiency (Eglob/mean(null model Eglob)) 

%    measures of segregation
%      graph_measures.catcode_[X].[ROI]_weighted_node_cc_[X]:  weighted nodal clustering coefficient (elec vector)
%      graph_measures.catcode_[X].[ROI]_weighted_cc_[X]:       mean weighted clustering coefficient across nodes (scalar)
%      graph_measures.catcode_[X].[ROI]_node_Eloc_[X]:              local efficiency (vector; analogous to clustering coeff)
%
%    measures of centrality
%      graph_measures.catcode_[X].[ROI]_node_betweenness_centrality_[X]:  node betweenness centrality (elec vector)
%
%    basic measures
%      graph_measures.catcode_[X].[ROI]_node_strength_[X]:   node strength (elec vector)
%%
%    small-world measures:
%      graph_measures.catcode_[X].[ROI]_gamma_[X]:   normalized clustering coefficient (weighted_cc/mean(null_weighted_cc))
%      graph_measures.catcode_[X].[ROI]_lambda_[X]:   normalized charactestic path length (char_path_length/mean(char_path_length))
%      graph_measures.catcode_[X].[ROI]_sigma_[X]:    small-world = gamma/lambda
%      graph_measures.catcode_[X].[ROI]_sigma_z_[X]:  z-value of real sigma compared to permuted sigma value (MXC 2014)
%
% Jeremy Harper, UMN, 07/11/15


% load erptfd_assoc_matrix
load(erptfd_assoc_matrix_name);

% parse components
components = [];
for j = 1:length(comps_defs_input(:,1)),
  [comps_defs(j).name, ...
   comps_defs(j).timesbin, ...
   comps_defs(j).timeebin, ...
   comps_defs(j).freqsbin, ...
   comps_defs(j).freqebin] = strread(char(comps_defs_input(j,:)),'%s%f%f%f%f');
end

% initalize measures structure
for aa = 1:length(comps_defs),
  name    = char(comps_defs(aa).name); 
  for bb = 1:length(erptfd_assoc_matrix(1).catcodes),
    if isnumeric(erptfd_assoc_matrix(1).catcodes(bb)),
      cur_catcode = num2str(erptfd_assoc_matrix(1).catcodes(bb));
    elseif ~isnumeric(erptfd_assoc_matrix(1).catcodes(bb)),
      cur_catcode = char(erptfd_assoc_matrix(1).catcodes(bb));
    end

    graph_measures.components.([name '_threshold_' cur_catcode]) = nan(length(erptfd_assoc_matrix),1);
    graph_measures.components.([name '_proportion_edges_retained_' cur_catcode]) = nan(length(erptfd_assoc_matrix),1);

    if isequal(preproc_options.binarize,1),

      graph_measures.components.([name '_bin_char_path_length_' cur_catcode]) = nan(length(erptfd_assoc_matrix),1);
      graph_measures.components.([name '_bin_node_betweenness_centrality_' cur_catcode]) = nan(length(erptfd_assoc_matrix(1).elec),length(erptfd_assoc_matrix));
      graph_measures.components.([name '_bin_node_cc_' cur_catcode]) = nan(length(erptfd_assoc_matrix(1).elec),length(erptfd_assoc_matrix));
      graph_measures.components.([name '_bin_cc_' cur_catcode]) = nan(length(erptfd_assoc_matrix),1);
      graph_measures.components.([name '_bin_node_degree_' cur_catcode]) = nan(length(erptfd_assoc_matrix(1).elec),length(erptfd_assoc_matrix));
      graph_measures.components.([name '_bin_Eglob_' cur_catcode]) = nan(length(erptfd_assoc_matrix),1);
      graph_measures.components.([name '_bin_node_Eloc_' cur_catcode]) = nan(length(erptfd_assoc_matrix(1).elec),length(erptfd_assoc_matrix));
      graph_measures.components.([name '_bin_gamma_' cur_catcode])  = nan(length(erptfd_assoc_matrix),1);
      graph_measures.components.([name '_bin_lambda_' cur_catcode]) = nan(length(erptfd_assoc_matrix),1);
      graph_measures.components.([name '_bin_nEglob_' cur_catcode]) = nan(length(erptfd_assoc_matrix),1);
      graph_measures.components.([name '_bin_sigma_' cur_catcode]) = nan(length(erptfd_assoc_matrix),1);
      graph_measures.components.([name '_bin_sigma_z_' cur_catcode]) = nan(length(erptfd_assoc_matrix),1);

    elseif isequal(preproc_options.binarize,0),

      graph_measures.components.([name '_wei_char_path_length_' cur_catcode]) = nan(length(erptfd_assoc_matrix),1);
      graph_measures.components.([name '_wei_node_betweenness_centrality_' cur_catcode]) = nan(length(erptfd_assoc_matrix(1).elec),length(erptfd_assoc_matrix));
      graph_measures.components.([name '_wei_node_cc_' cur_catcode]) = nan(length(erptfd_assoc_matrix(1).elec),length(erptfd_assoc_matrix));
      graph_measures.components.([name '_wei_cc_' cur_catcode]) = nan(length(erptfd_assoc_matrix),1);
      graph_measures.components.([name '_wei_node_strength_' cur_catcode]) = nan(length(erptfd_assoc_matrix(1).elec),length(erptfd_assoc_matrix));
      graph_measures.components.([name '_wei_Eglob_' cur_catcode]) = nan(length(erptfd_assoc_matrix),1);
      graph_measures.components.([name '_wei_node_Eloc_' cur_catcode]) = nan(length(erptfd_assoc_matrix(1).elec),length(erptfd_assoc_matrix));
      graph_measures.components.([name '_wei_gamma_' cur_catcode])  = nan(length(erptfd_assoc_matrix),1);
      graph_measures.components.([name '_wei_lambda_' cur_catcode]) = nan(length(erptfd_assoc_matrix),1);
      graph_measures.components.([name '_wei_nEglob_' cur_catcode]) = nan(length(erptfd_assoc_matrix),1);
      graph_measures.components.([name '_wei_sigma_' cur_catcode]) = nan(length(erptfd_assoc_matrix),1);
      graph_measures.components.([name '_wei_sigma_z_' cur_catcode]) = nan(length(erptfd_assoc_matrix),1);

    end
  end
end


for main_loop_counter = 1:length(erptfd_assoc_matrix),

  disp(['  Current subname: ' char(erptfd_assoc_matrix(main_loop_counter).subname)]); disp(blanks(1)');

  cur_erptfd_assoc_matrix = erptfd_assoc_matrix(main_loop_counter); % get cur sub data

  % fill fields 
  graph_measures.subnames(main_loop_counter,:)  = cur_erptfd_assoc_matrix.subname;
  graph_measures.subnum(main_loop_counter,:)    = cur_erptfd_assoc_matrix.subnum;
  graph_measures.catcodes(main_loop_counter,:)  = cur_erptfd_assoc_matrix.catcodes;
  graph_measures.elec(:,main_loop_counter)      = [cur_erptfd_assoc_matrix.elec]';

  % make symmetric grand average association matrix for use with small-world code below
  for ii = 1:length(cur_erptfd_assoc_matrix.catcodes),

    % get current sbj/catcode-specific network matrix
    if isnumeric(cur_erptfd_assoc_matrix.catcodes(ii)),
      cur_assoc_matrix = cur_erptfd_assoc_matrix.(['assoc_matrix_' num2str(cur_erptfd_assoc_matrix.catcodes(ii))]);
    elseif ~isnumeric(cur_erptfd_assoc_matrix.catcodes(ii)),
      cur_assoc_matrix = cur_erptfd_assoc_matrix.(['assoc_matrix_' char(cur_erptfd_assoc_matrix.catcodes(ii))]);
    end

    % make symmetric
    for fi=1:size(cur_assoc_matrix,3)
      for ti=1:size(cur_assoc_matrix,4)
        temp = squeeze(cur_assoc_matrix(:,:,fi,ti));
        temp = temp + triu(temp)'; 
        cur_assoc_matrix(:,:,fi,ti) = temp;
      end
    end
    tempmtrx(:,:,:,:,ii) = cur_assoc_matrix;
  end
  clear tempmtr cur_assoc_matrix temp ti fi ii;
  grand_avg_assoc_matrix = mean(tempmtrx,5); % calculate grand average TF association matrix

  disp(['    Computing grand average association matrix and null networks']); disp(blanks(1)');

  for comp_counter = 1:length(comps_defs),

    % current component defs
    name    = char(comps_defs(comp_counter).name); 
    c_timesbin = comps_defs(comp_counter).timesbin;
    c_timeebin = comps_defs(comp_counter).timeebin;
    c_freqsbin = comps_defs(comp_counter).freqsbin;
    c_freqebin = comps_defs(comp_counter).freqebin;

    % compute mean assoc_matrix for cTF ROI
    mean_grand_cur_assoc_matrix = mean(grand_avg_assoc_matrix(:,:,:,cur_erptfd_assoc_matrix.tbin+c_timesbin:cur_erptfd_assoc_matrix.tbin+c_timeebin),4);
    mean_grand_cur_assoc_matrix = mean(mean_grand_cur_assoc_matrix(:,:,c_freqsbin:c_freqebin),3);

    % threshold?
    if isequal(preproc_options.threshold,1) && isequal(preproc_options.thresh_type, '1SD>median'), % 1 SD >= median, from MXC (2011, 2012, 2014)
      tempsynch = nonzeros(mean_grand_cur_assoc_matrix); % vector of all nonzero PS vals
      thresh(comp_counter) = median(tempsynch) + std(tempsynch); % calc thresh 1 SD > median
     %mean_grand_cur_assoc_matrix(mean_grand_cur_assoc_matrix<thresh(comp_counter)) = 0; % set to 0 all values under threshold
      mean_grand_cur_assoc_matrix = threshold_absolute(mean_grand_cur_assoc_matrix, thresh(comp_counter)); % set to 0 all values under threshold
      disp(['    1 SD > Median thresholding performed']); disp(blanks(1)');
    elseif isequal(preproc_options.threshold, 1) && isequal(preproc_options.thresh_type, 'proportion'),
      thresh(comp_counter) = preproc_options.thresh_val;
      mean_grand_cur_assoc_matrix = threshold_proportional(mean_grand_cur_assoc_matrix, thresh(comp_counter)); % proportion of strongest weights
      disp(['    Proportion thresholding performed']); disp(blanks(1)');
    elseif isequal(preproc_options.threshold,0),
      disp(['    No thresholding performed']); disp(blanks(1)');
      thresh(comp_counter) = -999;
    end

    % binarize?
    if isequal(preproc_options.binarize,1),
      mean_grand_cur_assoc_matrix = weight_conversion(mean_grand_cur_assoc_matrix, 'binarize');
      disp(['    Binarization performed']); disp(blanks(1)');
    end

    disp(['    Begin computing randomized networks and null graph measures for component ' name]); disp(blanks(1)');
    % generate null models and compute null CC/path/Eglob for normalizing CC/char path length/Eglob
    for jj = 1:num_rand_null_networks,


      if isequal(preproc_options.threshold,0) & isequal(preproc_options.binarize,0),
        [rand_assoc_matrix] = null_model_und_sign(mean_grand_cur_assoc_matrix,5,1); % % randomize undirected weighted network while preserving degree-, weight- and strength-distribution; Rubinov and Sporns (2011) Neuroimage, says OK to use for [0,1] matrices
        %shuffle edge weights (Stam 2009, Bolanos), only good for full (non-thresh) matrices
        %n=length(mean_grand_cur_assoc_matrix);
        %mat=triu(mean_grand_cur_assoc_matrix,1);
        %loc=find(mat>0); 
        %y=randperm(length(loc));
        %z=loc(y);
        %rand_assoc_matrix=zeros(n);
        %rand_assoc_matrix(z)=mat(loc);
        %rand_assoc_matrix=rand_assoc_matrix+triu(rand_assoc_matrix)';
      elseif isequal(preproc_options.threshold,1),
        [rand_assoc_matrix] = randmio_und_connected(mean_grand_cur_assoc_matrix,1); % randomize undir network, preserve degree dist and maintain connectedness
      end

      if isequal(preproc_options.binarize,0), % weighted

        [Crand] = clustering_coef_wu(rand_assoc_matrix); % clustering coefficient (elec x 1)
        meanCrand = mean(Crand); % mean weighted clustering coeff (scalar)
        rand_clus_coef(jj,comp_counter) = meanCrand;

        [Lrand] = weight_conversion(rand_assoc_matrix, 'lengths'); % convert to connection lengths matrix
        [Drand Brand] = distance_wei(Lrand); % distance matrix
        [char_path_length_rand] = charpath(Drand); % char path length
        rand_char_path_length(jj,comp_counter) = char_path_length_rand;

        Eglob_rand = efficiency_wei(rand_assoc_matrix);  % glob eff (scalar); inversely related to lambda (char path length)
        rand_Eglob(jj,comp_counter) = Eglob_rand;

      elseif isequal(preproc_options.binarize,1), % binary

        [C] = clustering_coef_bu(rand_assoc_matrix); % clustering coefficient (elec x 1)
        meanCrand = mean(C); % mean binary clustering coeff (scalar)
        rand_clus_coef(jj,comp_counter) = meanCrand;

        [D] = distance_bin(rand_assoc_matrix); % distance matrix
        [char_path_length_rand] = charpath(D); % lambda = char path length (scalar)
        rand_char_path_length(jj,comp_counter) = char_path_length_rand;

        Eglob_rand = efficiency_bin(rand_assoc_matrix);  % glob eff (scalar); inversely related to lambda (char path length)
        rand_Eglob(jj,comp_counter) = Eglob_rand;

      end

    end

    disp(['    Finished computing randomized networks and null graph measures']); disp(blanks(1)');

    % compute permuted small-world-network-ness, from MXC 2014
    for permi=1:num_rand_null_networks
      whichnetworks2use = randsample(1:num_rand_null_networks,2);
      swn_permutations(permi,comp_counter)  = (rand_clus_coef(whichnetworks2use(1),comp_counter)/rand_clus_coef(whichnetworks2use(2),comp_counter)) / (rand_char_path_length(whichnetworks2use(1),comp_counter)/rand_char_path_length(whichnetworks2use(2),comp_counter));
    end

  end

  % display thresholds
  if isequal(preproc_options.threshold,1),
    disp(['    Thresholds used: ']);
    disp([thresh]'); disp(blanks(1)');
  end

  for catcode_loop_counter = 1:length(cur_erptfd_assoc_matrix.catcodes),

    % get current sbj/catcode-specific network matrix
    if isnumeric(cur_erptfd_assoc_matrix.catcodes(catcode_loop_counter)),
      cur_assoc_matrix = cur_erptfd_assoc_matrix.(['assoc_matrix_' num2str(cur_erptfd_assoc_matrix.catcodes(catcode_loop_counter))]);
      cur_catcode = num2str(cur_erptfd_assoc_matrix.catcodes(catcode_loop_counter));
    elseif ~isnumeric(cur_erptfd_assoc_matrix.catcodes(catcode_loop_counter)),
      cur_assoc_matrix = cur_erptfd_assoc_matrix.(['assoc_matrix_' char(cur_erptfd_assoc_matrix.catcodes(catcode_loop_counter))]);
      cur_catcode = char(cur_erptfd_assoc_matrix.catcodes(catcode_loop_counter));
    end

    disp(['    Begin computing graph measures for condition: ' cur_catcode]); disp(blanks(1)');

    % make symmetric
    for fi=1:size(cur_assoc_matrix,3)
      for ti=1:size(cur_assoc_matrix,4)
        temp = squeeze(cur_assoc_matrix(:,:,fi,ti));
        temp = temp + triu(temp)';  % make symmetric matrix
        cur_assoc_matrix(:,:,fi,ti) = temp;
      end
    end

    for comp_counter = 1:length(comps_defs),

      name    = char(comps_defs(comp_counter).name); 
      c_timesbin = comps_defs(comp_counter).timesbin;
      c_timeebin = comps_defs(comp_counter).timeebin;
      c_freqsbin = comps_defs(comp_counter).freqsbin;
      c_freqebin = comps_defs(comp_counter).freqebin;

      % compute mean assoc_matrix across specified TF ROI
      mean_cur_assoc_matrix = mean(cur_assoc_matrix(:,:,:,cur_erptfd_assoc_matrix.tbin+c_timesbin:cur_erptfd_assoc_matrix.tbin+c_timeebin),4);
      mean_cur_assoc_matrix = mean(mean_cur_assoc_matrix(:,:,c_freqsbin:c_freqebin),3);

    % threshold?
      if isequal(preproc_options.threshold, 1) && isequal(preproc_options.thresh_type, '1SD>median'), % 1 SD >= median, from MXC (2011, 2012, 2014)
       %mean_cur_assoc_matrix(mean_cur_assoc_matrix<thresh(comp_counter))=0; % set to 0 all values under threshold
        mean_cur_assoc_matrix = threshold_absolute(mean_cur_assoc_matrix, thresh(comp_counter)); % set to 0 all values under threshold
        disp(['      Proportion of edges retained for component ' name ': ' num2str(length(nonzeros(mean_cur_assoc_matrix))/((size(mean_cur_assoc_matrix,1)*size(mean_cur_assoc_matrix,2))-size(mean_cur_assoc_matrix,1)))]); disp(blanks(1)');
      elseif isequal(preproc_options.threshold, 1) && isequal(preproc_options.thresh_type, 'proportion'),
        mean_cur_assoc_matrix = threshold_proportional(mean_cur_assoc_matrix, thresh(comp_counter)); % proportion of strongest weights
        disp(['      Proportion of edges retained for component ' name ': ' num2str(length(nonzeros(mean_cur_assoc_matrix))/((size(mean_cur_assoc_matrix,1)*size(mean_cur_assoc_matrix,2))-size(mean_cur_assoc_matrix,1)))]); disp(blanks(1)');
      end

      % binarize?
      if isequal(preproc_options.binarize,1),
        mean_cur_assoc_matrix = weight_conversion(mean_cur_assoc_matrix, 'binarize');
      end

      graph_measures.components.([name '_threshold_' cur_catcode])(main_loop_counter) = thresh(comp_counter);
      graph_measures.components.([name '_proportion_edges_retained_' cur_catcode])(main_loop_counter) = [length(nonzeros(mean_cur_assoc_matrix))/((size(mean_cur_assoc_matrix,1)*size(mean_cur_assoc_matrix,2))-size(mean_cur_assoc_matrix,1))];

      % eval for binary or weighted
      if isequal(preproc_options.binarize,0), % weighted

        disp(['      Computing weighted measures']); disp(blanks(1)');

        % path-based measures
        [L] = weight_conversion(mean_cur_assoc_matrix, 'lengths'); % convert to connection lengths matrix
        [D B] = distance_wei(L); % distance matrix
        [char_path_length] = charpath(D); % lambda = char path length (scalar)
        [BC] = betweenness_wei(L); % betweenness centrality (elec x 1)
        graph_measures.components.([name '_wei_char_path_length_' cur_catcode])(main_loop_counter) = char_path_length;
        graph_measures.components.([name '_wei_node_betweenness_centrality_' cur_catcode])(:,main_loop_counter) = BC;

        % node clustering and strength
        [C] = clustering_coef_wu(mean_cur_assoc_matrix); % clustering coefficient (elec x 1)
        meanC = mean(C); % mean weighted clustering coeff (scalar)
        [str] = strengths_und(mean_cur_assoc_matrix); % component node strength (elec x 1)
        graph_measures.components.([name '_wei_node_cc_' cur_catcode])(:,main_loop_counter) = C;
        graph_measures.components.([name '_wei_cc_' cur_catcode])(main_loop_counter) = meanC;
        graph_measures.components.([name '_wei_node_strength_' cur_catcode])(:,main_loop_counter) = str';

        % efficiency measures
        Eglob = efficiency_wei(mean_cur_assoc_matrix);  % glob eff (scalar); inversely related to lambda (char path length)
        Eloc = efficiency_wei(mean_cur_assoc_matrix,1); % loc eff (elec x 1); analogous to clustering coeff
        graph_measures.components.([name '_wei_Eglob_' cur_catcode])(main_loop_counter) = Eglob;
        graph_measures.components.([name '_wei_node_Eloc_' cur_catcode])(:,main_loop_counter) = Eloc;

        % weighted CC, char path length, and Eglob normalized against random matrices
        gamma = (meanC/mean(rand_clus_coef(:,comp_counter)));
        lambda = (char_path_length/mean(rand_char_path_length(:,comp_counter)));
        nEglob = (Eglob/mean(rand_Eglob(:,comp_counter)));
        graph_measures.components.([name '_wei_gamma_' cur_catcode])(main_loop_counter)  = gamma;
        graph_measures.components.([name '_wei_lambda_' cur_catcode])(main_loop_counter) = lambda;
        graph_measures.components.([name '_wei_nEglob_' cur_catcode])(main_loop_counter) = nEglob;

        % compute small-world (sigma)
        swn_real       = gamma/lambda;
        graph_measures.components.([name '_wei_sigma_' cur_catcode])(main_loop_counter) = swn_real;

        % z-value of real sigma compared to permuted sigma value
        swn_z = (swn_real-mean(swn_permutations(:,comp_counter)))/std(swn_permutations(:,comp_counter));
        graph_measures.components.([name '_wei_sigma_z_' cur_catcode])(main_loop_counter) = swn_z;

      elseif isequal(preproc_options.binarize,1), % binary

        disp(['      Computing binary measures']); disp(blanks(1)');

        % path-based measures
        [D] = distance_bin(mean_cur_assoc_matrix); % distance matrix
        [char_path_length] = charpath(D); % lambda = char path length (scalar)
        [BC] = betweenness_bin(mean_cur_assoc_matrix); % betweenness centrality (elec x 1)
        graph_measures.components.([name '_bin_char_path_length_' cur_catcode])(main_loop_counter) = char_path_length;
        graph_measures.components.([name '_binnode_betweenness_centrality_' cur_catcode])(:,main_loop_counter) = BC';

        % node clustering and strength
        [C] = clustering_coef_bu(mean_cur_assoc_matrix); % clustering coefficient (elec x 1)
        meanC = mean(C); % mean binary clustering coeff (scalar)
        [id, od, deg] = degrees_dir(mean_cur_assoc_matrix); % component node degree (elec x 1)
        graph_measures.components.([name '_bin_node_cc_' cur_catcode])(:,main_loop_counter) = C;
        graph_measures.components.([name '_bin_cc_' cur_catcode])(main_loop_counter) = meanC;
        graph_measures.components.([name '_bin_node_degree_' cur_catcode])(:,main_loop_counter) = deg';

        % efficiency measures
        Eglob = efficiency_bin(mean_cur_assoc_matrix);  % glob eff (scalar); inversely related to lambda (char path length)
        Eloc = efficiency_bin(mean_cur_assoc_matrix,1); % loc eff (elec x 1); analogous to clustering coeff
        graph_measures.components.([name '_bin_Eglob_' cur_catcode])(main_loop_counter) = Eglob;
        graph_measures.components.([name '_bin_node_Eloc_' cur_catcode])(:,main_loop_counter) = Eloc;

        % binary CC, char path length, and Eglob normalized against random matrices
        gamma = (meanC/mean(rand_clus_coef(:,comp_counter)));
        lambda = (char_path_length/mean(rand_char_path_length(:,comp_counter)));
        nEglob = (Eglob/mean(rand_Eglob(:,comp_counter)));
        graph_measures.components.([name '_bin_gamma_' cur_catcode])(main_loop_counter)  = gamma;
        graph_measures.components.([name '_bin_lambda_' cur_catcode])(main_loop_counter) = lambda;
        graph_measures.components.([name '_bin_nEglob_' cur_catcode])(main_loop_counter) = nEglob;

        % compute small-world (sigma)
        swn_real       = gamma/lambda;
        graph_measures.components.([name '_bin_sigma_' cur_catcode])(main_loop_counter) = swn_real;

        % z-value of real sigma compared to permuted sigma value
        swn_z = (swn_real-mean(swn_permutations(:,comp_counter)))/std(swn_permutations(:,comp_counter));
        graph_measures.components.([name '_bin_sigma_z_' cur_catcode])(main_loop_counter) = swn_z;

      end % end weighted/binary loop

    end % end comp loop

    disp(['    Finshed computing graph measures for condition: ' cur_catcode]); disp(blanks(1)');

  end % end catcode loop

  disp(['  Finished computing graph measures for subname: ' char(erptfd_assoc_matrix(main_loop_counter).subname)]); disp(blanks(1)');

end % end main loop

graph_measures.elecnames = cellstr([cur_erptfd_assoc_matrix.elecnames]);
graph_measures.tbin      = cur_erptfd_assoc_matrix.tbin ;

disp(['  Saving out structure as : ' outfilename]);

save(outfilename,'graph_measures');


% BEGIN global measure tab delim file export

% make subname 
subname = cell(length(graph_measures.subnum),1);
for j = 1:length(subname),
  subname(j,:) = cellstr(graph_measures.subnames(graph_measures.subnum(j),:));
end
subname = char(subname);

data_out(1).var = {'subname'              }; data_out(1).name = {'subname' } ;
data_out(2).var = {'graph_measures.subnum'}; data_out(2).name = {'subnum'  } ;

compnames = fieldnames(graph_measures.components);
result_scalar = cellfun(@isempty,(regexp(compnames','node'))); % global measures (scalars)
compnames = compnames(result_scalar);

startnum = length(data_out);
comps_startnum = startnum+1; 
% add components 
for cur_cs = 1:length(compnames),
  data_out(startnum+cur_cs).var  = eval(['{''graph_measures.components.' char(compnames(cur_cs)) '''};' ]);
  data_out(startnum+cur_cs).name = {char(compnames(cur_cs))};
end

% build data matrix 
comps_mat = zeros(length(graph_measures.subnum),length(compnames)); 
for cur_cs = 1:length(compnames),
  comps_mat(:,cur_cs) =  eval([' graph_measures.components.' char(compnames(cur_cs)) ';' ]); 
end

% determine type of data in each variable (e.g. char vs numeric) 
for j=1:length(data_out),
  eval(['cur_var = ' char(data_out(j).var) ';']);
  if      ischar(cur_var),
    data_out(j).type = {'string'};
  elseif  isnumeric(cur_var),
      data_out(j).type = {'number'};
  end
end

% open file for writing 
fid=fopen([outfilename '.dat'],'w');

% write variable names as first line
for j=1:length(data_out)-1,
  fprintf(fid,'%s\t',char(data_out(j).name));
end
fprintf(fid,'%s'  ,char(data_out(j+1).name));
fprintf(fid,'\n');

% write data loop - row by row
for q=1:length(graph_measures.subnum)
  % write header vars
  for j=1:comps_startnum-1, % length(data_out)-1,
     eval(['cur_var = ' char(data_out(j).var) ';']);

    switch char(data_out(j).type),
      case 'string',  fprintf(fid,'%s\t',   cur_var(q,:));
      case 'integer', fprintf(fid,'%d\t',   cur_var(q,:));
      case 'number',  fprintf(fid,'%0.8g\t',cur_var(q,:)); end
    end

    % write comps matrix
    if (length(data_out)-(comps_startnum-1)) > 1,
      fprintf(fid,'%0.8g\t',comps_mat(q,1:end-1));
    end
      fprintf(fid,'%0.8g\n',comps_mat(q,end));

  end

fclose(fid);

% END global measure tab delim file export


% BEGIN local measure export

% make subname 
subname = cell(length(graph_measures.subnum),1);
for j = 1:length(subname),
  subname(j,:) = cellstr(graph_measures.subnames(graph_measures.subnum(j),:));
end
subname = char(subname);

% make elecname 
elecname = cell(length(graph_measures.elec),1);
for j = 1:length(elecname),
  elecname(j,:) = cellstr(graph_measures.elecnames(graph_measures.elec(j),:));
end
elecname = char(elecname);

local_network.subname  = repmat(blanks(size(subname,2)),[length(graph_measures.elec)*size(subname,1)],1);
local_network.subnum   = repmat(nan([length(graph_measures.elec)*size(subname,1)],1),1);
local_network.elec     = repmat(nan([length(graph_measures.elec)*size(subname,1)],1),1);
local_network.elecname = repmat(blanks(size(elecname,2)),[length(graph_measures.elec)*size(subname,1)],1);

for cur_sub = 1:size(subname,1)
  cur_idx = [(((length(graph_measures.elec)*(cur_sub-1))+1):(length(graph_measures.elec)*cur_sub))]';
  local_network.subname(cur_idx,:)  = repmat(subname(cur_sub,:),length(graph_measures.elec),1);
  local_network.subnum(cur_idx,:)   = repmat(graph_measures.subnum(cur_sub,:),length(graph_measures.elec),1);

  local_network.elecname(cur_idx,:) = elecname;
  local_network.elec(cur_idx,:)     = graph_measures.elec(:,cur_sub);
end

compnames = fieldnames(graph_measures.components);
result_vector = ~cellfun(@isempty,(regexp(compnames','node'))); % local measures (vectors)
compnames = compnames(result_vector);

% build data matrix 
for cur_cs = 1:length(compnames),
  temp_comp = eval([' graph_measures.components.' char(compnames(cur_cs)) ';' ]);
  local_network.components.([char(compnames(cur_cs))]) = reshape(temp_comp,length(graph_measures.elec)*length(graph_measures.subnum),1);
end

save([outfilename '_local_network'],'local_network');

% END local network export