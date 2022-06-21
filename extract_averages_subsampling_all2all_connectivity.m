function [erptfd_assoc_matrix] = extract_averages_subsampling_all2all_connectivity(innamebeg,subnames,innameend,outname,rsrate,domain,blsms,blems,startms,endms,P1,P2,OPTIONS,AT,catcodes2extract,elecs2extract,verbose),

%  function [erp, subsampling] = extract_averages_subsampling_all2all_connectivity(innamebeg,subnames,innameend,outname,rsrate,domain,blsms,blems,startms,endms,P1,P2,OPTIONS,AT,catcodes2extract,elecs2extract,verbose),
%
%  computes the interchannel phase-locking value (PLV) for every bivariate electrode pair (all2all connectivity)
% 
%  this code computes the requested time-frequency representation, and then computes the bivariate phase synchrony
%  value for every possible electrode pair for each point in time-frequency space. the function returns a structure (erptfd_assoc_matrix)
%  which contains a matrix (assoc_matrix_CATCODE) of all2all PLVs for each requested catcode. it uses the subsampling method to compute
%  the condition averages; please refer to the extract_averages_subsampling documentation for the required arguments (nothing additional is
%  required within this function).
%
%  erptfd_assoc_matrix can be passed to compute_graph_theoretic_measures to calculate various network statistics based on graph theory
%    see compute_graph_theoretic_measures for more information
%
%  PLV is calculated by comparing the phase values between two electrodes within each trial for each TF point (Lachaux et al., 1999, HBM)
%  these PLVs are then averaged across trials to calculate the phase consistency between two signals/electrodes
%
%  Jeremy Harper, UMN, 07/11/15

  if verbose >=1, disp([mfilename ': START']); 
     if ~isempty(outname), disp(['  Output file -- ' outname ]); end 
  end 

  % set timers 
  sub_clock_total= clock;

  disp(blanks(1)');

  % check if an external static_sets has been specified
  if (~isfield(OPTIONS.subsampling, 'static_sets')) || ...    % static_sets exists?
     (isfield(OPTIONS.subsampling, 'static_sets') && isempty(OPTIONS.subsampling.static_sets)) || ...    % static_sets empty?
     (isfield(OPTIONS.subsampling, 'static_sets') && isnumeric(OPTIONS.subsampling.static_sets) && isequal(OPTIONS.subsampling.static_sets,0)),    % static_sets a var but empty?
    subsampling.static_sets.subname        = [];
    subsampling.static_sets.sweeps         = [];
    OPTIONS.subsampling.static_sets        = 0;
    disp(['     MESSAGE: ' mfilename ' - creating new static_sets']); disp(blanks(1)'); 
  elseif isfield(OPTIONS.subsampling, 'static_sets'),
    if isfield(OPTIONS.subsampling, 'static_sets') & isstruct(OPTIONS.subsampling.static_sets), % if struct, unpack
      disp(['     MESSAGE: ' mfilename ' - using predefined subsampling parameters']);
      subsampling = OPTIONS.subsampling.static_sets;
      if ~isequal(subsampling.catcodes, catcodes2extract) | ~isequal(subsampling.num_of_subsamples, OPTIONS.subsampling.num_of_subsamples) | ~isequal(subsampling.subsample_length, OPTIONS.subsampling.subsample_length),
        disp(['ERROR: catcodes|num_of_subsamples|subsample_length in external file and OPTIONS.subsampling not equal; all must be equal. Exiting function.']); disp(blanks(1)');
        return;
      end
    elseif isfield(OPTIONS.subsampling, 'static_sets') & isstr(OPTIONS.subsampling.static_sets), % if string, load
      if exist([OPTIONS.subsampling.static_sets '.mat'],'file'),
        load(OPTIONS.subsampling.static_sets);
        disp(['     MESSAGE: ' mfilename ' - using predefined subsampling parameters']);
        if ~isequal(subsampling.catcodes, catcodes2extract) | ~isequal(subsampling.num_of_subsamples, OPTIONS.subsampling.num_of_subsamples) | ~isequal(subsampling.subsample_length, OPTIONS.subsampling.subsample_length),
          disp(['ERROR: catcodes|num_of_subsamples|subsample_length in external file and OPTIONS.subsampling not equal; all must be equal. Exiting function.']); disp(blanks(1)');
          return;
        end
      elseif ~exist([OPTIONS.subsampling.static_sets '.mat'],'file'),
        disp('ERROR: file specified in OPTIONS.subsampling.static_sets does not exist. Exiting function.'); disp(blanks(1)');
        return;
      end
    end
  end

  % check is an external predetermined_lengths has been specified % ADDED JH 081913
  if isfield(OPTIONS.subsampling, 'predetermined_lengths'),
    if isfield(OPTIONS.subsampling, 'predetermined_lengths') & isstruct(OPTIONS.subsampling.predetermined_lengths), % if struct, unpack
      disp(['     MESSAGE: ' mfilename ' - using predetermined lengths parameters']);
      predetermined_lengths = OPTIONS.subsampling.predetermined_lengths;
      if ~isequal({predetermined_lengths.subnames}', subnames),
        disp(['ERROR: subnames in external file and OPTIONS.predetermined_lengths not equal; all must be equal. Exiting function.']); disp(blanks(1)');
        return;
      end
    elseif isfield(OPTIONS.subsampling, 'predetermined_lengths') & isstr(OPTIONS.subsampling.predetermined_lengths), % if string, load
      if exist([OPTIONS.subsampling.predetermined_lengths '.mat'],'file'),
        load(OPTIONS.subsampling.predetermined_lengths);
        disp(['     MESSAGE: ' mfilename ' - using predetermined lengths parameters']);
        if ~isequal({predetermined_lengths.subnames}', subnames),
          disp(['ERROR: subnames in external file and OPTIONS.predetermined_lengths not equal; all must be equal. Exiting function.']); disp(blanks(1)');
          return;
        end
      elseif ~exist([OPTIONS.subsampling.predetermined_lengths '.mat'],'file'),
        disp('ERROR: file specified in OPTIONS.subsampling.predetermined_lengths does not exist. Exiting function.'); disp(blanks(1)');
        return;
      end
    end
  end

  % check for resampling parameters
  if (~isfield(OPTIONS,'subsampling')) ...   
     || ( ~isfield(OPTIONS.subsampling,'method') ...   
     | ~isfield(OPTIONS.subsampling,'subsample_length') ...   
     | ~isfield(OPTIONS.subsampling,'num_of_subsamples')),
    disp('ERROR: must specify resampling method, subsample_length, and num_of_subsamples'); 
    return;
  end

  if isnumeric(OPTIONS.subsampling.subsample_length),
    if length(OPTIONS.subsampling.subsample_length) > 1,
      if ~isequal(length(catcodes2extract), length(OPTIONS.subsampling.subsample_length));
        disp(['ERROR: length of catcodes and subsample_length not equal. Exiting function.']); disp(blanks(1)');
        return;
       end
    end
  end
  if length(OPTIONS.subsampling.num_of_subsamples) > 1,
    if ~isequal(length(catcodes2extract), length(OPTIONS.subsampling.num_of_subsamples));
      disp(['ERROR: length of catcodes and num_of_subsamples not equal. Exiting function.']); disp(blanks(1)');
      return;
     end
  end

  % if no boot_samples, set to default (0)
  if (~isfield(OPTIONS.subsampling, 'boot_samples')) || (isfield(OPTIONS.subsampling, 'boot_samples') && isempty(OPTIONS.subsampling.boot_samples)),
    OPTIONS.subsampling.boot_samples = 0;
  end

  % if no erp_averging, set to default (1)
  if (~isfield(OPTIONS.subsampling, 'erp_averaging')) || (isfield(OPTIONS.subsampling, 'erp_averaging') && isempty(OPTIONS.subsampling.erp_averaging)),
    OPTIONS.subsampling.erp_averaging = 1;
  end

  % evaluate for bootstrapping
  if OPTIONS.subsampling.boot_samples == 0,
    OPTIONS_BOOT.boot_samples = 0; % ADDED 051613 by JH .added OPTIONS_BOOT.boot_samples for aggregate_amnd_bootstrap_erp function which now needs OPTIONS.boot_samples parameter
    disp('     MESSAGE: no bootstrapping will be performed'); disp(blanks(1)');
  elseif OPTIONS.subsampling.boot_samples > 0,
    disp(['     MESSAGE: ' num2str(OPTIONS.subsampling.boot_samples) ' bootstrapped subsample(s) requested']); disp(blanks(1)');
    OPTIONS_BOOT.boot_samples = OPTIONS.subsampling.boot_samples;  % ADDED 051613 by JH .added OPTIONS_BOOT.boot_samples for aggregate_amnd_bootstrap_erp function which now needs OPTIONS.boot_samples parameter
  end

  % base vars 

  if exist('rsrate'          ,'var')==0, rsrate          =       0;      end
  if exist('P1'              ,'var')==0, P1              =       0;      end
  if exist('P2'              ,'var')==0, P2              =       0;      end
  if exist('XX'              ,'var')==0, XX              =       0;      end
  if exist('AT'              ,'var')==0, AT              =       'NONE'; end
  if exist('catcodes2extract','var')==0, catcodes2extract=       'ALL';  end
  if exist('elecs2extract'   ,'var')==0, elecs2extract   =       'ALL';  end
  if exist('verbose'         ,'var')==0, verbose         =       0;      end

  % handle catcodes == ALL. ADDED FROM EXTRACT_AVERAGES JH/EB 090513
    % evaluate catcodes2extract
  create_catcodes = extract_base_evaluate_2extract(catcodes2extract,'catcodes2extract');

  % determine single or multi subject inputs -- loads file if file-multi  
  extract_base_evaluate_singlemulti;

  % if cache_subject_files is wanted, make temp directory
  if isfield(OPTIONS.subsampling, 'cache_subject_files') & isequal(OPTIONS.subsampling.cache_subject_files, 1),
    mkdir('./','temp_sbjs');
    disp(['     MESSAGE: ' mfilename ' - Creating temp_sbjs directory ']); disp(blanks(1)');
    if ~isfield(OPTIONS.subsampling, 'temp_subject_suffix') || isempty(OPTIONS.subsampling.temp_subject_suffix),
      OPTIONS.subsampling.temp_subject_suffix = '_resampled';
    end
  elseif ~isfield(OPTIONS.subsampling, 'cache_subject_files') || isempty(OPTIONS.subsampling.cache_subject_files),
    OPTIONS.subsampling.cache_subject_files = 0;
  end

  % print catcodes
  disp(['     MESSAGE: ' mfilename ' - catcodes.name & text = ']); disp(blanks(1)');
  for zz = 1:length(catcodes2extract), disp(catcodes2extract(zz)); end

  % print subsample parameters
  disp(['     MESSAGE: ' mfilename ' - subsampling parameters = ']); disp(blanks(1)');
  disp(OPTIONS.subsampling);

  % print averaging parameter
  if isequal(OPTIONS.subsampling.erp_averaging, 1),
    disp(['     MESSAGE: ' mfilename ' - erp averaging performed']); disp(blanks(1)');
  elseif isequal(OPTIONS.subsampling.erp_averaging, 0),
    disp(['     MESSAGE: ' mfilename ' - no erp averaging']);  disp(blanks(1)');
  end

  % set subnum counter (used in extract_base_loadprep_data) 
  subnum = 0;
  % main loop
  for main_loop_counter=1:length(subnames(:,1)), 

    % set timers 
    sub_clock_current_subject = clock; 

    % load and prep data
    extract_base_loadprep_data; 

    % generate catcodes if not given (default ALL ttypes in file)  ADDED FROM EXTRACT_AVERAGES JH/EB 090513
    if create_catcodes>=1,
      clear ttypes catcodes catcodes2extract;  
      ttypes=unique(erp.ttype);
      for jj=1:length(ttypes), 
        cur_ttype = ttypes(jj);  
        eval(['catcodes2extract(' num2str(jj) ').name = ' num2str(cur_ttype) ';' ]);  
        eval(['catcodes2extract(' num2str(jj) ').text = ''erp.ttype==' num2str(cur_ttype) ''';' ]);
      end
    end

    erp_org = erp;
    clear erp;

    % loads minimums if user-defined, calculates if not
    cur_idx = strmatch(char(subnames(main_loop_counter,:)),{subsampling.static_sets.subname}','exact');
    catcode_minimums = [subsampling.static_sets(cur_idx).minimum]; 

    % begin catcode loop
    catcodes_loop_success = 1;  % ADDED JH 01.28 3:22
    for count_catcode = 1:length(catcodes2extract),

      if isnumeric(catcodes2extract(count_catcode).name),
        disp(['     MESSAGE: begin catcode: ' num2str(catcodes2extract(count_catcode).name) ]);
      else, 
        disp(['     MESSAGE: begin catcode: ' catcodes2extract(count_catcode).name]);
      end

      % define current subsample parameters mod JH 051613
      if ~isnumeric(OPTIONS.subsampling.subsample_length),  % edited JH 051613 to evaluate for string or numeric
        if isequal(OPTIONS.subsampling.subsample_length, 'LCD'),  % ADDED JH 051613, 
          cur_subsample_length = min(catcode_minimums);
        elseif isequal(OPTIONS.subsampling.subsample_length, 'predetermined'),  % ADDED JH 081913 to set length based on predetermined values individually for subjects.
          cur_idx = strmatch(char(subnames(main_loop_counter,:)),{predetermined_lengths.subnames}','exact');
          if length(predetermined_lengths(cur_idx).length) == 1,    % ADDED JH 090513 to allow for a row vector input of predetermined_lengths to apply individually across catcodes. if row vector, it must match length of catcodes
            cur_subsample_length = predetermined_lengths(cur_idx).length;
          elseif length(predetermined_lengths(cur_idx).length) > 1,
            if ~isequal(length(catcodes2extract), length(predetermined_lengths(cur_idx).length)),
              disp(['ERROR: length of catcodes and subsample_length not equal. Exiting function.']); disp(blanks(1)');
              return;
            end
            cur_subsample_length = predetermined_lengths(cur_idx).length(count_catcode);
          end
        end
      elseif isnumeric(OPTIONS.subsampling.subsample_length), 
        if length(OPTIONS.subsampling.subsample_length) == 1,
          cur_subsample_length  = OPTIONS.subsampling.subsample_length;  % sets cur_subsample_length
        elseif length(OPTIONS.subsampling.subsample_length) > 1,
          cur_subsample_length  = OPTIONS.subsampling.subsample_length(count_catcode);  % sets cur_subsample_length for specific catcode
        end
      end

      if length(OPTIONS.subsampling.num_of_subsamples) == 1,
        cur_num_of_subsamples = OPTIONS.subsampling.num_of_subsamples;
      elseif length(OPTIONS.subsampling.num_of_subsamples) > 1,
        cur_num_of_subsamples = OPTIONS.subsampling.num_of_subsamples(count_catcode);
      end

      disp(['     MESSAGE: cur_subsample_length = ' num2str(cur_subsample_length) ', cur_num_of_subsamples = ' num2str(cur_num_of_subsamples)]);

      % assign erp_temp with sweeps from current catcode 
      erp_temp = reduce_erp(erp_org,catcodes2extract(count_catcode).text);

      % test if there are enough sweeps for the cur_subsample_length, loop if not to forgo current catcode mod JH 051713; evaluates all catcode mins at once (saves time)
      if any(cur_subsample_length > catcode_minimums) | ... % if there aren't enough trials, remove current subject EDITED JH 051613 from >= to > for LCD method (won't subsample the lowest catcode since length = sweep count, but will others)
         ((isfield(OPTIONS.subsampling, 'min_trials') && ~isempty(OPTIONS.subsampling.min_trials)) && any(OPTIONS.subsampling.min_trials > catcode_minimums)), % ADDED JH 051613, evals for min_trials (e.g., if cur_subsample_length = 4, minimum sweep count across catcodes = 5, and min_trials = 6, this subject will be thrown out now [wouldn't in the past since min sweep count > cur_length]
        disp([' ']);  disp(['     WARNING: insufficient sweeps/trials for subsampling with current catcode -- skipping subject ']); disp([' ']);  
        catcodes_loop_success = 0;  % ADDED JH 01.28 3:22
       %subsampling.static_sets(subnum) = []; % ADDED JH 051713. Removes information for rejected subjects in subsampling (did not do this in the old code).
        if count_catcode > 1,% ADDED JH 021214. If count_catcode = 1, subsampling.static_sets(subnum) does not exist for current subject, so this removal is not needed and causes an index out-of-bounds error. If count_catocde is greater than 1, then subsampling.static_sets(subnum) does exist for current subject, and needs to be removed.
          subsampling.static_sets(subnum) = []; % ADDED JH 051713. Removes information for rejected subjects in subsampling (did not do this in the old code).
        end
        subnum = subnum - 1; 
        break % loop
      end  

      % define subsampling method parameters 
      % strip out current subject/catcode sweeps from static_sets
      cur_idx = strmatch(char(subnames(main_loop_counter,:)),{subsampling.static_sets.subname}','exact');  % ADDED JH 01.28 3:22 from subnum to main loop counter to index subnames properly; JH 02.24 added 'exact' to make sure correct subname was indexed
      disp(['     MESSAGE: loading predefined sweeps for subject: ' subsampling.static_sets(cur_idx).subname]);
      erp_temp.temp_resample.sweeps(:,:) = getfield(subsampling.static_sets(cur_idx).sweeps,['sweeps_' num2str(count_catcode)]);

     % run subsamling on the data (or not, if requested not to actually run it) 
      if isfield(OPTIONS.subsampling, 'erp_averaging') & isequal(OPTIONS.subsampling.erp_averaging, 1),

        % create index matrix for sweeps in each sample (0 = not included, 1 = included)
        erp_temp.temp_resample.trial_matrix = zeros(length(erp_temp.sweep),cur_num_of_subsamples);

        for cc = 1:cur_num_of_subsamples,

          % create new catcodes based on random_sampling
          catcodes_resamp(cc).name = cc; 
          catcodes_resamp(cc).text = ['(' catcodes2extract(count_catcode).text ')&erp.temp_resample.trial_matrix(:,' num2str(cc) ')==1'];

          % get current rand_sample sweeps and create cur_sweeps vector
          cur_sweeps2use = erp_temp.temp_resample.sweeps(:,cc);
          cur_sweeps     = zeros(size(erp_temp.sweep));

          for dd = 1:cur_subsample_length, % size of subsample
            cur_sweeps(erp_temp.sweep==[cur_sweeps2use(dd)])=1;
          end

          % populate matrix with cur_sweeps
          erp_temp.temp_resample.trial_matrix(:,cc) = cur_sweeps;

        end
      end

      % get TFR
      erp_temp.data = erp_temp.data(:,startbin:endbin); 
      erp_temp = generate_TFDs(erp_temp,domain,verbose);

      % transpose 3D data into 4D (elec x Hz x time x trl)
      for chani = 1:length(unique(erp_temp.elec)),
        cur_elec = erp_temp.elec==chani;
        erp_temp.bigmtx(chani,:,:,:) = shiftdim(erp_temp.data(cur_elec,:,:),1); % elec x Hz x time x trl
      end

      % initialize association matrix (elec x elec x Hz x time x subsample_iterations
      tfr_all2all = zeros(size(erp_temp.bigmtx,1),size(erp_temp.bigmtx,1),size(erp_temp.bigmtx,2),size(erp_temp.bigmtx,3),size(catcodes_resamp,2)); % elec x elec x Hz x time x num iterations

      % extract cur sweeps
      sweeps =  unique(erp_temp.sweep);imagesc

      % compute all2all connectivity (self-self connections are left as 0)
      disp(['     Begin computing all2all connectivity']);
      for chani=1:length(unique(erp_temp.elec)),
        for chanj=chani+1:length(unique(erp_temp.elec)),
          xsd = squeeze(erp_temp.bigmtx(chani,:,:,:).*conj(erp_temp.bigmtx(chanj,:,:,:))); % hz x time x trial
          for ii = 1:size(erp_temp.temp_resample.sweeps,2),  % use resampled sweeps to index all2all and compute resampled PLV
            sub_sweeps = erp_temp.temp_resample.sweeps(:,ii); % subsample sweep index
            [junk,sweepidx]=ismember(sub_sweeps,sweeps); % index of sweeps to average for current subsample
            temp_xsd = xsd(:,:,sweepidx);
            tfr_all2all(chani,chanj,:,:,ii) = abs(mean(exp(1i*angle(temp_xsd)),3)); % ispc_trials, mean PLV
          end
        end
      end

      %  average resampled PLVs
      if isnumeric(catcodes2extract(count_catcode).text),
        erp_all_subs(main_loop_counter).(['assoc_matrix_' num2str(catcodes2extract(count_catcode).name)]) = mean(tfr_all2all,5); % elec x elec x Hz x time
      elseif ~isnumeric(catcodes2extract(count_catcode).text),
        erp_all_subs(main_loop_counter).(['assoc_matrix_' catcodes2extract(count_catcode).name]) = mean(tfr_all2all,5); % elec x elec x Hz x time
      end
      disp(['     Finished computing all2all connectivity']); disp(blanks(1)');

      % get tbin (stim/resp onset)
      tbin = erp_temp.tbin;

      clear erp_temp tfr_all2all xsd temp_xsd catcodes_resamp cur_sweeps2use cur_sweeps sub_sweeps sweeps sweepidx;

    end % end catcode loop

    erp_all_subs(main_loop_counter).subname   = [];
    erp_all_subs(main_loop_counter).subname   = subnames(main_loop_counter,:);
    erp_all_subs(main_loop_counter).subnum    = main_loop_counter; 
    erp_all_subs(main_loop_counter).catcodes  = {catcodes2extract.name};
    erp_all_subs(main_loop_counter).elec      = unique(erp_org.elec);
    erp_all_subs(main_loop_counter).elecnames = erp_org.elecnames(unique(erp_org.elec),:);
    erp_all_subs(main_loop_counter).tbin      = tbin;
    erp_all_subs(main_loop_counter).rsrate    = rsrate;
    erp_all_subs(main_loop_counter).TF_params = [domain];

    clear erp_org catcode_minimums tbin;

    disp(['     Current subject processing time (secs):	' num2str(etime(clock,sub_clock_current_subject)) ]); 
    disp(['     Total processing time (secs): 		' num2str(etime(clock,sub_clock_total)) ]);
    disp(blanks(1)');

  end % end main loop

  erptfd_assoc_matrix           = erp_all_subs;

 % save outfile 
  if isempty(outname) == 0,
    if verbose >= 1, disp(['Writing file: ' outname ]); end
    save(outname,'erptfd_assoc_matrix');
  end
