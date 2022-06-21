function [raw, prctile95rndexplained, medianrndexplained, meanrndexplained] = parallel_analysis(PAparams);

%function [raw, prctile95rndexplained, medianrndexplained, meanrndexplained] = parallel_analysis(PAparams);
%
% This function loads time series/time-frequency structures, 
% structures them for SVD, and then performs parallel analysis 
% using column-wise permutation to empirically determine the number
% of components to retain
%
% example: 
%
%  PAparams.dataset_name   = 'adbf0_gng_stim_avg'; % target dataset
%  PAparams.rs             =  32;                  % resampling rate
%  PAparams.runtype        =  'pcatfd';            % pca or pcatfd
%  PAparams.verbose        =  1;                   % verbose (0=no, 1=yes)
%  PAparams.domain         = 'tfd';                % time or tfd
%  PAparams.num_iterations =  500;                 % number of iterations for parallel analysis
%
%  %runtype specific variables (these need to match the parameters from the target dataset)
%  switch PAparams.runtype
%  case {'pcatfd'}
%    PAparams.timebinss    = 32;
%    PAparams.startbin     =  0;
%    PAparams.endbin       = 32;
%    PAparams.fqbins       = 32;
%    PAparams.fqstartbin   = 1;
%    PAparams.fqendbin     = 25;
%    PAparams.fqamp01      = 'F';
%    PAparams.dmx          = 'acov';
%    PAparams.rot          = 'vmx';
%    PAparams.fa           = [];
%    PAparams.TFDbaseline.start =     -200;
%    PAparams.TFDbaseline.end   =       -1;
%  case {'pca'}
%    PAparams.startbin     = 0;
%    PAparams.endbin       = 64;
%    PAparams.dmx          = 'acov';
%    PAparams.rot          = 'vmx';
%    PAparams.fa           = [];
%  end
%
%  parallel_analysis_pca(PAparams);

  % load erp|erptfd structure (must have them computed from pcatfd_startup before continuing)

  IDvars = PAparams; clear PAparams;
  base_IDvars_unpack;

  if isequal(domain,'time'),
    erpname =  [dataset_name '_' (num2str(rs))];
    load(erpname);
    disp(['  Loading: ' erpname]);
  elseif isequal(domain,'tfd'),
    erpname =  [dataset_name '_' (num2str(rs))];
    load(erpname);
    disp(['  Loading: ' erpname]);
    erptfdname =  [dataset_name '_' (num2str(rs)) '_t' (num2str(timebinss)) 'f' (num2str(fqbins))];
    load(erptfdname);
    disp(['  Loading: ' erptfdname]);
  end

  % structure for SVD; below is from base_loaddata (cannot call function, as it has extra pieces we don't need)

  if isequal(domain,'time'),

    % define erp bins if undefined 
    if startbin==0 & endbin==0,
      startbin   = 1; endbin   = length(erp.data(1,:))-erp.tbin;
    end
  % display N in set 
    if verbose > 0, [disptrials,disptimebins] = size(erp.data); disp(['    ERP Data Matrix Size: ' num2str(disptrials) ' (rows/trials) X ' num2str(disptimebins) ' (bins)'  ]); clear disprows disptrials; end % TFDs 

  elseif isequal(domain,'tfd'),
  
    % if erptfd data is from newer structured erptfd (>=0.9.5), erptfd = erptfd.data from erptfd structure 
    if isstruct(erptfd),
      erptfd = erptfd.data;
    end
    % if erptfd data is from neweer structured erptfd (>=0.9.5), shift dimensions to old style for rest of toolbox processing 
    %   (This was also caused by a bug in 0.9.5 for a while where newer structured erptfd.data was in incorrect old-style toolbox erptfd matrix format. 
    %    Thus, check if format is old -- regarless of whether it came from structured or not -- and shift dimensions if needed.)   
    if ~isequal(size(erptfd,3),size(erp.elec,1)),
      erptfd = shiftdim(erptfd,1);
    end

    % prep tfds for current RunID
    base_loaddata_prepvars_from_tfd;
    SETvars.TFDbaseline = IDvars.TFDbaseline;
    base_loaddata_prep4run_tfd;

    % display N in set from base_loaddata code
    if verbose > 0, 
      [dispfreqbins,disptimebins,disptrials] = size(erptfd); 
      disp(['    TFD Data Matrix Size: ' ... 
            num2str(disptrials) ' (rows/trials) X ' ...  
             num2str(dispfreqbins*disptimebins) ... 
           ' (TFD: ' num2str(dispfreqbins) ' freq X ' num2str(disptimebins) ' time)' ]); 
      clear dispfreqbins disptimebins disptrials 
    end

  end

  % Parse variables that require values stored in data file, but AFTER processing for run  

  % define conversion factors 
  switch domain
    case 'time', 
      unit2binfactor = erp.samplerate/1000;
      SETvars.bin2unitfactor= 1000/erp.samplerate;  
    case 'freq', 
      unit2binfactor = length(erp.data(1,:))/fix(erp.samplerate/2); 
      SETvars.bin2unitfactor= fix(erp.samplerate/2)/length(erp.data(1,:)); 
  end 

  % end loaddata code

  % start pcatfd_get_components (only a portion)

  if isequal(domain,'time'),

    erpstart = startbin+erp.tbin; %round(erp.tbin*(rs/erp.samplerate));
    erpend   = endbin  +erp.tbin; %round(erp.tbin*(rs/erp.samplerate));
    disp(['']);
    disp(['Size of PCA matrix -- Rows:  ' num2str(length(erp.elec)) ...
          '     Cols:   ' num2str(length(erp.data(1,erpstart:erpend))) ]);
    [org.P, org.LATENT,org.EXPLAINED,fa,retval] = base_pcasvd(IDvars,SETvars,erp.data(:,erpstart:erpend));
    if retval == 0, disp('ERROR: running decomposition'); components=0; return; end 

    orgdata = erp.data(:,erpstart:erpend);

  elseif isequal(domain,'tfd'),

    % transform TFD surfaces to vectors for PCA  
    erptfdlong = base_mat2long(erptfd,fqendbin-fqstartbin+1);

    % run PCA 
    disp(['']);
    disp(['Size of PCA matrix -- Rows:  ' num2str(length(erp.elec)) ...
          ' Cols:       ' num2str(length(erptfdlong(1,:)))  ]);

    [org.P, org.LATENT,org.EXPLAINED,fa,retval] = base_pcasvd(IDvars,SETvars,erptfdlong);
    if retval == 0, disp('ERROR: running decomposition'); components=0; return; end


  end

  disp([' Begin Parallel Analysis']); disp(blanks(1)');

  if isequal(IDvars.dmx,'acov'),

    disp([' Association Covariance matrix decomposition']); 
    disp([' Permutate (column-wise) raw data matrix']); disp(blanks(1)');

  else,

    disp([' ERROR: Other decomposition types not yet implemented, acov is preferred anyways...']); disp(blanks(1)'); return;

  end

  for dd =  1:IDvars.num_iterations,

    [ncases,nvars] = size(erptfdlong);

    permmatrix = erptfdlong;
    for lupec = 1:nvars;
      for luper = 1:(ncases -1);
        k = fix( (ncases - luper + 1) * rand(1) + 1 )  + luper - 1;
        d = permmatrix(luper,lupec);
        permmatrix(luper,lupec) = permmatrix(k,lupec);
        permmatrix(k,lupec) = d;
      end
    end

    [m,n] = size(permmatrix);
    mdata         = mean(permmatrix);                     % mean
    data          = (permmatrix - mdata(ones(m,1),:));    % center
    data          = (permmatrix' * permmatrix) / (m-1);   % SSCP 
    [u,latent,pc] = svd(data,0);                          % decompose 
    randLATENT = diag(latent); 
    randEXPLAINED = 100 * randLATENT / sum(randLATENT); 

    rndSVsmatrix(:,dd)         = randLATENT;    % matrix of SV
    rndexplainedmatrix(:,dd)   = randEXPLAINED; % matrix of percentage explained

  end

  prctile95rndexplained = prctile(rndexplainedmatrix,95,2);
  medianrndexplained    = median(rndexplainedmatrix,2);
  meanrndexplained      = mean(rndexplainedmatrix,2);

  meanrndSVs      = mean(rndSVsmatrix,2);


  plot(org.EXPLAINED,'ob'); hold on;
  plot(medianrndexplained,'ok');
  plot(prctile95rndexplained,'or');
  title(' Variance Explained by PC');
  legend('actual','median','95prctile');

  disp([' First 15 %Explained Values']);
  disp(['   Actual MedianRnd 95Prctilernd MeanRnd']); disp(blanks(1)');
  disp([org.EXPLAINED(1:15) medianrndexplained(1:15) prctile95rndexplained(1:15) meanrndexplained(1:15) (1:15)']);


  % Alternate method using a random multivariate normal surrogate dataset to perform parallel analysis
  % Permutation is preferred given better statistical properties, but this is kept in for posterity
  %% begin creating rand mv norm vector with mu/sigma from raw data for PA
  %[m,n]         = size(orgdata);
  %morgdata      = mean(orgdata); % mean
  %stdorgdata    =  std(orgdata);  % std
  %for dd =  1:IDvars.num_iterations,
  %  newrandmatrix = zeros(m,n);
  %  if isequal(IDvars.dmx,'acov'),
  %    for ii = 1:n
  %     newrandmatrix(:,ii) = stdorgdata(ii).*randn(m,1)+morgdata(ii); % creates a matrix of randnorm variables with means/stddvs matching the actual data (used for aCOV matrix decomp)
  %    end
  %  elseif isequal(IDvars.dmx,'acor'),
  %    newrandmatrix = randn(m,n); % creates a matrix of randnorm variables with means/stddvs matching the actual data (used for aCOR matrix decomp)
  %  end
  %  evalc('[randPC, randLATENT, randEXPLAINED,fa,retval] = base_pcasvd(IDvars,SETvars,newrandmatrix)'); % evalc suppressess display text in window
  %  rndSVsmatrix(:,dd)         = randLATENT; % matrix of SV
  %  rndexplainedmatrix(:,dd)   = randEXPLAINED; % running avg of percentage explained
  %end
