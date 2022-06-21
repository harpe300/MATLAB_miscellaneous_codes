
function [erp] = subtract_averaged_ERP(erp,average_erp_filename,catcodes2extract);

% [erp] = subtract_averaged_ERP(erp,average_erp_filename,catcodes2extract);
% 
%  Should be placed in the preprocessing stream
%
%  Loads averaged ERP, and for each subject/electrode/catcode, 
%    this code will subtract the averaged ERP from the trial-level signals
%
%  The calculation of total power TF on the ERP-subtracted trial-level signals is the non-phase-locked power
%
%  required parameters: 
%  
%    erp                    -- trial-level time domain erp structure
%    average_erp_filename   -- filename of averaged ERP structure
%                         NOTE: subnames, catcodes, rsrate, and trial length must match across the ERP structure and NPL script 
%    catcodes2extract       -- catcode structure
%                         NOTE: must match catcodes used in the averaged ERP script
%
%  Refs: Kalcher and Pfurtscheller, 1995; Cohen and Donner, 2013
%
%  Jeremy Harper, UMN, 11/01/2016

  avg_erp = load(average_erp_filename); avg_erp = avg_erp.erp;
  disp(['  Loading ERP file: ' average_erp_filename]);

  cur_subname = erp.subs.name;

  % logical vector of the current subname from the large averaged ERP structure
  reduce_vec = avg_erp.subnum==strmatch(cellstr(cur_subname),avg_erp.subs.name);

  % reduce the large averaged ERP structure to the current subname
  avg_erp = reduce_erp(avg_erp,reduce_vec);
  disp(['  Reducing averaged ERP to subname: ' char(avg_erp.subs.name(unique(avg_erp.subnum)))]);

  if size(avg_erp.data,2) ~= size(erp.data,2),
    disp(['  WARNING: Length of erp.data and avg_erp.data unequal, aborting']);
    disp(['           erp.data and avg_erp.data must have same number of timepoints']);
    return;
  end

  for count_catcode = 1:length(catcodes2extract),

    % averaged ERP for the current catcode
    temp_avg_erp = reduce_erp(avg_erp,avg_erp.stim.catcodes==catcodes2extract(count_catcode).name);

    % logical vector of the current catcodes ttypes for the single-trial ERP
    eval(['temp_erp_ttype_reduce_vec = ' char(catcodes2extract(count_catcode).text) ';']); 

    % subtract the conditioned-average ERP from the single-trials by electrode
    erp.data(temp_erp_ttype_reduce_vec,:) = erp.data(temp_erp_ttype_reduce_vec,:) - repmat(temp_avg_erp.data,length(unique(erp.sweep(temp_erp_ttype_reduce_vec,:))),1);

    clear temp_avg_erp temp_erp_ttype_reduce_vec

  end

  disp(['  Subtraction of averaged ERP from trial-level signals complete. ']);

  clear avg_erp cur_subname
