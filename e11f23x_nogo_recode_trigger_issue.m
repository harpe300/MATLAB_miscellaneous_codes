% [OUTEEG] = nogo_recode_trigger_issue(INEEG, infilebeg, infileend) - recodes EEG.event.type trigger codes based on E-Prime edat record.
%
% Usage:
%   OUTEEG = nogo_recode_trigger_problem(INEEG);
%
% Input:
%   INEEG     - input EEG data structure
%   infilebeg - file prefix or path
%   infileend - file suffix
%
% Output:
%   OUTEEG - modified EEG data structure with corrected triggers
%
%
% Comments:
%   A problem was identified with the trial_list coding for the go/nogo task. This was fixed by JBH and SMM on 05/12/15. 
%   The correct response/end_of_trial triggers for Blocks 2 and 3 were not being read from the corresponding trial_lists.
%   Instead, E-Prime was recycling the correct response trigger vector from Block 1 for those blocks. 
%   Because of this, the correct response triggers for Blocks 2 and 3 were incorrect.
%   The program still coded overall correct/incorrect in the right way, but the stimulus/response trigger mapping was incorrect (e.g., [1, 13]  instead of [1, 11]).
%
%   This function recodes the triggers in EEG.event.type to accurately match the trial list coding.
%   The recoding is done for the correct response triggers. 
%   For example, all correct responses (11-14) following Go1 are recoded as 11s.
%   The result is checked against the task record produced by E-Prime (edat) to ensure that the recoding was performed correctly.
%
%   This should be used for the ES-FU4 (e11f23) and MF-FU6 (c11f34) data. AdBrain and ES-FU3 (e11f20) used a different E-Prime program and are OK.
%
% Jeremy Harper, UMN, 05/14/15
%
% CHANGELOG: MATLAB updated the dataset function, so it now properly handles underscores. changed Stim0x2ERESP to Stim_RESP in code. - JH 05.23.18 

function [EEG] = nogo_recode_trigger_issue(EEG, infilebeg, infileend)

% START

disp([' MESSAGE: begin fixing trigger issue ']); disp(blanks(1)');

subname = EEG.subject;  % get subID/subname

trials = dataset('File',[infilebeg subname infileend],'Delimiter','\t'); % read in edat data
trials.Stim_RESP(find(ismember(trials.Stim_RESP,'""')))={'0'}; % recode no resp ("") to 0
trials.RESP = str2num(char(trials.Stim_RESP)); % create numeric response vector for indexing below
%trials.Stim0x2ERESP(find(ismember(trials.Stim0x2ERESP,'""')))={'0'}; % recode no resp ("") to 0 % edited 05.23.18 - JH
%trials.RESP = str2num(char(trials.Stim0x2ERESP)); % create numeric response vector for indexing below % edited 05.23.18 - JH

  % create accuracy vector 0/1 (the one from edat isn't correct )
trials.ACC = nan(length(trials),1);
trials.ACC((trials.Nogo_t==1|trials.Nogo_t==2)&(trials.RESP==1|trials.RESP==11|trials.RESP==2|trials.RESP==22|trials.RESP==12|trials.RESP==21)) = 1; % correct go
trials.ACC((trials.Nogo_t==3|trials.Nogo_t==4)&(trials.RESP==0)) = 1; % correct nogo
trials.ACC((trials.Nogo_t==1|trials.Nogo_t==2)&(trials.RESP==0)) = 0; % incorrect go
trials.ACC((trials.Nogo_t==3|trials.Nogo_t==4)&(trials.RESP==1|trials.RESP==11|trials.RESP==2|trials.RESP==22|trials.RESP==12|trials.RESP==21)) = 0; % incorrect nogo

  % create response/accuracy trigger vector (accurate/corrected trigger mapping)
trials.ACC_t = nan(length(trials),1);
trials.ACC_t(trials.Nogo_t==1&trials.ACC==1) = 11; % correct go1
trials.ACC_t(trials.Nogo_t==1&trials.ACC==0) = 21; % incorrect go1
trials.ACC_t(trials.Nogo_t==2&trials.ACC==1) = 12; % correct go2
trials.ACC_t(trials.Nogo_t==2&trials.ACC==0) = 22; % incorrect go2
trials.ACC_t(trials.Nogo_t==3&trials.ACC==1) = 13; % correct nogo1
trials.ACC_t(trials.Nogo_t==3&trials.ACC==0) = 23; % incorrect nogo1
trials.ACC_t(trials.Nogo_t==4&trials.ACC==1) = 14; % correct nogo2
trials.ACC_t(trials.Nogo_t==4&trials.ACC==0) = 24; % incorrect nogo2

 % create response/accuracy trigger vector from the old inaccurate triggers sent by E-Prime
trials.ACC_t_org = nan(length(trials),1);
for jj = 1:length(trials),
  if trials.ACC(jj) == 1,
    trials.ACC_t_org(jj) = trials.CorrectResponse1_t(jj);
  elseif trials.ACC(jj) == 0,
    trials.ACC_t_org(jj) = trials.IncorrectResponse1_t(jj);
  end
end

%[trials.Nogo_t(1:20) trials.RESP(1:20) trials.ACC(1:20) trials.ACC_t(1:20) trials.ACC_t_org(1:20) trials.CorrectResponse1_t(1:20) trials.CorrectResponse2_t(1:20) trials.IncorrectResponse1_t(1:20) trials.IncorrectResponse2_t(1:20) trials.Block(1:20) trials.Trial(1:20)] % print out first 20 trials to visually inspect
%[trials.Nogo_t((end-20:end)) trials.RESP((end-20:end)) trials.ACC((end-20:end)) trials.ACC_t((end-20:end)) trials.ACC_t_org((end-20:end)) trials.CorrectResponse1_t((end-20:end)) trials.CorrectResponse2_t((end-20:end)) trials.IncorrectResponse1_t((end-20:end)) trials.IncorrectResponse2_t((end-20:end)) trials.Block((end-20:end)) trials.Trial((end-20:end))] % print out last 20 trials to visually inspect

 % create vectors interleaving stim/resp (this is used to compare with EEG.event.type)
stim_t = trials.Nogo_t(:)';          % row vector of stim trigs
org_resp_t = trials.ACC_t_org(:)';   % row vector of original resp trigs (wrong coding in E-Prime)
upd_resp_t = trials.ACC_t(:)';       % row vector of recoded resp trigs (right coding)
org_temp_vec = [stim_t;org_resp_t];  % stack row vectors
upd_temp_vec = [stim_t;upd_resp_t];  % stack row vectors
org_stimresp_vec = org_temp_vec(:);  % interleave (stim-resp-stim-resp)
upd_stimresp_vec = upd_temp_vec(:);  % interleave (stim-resp-stim-resp)

org_trigs = [{EEG.event.type}]';  % extract event triggers from EEG
org_trigs = str2num(char(org_trigs));
[task_trig_idx] = find(org_trigs~=-1&org_trigs~=-9&org_trigs~=170); % index of task-triggers, added 062015 JH
org_trigs = org_trigs(org_trigs~=-1&org_trigs~=-9&org_trigs~=170); % convert to numeric vector, remove -1s (1s markers), -9s (artifacts), and 170 (beg trig)

  % if there are missing triggers in EEG, find and remove from edat trial structure
  % system can miss response triggers for some reason (I've seen this before on different systems with this task - JH)
if length(org_trigs)<length(org_stimresp_vec),
  num_diff = (length(org_stimresp_vec)-length(org_trigs)); 
  disp([' WARNING: ' num2str(num_diff) ' triggers are missing between edat record and EEG']); disp(blanks(1)');
  org_trigs(end+1:end+num_diff) = 0; % append 0s to org_trigs vector to match lengths
  for kk = 1:num_diff,
    diff_idx = find([org_trigs-org_stimresp_vec]~=0,1);  % find first mismatch
    org_stimresp_vec(diff_idx) = []; % remove extra trigger
    upd_stimresp_vec(diff_idx) = []; % remove extra trigger
    org_trigs(end) = []; % remove trailing 0 
  end
end
if length(org_trigs)>length(org_stimresp_vec),  % added 060815 by JH to account for more triggers in EEG than full lines in edat (e.g., if e-prime/recording stopped mid-trial)
  num_diff = (length(org_trigs)-length(org_stimresp_vec));
  disp([' WARNING: ' num2str(num_diff) ' triggers are missing between edat record and EEG']); disp(blanks(1)');
  org_stimresp_vec(end+1:end+num_diff) = 0; % append 0s to org_trigs vector to match lengths
  for kk = 1:num_diff,
    diff_idx = find([org_stimresp_vec-org_trigs]~=0,1);  % find first mismatch
    org_trigs(diff_idx) = []; % remove extra trigger
   %EEG = pop_selectevent( EEG, 'omitevent',diff_idx,'deleteevents','on'); %remove extra trigger from EEG, calls eeg_checkset which resorts by latency and messes up EEG.event ordering
    EEG.event(diff_idx) = []; % remove extra trigger
    org_stimresp_vec(end) = []; % remove trailing 0 
  end
end

  % verify that the event triggers from EEG match the original from E-Prime
if ~isequal(org_trigs, org_stimresp_vec),
  disp([' WARNING: event triggers from EEG and edat do not match. aborting.']);
  EEG.event(:) = [];
  return;
end 

%[trials.Nogo_t trials.stim_resp trials.ACC trials.ACC_t org_trigs(1:2:end) org_trigs(2:2:end) ]

 % construct structure of new trigger mappings
clear REC;
REC(1).label   = {'11'}; REC(1).template   = {'1' '11'}; REC(1).position   = 2; % go
REC(2).label   = {'11'}; REC(2).template   = {'1' '12'}; REC(2).position   = 2;
REC(3).label   = {'11'}; REC(3).template   = {'1' '13'}; REC(3).position   = 2; 
REC(4).label   = {'11'}; REC(4).template   = {'1' '14'}; REC(4).position   = 2;
REC(5).label   = {'12'}; REC(5).template   = {'2' '12'}; REC(5).position   = 2; % go
REC(6).label   = {'12'}; REC(6).template   = {'2' '11'}; REC(6).position   = 2;
REC(7).label   = {'12'}; REC(7).template   = {'2' '13'}; REC(7).position   = 2; 
REC(8).label   = {'12'}; REC(8).template   = {'2' '14'}; REC(8).position   = 2;
REC(9).label   = {'13'}; REC(9).template   = {'3' '13'}; REC(9).position   = 2; % nogo
REC(10).label  = {'13'}; REC(10).template  = {'3' '11'}; REC(10).position  = 2;
REC(11).label  = {'13'}; REC(11).template  = {'3' '12'}; REC(11).position  = 2; 
REC(12).label  = {'13'}; REC(12).template  = {'3' '14'}; REC(12).position  = 2;
REC(13).label  = {'14'}; REC(13).template  = {'4' '14'}; REC(13).position  = 2; % nogo
REC(14).label  = {'14'}; REC(14).template  = {'4' '11'}; REC(14).position  = 2;
REC(15).label  = {'14'}; REC(15).template  = {'4' '12'}; REC(15).position  = 2; 
REC(16).label  = {'14'}; REC(16).template  = {'4' '13'}; REC(16).position  = 2;

  % recode correct respose triggers, edited by JH on 062015 from original to index triggers based on position in EEG.event.type
for ll = 1:length(REC),
  cur_rec = REC(ll).template; cur_rec_pos = REC(ll).position-1; cur_rec_label = REC(ll).label;  % get current recode information
  REC_idx = strfind(str2num(char({EEG.event(task_trig_idx).type}))',str2num(char(cur_rec))');  % index event triggers matching current template
  cur_task_trig_idx = task_trig_idx(REC_idx+cur_rec_pos); % get index for task-trigger to recode
  for mm = 1:length(REC_idx),
    EEG.event(cur_task_trig_idx(mm)).type = char(cur_rec_label);  % recode correct_response trigger with accurate trigger mapping 
  end
end
clear REC;

new_trigs = [{EEG.event.type}]';
new_trigs = str2num(char(new_trigs));
[new_task_trig_idx] = find(new_trigs~=-1&new_trigs~=-9&new_trigs~=170); % index of new task-triggers
new_trigs = new_trigs(new_trigs~=-1&new_trigs~=-9&new_trigs~=170); 

  % verify that the updated event triggers from EEG match the right/correct trial coding from edat
if ~isequal(new_trigs, upd_stimresp_vec),
  disp([' WARNING: event triggers from EEG and edat do not match. aborting.']);
  EEG.event(:) = [];
  return;
end 
  % verify that the updated event trigger index from EEG match the right/correct trial index from edat
if ~isequal(new_task_trig_idx,task_trig_idx),
  disp([' WARNING: event trigger index from EEG and edat do not match. aborting.']);
  EEG.event(:) = [];
  return;
end 

disp([' MESSAGE: finished fixing trigger issue ']); disp(blanks(1)');

% END
