


  % trigger coding
  % 4 = stim-locked congruent   (S 1)
  % 6 = stim-locked incongruent (S 2)
  % 7 = resp-locked correct     (S 3)
  % 8 = resp-locked correct     (S 3)
  % 9 = resp-locked error       (S 9)
  % 5 = resp-locked error       (S10)

addpath /labs/mctfr-psyphys/data/public/labuser/eeg_commander
addpath /labs/mctfr-psyphys/data/public/jharper/components.care/eeglab14_1_1b
eeglab;
eeg_commander_startup;


subnames = {
'CARE001'
'CARE002'
'CARE003'
};


for jj = 1:length(subnames),


  % read in data
  EEG = pop_biosig(['./data/' char(subnames(jj)) '_flanker.vhdr'], 'blockepoch','off');
  disp(['Loading : ' subnames(jj)]);


  % remove initial triggers sent at the beginning of the file (all have the same latency, which causes problems down the line)
  EEG = pop_selectevent(EEG, 'type', [4 5 6 7 8 9], 'deleteevents', 'on');


  % convert events to numeric
  for ii=1:length(EEG.event),
     EEG.event(ii).type = num2str(EEG.event(ii).type);
  end;
  for hh=1:length(EEG.urevent),
     EEG.urevent(hh).type = num2str(EEG.urevent(hh).type);
  end;


  % define some vars
  EEG.filename = [char(subnames(jj)) '_flanker.vhdr'];
  EEG.subject  = char(subnames(jj));


  %[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'gui','off'); 
  %eeglab redraw

  % uppercase elec labels
  for aa = 1:length({EEG.chanlocs.labels}),
    EEG.chanlocs(aa).labels = strtrim(upper(EEG.chanlocs(aa).labels));
    if isequal(EEG.chanlocs(aa).labels, 'M1'),
      EEG.chanlocs(aa).labels = 'Ref1';
    elseif isequal(EEG.chanlocs(aa).labels, 'M2'),
      EEG.chanlocs(aa).labels = 'Ref2';
    elseif isequal(EEG.chanlocs(aa).labels, 'VEO+'),
      EEG.chanlocs(aa).labels = 'Veog+';
    elseif isequal(EEG.chanlocs(aa).labels, 'VEO-'),
      EEG.chanlocs(aa).labels = 'Veog-';
    elseif isequal(EEG.chanlocs(aa).labels, 'HEOR'),
      EEG.chanlocs(aa).labels = 'Heogr';
    elseif isequal(EEG.chanlocs(aa).labels, 'HEOL'),
      EEG.chanlocs(aa).labels = 'Heogl';
    end
  end


  % event recode
  REC( 1).label  = {'14'}; REC( 1).template  = {'4'    '7'}; REC( 1).position  = 1; % stim: correct   (commission)
  REC( 2).label  = {'14'}; REC( 2).template  = {'4'    '8'}; REC( 2).position  = 1;
  REC( 3).label  = {'16'}; REC( 3).template  = {'6'    '7'}; REC( 3).position  = 1;
  REC( 4).label  = {'16'}; REC( 4).template  = {'6'    '8'}; REC( 4).position  = 1;

  REC( 5).label  = {'17'}; REC( 5).template  = {'14'   '7'}; REC( 5).position  = 2; % resp: correct   (commission)
  REC( 6).label  = {'18'}; REC( 6).template  = {'14'   '8'}; REC( 6).position  = 2;
  REC( 7).label  = {'17'}; REC( 7).template  = {'16'   '7'}; REC( 7).position  = 2;
  REC( 8).label  = {'18'}; REC( 8).template  = {'16'   '8'}; REC( 8).position  = 2;

  REC( 9).label  = {'24'}; REC( 9).template  = {'4'    '5'}; REC( 9).position  = 1; % stim: incorrect   (commission)
  REC(10).label  = {'24'}; REC(10).template  = {'4'    '9'}; REC(10).position  = 1;
  REC(11).label  = {'26'}; REC(11).template  = {'6'    '5'}; REC(11).position  = 1;
  REC(12).label  = {'26'}; REC(12).template  = {'6'    '9'}; REC(12).position  = 1;

  REC(13).label  = {'25'}; REC(13).template  = {'24'   '5'}; REC(13).position  = 2; % resp: incorrect   (commission)
  REC(14).label  = {'29'}; REC(14).template  = {'24'   '9'}; REC(14).position  = 2;
  REC(15).label  = {'25'}; REC(15).template  = {'26'   '5'}; REC(15).position  = 2;
  REC(16).label  = {'29'}; REC(16).template  = {'26'   '9'}; REC(16).position  = 2;

  for bb = 1:length(REC),
    cur_rec = REC(bb).template; cur_rec_pos = REC(bb).position-1; cur_rec_label = REC(bb).label;
    REC_idx = strfind(str2num(char({EEG.event.type}))',str2num(char(cur_rec))');
    for kk = 1:length(REC_idx+cur_rec_pos),
      EEG.event(REC_idx(kk)+cur_rec_pos).type = char(cur_rec_label);
    end
  end
  clear REC cur_rec REC_idx;


  %EEG = pop_select(EEG,'time',[2 100]);


  % import chanlocs
  EEG=pop_chanedit(EEG, 'lookup','/labs/mctfr-psyphys/data/public/jharper/components.care/Standard-10-20-Cap34_care.ced');


  % downsample to 256 Hz
  EEG = pop_resample (EEG, 256);


  % high-pass 0.1 Hz with KAISER
  EEG = eeg_hpfilter(EEG, .1, 1, 0);


  % cfg import
  process.name     = 'temp_CARE';
  process.logfile  =  0;
  process.opts.EEG =  EEG;
  process.loaddir  = 'EEG';
  process.abortif  = 'EEG.nbchan<2';
  run(['/labs/mctfr-psyphys/data/public/jharper/components.care/cfg_import_general.m']);
  EEG = proc_commander(process, config);
  clear process config;


  % run recode to get rid of -1 triggers (which demarcate 1-s epochs used for import processing)
  clear REC
  REC(1).label = {'14'}; REC(1).template  = {'14'}; REC(1).position  = 1;  % stim:   correct   (commission)
  REC(2).label = {'16'}; REC(2).template  = {'16'}; REC(2).position  = 1;
  REC(3).label = {'17'}; REC(3).template  = {'17'}; REC(3).position  = 1;  % resp:   correct   (commission)
  REC(4).label = {'18'}; REC(4).template  = {'18'}; REC(4).position  = 1;
  REC(5).label = {'24'}; REC(5).template  = {'24'}; REC(5).position  = 1;  % stim: incorrect   (commission)
  REC(6).label = {'26'}; REC(6).template  = {'26'}; REC(6).position  = 1;
  REC(7).label = {'25'}; REC(7).template  = {'25'}; REC(7).position  = 1;  % resp: incorrect   (commission)
  REC(8).label = {'29'}; REC(8).template  = {'29'}; REC(8).position  = 1; 
  EEG = eeg_recode(EEG, REC, 1);


  % excise artifact from continuous data
  disp(['Removing these timepoints from cont data: ']); disp(num2str(EEG.reject.rejmanualtimes/EEG.srate));
  EEG = eeg_eegrej(EEG, EEG.reject.rejmanualtimes);


  % select and remove ocular ICs
  evalc(['pop_topoplot(EEG,0, [1:length(EEG.icawinv)] ,''ICs'',[6 7] ,0,''electrodes'',''on'');']);

  disp(['What IC timecourses do you want to view? Input as (e.g.) ICview = [];']);
  keyboard;
  tmpdata = eeg_getdatact(EEG, 'component', [1:size(EEG.icaweights,1)]);
  eegplot(tmpdata(ICview,:), 'srate', EEG.srate); clear tmpdata ICview;

  disp(['Enter the ICs you want to reject. eg., ICrej = [1 3 7];']);
  disp(['NOTE: You must follow the format of ICrej = [];']);
  keyboard;
  disp(['NOTE: Make sure you specified which ICs to remove!']);
  keyboard;

  disp(['ICs selected: ' num2str(ICrej)]);
  EEG = pop_subcomp(EEG, ICrej, 0); clear ICrej;


  % cfg export
  process.name     = 'temp_CARE';
  process.logfile  =  0;
  process.opts.EEG =  EEG;
  process.loaddir  = 'EEG';
  process.abortif  = 'EEG.nbchan<2';

  % -- epoch (note that this baseline corrects for the entire epoch)
  config(1).epoch.events =  {   '14'  '16'  '17'  '18'  '24'  '25'  '26'  '29' }; % ('all','ALL',[]) or cellstr
  config(1).epoch.winsec =       [-2    2]; % in sec
  config(1).epoch.rmbase =              []; % in msec, or 'no'/'NO', or []
  config(1).epoch.rdcevt =               1; % reduces out all other events

  % -- remove blink channels
  config(2).userinput = 'EEG = pop_select(EEG,''nochannel'',{''Veog+'',''Veog-'',''Heogr'',''Heogl''});';

  % -- artifact elec-sweeps
  config(3).artifact.type         =     'matrix'; % see artifact_channels, thresholder
  config(3).artifact.datafilt     =           0 ; % use only non-tagged epochs/channels
  config(3).artifact.pthreshchan  =         .75 ; 
  config(3).artifact.minchans     =           1 ;
  config(3).artifact.rejchan      =           1 ;
  config(3).artifact.rejtrial     =           1 ;
  config(3).artifact.opts(1).meas =           {    'minvar'}; % salt-bridging, sample-wide
  config(3).artifact.opts(1).stat =                { 'abs'};
  config(3).artifact.opts(1).crit =                {  [2] };
  config(3).artifact.opts(1).joint=                      0 ;
  config(3).artifact.opts(2).meas =           {'Vdist-min'}; % salt-bridging, sample-wide
  config(3).artifact.opts(2).stat =                { 'abs'};
  config(3).artifact.opts(2).crit =                {  [2] };
  config(3).artifact.opts(2).joint=                      0 ;

  config(4).artifact.type         =     'matrix'; % see artifact_channels, thresholder
  config(4).artifact.datafilt     =           0 ; % use only non-tagged epochs/channels
  config(4).artifact.pthreshchan  =         .75 ;
  config(4).artifact.minchans     =           1 ;
  config(4).artifact.rejchan      =           1 ;
  config(4).artifact.rejtrial     =           1 ;
  config(4).artifact.opts(1).meas = {'Vdist-nearest'  'var'}; % large deviation 
  config(4).artifact.opts(1).stat = {          'abs' 'nMAD'};
  config(4).artifact.opts(1).crit = {          [100]    [4]};
  config(4).artifact.opts(1).joint=                      1 ;
  config(4).artifact.opts(2).meas = {      'fqvar' 'range'}; % muscle 
  config(4).artifact.opts(2).stat = {       'nMAD'   'abs'};
  config(4).artifact.opts(2).crit = {          [4]   [200]};
  config(4).artifact.opts(2).joint=                      1 ;

  % run!
  EEG = proc_commander(process, config);
  clear process config; 


  % mark-up trials to exclude, if no trials, skip and do not save (reject)
  minchans = 26;
  at_trials= double(sum(EEG.reject.rejmanualE==0)<minchans);
  if ~isempty(sum(at_trials))&&sum(at_trials)<EEG.trials,
    t_rejmanualE = EEG.reject.rejmanualE;
    t_rejmanual  = EEG.reject.rejmanual ;
    t_rejmanualE(:,find(at_trials)) = '';
    t_rejmanual(   find(at_trials)) = '';
    t_rejmanualE(find(ismember({EEG.chanlocs.labels},{'Veog+','Veog-','Heogr','Heogl'})),:) = '';
    EEG = pop_select(EEG, 'notrial', find(at_trials),'nochannel', {'Veog+','Veog-','Heogr','Heogl'}); %drops trials/EOGs
    EEG.reject.rejmanualE = t_rejmanualE; EEG.reject.rejmanual = t_rejmanual;
  elseif sum(at_trials)==EEG.trials,
    t_rejmanualE = EEG.reject.rejmanualE;
    t_rejmanual  = EEG.reject.rejmanual ;
    t_rejmanualE(find(ismember({EEG.chanlocs.labels},{'Veog+','Veog-','Heogr','Heogl'})),:) = '';
    EEG = pop_select(EEG, 'nochannel', {'Veog+','Veog-','Heogr','Heogl'}); %drops EOGs, if no trials are to be dropped
    disp([' ERROR: cur_sub has less than 26 clean elecs (<75%). Has ' num2str(size(EEG.chanlocs,2)) ' elecs. rejecting.']);
    return
  end
  clear minchans at_trials t_rejmanualE t_rejmanual;


  % interpolate to full-montage
  process.loaddir = 'EEG';
  process.opts.EEG = EEG ;
  config(1).interp.type      =                         'full'; % or 'artifact'
  config(1).interp.montage   = 'Standard-10-20-Cap34_care.ced'; % see readlocs
  config(1).interp.interpmaj =                            0  ; % if >50% bad, interp all
  config(1).interp.rejthresh =                           25  ; % lower-bound (pct)
  config(1).interp.rejreject =                            1  ; % reject non-interpolated
  EEG = proc_commander(process,config);
  clear process config;


  EEG = pop_saveset(EEG, 'filename', [EEG.subject '_flanker_cnt_epoch'], 'filepath', '/labs/mctfr-psyphys/data/public/jharper/components.care/data/exported/', 'savemode', 'onefile');


  clear EEG process config;
  close all;

  disp(blanks(2)');

end

