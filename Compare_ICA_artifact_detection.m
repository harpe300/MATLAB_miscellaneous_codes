

% This script was written to visually and empirically compare performance of two automated methods for detecting and removing artifacts from high-dimensional data using ICA.

e11f23x_nogo_subnames;
subnames = char(subnames);

cd ./ICA_figures

diary ./ICA_plots_OLDvNEWBlinkHotrSep.log

eeglab

dirlist = dir('/home/harpe300/Desktop/mctfr-psyphys/data/main/biosemi/nogo/e11f23-x/*cnt.set');

subnames = deblank(char({dirlist.name}));

config.epoch.events =         {'G','N'}; % ('all','ALL',[]) or cellstr
config.epoch.winsec =       [-2    2]; % in sec
config.epoch.rmbase =              []; % in msec, or 'no'/'NO', or []
config.epoch.rdcevt =               1; % reduces out all other events

eogics.subname = zeros(length(subnames),8);
eogics.ics_old = zeros(length(subnames),1);
eogics.ics_new = zeros(length(subnames),1);
eogics.ics_old_flag = zeros(length(subnames),1);
eogics.ics_new_flag = zeros(length(subnames),1);

eogics = struct('subname', cell(1, length(subnames)), 'ics_old', cell(1, length(subnames)), 'ics_new', cell(1, length(subnames)), 'ics_old_flag', cell(1, length(subnames)), 'ics_new_flag', cell(1, length(subnames)));

for jj = 1:length(subnames),

  EEG = pop_loadset([deblank(subnames(jj,:))]);
  disp(deblank(subnames(jj,:)));

  eogics(jj).subname = EEG.subject;

  recode_nogo;
  EEG = eeg_eegrej(EEG, EEG.reject.rejmanualtimes);

% this chunk is adapted from eeg_commander (https://github.com/sjburwell/eeg_commander)
% ------------------------- Vertical Ocular Movement ----------------------
stat = 'EM';
crit = 1;

%threshold
THRESHOLD(1).stat = stat; THRESHOLD(1).crit = crit;
THRESHOLD(2).stat = stat; THRESHOLD(2).crit = crit;

%spatial component
SPATIAL    = 'blink_ICWINV67.mat';

%temporal component
CRIT_ELECS = [];
if ~isempty(strmatch('Veog+',{EEG.chanlocs.labels}))&&~isempty(strmatch('Veog-',{EEG.chanlocs.labels}))
    CRIT_ELECS = {'Veog+'; 'Veog-'};
elseif ~isempty(strmatch('Veog+',{EEG.chanlocs.labels}))&&isempty(strmatch('Veog-',{EEG.chanlocs.labels}))
    CRIT_ELECS = {'Veog+'};
elseif isempty(strmatch('Veog+',{EEG.chanlocs.labels}))&&~isempty(strmatch('Veog-',{EEG.chanlocs.labels}))&&~isempty(strmatch('FP2',{EEG.chanlocs.labels})),
    CRIT_ELECS = { 'FP2';  'Veog-'};
elseif isempty(strmatch('Veog+',{EEG.chanlocs.labels}))&&~isempty(strmatch('Veog-',{EEG.chanlocs.labels}))&&~isempty(strmatch('FP1',{EEG.chanlocs.labels})),
    CRIT_ELECS = { 'FP1';  'Veog-'};
elseif isempty(strmatch('Veog+',{EEG.chanlocs.labels}))&&~isempty(strmatch('Veog-',{EEG.chanlocs.labels}))&&~isempty(strmatch('FPZ',{EEG.chanlocs.labels})),
    CRIT_ELECS = { 'FPZ';  'Veog-'};
elseif isempty(strmatch('Veog+',{EEG.chanlocs.labels}))&& isempty(strmatch('Veog-',{EEG.chanlocs.labels}))&&~isempty(strmatch('FP2',{EEG.chanlocs.labels})),
    CRIT_ELECS = { 'FP2' };
elseif isempty(strmatch('Veog+',{EEG.chanlocs.labels}))&& isempty(strmatch('Veog-',{EEG.chanlocs.labels}))&&~isempty(strmatch('FP1',{EEG.chanlocs.labels})),
    CRIT_ELECS = { 'FP1' };
elseif isempty(strmatch('Veog+',{EEG.chanlocs.labels}))&& isempty(strmatch('FP1',{EEG.chanlocs.labels}))  && isempty(strmatch('FP2',{EEG.chanlocs.labels})),
    if ~isempty(strmatch('FPZ',{EEG.chanlocs.labels})),
        CRIT_ELECS = { 'FPZ' }; 
    elseif ~isempty(strmatch('AF4',{EEG.chanlocs.labels})),
        CRIT_ELECS = { 'AF4' };
    elseif ~isempty(strmatch('AF3',{EEG.chanlocs.labels})),
        CRIT_ELECS = { 'AF3' };
    elseif ~isempty(strmatch('AFZ',{EEG.chanlocs.labels})),
        CRIT_ELECS = { 'AFZ' };
    else,
        disp(['   eeg_icartifact_eyes; No criterion channel found for ' EEG.filename ', aborting blink correction...']);
        ICs = [];
        return
    end
end
[ICs_blink, blink_corrs] = ic_temporalspatial(EEG, SPATIAL, CRIT_ELECS, THRESHOLD, 1);

% ------------------------- Horizontal Ocular Movement ----------------------
%threshold
THRESHOLD(1).stat = stat; THRESHOLD(1).crit = crit;
THRESHOLD(2).stat = stat; THRESHOLD(2).crit = crit;

%spatial component
SPATIAL    = 'hem_ICWINV67.mat';

%temporal component
CRIT_ELECS = [];
if     ~isempty(strmatch('Heogr',{EEG.chanlocs.labels})) && ~isempty(strmatch('Heogl',{EEG.chanlocs.labels}))
    CRIT_ELECS = {'Heogr'; 'Heogl'};
elseif ~isempty(strmatch('Heogr',{EEG.chanlocs.labels})) &&  isempty(strmatch('Heogl',{EEG.chanlocs.labels}))
    CRIT_ELECS = {'Heogr'};
elseif ~isempty(strmatch('Heogl',{EEG.chanlocs.labels})) &&  isempty(strmatch('Heogr',{EEG.chanlocs.labels}))
    CRIT_ELECS = {'Heogl'};
elseif ~isempty(strmatch(   'F8',{EEG.chanlocs.labels})) && ~isempty(strmatch(   'F7',{EEG.chanlocs.labels}))
    CRIT_ELECS = {   'F8';    'F7'};
else
    disp(['   eeg_icartifact_eyes; No criterion channel found for ' EEG.filename ', aborting horizontal eye movement correction...']);
    ICs = [];
    return
end
[ICs_hem, hem_corrs] = ic_temporalspatial(EEG, SPATIAL, CRIT_ELECS, THRESHOLD, 1);


% ------------------------- Blink and Horizontal Separately, Joint Temporal+Spatial. ----------------------
all_ics = [1:size(EEG.icawinv,2)];
if ~isempty(unique([ICs_blink, ICs_hem]));
   
    disp([EEG.filename '; threshold(s): [' num2str([THRESHOLD.crit]) ']'])

    ICs_blink_hem = unique([ICs_blink, ICs_hem]);
    disp(['Spatial-temporal ICs matched to blinks   : ' num2str(ICs_blink)]);
    disp(['Spatial-temporal ICs matched to horizonal: ' num2str(ICs_hem)]);
    disp(['Spatial-temporal ICs matched: ' num2str(ICs_blink_hem)]);
end

eogics(jj).ics_new = ICs_blink_hem;


% ------------------------- Joint. Joint Blink+Horizontal, Temporal and Spatial separately ----------------------
eyespace  = abs(atanh([blink_corrs(:,1); hem_corrs(:,1)]));
[spacemeas,spacecrit] = thresholder(stat, eyespace, crit);

eyetime   = abs(atanh([blink_corrs(:,2); hem_corrs(:,2)])); 
[timemeas ,timecrit ] = thresholder(stat, eyetime , crit);

all_ics = [1:size(EEG.icawinv,2), 1:size(EEG.icawinv,2)];
if ~isempty(intersect(find(spacemeas>spacecrit),find(timemeas>timecrit)));
    disp([EEG.filename '; threshold(s): [' num2str([THRESHOLD.crit]) ']'])
    ICs = all_ics(intersect(find(spacemeas>spacecrit),find(timemeas>timecrit)));
    ICs = unique(ICs);
    disp(['Spatial-temporal ICs matched: ' num2str(ICs)]);
else
    [spacemeas,spacecrit] = thresholder(stat, eyespace, crit/2);
    [timemeas ,timecrit ] = thresholder(stat, eyetime , crit*2);
    if ~isempty(intersect(find(spacemeas>spacecrit),find(timemeas>timecrit)));
        disp([EEG.filename '; threshold(s): [' num2str([crit/2 crit*2]) ']'])
        ICs = all_ics(intersect(find(spacemeas>spacecrit),find(timemeas>timecrit)));
        ICs = unique(ICs);
        disp(['Spatial-temporal ICs matched: ' num2str(ICs)]);
    else,
        [spacemeas,spacecrit] = thresholder(stat, eyespace, crit*2);
        [timemeas ,timecrit ] = thresholder(stat, eyetime , crit/2);
        if ~isempty(intersect(find(spacemeas>spacecrit),find(timemeas>timecrit)));
            disp([EEG.filename '; threshold(s): [' num2str([crit*2 crit/2]) ']'])
            ICs = all_ics(intersect(find(spacemeas>spacecrit),find(timemeas>timecrit)));
            ICs = unique(ICs);
            disp(['Spatial-temporal ICs matched: ' num2str(ICs)]);
        else,
            disp([EEG.filename '; no temporal-spatial ICs identified, skipping IC-removal.']);
            ICs = [];
        end
    end
end
% end eeg_commander chunk

  eogics(jj).ics_old = ICs;


  blink_crit        = 60; % detection threshold for target in time series
  proportion_trials =.25; % threshold for proportion of targets across time series

  EEGsubcompOLD = pop_subcomp(EEG, ICs);
  EEGsubcompNEW = pop_subcomp(EEG, ICs_blink_hem);

  EEGsubcompOLD = proc_epoch(EEGsubcompOLD, config.epoch);
  EEGsubcompNEW = proc_epoch(EEGsubcompNEW, config.epoch);

  % Test if data exceeds proportion threshold using old method
  numblinkL  = max(EEGsubcompOLD.data(strmatch('FP1',{EEGsubcompOLD.chanlocs.labels}),:,:))>blink_crit;
  numblinkR  = max(EEGsubcompOLD.data(strmatch('FP2',{EEGsubcompOLD.chanlocs.labels}),:,:))>blink_crit;

  trialblink = length(find(numblinkL&numblinkR));
  if trialblink>(proportion_trials*EEGsubcompOLD.trials), 
    disp(['   check_ocular_correction; (' EEGsubcompOLD.filename ') evidence that blinks exist in OLD IC data.']);
    eogics(jj).ics_old_flag = 1;
  else, 
    eogics(jj).ics_old_flag = 0;
  end


  % Test if data exceeds proportion threshold using new method
  numblinkL  = max(EEGsubcompNEW.data(strmatch('FP1',{EEGsubcompNEW.chanlocs.labels}),:,:))>blink_crit;
  numblinkR  = max(EEGsubcompNEW.data(strmatch('FP2',{EEGsubcompNEW.chanlocs.labels}),:,:))>blink_crit;

  trialblink = length(find(numblinkL&numblinkR));
  if trialblink>(proportion_trials*EEGsubcompNEW.trials), 
    disp(['   check_ocular_correction; (' EEGsubcompNEW.filename ') evidence that blinks exist in NEW IC data.']);
    eogics(jj).ics_new_flag = 1;
  else, 
    eogics(jj).ics_new_flag = 0;
  end

  % Print ICA component numbers
  disp(blanks(1));
  eogics(jj)
  disp(blanks(1));


  % Print eps file for visual inspection and posterity
  epsname = ['./' subnames(jj,1:13) '_plots_OldvNewBlinkHorSep'];

  %pop_eegplot( EEG, 1, 1, 1); title('EEG');
  eegplot(EEG.data,'eloc_file',EEG.chanlocs,'winlength',20,'spacing',50); title('EEG');
  orient landscape; print(epsname,'-dpsc', '-fillpage'); close all;

  %pop_eegplot( EEG, 0, 1, 1); title('ICs');
  tmpdata = eeg_getdatact(EEG, 'component', [1:size(EEG.icaweights,1)]);
  eegplot(tmpdata,'winlength',20); title('IC activations');
  orient landscape; print(epsname,'-dpsc','-append', '-fillpage'); close all;

   %evalc('pop_topoplot(EEG,0, [1:size(EEG.icaweights,1)],''All ICs'',0,''electrodes'',''on'')');
   evalc('pop_topoplot(EEG,0, [1:16],''ICs 1:16'',0,''electrodes'',''on'')');
   orient portrait; print(epsname,'-dpsc','-append', '-fillpage'); close all;

  if length(ICs) == 1,
    evalc('pop_topoplot(EEG,0, [ICs ICs] ,''OLD Tagged ICs'',0,''electrodes'',''on'')');
    orient portrait; print(epsname,'-dpsc','-append'); close all;
  else,
    evalc('pop_topoplot(EEG,0, [ICs] ,''OLD Tagged ICs'',0,''electrodes'',''on'')');
    orient portrait; print(epsname,'-dpsc','-append'); close all;
  end

  eegplot(tmpdata(ICs,:),'winlength',20); title('OLD Tagged IC activations');
  orient landscape; print(epsname,'-dpsc','-append', '-fillpage'); close all;

  if length(ICs_blink_hem) == 1,
    evalc('pop_topoplot(EEG,0, [ICs_blink_hem ICs_blink_hem] ,''NEW Tagged ICs'',0,''electrodes'',''on'')');
    orient portrait; print(epsname,'-dpsc','-append'); close all;
  else,
    evalc('pop_topoplot(EEG,0, [ICs_blink_hem] ,''NEW Tagged ICs'',0,''electrodes'',''on'')');
    orient portrait; print(epsname,'-dpsc','-append'); close all;
  end

  eegplot(tmpdata(ICs_blink_hem,:),'winlength',20); title('NEW Tagged IC activations');
  orient landscape; print(epsname,'-dpsc','-append', '-fillpage'); close all;

  all_ICs = unique([ICs, ICs_blink_hem]);
  for aa = 1:length(all_ICs),
    evalc('pop_prop( EEG, 0, all_ICs(aa), NaN, {''freqrange'' [2 50] })');
    orient portrait; print(epsname,'-dpsc','-append'); close all;
  end

  %EEGsubcompOLD = pop_subcomp(EEG, ICs);
  %pop_eegplot( EEGsubcomp, 1, 1, 1); title('EEG minus OLD tagged ICs');
  eegplot(EEGsubcompOLD.data,'eloc_file',EEGsubcompOLD.chanlocs,'winlength',20,'spacing',50); title('EEG minus tagged ICs'); 
  orient landscape; print(epsname,'-dpsc','-append', '-fillpage');

  %EEGsubcompNEW = pop_subcomp(EEG, ICs_blink_hem);
  %pop_eegplot( EEGsubcomp, 1, 1, 1); title('EEG minus NEW tagged ICs');
  eegplot(EEGsubcompNEW.data,'eloc_file',EEGsubcompNEW.chanlocs,'winlength',20,'spacing',50); title('EEG minus tagged ICs'); 
  orient landscape; print(epsname,'-dpsc','-append', '-fillpage');


  clearvars  -except subnames jj eogics config

  close all

end

