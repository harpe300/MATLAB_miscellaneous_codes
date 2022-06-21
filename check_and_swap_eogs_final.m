% Script to empirically check for evidence that the EOG channels were incorrectly placed during recording

diary check_and_swap_eogs_final.log

for jj = 1:length(subnames),

  [junk,EEG] = evalc('pop_loadset([deblank(subnames(jj,:))])');
  disp(deblank(subnames(jj,:)));

  eogics(jj).subname = EEG.subject;

if ~isempty(strmatch('Veog-',{EEG.chanlocs.labels})) | ~isempty(strmatch('Veog+',{EEG.chanlocs.labels})),

% make virtual blink/horizontal eye movement leads, lowpass filter        
pts    = floor(EEG.pnts/6):floor(EEG.pnts-(EEG.pnts/6));
if length(find(ismember({EEG.chanlocs.labels},{'FP1' 'FPZ' 'FP2' 'AF3' 'AFZ' 'AF4'}))) == 1,
  blink  = (EEG.data(find(ismember({EEG.chanlocs.labels},{'FP1' 'FPZ' 'FP2' 'AF3' 'AFZ' 'AF4'})),pts));
elseif isempty(find(ismember({EEG.chanlocs.labels},{'FP1' 'FPZ' 'FP2' 'AF3' 'AFZ' 'AF4'}))),
  disp(['ATTN: No frontal elecs available, skipping...']);
  eogics(jj).blink = -999
  continue
else,
  blink  = mean(EEG.data(find(ismember({EEG.chanlocs.labels},{'FP1' 'FPZ' 'FP2' 'AF3' 'AFZ' 'AF4'})),pts));
end

%hem    = EEG.data(strmatch('F8',{EEG.chanlocs.labels}),pts) - ...
%         EEG.data(strmatch('F7',{EEG.chanlocs.labels}),pts);
[results,blinkfilt] = evalc(['eegfilt(blink,EEG.srate,0,6)']);
%hem_filt  = eegfilt(hem ,EEG.srate,0,6);
blinkpts  =  abs(zscore(blinkfilt))>3;  % get indices of blinks with Z > 3
%hempts    =  abs(zscore(hem_filt))>3;
mad_blink = ( blinkfilt-median(blinkfilt) ) ./ (1.4826*mad(blinkfilt,1));  % calc Z-score but with MAD in place
%mad_hem   = (  hem_filt-median( hem_filt) ) ./ (1.4826*mad( hem_filt,1));
blinkpts  = (3<   ( mad_blink )); %&(   ( mad_blink )<5);
%hempts    = (3<abs( mad_hem   )); %&(abs( mad_hem   )<5);
%eyepts    =  unique([find(blinkpts) find(hempts)]);
eyepts    =  unique([find(blinkpts)]);
blink     =  blinkfilt(eyepts);
%hem       =  hem_filt(eyepts);


if ~isempty(strmatch('Veog+',{EEG.chanlocs.labels})) 
  Vplus   = EEG.data(strmatch('Veog+',{EEG.chanlocs.labels}),:);

  [results,Vpfilt]    = evalc('eegfilt(Vplus(pts),EEG.srate,0,6)');

  Vpfilt10 = prctile(Vpfilt',10);
  Vpfilt90 = prctile(Vpfilt',90);

  Vpfilt    = Vpfilt(eyepts);  % for corr with virtual elecs, keeps only data that corresponds with blinks in virtual eyes defined above

end


% grab EOG data, do the same as above
if ~isempty(strmatch('Veog-',{EEG.chanlocs.labels})),
  Vminus    = EEG.data(strmatch('Veog-',{EEG.chanlocs.labels}),:);

  [results,Vmfilt]    = evalc('eegfilt(Vminus(pts),EEG.srate,0,6)');

  Vmfilt10 = prctile(Vmfilt',10);
  Vmfilt90 = prctile(Vmfilt',90);

  Vmfilt    = Vmfilt(eyepts);  % for corr with virtual elecs, keeps only data that corresponds with blinks in virtual eyes defined above


end 

eogics(jj).blink = -999;

if ~isempty(strmatch('Veog-',{EEG.chanlocs.labels})) && ~isempty(strmatch('Veog+',{EEG.chanlocs.labels})),
  disp(['Corr with frontal blinks: Veog+ = ' num2str(corr(blink',Vpfilt')) ' | Veog- = ' num2str(corr(blink',Vmfilt'))]);
  disp(['10th/90th Perc: Veog+ = [' num2str(Vpfilt10) ', ' num2str(Vpfilt90)  '] | Veog- = [' num2str(Vmfilt10) ', ' num2str(Vmfilt90) ']']);

  if corr(blink',Vpfilt') < corr(blink',Vmfilt') && (Vpfilt90 < Vmfilt90),
    disp(['   switch_eogs (' EEG.filename '); Evidence that Veog- and Veog+ may be wrongly assigned in BDF']);
    eogics(jj).blink = 1;
  else,
    eogics(jj).blink = 0;
  end

else,
  disp(['  Cannot compare Veog+ and Veog- because ID is lacking both channels']);

end

else,
  disp(['  Cannot compare Veog+ and Veog- because ID is lacking both channels']);
  eogics(jj).blink = -999;

end


  clearvars  -except subnames jj eogics config

disp(blanks(4));


end

diary off
