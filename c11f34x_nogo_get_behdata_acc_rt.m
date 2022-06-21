%
%
% Reads, processes, and scores E-Prime data for each requested file
% Creates and exports a datafile of scored behavioral data 
%
% Input:
%   subname   - input subname
%   infilebeg - file prefix or path
%   infileend - file suffix
%
% Output:
%   respdata - vector of response hand used (1,2)
%
%
% Jeremy Harper, UMN, 10/25/21

  diary c11f34x_nogo_behdata_acc_rt.log;

  c11f34x_nogo_subnames; subnames = char(subnames);

  for aa = 1:length(subnames),

    respdata(aa).subname = subnames(aa,:); 
    subname = subnames(aa,:); disp(subname); disp(blanks(1)');

    trials = dataset('File',['/labs//mctfr-psyphys/data/public/labuser/eeg.toolbox-v0.5/scripts/c11f34x/nogo/eprime_txt_files/' subname '_nogo.txt.tab'],'Delimiter','\t'); % read in edat data
    trials.Stim_RESP(find(ismember(trials.Stim_RESP,'""')))={'0'}; % recode no resp ("") to 0
    trials.RESP = str2num(char(trials.Stim_RESP)); % create numeric response vector for indexing below

    if isempty(strmatch('Nogo_t',get(trials,'VarNames'))),

      % create accuracy vector 0/1 (the one from edat isn't correct )
      trials.ACC = nan(length(trials),1);
      trials.ACC((trials.Flanker_t==1|trials.Flanker_t==2)&(trials.RESP==1|trials.RESP==11|trials.RESP==2|trials.RESP==22|trials.RESP==26|trials.RESP==12|trials.RESP==21)) = 1; % correct go
      trials.ACC((trials.Flanker_t==3|trials.Flanker_t==4)&(trials.RESP==0)) = 1; % correct nogo
      trials.ACC((trials.Flanker_t==1|trials.Flanker_t==2)&(trials.RESP==0)) = 0; % incorrect go
      trials.ACC((trials.Flanker_t==3|trials.Flanker_t==4)&(trials.RESP==1|trials.RESP==11|trials.RESP==2|trials.RESP==22|trials.RESP==26|trials.RESP==12|trials.RESP==21)) = 0; % incorrect nogo

      % create response/accuracy trigger vector (accurate/corrected trigger mapping)
      trials.ACC_t = nan(length(trials),1);
      trials.ACC_t(trials.Flanker_t==1&trials.ACC==1) = 11; % correct go1
      trials.ACC_t(trials.Flanker_t==1&trials.ACC==0) = 21; % incorrect go1
      trials.ACC_t(trials.Flanker_t==2&trials.ACC==1) = 12; % correct go2
      trials.ACC_t(trials.Flanker_t==2&trials.ACC==0) = 22; % incorrect go2

      trials.ACC_t(trials.Flanker_t==3&trials.ACC==1) = 13; % correct nogo1
      trials.ACC_t(trials.Flanker_t==3&trials.ACC==0) = 23; % incorrect nogo1
      trials.ACC_t(trials.Flanker_t==4&trials.ACC==1) = 14; % correct nogo2
      trials.ACC_t(trials.Flanker_t==4&trials.ACC==0) = 24; % incorrect nogo2

    elseif ~isempty(strmatch('Nogo_t',get(trials,'VarNames'))),
      % create accuracy vector 0/1 (the one from edat isn't correct )
      trials.ACC = nan(length(trials),1);
      trials.ACC((trials.Nogo_t==1|trials.Nogo_t==2)&(trials.RESP==1|trials.RESP==11|trials.RESP==2|trials.RESP==22|trials.RESP==26|trials.RESP==12|trials.RESP==21)) = 1; % correct go
      trials.ACC((trials.Nogo_t==3|trials.Nogo_t==4)&(trials.RESP==0)) = 1; % correct nogo
      trials.ACC((trials.Nogo_t==1|trials.Nogo_t==2)&(trials.RESP==0)) = 0; % incorrect go
      trials.ACC((trials.Nogo_t==3|trials.Nogo_t==4)&(trials.RESP==1|trials.RESP==11|trials.RESP==2|trials.RESP==22|trials.RESP==26|trials.RESP==12|trials.RESP==21)) = 0; % incorrect nogo

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

    end

    acc = [trials.ACC_t];

    block = [trials.Block];


    respdata(aa).block_1_acc_total_per = -999;
    respdata(aa).block_1_acc_nogo_per = -999;
    respdata(aa).block_1_acc_go_per = -999;

    respdata(aa).block_2_acc_total_per = -999;
    respdata(aa).block_2_acc_nogo_per = -999;
    respdata(aa).block_2_acc_go_per = -999;

    respdata(aa).block_3_acc_total_per = -999;
    respdata(aa).block_3_acc_nogo_per = -999;
    respdata(aa).block_3_acc_go_per = -999;

    if length(block(block==1)) > 0,
      respdata(aa).block_1_acc_total_per = 100*(length(acc((acc==11|acc==12|acc==13|acc==14)&block==1))/length(acc((acc==11|acc==12|acc==13|acc==14|acc==21|acc==22|acc==23|acc==24)&block==1)));
      respdata(aa).block_1_acc_nogo_per = 100*(length(acc((acc==13|acc==14)&block==1))/length(acc((acc==13|acc==14|acc==23|acc==24)&block==1)));
      respdata(aa).block_1_acc_go_per = 100*(length(acc((acc==11|acc==12)&block==1))/length(acc((acc==11|acc==12|acc==21|acc==22)&block==1)));
    end
    if length(block(block==2)) > 0,
      respdata(aa).block_2_acc_total_per  = 100*(length(acc((acc==11|acc==12|acc==13|acc==14)&block==2))/length(acc((acc==11|acc==12|acc==13|acc==14|acc==21|acc==22|acc==23|acc==24)&block==2)));
      respdata(aa).block_2_acc_nogo_per  = 100*(length(acc((acc==13|acc==14)&block==2))/length(acc((acc==13|acc==14|acc==23|acc==24)&block==2)));
      respdata(aa).block_2_acc_go_per  = 100*(length(acc((acc==11|acc==12)&block==2))/length(acc((acc==11|acc==12|acc==21|acc==22)&block==2)));
    end
    if length(block(block==3)) > 0,
      respdata(aa).block_3_acc_total_per = 100*(length(acc((acc==11|acc==12|acc==13|acc==14)&block==3))/length(acc((acc==11|acc==12|acc==13|acc==14|acc==21|acc==22|acc==23|acc==24)&block==3)));
      respdata(aa).block_3_acc_nogo_per = 100*(length(acc((acc==13|acc==14)&block==3))/length(acc((acc==13|acc==14|acc==23|acc==24)&block==3)));
      respdata(aa).block_3_acc_go_per = 100*(length(acc((acc==11|acc==12)&block==3))/length(acc((acc==11|acc==12|acc==21|acc==22)&block==3)));
    end

    respdata(aa).Go_c_per   = -999;
    respdata(aa).Nogo_c_per = -999;
    respdata(aa).Go_e_per   = -999;
    respdata(aa).Nogo_e_per = -999;
    respdata(aa).Go_c   = -999;
    respdata(aa).Nogo_c = -999;
    respdata(aa).Go_e   = -999;
    respdata(aa).Nogo_e = -999;

    respdata(aa).Go_c_per   = 100*(length(acc(acc==11|acc==12))/length(acc(acc==11|acc==12|acc==21|acc==22)));
    respdata(aa).Nogo_c_per = 100*(length(acc(acc==13|acc==14))/length(acc(acc==13|acc==14|acc==23|acc==24)));
    respdata(aa).Go_e_per   = 100*(length(acc(acc==21|acc==22))/length(acc(acc==11|acc==12|acc==21|acc==22)));
    respdata(aa).Nogo_e_per = 100*(length(acc(acc==23|acc==24))/length(acc(acc==13|acc==14|acc==23|acc==24)));
    respdata(aa).Go_c   = length(acc(acc==11|acc==12));
    respdata(aa).Nogo_c = length(acc(acc==13|acc==14));
    respdata(aa).Go_e   = length(acc(acc==21|acc==22));
    respdata(aa).Nogo_e = length(acc(acc==23|acc==24));

    respcode = trials.RESP; % create numeric response vector for indexing below
    respcode(respcode==0) = [];
    if length(find(unique(respcode)==1|unique(respcode)==2|unique(respcode)==11|unique(respcode)==12|unique(respcode)==16|unique(respcode)==21|unique(respcode)==22|unique(respcode)==26))>1,
      if length(unique(respcode)) == 2 && (isequal(unique(respcode),[1;11]) | isequal(unique(respcode),[2;22])),
        if unique(respcode)== [1;11],
          respcode = 1;
        elseif unique(respcode)== [2;22],
          respcode = 2;
        end
      else,
        disp(['WARNING: More than one response hand used. Evaluate and decide.']); 
        disp(unique(respcode)');
        keyboard;
      end
    end

    respdata(aa).respcode = unique(respcode);

  % d-prime with loglinear approach to correct for FA or HIT of 0 or 1

    respdata(aa).d_prime = -999;

    hitrate    = (0.5+respdata(aa).Go_c)/(1+(respdata(aa).Go_c+respdata(aa).Go_e));
    falsealarm = (0.5+respdata(aa).Nogo_e)/(1+(respdata(aa).Nogo_c+respdata(aa).Nogo_e));
    respdata(aa).d_prime = norminv(hitrate)-norminv(falsealarm);

    rt = [trials.Stim_RT];

    respdata(aa).Go_c_rt       = -999;
    respdata(aa).Go_c_rt_rtv   = -999;
    respdata(aa).Go_c_rt_icv   = -999;

    respdata(aa).Nogo_e_rt     = -999;
    respdata(aa).Nogo_e_rt_rtv = -999;
    respdata(aa).Nogo_e_rt_icv = -999;

    respdata(aa).Go_c_rt       = mean(rt(acc==11|acc==12));
    respdata(aa).Go_c_rt_rtv   = std(rt((acc==11|acc==12)));
    respdata(aa).Go_c_rt_icv   = (std(rt((acc==11|acc==12))))/(mean(rt((acc==11|acc==12))));

    if respdata(aa).Nogo_e_per > 0,
      respdata(aa).Nogo_e_rt     = mean(rt(acc==23|acc==24));
      respdata(aa).Nogo_e_rt_rtv = std(rt((acc==23|acc==24)));
      respdata(aa).Nogo_e_rt_icv = (std(rt((acc==23|acc==24))))/(mean(rt((acc==23|acc==24))));
    end

    disp(respdata(aa)); disp(blanks(1));

    clear subname trials rt acc block falsealarm hitrate;

  end

subname = str2num([(char({respdata.subname}))]);
respcode = [respdata.respcode]';
block_1_acc_total_per = [respdata.block_1_acc_total_per]';
block_1_acc_nogo_per = [respdata.block_1_acc_nogo_per]';
block_1_acc_go_per = [respdata.block_1_acc_go_per]';
block_2_acc_total_per = [respdata.block_2_acc_total_per]';
block_2_acc_nogo_per = [respdata.block_2_acc_nogo_per]';
block_2_acc_go_per = [respdata.block_2_acc_go_per]';
block_3_acc_total_per = [respdata.block_3_acc_total_per]';
block_3_acc_nogo_per = [respdata.block_3_acc_nogo_per]';
block_3_acc_go_per = [respdata.block_3_acc_go_per]';
Go_c_per = [respdata.Go_c_per]';
Nogo_c_per = [respdata.Nogo_c_per]';
Go_e_per = [respdata.Go_e_per]';
Nogo_e_per = [respdata.Nogo_e_per]';
Go_c = [respdata.Go_c]';
Nogo_c = [respdata.Nogo_c]';
Go_e = [respdata.Go_e]';
Nogo_e = [respdata.Nogo_e]';
d_prime = [respdata.d_prime]';
Go_c_rt = [respdata.Go_c_rt]';
Nogo_e_rt = [respdata.Nogo_e_rt ]';
Go_c_rt = [respdata.Go_c_rt]';
Go_c_rt_rtv = [respdata.Go_c_rt_rtv]';
Go_c_rt_icv = [respdata.Go_c_rt_icv]';
Nogo_e_rt = [respdata.Nogo_e_rt]';
Nogo_e_rt_rtv = [respdata.Nogo_e_rt_rtv]';
Nogo_e_rt_icv = [respdata.Nogo_e_rt_icv]';

fid = fopen('c11f34x_nogo_behdata_acc_rt.txt','wt');
fprintf(fid, repmat('%s\t',1,26), 'subname', 'respcode', 'block_1_acc_total_per', 'block_1_acc_nogo_per','block_1_acc_go_per','block_2_acc_total_per','block_2_acc_nogo_per','block_2_acc_go_per','block_3_acc_total_per','block_3_acc_nogo_per','block_3_acc_go_per','Go_c_per','Nogo_c_per','Go_e_per','Nogo_e_per','Go_c','Nogo_c','Go_e','Nogo_e','d_prime','Go_c_rt','Go_c_rt_rtv','Go_c_rt_icv','Nogo_e_rt','Nogo_e_rt_rtv','Nogo_e_rt_icv'); 
fprintf(fid,'\n');
fclose(fid);
% have a matrix M(N,3) ready to go
dlmwrite('c11f34x_nogo_behdata_acc_rt.txt',[subname respcode block_1_acc_total_per block_1_acc_nogo_per block_1_acc_go_per block_2_acc_total_per block_2_acc_nogo_per block_2_acc_go_per block_3_acc_total_per block_3_acc_nogo_per block_3_acc_go_per Go_c_per Nogo_c_per Go_e_per Nogo_e_per Go_c Nogo_c Go_e Nogo_e d_prime Go_c_rt Go_c_rt_rtv Go_c_rt_icv Nogo_e_rt Nogo_e_rt_rtv Nogo_e_rt_icv],'delimiter','\t','precision',8, '-append')

diary off;

% END
