function [erp] = erp_recode_ttypes(erp,REC),

%function [erp] = erp_recode_ttypes(erp,REC),
%
%  Recodes ttypes using a static set of numeric values (REC.newtypes)
%  that are linked to erp.ttype through erp.stimkeys.ttype.
%  This function is mainly used to ensure that ttypes are comparable
%  across subjects (e.g., ttypes returned by eeglab2ptb may vary on a 
%  subject-by-subject basis). Must be used on a subject-level file.
%
%  erp - subject-level erp structure
%  REC - structured variable containing
%        .newtype  - contains the desired numeric ttype
%        .template - contains the corresponding character/value from
%                    erp.stimkeys.ttype
%
%  Example:
%  REC(1).newttype  = [ 1]; REC(1).template  = {'C'};
%  REC(2).newttype  = [ 2]; REC(2).template  = {'E'};
%  erp = erp_recode_ttypes(erp,REC);
%
%  Jeremy Harper, UMN, 121514
%


for aa = 1:length(REC),
  match_idx = find(ismember(erp.stimkeys.ttype, REC(aa).template));
  if isempty(match_idx),
    REC(aa).oldttype = -999;
  else,
    REC(aa).oldttype = find(ismember(erp.stimkeys.ttype, REC(aa).template));
  end
end

newttypevec = zeros(size(erp.ttype));

for bb = 1:length(erp.ttype),
  cur_idx = find([REC.oldttype]==erp.ttype(bb));
  newttypevec(bb) = [REC(cur_idx).newttype];
end

erp.ttype = newttypevec;
erp.stimkeys.ttype = cell([REC.template]);
