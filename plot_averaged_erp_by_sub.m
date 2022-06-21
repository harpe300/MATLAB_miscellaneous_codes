
% small script to step through each file and visually inspect the data at multiple points to check for gross artifacts

load sbjavgs_e11f23x_gng_stim_s25lLCD_erp_55Hz_256; erp55 = erp; clear erp
load sbjavgs_e11f23x_gng_stim_s25lLCD_erp_256;

subnums = unique(erp.subnum);

for aa  = 1:length(unique(erp.subnum)),

  cur_sub = subnums(aa);

  disp(['Cur subname = ' erp.subs.name(cur_sub,:)]);

  subplot(3,2,1); 
  plot(erp.data(erp.subnum==cur_sub&erp.elec==22&erp.stim.catcodes=='N',erp.tbin-128:erp.tbin+256),'r'); hold on;
  plot(erp.data(erp.subnum==cur_sub&erp.elec==22&erp.stim.catcodes=='G',erp.tbin-128:erp.tbin+256),'b'); hold off;
  title('FCZ');
  subplot(3,2,2); 
  plot(erp55.data(erp55.subnum==cur_sub&erp55.elec==22&erp55.stim.catcodes=='N',erp55.tbin-128:erp55.tbin+256),'r'); hold on;
  plot(erp55.data(erp55.subnum==cur_sub&erp55.elec==22&erp55.stim.catcodes=='G',erp55.tbin-128:erp55.tbin+256),'b'); hold off;
  title('FCZ-55hz');

  subplot(3,2,3);
  plot(erp.data(erp.subnum==cur_sub&erp.elec==31&erp.stim.catcodes=='N',erp.tbin-128:erp.tbin+256),'r'); hold on;
  plot(erp.data(erp.subnum==cur_sub&erp.elec==31&erp.stim.catcodes=='G',erp.tbin-128:erp.tbin+256),'b'); hold off;
  title('CZ');
  subplot(3,2,4);
  plot(erp55.data(erp55.subnum==cur_sub&erp55.elec==31&erp55.stim.catcodes=='N',erp55.tbin-128:erp55.tbin+256),'r'); hold on;
  plot(erp55.data(erp55.subnum==cur_sub&erp55.elec==31&erp55.stim.catcodes=='G',erp55.tbin-128:erp55.tbin+256),'b'); hold off;
  title('CZ-55hz');

  subplot(3,2,5);
  plot(erp.data(erp.subnum==cur_sub&erp.elec==49&erp.stim.catcodes=='N',erp.tbin-128:erp.tbin+256),'r'); hold on;
  plot(erp.data(erp.subnum==cur_sub&erp.elec==49&erp.stim.catcodes=='G',erp.tbin-128:erp.tbin+256),''); hold off;
  title('PZ');
  subplot(3,2,6);
  plot(erp55.data(erp55.subnum==cur_sub&erp55.elec==49&erp55.stim.catcodes=='N',erp55.tbin-128:erp55.tbin+256),'r'); hold on;
  plot(erp55.data(erp55.subnum==cur_sub&erp55.elec==49&erp55.stim.catcodes=='G',erp55.tbin-128:erp55.tbin+256),''); hold off;
  title('PZ-55hz');

  keyboard;

  clear cur_sub;

end

