%% WADCP PACKAGE 0.1 -- Matlab routines for ADCP data post-processing
%%
% Routine used to calibrate ADCP backscatter with moored turbidimeter data
% 
% Input turbidity data : txt file with columns as follow: [Year] [Month] [Day] [Hour] [Minute] [Turbidity g/l]  

function cal_turbidity_plot

%global icolor filein_path path_routines
global IV xrange yrange 
%global t_calib0,turb_int0,pval0,t_calib,turb_10log,turb_alt,ivcal,nens_calib,coef,er2

display('CALIBRATION RESULTS 10log10(M)=aIV+b')
display(['a= ' num2str(ivcal(1))])
display(['b= ' num2str(ivcal(2))])
display(['number of data used N= ' num2str(nens_calib)])
display(['correlation coefficient R= ' num2str(coef(1,2))])
display(['root mean square error RMS= ' num2str(er2) ' %'])

% FIGURES TO CONTROL RESULTS OF CALIBRATION
%--------------------------------------------------------------------------

% figure
% hold on
% plot(IV(t_calib,turb_alt),pval,'r')
% plot(IV(t_calib,turb_alt),turb_10log,'.','Color',icolor(:,2))
% 
% figure; plot(xrange(t_calib),IV(turb_alt,t_calib));
% figure;plot(xrange(t_calib),10.^(pval/10));

%%%time series of ADCP and TBD en 10LOG used for the calibration
%%% & time series of concentrations over all the data period
figure; 
%plot(xrange(t_calib),turb_10log,'r',xrange(t_calib),pval,'b');ylim([0 max(turb_10log)]);
plot(xrange(t_calib0),turb_int0,'r',xrange(t_calib0),10.^(pval0/10),'b'); ylim([0 max(turb_int0)]);
legend('TBD','ADCP');
set(gca,'Xlim',[fix(xrange(t_calib(1))) fix(xrange(t_calib(end))+1)],'Xtick',[fix(xrange(t_calib(1))):1:fix(xrange(t_calib(end))+1)]);
datetick('x',19,'keeplimits','keepticks');
xlabel(num2str(max(datevec(xrange(1)))));

%%%% scatter plot IV/10log10(M):
gam_iv=[fix(min(IV(turb_alt,t_calib))) fix(max(IV(turb_alt,t_calib)))];
figure;
plot(IV(turb_alt,t_calib0),10*log10(turb_int0),'.','Color',[.5 .5 .5],'MarkerSize',6);
ylabel('10 * log_{10}( TBD [g.l^{-1}] )','FontSize',7);xlabel('IV [dB_{/1m^3}]','FontSize',7);
hold on; plot(IV(turb_alt,t_calib),turb_10log,'.r','MarkerSize',6); hold off; 
%axis([-60 -30 4 22]);set(gca,'XTick',[-60:5:-20],'FontSize',7);
hold on; plot([gam_iv(1):1:gam_iv(2)],ivcal(1).*[gam_iv(1):1:gam_iv(2)]+ivcal(2),'k'); hold off;
text(gam_iv(1)+2,max(turb_10log)-2,strcat(['y=' num2str(ivcal(1)) ' x+' num2str(ivcal(2)) '   R=' num2str(coef(1,2))]),'Color','k','FontSize',7);

%%
