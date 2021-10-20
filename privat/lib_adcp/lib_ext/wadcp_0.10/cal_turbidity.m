%% WADCP PACKAGE 0.3 -- Matlab routines for ADCP data post-processing
%%
% Routine used to calibrate ADCP backscatter with moored turbidimeter data
% 
% Input turbidity data : txt file with columns as follow: [Year] [Month] [Day] [Hour] [Minute] [Second] [Turbidity g/l]  

function [bical,pval,turb_alt_ind,t_calib0,turb_int0,nens_calib,t_calib,turb_int]=cal_turbidity


global icolor filein_path path_routines fileout_path mode_exp fileturbtemp turb_alt fidlog
global BI xrange yrange filein tminproc tmaxproc spm_threshold db_threshold

display('ADCP backscatter calibration in progress...')
cd(filein_path)
if mode_exp==0
[fileturbtemp,path1]=uigetfile({'*.txt';'*.*'},'Open the .txt file turbidity data....');
cd(path1);
end

file_turb=load(fileturbtemp);
cd(path_routines);

%%
if mode_exp==0
prompt = {'Altitude of the turbidimeter above bottom [m]','SPM thresold (g/l)','dB threshold (dB)'};
dlg_title = 'Turbidity calibration';
num_lines = 1;
def = {'1','0.000','-90'};
tmp = inputdlg(prompt,dlg_title,num_lines,def);
turb_alt=eval(char(tmp{1}));
spm_threshold=eval(char(tmp{2}));
db_threshold=eval(char(tmp{3}));
end

turb_alt_ind=max(find(turb_alt>=yrange));
%%%%turb_alt=turb_alt(1); % or -1 if the cell just below the turbidity meter
%!! turb_alt is now an indice of vector yrange, and no more the real level from the bottom of the turbidimeter !!

% le croisic first cell 
% turb_alt_ind=1;

if isempty(turb_alt_ind)==1 
    display('Problem with position of turbidimeter')
    display('no correspondance with ADCP data, first cell choosen')

    turb_alt_ind=1;
end

time_turb=datenum(file_turb(:,1),file_turb(:,2),file_turb(:,3),file_turb(:,4),file_turb(:,5),file_turb(:,6));
turb=file_turb(:,7);


% Take all the available data of TBD :
t_calib0=find(xrange>min(time_turb) & xrange<max(time_turb));

BItmp0=BI(turb_alt_ind,t_calib0);

t_calib0(isnan(BItmp0))=[];
BItmp0(isnan(BItmp0))=[];
turb_int0=interp1(time_turb,turb,xrange(t_calib0),'linear');

% Ask user-defined dates to consider for the calibration BI/tbd :

t_min=min(time_turb);
t_max=max(time_turb);
if mode_exp==0
prompt={'TBD tmin', 'TBD tmax'};
dlg_title='ADCP CALIBRATION WITH TBD DATA - time limits';
num_lines=1;
def={datestr(t_min),datestr(t_max)};
input1=inputdlg(prompt,dlg_title,num_lines,def);
tinf=datenum(input1(1));
tsup=datenum(input1(2));
else
tinf=tminproc;
tsup=tmaxproc;
end

% defined time calibration period :
%t_calib=xrange(xmin:xmax);
t_calib=find(xrange>tinf & xrange<tsup);

BItmp=BI(turb_alt_ind,t_calib);

t_calib(isnan(BItmp))=[];
%interpolate TBD data on ADCP time axis
turb_int=interp1(time_turb,turb,xrange(t_calib),'linear');

% turb_10log(isnan(BItmp))=[];
BItmp(isnan(BItmp))=[];

figure
plot(turb_int0)
% processing for turbidity values above an SPM concentration threshold
% (g/l)


BItmp(turb_int<spm_threshold)=[];
t_calib(turb_int<spm_threshold)=[];
turb_int(turb_int<spm_threshold)=[];

turb_int(BItmp<db_threshold)=[];
t_calib(BItmp<db_threshold)=[];
BItmp(BItmp<db_threshold)=[];

turb_10log=10*log10(turb_int);
nens_calib=length(BItmp);
%calculate coefficient calibration 10log(Mtbd)=a*BI+b
bical=polyfit(BItmp,turb_10log',1);


%ADCP values (10log(Madcp))
pval=polyval(bical,BItmp);   %over calibration period
pval0=polyval(bical,BItmp0); %over all available TBD period


%correlation matrix
coef=corrcoef(pval,turb_10log');
% relative error 
er=abs(10.^(pval0/10)-turb_int0');
er1=sum(er./turb_int0')/nens_calib;
er2=sqrt(sum(er.^2)/sum(turb_int0.^2));  %=17.8%


% %recherche minimisation somme des ecarts au carre (cf fonction crpi.m):
% %OK idem anterior
% mes=turb_10log;
% BI_cal=BI(turb_alt,t_calib)';
% v0=[0.54 24];       %parametres initiaux 
% [vf,FVAL,EXITFLAG]=fminsearch('crpi',v0)  
% yy=vf(1).*BI+vf(2); %MES estimees
% coef=corrcoef(yy,mes)    %matrice de correlation 


%%
% FIGURES TO CONTROL RESULTS OF CALIBRATION
%--------------------------------------------------------------------------

%cal_turbidity_plot

    fprintf(fidlog,'Turbidity sensor file : %s\n',fileturbtemp);    
    fprintf(fidlog,'Turbidity sensor height above bed : %2.2f\n', turb_alt);   
    fprintf(fidlog,'Calibration period : from %s to %s\n', datestr(min(xrange(t_calib)),1), datestr(max(xrange(t_calib)),1));      
    fprintf(fidlog,'\n');
    fprintf(fidlog,'CALIBRATION RESULTS 10log10(M)=aBI+b\n');
    fprintf(fidlog,'a= %f / b= %f\n',[bical(1) bical(2)]);
     fprintf(fidlog,'\n');
    %   fprintf(fidlog,'number of data used N= %f\n', num2str(nens_calib));
 %   fprintf(fidlog,'correlation coefficient R= %f\n', num2str(coef(1,2)));
 %   fprintf(fidlog,'root mean square error RMS= %f \n', num2str(er2*100));
    
    
display('CALIBRATION RESULTS 10log10(M)=aBI+b')
display(['a= ' num2str(bical(1))])
display(['b= ' num2str(bical(2))])
display(['number of data used N= ' num2str(nens_calib)])
display(['correlation coefficient R= ' num2str(coef(1,2))])
display(['root mean square error RMS= ' num2str(er2*100) ' %'])

%%%time series of ADCP and TBD concentrations 
%%% over all the data period (blue) & used for the calibration (red)
figure; 
subplot(1,2,1)
plot(xrange(t_calib0),turb_int0,'k',xrange(t_calib0),10.^(pval0/10),'b'); ylim([0 max(turb_int0)]);
legend('TBD','ADCP');
%hold on; plot(xrange(t_calib),turb_int,'r');hold off;
set(gca,'Xlim',[fix(xrange(t_calib(1))) fix(xrange(t_calib(end))+1)],'Xtick',[fix(xrange(t_calib(1))):1:fix(xrange(t_calib(end))+1)]);
datetick('x',19,'keeplimits','keepticks');
xlabel(num2str(max(datevec(xrange(1)))));


%%%% scatter plot BI/10log10(M):
gam_bi=[fix(min(BI(turb_alt_ind,t_calib))) fix(max(BI(turb_alt_ind,t_calib)))];
subplot(1,2,2)
plot(BI(turb_alt_ind,t_calib0),10*log10(turb_int0),'.','Color',[.5 .5 .5],'MarkerSize',6);
ylabel('10 * log_{10}( TBD [g.l^{-1}] )','FontSize',7);xlabel('BI [dB_{/1m^3}]','FontSize',7);
hold on; plot(BI(turb_alt_ind,t_calib),turb_10log,'.r','MarkerSize',6); hold off; 
%axis([-60 -30 4 22]);set(gca,'XTick',[-60:5:-20],'FontSize',7);
hold on; plot([gam_bi(1):1:gam_bi(2)],bical(1).*[gam_bi(1):1:gam_bi(2)]+bical(2),'k'); hold off;
text(gam_bi(1)+2,max(turb_10log)-2,strcat(['y=' num2str(bical(1)) ' x+' num2str(bical(2)) '   R=' num2str(coef(1,2))]),'Color','k','FontSize',7);

out_cal(:,1)=BI(turb_alt_ind,t_calib0);
out_cal(:,2)=10*log10(turb_int0);
out_cal(:,4:9)=datevec(xrange(t_calib0));
out_cal(:,3)=turb_int0;
cd(filein_path)
save([filein,'_calibration.txt'],'out_cal','-ASCII')
cd(path_routines);

if mode_exp==0
%%
% *Possibility to calibrate once again*   

button = questdlg('Calibration completed -- Do you want to continue or re-calibrate ?',... 
'Data processing','Continue','Calibrate again','Help','Continue'); 
if strcmp(button,'Continue') 
                                                    
elseif strcmp(button,'Calibrate again') 
    [bical,pval,turb_alt,t_calib0,turb_int0]=cal_turbidity;  
elseif strcmp(button,'Help')  
   disp('Sorry, no help available')   
end


end


