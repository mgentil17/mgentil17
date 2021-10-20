%% WADCP PACKAGE 0.3 -- Matlab routines for ADCP data post-processing
%%
% Routine used to calibrate ADCP backscatter with moored turbidimeter data
% 
% Input turbidity data : txt file with columns as follow: [Year] [Month] [Day] [Hour] [Minute] [depth] [distance to bottom] [Turbidity g/l]  

function [bical,pval,path1,file2,flag,BI_profile]=cal_turbidity_profiles

global icolor filein_path path_routines
global BI xrange yrange

display('ADCP backscatter calibration in progress...')
cd(filein_path)
[file2,path1]=uigetfile({'*.txt';'*.*'},'Open the .txt file turbidity data....','MultiSelect','on');

filetmp=file2;
file2(end)=filetmp(1);
for i=1:length(file2)-1
file2(i)=filetmp(i+1);
end
clear filetmp
nbf=length(file2);

turb_10log=[];
turb_adcp=[];
turb_flag=[];


for k=1:nbf
cd(path1);
file_turb=load(file2{k});

cd(path_routines);


%%
% only the first time/date
% search and rank CTD profiles

time_profile(k)=datenum(file_turb(1,1),file_turb(1,2),file_turb(1,3),file_turb(1,4),file_turb(1,5),file_turb(1,6));
turb_profile=file_turb(:,9)/1000; %in g/l
%depth_profile=file_turb(:,7);
depth_profile=file_turb(:,8);
% nota : depth must be bottom-referenced (distance to bottom)

% find ADP times around TBD profiles
tadcp_turb=find(xrange>time_profile(k));
tadcp_turb=tadcp_turb(1);
display(['PROFILE N°', num2str(k), '---------'])
display(['TBD profil :', datestr(time_profile(k))])
display(['ADCP time 1:', datestr(xrange(tadcp_turb-1))])
display(['ADCP time 2:', datestr(xrange(tadcp_turb))])
% BI data at the time just after of TBD profile :
%turb_adcp=[turb_adcp; BI(:,tadcp_turb)];
% interpolate BI data at the time of TBD profile :
BI_profile(:,k)=interp1(xrange(tadcp_turb-1:tadcp_turb),BI(:,tadcp_turb-1:tadcp_turb)',time_profile(k),'linear');
turb_adcp=[turb_adcp; BI_profile(:,k)];
%mean BI around time of TBD profile :
%turb_adcp=[turb_adcp; nanmean(BI(:,tadcp_turb-1:tadcp_turb),2)];


% interpolate every TBD profiles over water column
%turb_int=interp1(depth_profile,turb_profile,yrange,'linear');
% or mean over each ADCP bin size :
bin_size=yrange(2)-yrange(1); 
for j=1:length(yrange)
    %yrange(j) is at the center of the bin j
    indices=find(depth_profile>=(yrange(j)-bin_size/2)&depth_profile<=(yrange(j)+bin_size/2));
    turb_int(j)=nanmean(turb_profile(indices));
    clear indices
end
turb_10log=[turb_10log; 10*log10(turb_int)'];

figure(100)
subplot(1,2,1)
plot(turb_int,yrange,'o-r')
hold on;
plot(turb_profile,depth_profile,'.k')
hold off
ylim([0 max(depth_profile)]);
legend('mean','raw'); xlabel('TBD - M [g.l^{-1}]');
title(['PROFILE N°', num2str(k)])

subplot(1,2,2)
plot(BI(:,tadcp_turb-1),yrange,'.-k')
hold on
plot(BI(:,tadcp_turb),yrange,'.-b')
plot(BI_profile(:,k),yrange,'.-r')
hold off
ylim([0 max(depth_profile)]);
legend('t1','t2','t'); xlabel('ADCP - BI [dB_{/1m^3}]');

pause
end

%%
% Ask user-defined profiles to consider for the calibration BI/tbd:
cd(path1);
%create listing of files of charging profiles
fid=fopen('list_files.txt','w');

for k=1:nbf
%fprintf(fid,'%s\n',[num2str(round(datevec(time_profile(k)))),'  ', file2{k},' 01']);
fprintf(fid,'%s\n',[datestr(time_profile(k)),'  ', file2{k},' 01']);
end
fclose(fid);

%ask user to confirm flag
display('confirm flags of profiles to consider for the calibration')
display('01:yes 00:No')
display('don t forget to save the file before to press key...')
open list_files.txt
pause
%reload file of flags
fid=fopen('list_files.txt','r');
for k=1:nbf
%list_files(k)=fscanf(fid,'%s',[length(file2)]);
line=fgetl(fid);
flag(k)=str2num(line(:,end-2:end));
turb_flag=[turb_flag;flag(k)+(0*yrange)'];
end
fclose(fid);
%turb_flag(find(turb_flag==0))=NaN;
%turb_flag(isnan(turb_flag))=0;

cd(path_routines);

%%
%calibration BI/tbd

% without NaN values :
nan_flag(1:length(turb_adcp),1)=1;
nan_flag(isnan(turb_adcp.*turb_10log))=0;   %without NaN values
nan_flag(find(turb_10log.*turb_flag==0))=0; %without user-choosed profiles

turb_adcp2=turb_adcp;
turb_adcp2(find(nan_flag==0))=[];
turb_10log2=turb_10log;
turb_10log2(find(nan_flag==0))=[];

%calculate coefficient calibration 10log(Mtbd)=a*BI+b
bical=polyfit(turb_adcp2,turb_10log2,1);
%ADCP values (10log(Madcp))
pval=polyval(bical,turb_adcp2);   
%correlation matrix
coef=corrcoef(pval,turb_10log2);

nens_calib=length(turb_adcp2);
turb_tbd2=10.^(turb_10log2/10);
% relative error 
er=abs(10.^(pval/10)-turb_tbd2);
er1=sum(er./turb_tbd2)/nens_calib;
er2=sqrt(sum(er.^2)/sum(turb_tbd2.^2));  %=17.8%

display('CALIBRATION RESULTS 10log10(M)=aBI+b')
display(['a= ' num2str(bical(1))])
display(['b= ' num2str(bical(2))])
display(['number of data used N= ' num2str(nens_calib)])
display(['correlation coefficient R= ' num2str(coef(1,2))])
display(['root mean square error RMS= ' num2str(er2*100) ' %'])

%%
% FIGURES TO CONTROL RESULTS OF CALIBRATION
%--------------------------------------------------------------------------
gam_bi=[fix(min(turb_adcp2)) fix(max(turb_adcp2))];
figure
hold on;
plot(turb_adcp,turb_10log,'.','Color',[.5 .5 .5],'MarkerSize',6)
plot(turb_adcp2,turb_10log2,'.r','MarkerSize',6)
plot(turb_adcp2,pval,'b')
%plot([gam_bi(1):1:gam_bi(2)],bical(1).*[gam_bi(1):1:gam_bi(2)]+bical(2),'b'); 
hold off;
%xlim(gam_bi);ylim(gam_bi);
ylabel('10 * log_{10}( TBD [g.l^{-1}] )','FontSize',7);xlabel('BI [dB_{/1m^3}]','FontSize',7);
text(gam_bi(1)+2,max(turb_10log)-2,strcat(['y=' num2str(bical(1)) ' x+' num2str(bical(2)) '   R=' num2str(coef(1,2))]),'Color','k','FontSize',7);

% figure
% hold on
% plot(turb_adcp2,turb_10log2,'.','Color',icolor(:,2))
% hold off



%%
% *Possibility to calibrate once again*   

button = questdlg('Calibration completed -- Do you want to continue or re-calibrate ?',... 
'Data processing','Continue','Calibrate again','Help','Continue'); 
if strcmp(button,'Continue') 
                                                    
elseif strcmp(button,'Calibrate again') 
    [bical,pval,path1,file2,flag,BI_profile]=cal_turbidity_profiles;  
elseif strcmp(button,'Help')  
   disp('Sorry, no help available')   
end

