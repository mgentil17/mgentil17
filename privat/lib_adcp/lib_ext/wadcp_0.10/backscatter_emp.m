
%% WADCP PACKAGE 0.3 -- function backscatter_emp(turbi_calib_choice)
%%
%
% Backscatter signal processing toolbox
%
function backscatter_emp(turbi_calib_choice)

global icolor filein_path path_routines mode_exp
global BI xrange yrange time Depth_adcp Mgl

    
%
if turbi_calib_choice==1
   %calibration of Scatter Index with TBD time serie at one depth
   [ivcal,pcal,turb_alt,t_calib0,turb_int0,nens_calib,t_calib,turb_int]=cal_turbidity ;
elseif turbi_calib_choice==2
   %calibration of Scatter Index with TBD vertical profiles
   [ivcal,pcal,path1,file2,flag,BI_profile]=cal_turbidity_profiles ;
elseif turbi_calib_choice==3
if mode_exp==0
    prompt = {'ival(1)';'ivcal(2)'};
dlg_title = 'ADCP Calibratin coefficients ival(1)*IB+ivcal(2)';
num_lines = 1;
%def = {'0.42794';'2.8907'};
def = {'0.85';'67'};
tmp = inputdlg(prompt,dlg_title,num_lines,def);
ivcal(1)=eval(char(tmp{1}));
ivcal(2)=eval(char(tmp{2}));
end   
    
end

% ivcal(1)=0.85;
% ivcal(2)=-33;

%%
% founded relation applied to all ADCP data 

M_ivcal=10.^((ivcal(1).*BI+ivcal(2))/10); %en g/l !!!
Mgl=M_ivcal;

maxMgl=max(Mgl(1,:));

figure;%colormap(jet(12));
imagesc(xrange,yrange,Mgl*1000);
axis xy; caxis([0 maxMgl*1000]);hcb=colorbar;
set(hcb,'YTick',[0:10:maxMgl*1000]);
set(gca,'Xlim',[fix(xrange(1)) fix(xrange(end)+1)],'Xtick',[fix(xrange(1)):1:fix(xrange(end)+1)]);
datetick('x',19,'keeplimits','keepticks');
xlabel(num2str(max(datevec(xrange(1)))));
ylim([0 max(Depth_adcp)+1]); ylabel('height a.b. [m]');
title('SPM concentration [mg/l]');
hold on;
plot(time,Depth_adcp,'w');
hold off;

% %%
% % figures for validation of calibration
% 
% if turbi_calib_choice==1
%     
% figure;
% plot(xrange(t_calib0),turb_int0,'r',xrange,Mgl(turb_alt,:)','b'); 
% legend('TBD','ADCP');ylim([0 max(turb_int0)]);
% set(gca,'Xlim',[fix(xrange(1)) fix(xrange(end)+1)],'Xtick',[fix(xrange(1)):1:fix(xrange(end)+1)]);
% datetick('x',19,'keeplimits','keepticks');
% xlabel(num2str(max(datevec(xrange(1)))));
% 
% elseif turbi_calib_choice==2
%     
% nbf=length(file2);
% for k=1:nbf
% %reload TBD profiles
% cd(path1);
% file_turb=load(file2{k});
% time_profile(k)=datenum(file_turb(1,1),file_turb(1,2),file_turb(1,3),file_turb(1,4),file_turb(1,5),file_turb(1,6));
% turb_profile=file_turb(:,9)/1000; %in g/l
% %depth_profile=file_turb(:,7);
% depth_profile=file_turb(:,8);
% cd(path_routines);
% 
% %calculate ADCP concentration profiles at the same time
% pval=polyval(ivcal,BI_profile(:,k));
% Madcp_profile(:,k)=10.^(pval./10);
% 
% %figure of all profiles
% nbp=10;
% kk=fix(k/nbp); 
% switch k
%     case nbp, kk=0,
%     case 2*nbp, kk=1,
%     case 3*nbp, kk=2,
%     case 4*nbp, kk=3,
%     case 5*nbp, kk=4,
% end
%     
% figure(kk+110);
% subplot(2,5,k-kk*10)
% if flag(k)==1
%     color_p=icolor(:,5);
% else
%     color_p=icolor(:,7);
% end
% plot(turb_profile,depth_profile,'.','Color',color_p);
% hold on; 
% plot(Madcp_profile(:,k),yrange,'.','Color',icolor(:,2));
% hold off;
% ylim([0 max(depth_profile)]);
% date_profile=datestr(time_profile(k),0);
% title(['P',num2str(k),' ',date_profile]);
% ylabel('h a.b. [m]'), xlabel('M [g.l^{-1}]')
% end
% 
% end
% 
% 
% 
% % 
% % 
% % 
