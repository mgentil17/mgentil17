% WADCP PACKAGE 0.3 -- function process_adcp
%%
%
% Backscatter signal processing toolbox
%
%

function process_adcp

%%
% _deployment parameters_
global nbens theta WN WS WB IMM ens_int f a_t ouv Kc EC0 B SL0 R1 Lt binpos binpos_bt inst_mode
global time ve vn vu dir spd a1 a2 a3 a4 amean pitch roll heading c1 c2 c3 c4 cmean BI Mgl
global Batt Tb_adcp Depth_adcp Sb_adcp mode_exp beam_choice ext fidlog filter_opt

%%
% _processing parameters_ 
global  icolor inst_type inst_freq 

%%
% _file names and paths_
global filein_name fileout_name path_routines filein_path fileout_path ncout filein 

if mode_exp==0

[filein_name,filein_path]=uigetfile({'*.mat';'*.nc'},'Open formatted data file....');
[pathstr, name, ext, versn] = fileparts(filein_name);

end

if strcmp(ext, '.nc')
    nc_load
else
    mat_load
end
    cd(filein_path)
    fidlog=fopen([filein,'_proc.log'],'w');
    cd(path_routines)
    
    fprintf(fidlog,'*****************************************************\n');
    fprintf(fidlog,'******  WADCP 0.5 - Log file data processing ********\n');
    fprintf(fidlog,'*****************************************************\n');
    fprintf(fidlog,'\n');    
    fprintf(fidlog,'%s\n',datestr(date,1));
    fprintf(fidlog,'file : %s\n', filein);
    fprintf(fidlog,'\n');
    fprintf(fidlog,'Reminder Deployment information\n');
    fprintf(fidlog,'*****************************************************\n');
    fprintf(fidlog,'\n');
    fprintf(fidlog,'Distance from the bed : IMM : %1.2f\n', IMM);
    fprintf(fidlog,'Blanking distance : WB : %1.2f\n', WB);    
    fprintf(fidlog,'Number of cells : WN : %2.0f\n', WN);    
    fprintf(fidlog,'Cell size : WS : %2.2f\n', WS);        
    fprintf(fidlog,'Transducer diameter : a_t : %1.4f\n', a_t);
    fprintf(fidlog,'Transducer angle aperture : ouv : %1.4f\n', ouv);
    fprintf(fidlog,'Transducer angle : theta : %2.2f\n', theta);
    fprintf(fidlog,'\n');

    
%%
% *Filtering surface and in air data

if mode_exp==1
    if strcmp(filter_opt,'On')
        filter_adcp
    end
end

display('Filtering completed')

tmin=min(time);
tmax=max(time);
depthmax=Inf;
%depthmax=max(Depth_adcp);


figure
hold on
if inst_mode==0
imagesc(time,binpos,a1')
else
    imagesc(time,binpos_bt,a1')
end
set(gca,'XTick',(tmin:(tmax-tmin)/10:tmax),'XMinorTick','on','TickDir','out','XTickLabel',datestr((tmin:(tmax-tmin)/10:tmax),19),'FontSize',14,'Box','on');
shading flat
caxis([50 150])
colorbar
if inst_mode==0
plot(time,Depth_adcp,'w')
end
axis([tmin tmax 0 Inf])
hold off
title('Backscatter amplitude')
%%
% *changing processing parameters
   if mode_exp==0
prompt = {'Kc : dB/counts conversion coefficient (ranging from 0.4 and 0.45)';...
    'EC0 : Internal noise [count] i.e. Echo value in air (ranging from 40 to 65)';...
    'B : noise [dB/1microPa] (ranging from 70 to 96)';...
    'SL0: Source level [dB/1µPa/1m] (ranging from 214 to 218)'};
dlg_title = 'Summary of ADCP characteristics - can be changed by users';
num_lines = 1;
def = {num2str(Kc);num2str(EC0);num2str(B);num2str(SL0)};
tmp = inputdlg(prompt,dlg_title,num_lines,def);
Kc=eval(char(tmp{1}));
EC0=eval(char(tmp{2}));
B=eval(char(tmp{3}));
SL0=eval(char(tmp(4)));



beam_txt=char('beam 1','beam 2','beam 3','beam 4 [RDI only]','averaged beams');
[beam_choice,presult]=listdlg('PromptString','Beam choice for backscatter signal processing','ListString',beam_txt);
   end

if beam_choice==1
beam='beam 1';
elseif beam_choice==2
beam='beam 2';
elseif beam_choice==3
beam='beam 3';
elseif beam_choice==4
beam='beam 4 [RDI only]';
elseif beam_choice==5
beam='averaged beams';
end


%
%* Processing backscatter data

display('Processing ADCP backscatter data...') 

if strcmp(beam,'beam 1')
    backscatter(a1);
elseif strcmp(beam,'beam 2')
    backscatter(a2);    
elseif strcmp(beam,'beam 3')
    backscatter(a3);
elseif strcmp(beam,'beam 4 [RDI only]')
    backscatter(a4);    
elseif strcmp(beam,'averaged beams')
    backscatter(amean);    
end

display('Backscatter processing completed')

%plot_adcp



% plot_adcp


%%
% *Save processed data

if strcmp(ext, '.nc')
    nc_save_proc
else
    mat_save_proc
end



