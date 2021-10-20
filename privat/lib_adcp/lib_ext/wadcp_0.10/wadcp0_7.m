%% WADCP PACKAGE 0.7 -- Matlab routines for ADCP data post-processing
% *Introduction to WADCP0.7 : Evaluation of currents, waves and backscatter for moored and down-looking ADCP*
%
% _Written by C. Tessier and R. Verney, IFREMER |contact : caroline.tessier and romaric.verney@ifremer.fr | Download section : http://www.ifremer.fr/dyneco/_ 
%
% This suite of Matlab(R) routines is dedicated to the post-processing of
% ADP data [from Nortek, RDI and Sontek instruments].
%
% These routines were written to be independent of any Matlab Toolboxes. Morevover, they have been successfully executed from Matlab Release version 6.5.
% If users encounter problems when running please contact us! 
% 
% In order to correctly execute the program, all the routines must be saved in a identical folder, but not necessarily where the data are saved.
% Several routines are called successively from wadcpx_x.m and are detailed below.
%
% Latest update : August 25, 2009
%
% _What's new?_
% * Possibility to use CTD profiles to calibrate the backscatter index 
% * Integration of bottom track mode 
% * Change name of variable in english
% * Saving all variables in matlab structure in addition to netcdf
% * Using a namelist instead of the interactive menu, chose wadcp0_6(1) to
% enter this expert mode and modifiy wadcp_config.txt
% * save main information in a log file
% * Saving output in matlab format
% * Using an SPM threshold for BI calibration
%_Bugs to solve_ 
%
% 
%

%% Function wadcp0_7
% Main function that calls all the subroutines to process the data. 

function wadcp0_7(varargin)

%clear all
%clear

%%
% *Declaration of global variables* 
%
% Common library to all routines, even if variables are not  used in this
% routine
%


%%
% _deployment parameters_
global nbens theta WN WS WB IMM ens_int f a_t ouv Kc EC0 B SL0 R1 Lt binpos
global time ve vn vu dir spd a1 a2 a3 a4 amean pitch roll heading c1 c2 c3 c4 cmean
global Batt Tb_adcp Depth_adcp Sb_adcp Mgl

%%
% _processing parameters_ 
global  icolor inst_type inst_freq inst_mode mode_exp qst_mode beam_choice fchoice cellskip cor_threshold postproc_choice tminproc tmaxproc
global  varia_batt qst_nf backs_calib_choice fidlog spm_threshold filter_opt db_threshold

%%
% _file names and paths_
global filein_name fileout_name path_routines filein_path fileout_path ncout filein ext fileturbtemp turb_alt


warning off all

%%
% * User-defined colors for all plots
icolor(:,1)=[162 38 109]/255;  %purple
icolor(:,2)=[0 96 255]/255;    %blue
icolor(:,3)=[80 162 84]/255;   %green
icolor(:,4)=[192 119 38]/255;  %brown
icolor(:,5)=[217 23 222]/255;  %pink
icolor(:,6)=[167 211 247]/255; %celeste
icolor(:,7)=[0.5 0.5 0.5];     %grey

%%
% *Dialog box open data files*
%
% * As data formats are dependant on the instrument manufacturer, all data are first loaded and save in a single common format 

path_routines=cd;
fileout_name=NaN;
if isempty(varargin)
    mode_exp=0;
else
    
    mode_exp=varargin{1};
end

if mode_exp==0; % full procedure (or 1 = expert)

qst_proc=questdlg('Load raw data, process data or playback netcdf files?','Process choice','load raw data','process data','playback','load raw data');

if strcmp(qst_proc,'playback')
   load_nc
   
elseif strcmp(qst_proc,'load raw data')

    qst_mode=questdlg('Instrument mode (Moored or Survey)?','Process choice','Moored','Survey','Moored');
if strcmp(qst_mode,'Moored')
    inst_mode=0;
else % survey
    inst_mode=1;
end

  [filein_name,filein_path]=uigetfile({'*.mat;*.ctl;*.hdr;*.000;*.001;*.002;*.PD0','ADP raw Files (*.mat,*.ctl,*.hdr,*.000,*.001,*.002,*.PD0'},'Open data file....');
  filein=filein_name(1:end-4);
  
  
  
  if strcmp(filein_name(end-2:end),'mat')
    inst_type='RDI';
    open_adcp_rdi
  elseif strcmp(filein_name(end-2:end),'ctl')
    inst_type='Sontek';
    open_adcp_sontek
  elseif strcmp(filein_name(end-2:end),'hdr')
    inst_type='Nortek';
    open_hdr_adcp_nortek
    open_adcp_nortek_2
  elseif (strcmp(filein_name(end-2:end),'000') || strcmp(filein_name(end-2:end),'PD0') || strcmp(filein_name(end-2:end),'001') || strcmp(filein_name(end-2:end),'002'))
    inst_type='RDI';
    open_rdi_binary

  else
    display('Warning : unrecognized data format')
  end

%ask for ref time in ADCP file before to save
  t_lag=0;
  prompt = {'t_lag [hours] : time-lag of ADCP data (from TU reference) '};
  dlg_title = 'Time reference of ADCP - to be confirmed by user';
  num_lines = 1;
  def = {num2str(t_lag)};
  tmp = inputdlg(prompt,dlg_title,num_lines,def);
  t_lag=eval(char(tmp{1}));
%change time vector in UT before to save   
  time=time-t_lag/24;
  
    qst_output=questdlg('Output format?','Process choice','Netcdf','Matlab','Matlab');
  if strcmp(qst_output,'Netcdf')
%save data in NetCDF file (in UT ref time)
  nc_save
  else
  mat_save
  end
  
  display('raw data saved as netcdf data')
  display('restart wadcp for processing...')

elseif strcmp(qst_proc,'process data')

  process_adcp
%  plot_adcp

end
 
else %expert mode, read wadcp_config.txt
    

    fid=fopen('wadcp_config.txt');
    for i=1:8
        fgetl(fid);
    end
    
    qst_proc=sscanf(fgetl(fid),'%f %s',1);
    qst_mode=sscanf(fgetl(fid),'%f %s',1);
    filetemp=sscanf(fgetl(fid),'%s %s',1);
    t_lag=sscanf(fgetl(fid),'%f %s',1);
    qst_output=sscanf(fgetl(fid),'%s %s',1);
    inst_freq=sscanf(fgetl(fid),'%s %s',1);
    filter_opt=sscanf(fgetl(fid),'%s %s',1);
    fgetl(fid);
    
  IMM=sscanf(fgetl(fid),'%f %s',1);
  WB=sscanf(fgetl(fid),'%f %s',1);
  a_t=sscanf(fgetl(fid),'%f %s',1);
  theta=sscanf(fgetl(fid),'%f %s',1);
  ouv=sscanf(fgetl(fid),'%f %s',1);
  Kc=sscanf(fgetl(fid),'%f %s',1);
  EC0=sscanf(fgetl(fid),'%f %s',1);
  B=sscanf(fgetl(fid),'%f %s',1);
  SL0=sscanf(fgetl(fid),'%f %s',1);
  beam_choice=sscanf(fgetl(fid),'%f %s',1);
  varia_batt=sscanf(fgetl(fid),'%f %s',1);
  qst_nf=sscanf(fgetl(fid),'%s %s',1);
  backs_calib_choice=sscanf(fgetl(fid),'%f %s',1);
  fileturbtemp=sscanf(fgetl(fid),'%s %s',1);
  spm_threshold=sscanf(fgetl(fid),'%f %s',1);
  db_threshold=sscanf(fgetl(fid),'%f %s',1);
  turb_alt=sscanf(fgetl(fid),'%f %s',1);
  fchoice=sscanf(fgetl(fid),'%f %s',1);
  cellskip=sscanf(fgetl(fid),'%f %s',1);
  cor_threshold=sscanf(fgetl(fid),'%f %s',1);
  postproc_choice=sscanf(fgetl(fid),'%f %s',1);
  tminproc=sscanf(fgetl(fid),'%s %s',1);
  tmaxproc=sscanf(fgetl(fid),'%s %s',1);
  tminproc=datenum(tminproc);
  tmaxproc=datenum(tmaxproc);
  f=str2num(inst_freq);

      [filein_path,filein,ext] = fileparts(filetemp);
    filein_name=[filein,ext];
%     filein_path=strrep(filetemp,filein_name,'');
    fclose(fid)

if qst_proc==1

  
    
    if strcmp(ext,'.ctl')
        inst_type='Sontek';
        open_adcp_sontek
    elseif strcmp(ext,'.hdr')
        inst_type='Nortek';
        open_hdr_adcp_nortek
        open_adcp_nortek_2
    elseif strcmp(ext,'.mat')
        inst_type='RDI';
        open_adcp_rdi
    elseif (strcmp(ext,'.000') || strcmp(ext,'.PD0'))
        inst_type='RDI';
        open_rdi_binary
    else
        display('Warning : unrecognized data format')        
    end
  
    
      time=time-t_lag/24;
    
  if strcmp(qst_output,'Netcdf')
%save data in NetCDF file (in UT ref time)
  nc_save
  else
  mat_save
  %mpat_save_test2
  end
    
    
elseif qst_proc==2
      process_adcp
elseif qst_proc==3
    load_nc
end



end

