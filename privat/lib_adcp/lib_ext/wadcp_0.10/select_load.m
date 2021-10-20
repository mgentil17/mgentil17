function select_load

global time ve vn vu a1 a2 a3 a4 amean pitch roll heading c1 c2 c3 c4 cmean spd dir
global path_routines filein_path ncin fileout_name fileout_path filein_name
global Batt Tb_adcp Depth_adcp Sb_adcp
global nbens theta WN WS WB IMM ens_int f a_t ouv Kc EC0 B SL0 Lt binpos binpos_bt inst_mode

if isnan(fileout_name)
    
%cd(filein_path);
 [filein_name,filein_path]=uigetfile({'*.mat;*.nc','ADP formatted Files (*.mat,*nc)'},'Open data file....');
cd(filein_path)
[pathstr, name, ext, versn] = fileparts(filein_name);
if strcmp(ext,'nc')
    nc_load
else
    mat_load
end

else
    cd(fileout_path)
[pathstr, name, ext, versn] = fileparts(fileout_name);
if strcmp(ext,'nc')
    nc_load
else
    mat_load
end
    
end
cd(path_routines)


