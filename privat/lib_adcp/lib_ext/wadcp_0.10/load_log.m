function load_log

global fileout_name
global path_routines
global matsave

matsave=[];
path_routines=cd;
[filein_name,filein_path]=uigetfile('*.TXT','Open the .Dat file with the pre-processed data','MultiSelect','on');
cd(filein_path);

if iscell(filein_name)

filetmp=filein_name;
filein_name(end)=filetmp(1);


for i=1:length(filein_name)-1
filein_name(i)=filetmp(i+1);
end
clear filetmp

for k=1:length(filein_name)

mat=load(filein_name{k});

matsave=[matsave;mat];
clear mat filein
filein_name=filein_name(1);
end



else
    
   mat=load(filein_name);
   matsave=[matsave;mat];    
    
end


[fileout_name,fileout_path]=uiputfile('*.mat','Name of the matlab file where to save formatted raw data:',filein_name(1:end-4));

cd(fileout_path);
save(fileout_name,'matsave')

cd(path_routines)

save_log_mat
