function [ind_start_sec, ind_end_sec] = delim_glid_section_v2(glid)

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% User define glider section according to the deployment
% It's a manual part
% The user must make changes in this code
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


%==========================================================================
% Define variables
%==========================================================================
time = glid.Dac(:,1);
lat = glid.Dac(:,2);
dpth = glid.nav(:,end);

%==========================================================================
% Plot Glider Path
%==========================================================================
h = figure('units','normalized','outerposition',[0 0 1 1]);
plot(time,lat,'Linewidth',2); grid minor;
ax = gca;
ax.FontSize = 12;
set(ax,'TickDir','out');
datetick('x','dd/mm HH:MM','keepticks');
ylabel('latitude (°)','Fontweight','bold');
pause(5);

%==========================================================================
% Select Glider Sections
%==========================================================================
% Must be adjusted by the user
% and allocate outputs
disp('#### Manual delimitation of the sections to be studied: User defined it in delim_glid_section_v2.m from line 151 to 165 ####');
disp('#### For each selected section click twice on the figure (start and end point) #####');

% Indicates the number of sections do you want process
answer = inputdlg('Enter sections number to process:',...
             'Input', [1 50]);
user_val = str2num(answer{1});      

val_sec = [];
ind_sec = [];
for ii = 1:user_val
    % Select the sections
    shg;
    [x,~] = ginput(2);
    tmp = [x(1) x(2)];            
    val_sec = [val_sec;tmp];
    
    % Find index corresponding to time value of sections
    [~,id_start] = min(abs(time-tmp(1)));
    [~,id_end] = min(abs(time-tmp(2)));
    tmpp = [id_start id_end];
    ind_sec = [ind_sec;tmpp];
end

% % For the example toolbox
% f = msgbox({'Selection of only a section for the rest of processing (example)'; 'because very time consumming.'},'Warn','warn');
% ind_start_sec = ind_sec(2,1);
% ind_end_sec = ind_sec(2,2);

% For all deployments
ind_start_sec = ind_sec(:,1);       % first column ind_start_section
ind_end_sec = ind_sec(:,2);         % second column: ind_end_section

end
