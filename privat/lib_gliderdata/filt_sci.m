function [sci,opt,out_sci,out_opt] = filt_sci(sci,thresh)

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Cleaning glider science data
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

no_data = 1:length(sci.data(:,1));

%==========================================================================
% Variables
%==========================================================================
opt.sensors = [sci.sensors(1) sci.sensors(5:7)];
opt.data = [sci.data(:,1) sci.data(:,5:7)];

sci.sensors = sci.sensors(1:4);
sci.data = sci.data(:,1:4);

%==========================================================================
% Kill NaN
%==========================================================================
i_bad = find(isnan(sci.data(:,1)) | isnan(sci.data(:,2)) | ...
    isnan(sci.data(:,3)) | isnan(sci.data(:,4)));

fprintf('##### %d Remove NaN values of science data\n #####',length(i_bad));
sci.data(i_bad,:) = [];

ii_bad = find(isnan(opt.data(:,1)) | isnan(opt.data(:,2)) | ...
    isnan(opt.data(:,3)) | isnan(opt.data(:,4)));

fprintf('##### %d Remove NaN values of optical data\n #####',length(ii_bad));
opt.data(ii_bad,:) = [];

%==========================================================================
% Remove null values
%==========================================================================
iii_bad = find(sci.data(:,2)==0 | sci.data(:,3)==0 | sci.data(:,4)==0);
fprintf('##### %d Remove null values of science data\n #####',length(iii_bad));
sci.data(iii_bad,:) = [];

iv_bad = find(opt.data(:,2)==0 | opt.data(:,3)==0 | opt.data(:,4)==0);
fprintf('##### %d Remove null values of optical data\n #####',length(iv_bad));
opt.data(iv_bad,:) = [];

%==========================================================================
% Remove gap of pressure
%==========================================================================
% sci.data(:,3)=sci.data(:,3);
dP = diff(sci.data(:,3));              % difference of pressure
v_bad = find(dP <= thresh & dP >= -thresh);
v_bad = v_bad +1;

fprintf('##### %d bad scans of derivative pressure\n #####',length(v_bad));

% figure; plot(sci.data(:,1),-sci.data(:,3),'.b'); hold on; plot(sci.data(v_bad,1),-sci.data(v_bad,3),'or')

sci.data(v_bad,:) = [];
out_sci = length(i_bad) + length(iii_bad) + length(v_bad);
out_opt = length(ii_bad) + length(iv_bad);

end