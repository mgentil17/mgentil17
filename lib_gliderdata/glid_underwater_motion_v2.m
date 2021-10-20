function [i_GuwM] = glid_underwater_motion_v2(int_glider,int_glider_derived,param0,idx_time,thresh,options,in_comp,plot_out,dpath)

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
% Compute Glider underwater motion 
% Use the steady state flight model developped by Merckelbach et al 2009
% Optimized glider parameters
%
% version 2
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%


%==========================================================================
% PREPARING GLIDER MATRIX FOR THE OPTIMIZATION
%==========================================================================
glider.time = int_glider.nav(:,1);                                         % time
glider.lat = int_glider.nav(:,2);                                          % latitude
glider.lon = int_glider.nav(:,3);                                          % longitude
glider.pres = int_glider.sci(:,3);                                         % pressure
glider.dens = int_glider_derived.dens;                                     % density
glider.pitch = int_glider.nav(:,5);                                        % pitch
glider.balla = int_glider.nav(:,4);                                        % balla
glider.head = int_glider.nav(:,6);                                         % heading
glider.temp = int_glider.sci(:,2);                                         % temperature

%==========================================================================
% COMPUTE GLIDER VERTICAL VELOCITIES
%==========================================================================
% Convert time units
time = glider.time;
sinceEden = datenum(1970,1,2) - datenum(0000,1,1);
time = time - sinceEden;
time = time*86400;

wg_ctd = -diff(glider.pres)./diff(time);   
Wglider = wg_ctd;

%==========================================================================
% OPTIMIZATION TO RETRIEVE ug/wg
%==========================================================================
% Defining the global variables useful for the cost.m optimising function
global pres dens pitch balla Wglider ind temp

% Cutting in intervals corresponding to the parameter "period" for the optimisation
time_decoupage = nanmin(time):idx_time:nanmax(time);
fprintf('%d timesteps for time decoupage\n',length(time_decoupage));

% Initialized optimized variables
Wm_opt = nan(1,length(Wglider));                                           % w flight model 
Um_opt = nan(1,length(Wglider));                                           % U horizontal velocities flight model along the glider path
attm_opt = nan(1,length(Wglider));
Hm_opt = nan(1,length(Wglider));                                           % optimized horizontal model velocity
um_opt = nan(1,length(Wglider));                                           
vm_opt = nan(1,length(Wglider));  
Vg = nan(1,length(time_decoupage)-1);
eps = nan(1,length(time_decoupage)-1);
Cd0 = nan(1,length(time_decoupage)-1);
cost_opt = nan(1,length(time_decoupage)-1);

%--------------------------------------------------------------------------
% Computation of flight variables using initial glider's parameters 
if in_comp == 1
    pres = glider.pres; pres(pres<thresh ) = NaN; ind_NaN_pres = find(isnan(pres));
    dens = real(glider.dens); dens(ind_NaN_pres) = NaN;
    pitch = glider.pitch; pitch(ind_NaN_pres) = NaN;
    balla = glider.balla; balla(ind_NaN_pres) = NaN;
    temp = glider.temp; temp(ind_NaN_pres) = NaN;
    head = glider.head; head(ind_NaN_pres) = NaN;  
    [U0,Wmodel0,att0,Fg0,Fb0,Fl0,Fd0] = flight_model_v2(pres,dens,pitch,balla,temp,param0(1),param0(2),param0(3),param0(4));
end
%--------------------------------------------------------------------------

for i = 1:length(time_decoupage)-1    
    disp(['### days ' num2str(i) '/' num2str(length(time_decoupage)-1) ]);
    ind = find(time>time_decoupage(i) & time<=time_decoupage(i+1));

    % Set surface data to NaN
    pres = glider.pres(ind); pres(pres<thresh ) = NaN; ind_NaN_pres = find(isnan(pres));
    dens = real(glider.dens(ind)); dens(ind_NaN_pres) = NaN;
    pitch = glider.pitch(ind); pitch(ind_NaN_pres) = NaN;
    balla = glider.balla(ind); balla(ind_NaN_pres) = NaN;
    temp = glider.temp(ind); temp(ind_NaN_pres) = NaN;
    head = glider.head(ind); head(ind_NaN_pres) = NaN;    
    
    % Optimization of parameters (eps,Cd0,Vg) by minimazing cost fuction 'cost'
    [param_opt,fval,exitflag,output]= fminsearch('cost_v2',param0,options);
    
    % Computation of flight variables using optimized glider's parameters
    [U,Wmodelopt,att,Fg,Fb,Fl,Fd] = flight_model_v2(pres,dens,pitch,balla,temp,param_opt(1),param_opt(2),param_opt(3),param_opt(4));

    Wm_opt(ind) = Wmodelopt;                                               % optimized glider vertical velocity 
    Um_opt(ind) = U;                                                       % optimized glider horizontal velocity along a glider path
    attm_opt(ind) = att*(180/pi);                                          % optimized attack angle 
    Vg(i) = param_opt(1);                                                  % optimizedglider volume 
    eps(i) = param_opt(2);                                                 % optimizedglider hull compressibiliy 
    Cd0(i) = param_opt(3);                                                 % optimized glider parasite drag 
    mg(i) = param_opt(4);                                                  % optimized glider mass
    
    % Compute horizontal glider velocitites component
    hm_temp = U .* cosd(pitch+(att*(180/pi))); 
    Hm_opt(ind) = hm_temp;                                                 % optimized horizontal model velocity
    [um_opt(ind),vm_opt(ind)] = polar2uv(head,hm_temp);
    
    [cout] = cost_v2(param_opt);
    cost_opt(i) = cout;
    
    fprintf('%10s %10s %10s %10s \r \n','Vg','eps','Cd0', 'mg');
    result = [param_opt(1), param_opt(2), param_opt(3),param_opt(4)];
    fprintf('%10.5f %10.5e %10.5f %10.5f \n',result);
end

%==========================================================================
% COMPARISON OF VERTICAL VELOCITIES
%==========================================================================
RMSE_W_glider = sqrt(nanmean((Wglider - Wm_opt').^2));          % in m.s-1

% Colors
pink = rgb('DeepPink');
green = rgb('YellowGreen');
blue = rgb('SteelBlue');

if plot_out
    % Plot
    h = figure('units','normalized','outerposition',[0 0 1 1]);
    plot(glider.time(2:end),Wglider,'Color',green,'Linewidth',2);                 % Vertical glider velocity
    hold on; grid minor;
    plot(glider.time(2:end),Wm_opt','--','Color',pink,'Linewidth',2);             % Optimized vertical glider velocity
%     plot(glider.time(2:end),Wglider-Wm_opt','.','Color',blue,'MarkerSize',5);     % Differences between model and measurements glider velocity
    plot(glider.time(2:end),Wglider-Wm_opt','-','Color',blue,'Linewidth',2);     % Differences between model and measurements glider velocity
    ylabel('Speed [m s^-^1]','Fontweight','bold');
    ylim([-0.2 0.35]);
    title('Glider Vertical Velocities');
    legend('W_o_b_s','W_m_o_d','W_d_i_f_f');
    xlim([glider.time(2) glider.time(2000)]);
    ax = gca;
    ax.FontSize = 12;
    set(ax,'TickDir','out');
    datetick('x','mm/dd HH:MM','keepticks');
    cd(dpath);
    saveas(h,'Glider_vertical_velocities.jpeg');
end

%==========================================================================
% ALLOCATE OUTPUTS
%==========================================================================
i_GuwM.time = glider.time(2:end);                                          % due to the difference on vertical velocities
i_GuwM.pres = glider.pres(2:end);
i_GuwM.att = attm_opt;
i_GuwM.wg = Wm_opt;
i_GuwM.ug = um_opt;
i_GuwM.vg = vm_opt;
i_GuwM.rmse = RMSE_W_glider;
i_GuwM.Hg = Hm_opt;
i_GuwM.Ug = Um_opt;
i_GuwM.Vg = Vg.*1e-2;
i_GuwM.eps = eps.*1e-10;
i_GuwM.Cd0 = Cd0.*1e-2;
i_GuwM.mg = mg.*1e1;

%==========================================================================
% PLOT OPTIMIZED VELOCITIES
%==========================================================================
% Comparison of vertical and horizontal velocities with  initial and
% optimized parameters
if in_comp == 1
    h1 = figure('units','normalized','outerposition',[0 0 1 1]);
    
    subplot 211; 
    plot(glider.time,Wmodel0,'Color',green,'Linewidth',2);
    hold on; grid minor;
    plot(i_GuwM.time,i_GuwM.wg,'-','Color',pink,'Linewidth',2);
    legend('Initial_p_a_r_a_m', 'Optimized_p_a_r_a_m');
    title('Glider Vertical Velocities');
    ylabel('Speed [m s^-^1]','Fontweight','bold');
    ylim([-0.2 0.35]);
    xlim([glider.time(2) glider.time(2000)]);
    ax = gca;
    ax.FontSize = 12;
    set(ax,'TickDir','out');
    datetick('x','mm/dd HH:MM','keepticks');
    set(ax,'xticklabel',[]);    
    
    subplot 212; 
    plot(glider.time,U0,'Color',green,'Linewidth',2);
    hold on; grid minor;
    plot(i_GuwM.time,i_GuwM.Ug,'-','Color',pink,'Linewidth',2);
    legend('Initial_p_a_r_a_m', 'Optimized_p_a_r_a_m');
    title('Glider Horizontal Velocities');
    ylabel('Speed [m s^-^1]','Fontweight','bold');
    ylim([0.2 0.6]);
    xlim([glider.time(2) glider.time(2000)]);
    ax = gca;
    ax.FontSize = 12;
    set(ax,'TickDir','out');
    datetick('x','mm/dd HH:MM','keepticks');
    
    cd(dpath);
    saveas(h1,'Glider_FM_optimization.jpeg');
end

end