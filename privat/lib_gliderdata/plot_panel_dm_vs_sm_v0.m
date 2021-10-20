function plot_panel_dm_vs_sm_v0(dm_vel_filt,sm_vel_filt,DAC_cur,dpath)

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Comparison of  :
%
% Direct method - Difference in measured ADCP speeds and static speed
% of the glider determined by the flight model
%
% Shear method - The vertical derivatives of the current profiles
% (according to each component)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

dm = dm_vel_filt;
sm = sm_vel_filt;
dac = DAC_cur;
names = fieldnames(dm_vel_filt);                       % number of sections

for ii = 1:length(names)
    %======================================================================
    % DEFINE VARIABLES
    %======================================================================
    time = dm.(names{ii}).time;
    dpth = dm.(names{ii}).dpth;
  
    % Depth vecor
    max_dpth = (max(dpth));
    [~,I] = max(max_dpth);
    z = dpth(:,I);
    
    % DAC
    vx_mean = nanmean(dac.(names{ii}).vx);
    vy_mean = nanmean(dac.(names{ii}).vy);
    
    % Direct method
    u_dm = dm.(names{ii}).u;
    v_dm = dm.(names{ii}).v;
    ustd_dm = dm.(names{ii}).ustd;
    vstd_dm = dm.(names{ii}).vstd;
    
    % Shear method
    u_sm = sm.(names{ii}).u;
    v_sm = sm.(names{ii}).v;
    ustd_sm = dm.(names{ii}).ustd;
    vstd_sm = dm.(names{ii}).vstd;
    
    % Depth
    z_dm = dm.(names{ii}).dpth;
    z_sm = sm.(names{ii}).dpth;
    
    % Limits
    lim_min = -0.5;
    lim_max = 0.5;
    
    %======================================================================
    % Velocity interpolation on a regular depth grid
    %======================================================================
    % Define reference depth vector
    z_ref = z;
    
    % Remove NaN value
    ind_NaN = find(isnan(z_ref));
    z_ref(ind_NaN) = [];
    
    % Define size matrix output
    [row,~] = size(z_ref);
    [~,col] = size(u_dm);
    
    % Initialize matrix
    u_int_dm = NaN(row,col);
    v_int_dm = NaN(row,col);
    ustd_int_dm = NaN(row,col);
    vstd_int_dm = NaN(row,col);
    u_int_sm = NaN(row,col);
    v_int_sm = NaN(row,col);
    ustd_int_sm = NaN(row,col);
    vstd_int_sm = NaN(row,col);
    
    % Interpolation on refernce depth vector
    for i = 1:col
        % Direct method
        iok = find(~isnan(z_dm(:,i)));
        u_int_dm (:,i) = interp1(z_dm(iok,i),u_dm(iok,i),z_ref);
        v_int_dm (:,i) = interp1(z_dm(iok,i),v_dm(iok,i),z_ref);
        ustd_int_dm (:,i) = interp1(z_dm(iok,i),ustd_dm(iok,i),z_ref);
        vstd_int_dm (:,i) = interp1(z_dm(iok,i),vstd_dm(iok,i),z_ref);
        
        % Shear method
        iok = find(~isnan(z_sm(:,i)));
        u_int_sm (:,i) = interp1(z_sm(iok,i),u_sm(iok,i),z_ref);
        v_int_sm (:,i) = interp1(z_sm(iok,i),v_sm(iok,i),z_ref);
        ustd_int_sm (:,i) = interp1(z_sm(iok,i),ustd_sm(iok,i),z_ref);
        vstd_int_sm (:,i) = interp1(z_sm(iok,i),vstd_sm(iok,i),z_ref);
    end
    
    %======================================================================
    % Compute variables of data distribution
    %======================================================================
    % profiles number
    nbins = size(u_int_dm,2);
    
    % Compute the difference between the shear method and direct method
    u = u_int_dm - u_int_sm;
    v = v_int_dm - v_int_sm;
    
    % Compute the 5 and 95 percentiles
    u_perc_5 = ones(size(u_int_dm,2),1)*(prctile(u(:),5));
    u_perc_95 = ones(size(u_int_dm,2))*(prctile(u(:),95));
    
    v_perc_5 = ones(size(u_int_dm,2))*(prctile(v(:),5));
    v_perc_95 = ones(size(u_int_dm,2))*(prctile(v(:),95));
    
    %======================================================================
    % Plot Velocity sections 
    % from direct method vs shear method
    % and data distribution
    %======================================================================
    h = figure('units','normalized','outerposition',[0 0 1 1]);   
    
    % Plot velocity section
    iok = find(~isnan(z));
    
    subplot 421; 
    contourf(time,-z(iok),u_int_dm(iok,:)); caxis([lim_min lim_max]); ylim([-120 0]); colorbar;
    hold on;
    [M,c] = contour(time,-z(iok),u_int_dm(iok,:),'-k');    
    clabel(M,c,'Fontweight','bold');
    cmocean('balance');
    title('Direct method : u (m s^-^1)');
    ylabel('Depth (m)','Fontweight','bold'); 
    ax = gca;
    ax.FontSize = 12;
    set(ax,'TickDir','out');
    set(ax,'xticklabel',[]);
    
    subplot 422; 
    contourf(time,-z(iok),v_int_dm(iok,:)); caxis([lim_min lim_max]); ylim([-120 0]); colorbar;
    hold on;
    [M,c] = contour(time,-z(iok),v_int_dm(iok,:),'-k');    
    clabel(M,c,'Fontweight','bold');
    cmocean('balance');    
    title('Direct method : v (m s^-^1)');
    ax = gca;
    ax.FontSize = 12;
    set(ax,'TickDir','out');
    set(ax,'xticklabel',[]);
    
    subplot 423;
    contourf(time,-z(iok),u_int_sm(iok,:)); caxis([lim_min lim_max]); ylim([-120 0]); colorbar;
    hold on;
    [M,c] = contour(time,-z(iok),u_int_sm(iok,:),'-k');    
    clabel(M,c,'Fontweight','bold');
    cmocean('balance');    
    title('Shear method : u (m s^-^1)');
    ylabel('Depth (m)','Fontweight','bold');
    ax = gca;
    ax.FontSize = 12;
    set(ax,'TickDir','out');
    set(ax,'xticklabel',[]);
    
    subplot 424; 
    contourf(time,-z(iok),v_int_sm(iok,:)); caxis([lim_min lim_max]); ylim([-120 0]); colorbar;
    hold on;
    [M,c] = contour(time,-z(iok),v_int_sm(iok,:),'-k');    
    clabel(M,c,'Fontweight','bold');
    cmocean('balance');
    title('Shear method : v (m s^-^1)');
    ax = gca;
    ax.FontSize = 12;
    set(ax,'TickDir','out');
    set(ax,'xticklabel',[]);
    
    subplot 425; 
    contourf(time,-z(iok),u(iok,:)); caxis([lim_min lim_max]); ylim([-120 0]); colorbar;
    hold on;
    [M,c] = contour(time,-z(iok),u(iok,:),'-k');    
    clabel(M,c,'Fontweight','bold');
    cmocean('balance');    
    title('Difference from both methods : u (m s^-^1)')
    ylabel('Depth (m)','Fontweight','bold');
    ax = gca;
    ax.FontSize = 12;
    set(ax,'TickDir','out');
    datetick('x','HH:MM','keepticks');
    
    subplot 426; 
    contourf(time,-z(iok),v(iok,:)); caxis([lim_min lim_max]); ylim([-120 0]); colorbar;
    hold on;
    [M,c] = contour(time,-z(iok),v(iok,:),'-k');    
    clabel(M,c,'Fontweight','bold');
    cmocean('balance');
    title('Difference from both methods : v (m s^-^1)');
    ax = gca;
    ax.FontSize = 12;
    set(ax,'TickDir','out');
    datetick('x','HH:MM','keepticks');
    
    % Plot data distribution    
    subplot 427;
    m = histfit(u(:),nbins,'kernel');
    grid minor;
    x = get(m(2),'xdata');
    y = get(m(1),'ydata');
    set(gca, 'XLim', [min(x) max(x)]);
    set(gca, 'YLim', [min(y) max(y)]);
    hold on;
    plot(u_perc_5(1:length(y),1),y,'--k','Linewidth',1.5);    
    plot(u_perc_95(1:length(y),1),y,'--k','Linewidth',1.5);
    title('Kernel distribution : u(m.s^-^1)');
    xlabel('u (m s^-^1)','Fontweight','bold');
    set(gca,'YTick', []);
    ax = gca;
    ax.FontSize = 12;
    
    subplot 428;
    m = histfit(v(:),nbins,'kernel');
    grid minor;
    x = get(m(2),'xdata');
    y = get(m(1),'ydata');
    set(gca, 'XLim', [min(x) max(x)]);
    set(gca, 'YLim', [min(y) max(y)]);
    hold on;
    plot(v_perc_5(1:length(y),1),y,'--k','Linewidth',1.5);    
    plot(v_perc_95(1:length(y),1),y,'--k','Linewidth',1.5);
    title('Kernel distribution');
    xlabel('v (m s^-^1)','Fontweight','bold');
    set(gca,'YTick', []);
    ax = gca;
    ax.FontSize = 12; 
    
    % Save figure
    cd(dpath)
    saveas(h,'abs_vel_comparisons.jpeg');
    
end

end