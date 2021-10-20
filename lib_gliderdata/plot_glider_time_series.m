function plot_glider_time_series(X,Xbis,Y,Ymax,Z,lat,Col,n,Var,glid_l,deploy,dpath)

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Plot Glider Time-Series
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

names = fieldnames(X);

% Color
grey = [0.5 0.5 0.5];

for i = 1:length(names)
    
    h = figure;
    set(h,'Units','Normalized','Outerposition',[0 0 0.85 1]);
    
    % Define Variables
    X_plot = X.(names{i});
    Xbis_plot = Xbis.(names{i});
    Y_plot = Y.(names{i});
    Z_plot = Z.(names{i});
    lat_plot = lat.(names{i});
    
    % Limits
    Xmax = nanmean(Xbis_plot(:,end));
    Xmin = nanmean(Xbis_plot(:,1));
    
    % Temporal variable
    subplot(2,1,1);
    box on; hold on;
    pcolor(Xbis_plot,-Y_plot,Z_plot); shading interp;
    colormap(Col);
    set(gca,'Color',grey)
    cb = colorbar;
    if Var == 8
        cb.Ticks = [0:0.25:0.5];                   % for Richardson number
    end
    ylabel('Depth (m)', 'FontSize', 12, 'FontWeight', 'bold');
    SW = [Xbis_plot(1,10) -140];
    text(SW(1), SW(2), n, 'VerticalAlignment', 'bottom', ...
    'HorizontalAlignment', 'left', 'color', 'black', 'Fontsize', 16, ...
    'Fontweight','bold')
    ylim([-Ymax 0]);
    xlim([Xmin Xmax]);
    ax = gca;
    ax.FontSize = 12;
    set(ax,'TickDir','out');
    datetick('x','mm/dd HH:MM','keepticks');
    title([glid_l '/' deploy '/' names{i}]);
    
    % Distance
    subplot(2,1,2);
    box on; hold on; grid minor;
    scatter(nanmean(Xbis_plot),nanmean(lat_plot),[],nanmean(X_plot),'filled');
    colormap(Col);
    ylabel('Latitude (°)', 'FontSize', 12, 'FontWeight', 'bold');
    c = colorbar;
    xlim([Xmin Xmax]);
    text(SW(1), min(nanmean(lat_plot))+0.04, 'Distance from the coast [km]', 'VerticalAlignment', 'bottom', ...
     'HorizontalAlignment', 'left', 'color', 'black', 'Fontsize', 13, ...
     'Fontweight','bold');
    ax = gca;
    ax.FontSize = 12;
    set(ax,'TickDir','out');
    datetick('x','mm/dd HH:MM','keepticks');
    
    cd(dpath);
    saveas(h,[n ' - section' num2str(i) '.jpeg']);
end

end