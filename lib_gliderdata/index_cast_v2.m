function [cast] = index_cast_v2(depth,time,threshold,filt_order)
    %----------------------------------------------------------------------
    % Identifies upcasts, downcasts and surfcasts made by the glider
    %
    % Inputs
    %   - depth parameter
    %   - time
    %   - threshold (empirically defined according to depth data)
    %   - filter order (order of median)
    %
    % Outputs
    %   - cast
    %       *flag 1 : downcasts
    %       *flag 2 : upcasts
    %       *flag 3 : surfacings (surfcasts)
    %----------------------------------------------------------------------
    
    % Compute derivative depth 
    diff_depth = diff(depth);
% %     figure; plot(diff_depth,'-r');                                                              
% %     title('Derivative depth'); ylabel('spped of glider (m/s)');
    
    % Time
    time = time(2:end);                                                    % same size as diff depth
    
    % Smooth depth derive with median filter 
    diff_depth_smooth = medfilt1(diff_depth,filt_order,'omitnan');                     
% %     figure; plot(diff_depth_smooth,'o:r');
% %     title('Smooth derivative depth'); ylabel('speed of glider (m/s)');
    
    % Apply threshold on descent/ascent glider speed
    % to discriminate downcasts/upcasts/surfcasts.
    % Work on duplicated depth vector    
    depth_down = diff_depth;
    depth_down(diff_depth_smooth<threshold) = NaN;                  % downcasts
    
    depth_up = diff_depth;
    depth_up(diff_depth_smooth>-threshold) = NaN;                   % upcasts    
    
    depth_surf = diff_depth;
    depth_surf(diff_depth_smooth<=-threshold | diff_depth_smooth>=threshold) = NaN;     % surfcasts     
    
% %     % Plot filtered depths 
% %     % to check downcast/upcasts/surfacasts extraction
% %     figure; plot(-depth_up,'b.'); hold on; plot(-depth_down,'r-o'); plot(-depth_surf,'gx');
% %     ylabel('depth (m)'); title('Indexing of glider depths : upcasts (blue) / downcasts (red) / surfcasts (green)');
    
    % Flag of indices up/down/surfcasts 
    ind_up = find(~isnan(depth_up));
    ind_down = find(~isnan(depth_down));
    ind_surf = find(~isnan(depth_surf));
    
    % Plot downcasts/upcasts/surfcasts
    figure; plot(time(ind_down),-depth(ind_down),'or'); hold on; 
    plot(time(ind_up),-depth(ind_up),'ob'); plot(time(ind_surf),-depth(ind_surf),'+g');
    ylabel('depth (m)'); title( 'depth : downcasts (red) / upcasts (blue) / surfcasts (green)');
    
    % Flag up/down/surfcasts 
    flag_downcast(1:length(ind_down)) = 1;
    flag_upcast(1:length(ind_up)) = 2;
    flag_surfcast(1:length(ind_surf)) = 3;
    
    % Concatenate vectors
    casts(:,1) = vertcat(ind_down,ind_up,ind_surf);
    casts(:,2) = vertcat(flag_downcast',flag_upcast',flag_surfcast');
    casts= sortrows(casts,1);
    
    % Outputs
    cast.data = casts(:,2);
    cast.var = {'cast'};

end