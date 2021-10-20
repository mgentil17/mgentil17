function [cost] = cost_v2(param)
% Compute differences between the measurement and the model 
% on vertical velocities
% Inputs
%   Parameters(Vg,eps,Cd0)
%
% Outputs
%   Cost value : sum((Wglider-Wmodel).^2)/ Ndata

global pres dens pitch balla Wglider ind temp mg
[U,Wmodel,att,Fg,Fb,Fl,Fd] = flight_model_v2(pres,dens,pitch,balla,temp,param(1),param(2),param(3),param(4));
cost = nansum((Wglider(ind)-Wmodel).^2)/length(Wmodel);
end