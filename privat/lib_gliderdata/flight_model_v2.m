function [U,wg,att,Fg,Fb,Fl,Fd] = flight_model_v2(p,rho,pitch,dvbp,tt,Vg,eps,Cd0,mg)

% Steady State Model based on Merckelbach et al., 2009
%
% Inputs : (ajouter unités)
%   p : pressure
%   rho : density
%   pitch : pitch
%   dvbp : ballast pumped
%   tt : temperature
%   Vg : volume glider
%   eps : hull compressibility
%   Cd0 : parasite drag
%   mg : mass of glider
%
% Outputs
%   U : glider velocity through water along the glider path
%   wg : vertical velocity of glider
%   att : attack angle
%   Fg : gravity force
%   Fb : buoyancy force
%   Fl : lift force
%   Fd : drag force
%   ug : horizontal velocity of glider


% Constant parameters
g = 9.81;           % Acceleration of gravity
S = 0.1;            % Wing surface area (convention used in aerodynamics)
aw = 3.7;           % lift-slope coefficient for the wings
ah = 2.4;           % lift-slope coefficient for the hull
Cd1w = 0.78;        % drag coefficient of the wings
Cd1h = 2.1;         % drag coefficient of the hull

% Convert pitch angles from degrees to radians
pitch_rad = pitch .* (pi/180);

% Using  actual values of the optimized parameters 
% to compute the flight model
xVg = Vg*1e-2;
xeps = eps*1e-10;
xCd0 = Cd0*1e-2;
xmg = mg*1e1;

% Resolution of second degree equation to estimate attack angle
% Assignment of solutions according to up/downcasts
a = (Cd1w+Cd1h).*ones(size(pitch_rad));
b = -(aw+ah)*tan(pitch_rad);
c = xCd0.*ones(size(pitch_rad));
delta = b.^2 - 4*(a.*c);
% Set rare negative values (during change of attitude between up and downcast) to zero
delta(delta<0) = 0;    
att1 = (-b + sqrt(delta)) ./ (2*a);
att2 = (-b - sqrt(delta)) ./ (2*a);
att = sign(pitch_rad) .* min(abs(att1),abs(att2));
  
% Vertical forces
Fg = xmg*g;
Fb = g.*(rho).*(xVg.*(1 - xeps.*(p.*10000) + 7.05e-5.*(tt-13.2)) + dvbp./1000000);

% Speeds
U2 = 2*(Fb-Fg).*sin(pitch_rad+att)./((rho).*S.*(xCd0+(Cd1w+Cd1h).*(att.^2)));
% Set rare negative values (during change of attitude between up and downcast) to zero
U2(U2<0) = 0; 
U = sqrt(U2);

wg = U.*sin(pitch_rad+att);
% ug = U.*cos(pitch_rad+att);

% Drag and lift
Fd = 0.5.*rho.*(xCd0+(Cd1w+Cd1h).*(att.^2)).*S.*U.*U;
Fl = 0.5.*rho.*(ah+aw).*att.*S.*U.*U;                   % positive upward
end