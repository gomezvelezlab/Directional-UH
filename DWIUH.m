%   Coded by Gabriel Perez
%   Repository : Directional-UH
%   Email:   perezmesagj@ornl.gov
%	Last update: 06/07/2023,   MATLAB	2019b  version
%	IF  YOU	PUBLISH  WORK  BENEFITING  FROM  THIS  M-FILE,   PLEASE  CITE  IT AS:
%	Perez. G., et al., (2023)  The Directional Unit Hydrograph Model: 
%   Connecting Streamflow Response to Storm Dynamics, Journal of Hydrology
%   https://doi.org/10.1016/j.jhydrol.2023.130422

%% Direction Width Function Unit Hydrograph

function [Q_conv_storm,wf_conv_storm,wf_chs,t_final_storm,Geo_Rain,d_vector] = DWIUH(Table_river,Vc,Vs,D,th_storm,i_rain,Duration_rain,wf_time_Start,wf_time_End,t0,D_t,x_centroid,y_centroid)

%% Inputs:
% Table_river: Table describing the below HRUs characteristics. 
    % LinkID :      Unique indentifier
    % A_up :        Upstream area [km2]
    % x_c :         x-coordinate centroid of the link [decimal degrees]
    % y_c :         y-coordinate centroid of the link [decimal degrees]
    % A_h :         Hillslope area [km2]
    % L_out :       Length to the basin outlet [km]. This includes the local link length
    % L_h :         Link length [km]
% Vc:            Channel celerity [m/s] 
% Vs:            Storm speed [m/s] 
% D:             Dispersion Coefficient [m2/s]
% th_storm:      Storm direction [degrees]. Zero degrees is East. The angles are measured counter clockwise direction
% i_rain:        Effective Rainfall intensity [mm/hr]
% Duration_rain: Rainfall duration [s]
% wf_time_Start: Vector of the initial travel time to get to the outlet for each link [s]
% wf_end_Start: Vector of the final travel time to get to the outlet for each link [s]
% t0:            Initial time of storm event [s]
% D_t:           Delta time to compute convolution [s]
% x_centroid:    x-coordinate for the watershed centroid [degrees]
% y_centroid:    x-coordinate for the watershed centroid [degrees]

%% Outputs: 
% Q_conv_storm:     Time series of streamflow at the basin outlet [m3/s]
% wf_conv_storm:    Transfer function based on the directional width function and a Deltat
% wf_chs:           Directional Width function based on a D_x (Vc*D_t)
% t_final_storm:    Longest travel time in the system [s]
% Geo_rain:         Spatial coordinates of the rectangular storm at the initial time t0
% d_vector:         Distance related to each element in wf_chs

%% Set rainfall vector for the convolution
P_rain_conv = i_rain.*ones(size((t0:D_t:Duration_rain)))./(1000*60*60); % Hyetogram (Rainfall Depth) from [mm/hr] to [m/s]

%%  Set the initial location of the rectangular moving storm. 
% The storm need to start moving outside of the watershed domain

% Set the distance between the basin centroid and the initial location of
% the rectangular storm
% radius in degrees (1 degree ~ 100km rough approximation)
% Use the max distance basin bound length as proxy
r = max(abs(max(Table_river.x_c)-min(Table_river.x_c)), abs(max(Table_river.y_c)-min(Table_river.y_c)))/2; 
r = r + 0.5; % Add 0.5 degrees to make sure that the storm is outside the basin

% Compute coordinates of the circumscribed circle for plotting
th = 0:pi/50:2*pi; % Angles to plot [rad]
Geo_Rain.xunit = r * cos(th) + x_centroid;
Geo_Rain.yunit = r * sin(th) + y_centroid;

% Compute coordinates of initial storm front location
x_start_storm = r * cos(th_storm) + x_centroid;
y_start_storm = r * sin(th_storm) + y_centroid;

% Storm Width based on the rainfall duration and storm velocity. Conversion factor from m to degrees
Width_Storm = (Vs*Duration_rain)/(1000*100); 

% Storm lenght. This need to be larger than the watershed diameter
Lenght_storm = 2*r; % [degrees] (1 ~ 100km)

% Calculate the points of the "base" of the storm at the initial location
Geo_Rain.Pt1_y = y_start_storm+cos(pi-th_storm)*Lenght_storm/2;
Geo_Rain.Pt1_x = x_start_storm+sin(pi-th_storm)*Lenght_storm/2;
Geo_Rain.Pt2_y = y_start_storm-cos(pi-th_storm)*Lenght_storm/2;
Geo_Rain.Pt2_x = x_start_storm-sin(pi-th_storm)*Lenght_storm/2;

% Calculate slope and intercept of storm base
m_slope = (Geo_Rain.Pt2_y-Geo_Rain.Pt1_y)/(Geo_Rain.Pt2_x-Geo_Rain.Pt1_x);
b_intercept = Geo_Rain.Pt1_y-m_slope*Geo_Rain.Pt1_x;

% Calculate the point of the opposite side of the base
if th_storm<=pi && th_storm>=0
    Geo_Rain.Pt1_y_w = Geo_Rain.Pt1_y+Width_Storm*sin(pi/2-atan(m_slope));
    Geo_Rain.Pt1_x_w = Geo_Rain.Pt1_x-Width_Storm*cos(pi/2-atan(m_slope));
    Geo_Rain.Pt2_y_w = Geo_Rain.Pt2_y+Width_Storm*sin(pi/2-atan(m_slope));
    Geo_Rain.Pt2_x_w = Geo_Rain.Pt2_x-Width_Storm*cos(pi/2-atan(m_slope));
else
    Geo_Rain.Pt1_y_w = Geo_Rain.Pt1_y-Width_Storm*sin(pi/2-atan(m_slope));
    Geo_Rain.Pt1_x_w = Geo_Rain.Pt1_x+Width_Storm*cos(pi/2-atan(m_slope));
    Geo_Rain.Pt2_y_w = Geo_Rain.Pt2_y-Width_Storm*sin(pi/2-atan(m_slope));
    Geo_Rain.Pt2_x_w = Geo_Rain.Pt2_x+Width_Storm*cos(pi/2-atan(m_slope));
end

%% Calculate the distance of each link to the storm front
if isinf(m_slope) % This is special case for 0 and 180 degrees
    wf_storm = 100*abs(Table_river.x_c-Geo_Rain.Pt2_x); % [km] Rough approximation from degrees to km
else
    wf_storm = 100*abs(-m_slope*Table_river.x_c+1*Table_river.y_c-b_intercept)./sqrt(m_slope^2+1^2); % [km] Rough approximation from degrees to km
end

% Rescaled storm distance
wf_storm_time_rescaled = (Vc/Vs).*(1000*wf_storm);  % [m]

% Times
%wf_storm_time = wf_storm.*1000./Vs; % Time for the storm to arrive to each "link" [s]
wf_total_time_Start = wf_time_Start + wf_storm_time_rescaled./Vc; % Total time: Travel time of storm to the link + Flow throught the river network [s]
wf_total_time_End = wf_time_End + wf_storm_time_rescaled./Vc; % [s]
t_final_storm = nanmax(wf_total_time_End); % [s]

% Rescaled distances
wf_total_distance_Start = wf_total_time_Start*Vc; % [m]
wf_total_distance_End = wf_total_time_End*Vc; % [m]
d_final_storm = nanmax(wf_total_distance_End); % [m]
%Area_Delta=Area_h./((wf_total_time_End-wf_total_time_Start)); % Delta Area for each Link segment in D_t

D_x = Vc*D_t; % [m]
wf_chs = zeros(size(0:D_x:round(d_final_storm),2),1);
d_vector=0:D_x:round(d_final_storm);
for i=1:size(wf_total_distance_Start,1) % Number of Links
    X=[wf_total_distance_Start(i) wf_total_distance_End(i)];
    Cte_TTD=1/(wf_total_distance_End(i)-wf_total_distance_Start(i)); % Rectangle Height
    Y=[Cte_TTD Cte_TTD];
    Yi=interp1q(X',Y',d_vector').*(Table_river.A_h(i).*1000^2); Yi(isnan(Yi))=0; 
    wf_chs=Yi+wf_chs;
end
wf_chs = wf_chs./sum(Table_river.A_h.*1000^2); % Normalize []
nansum(wf_chs.*D_x) % Must be equal to 1

if D == 0 % Compute the transfer function. Kinematic case
    wf_conv_storm = wf_chs./sum(wf_chs.*D_t); % Transform  to time and normalize
    nansum(wf_conv_storm.*D_t) % Must be equal to 1
else % Dispersion
    Delta_Fix = 100; % [m] This add a delta distance in the vector to avoid innestability in the equation for wf_conv
    Distance_X = Delta_Fix:D_x:round(d_final_storm)+Delta_Fix; % [m]
    Distance_X = Distance_X';
    wf_conv_storm = zeros(size(0:D_t:round(t_final_storm),2),1);
    Pos = 1;
    for ii = 0:D_t:t_final_storm
        wf_conv_storm(Pos,1) = nansum(D_x.*((Distance_X.*wf_chs)./sqrt(4*pi*D*(ii^3))).*(exp(((Distance_X-Vc.*ii).^2)./(-4*D*ii))));        
        Pos = 1+Pos;
    end
    nansum(wf_conv_storm.*D_t) % Must be equal to 1 
end

%% Compute hydrograph
Q_conv_storm = 1000^2*sum(Table_river.A_h).*conv(wf_conv_storm,P_rain_conv).*D_t; % [m3/s]


end