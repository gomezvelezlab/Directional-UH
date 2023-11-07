%   Coded by Gabriel Perez
%   Repository : Directional-UH
%   Email:   perezmesagj@ornl.gov
%   Last update: 06/07/2023,   MATLAB   2019b  version
%   IF  YOU PUBLISH  WORK  BENEFITING  FROM  THIS  M-FILE,   PLEASE  CITE  IT AS:
%   Perez. G., et al., (2023)  The Directional Unit Hydrograph Model: 
%   Connecting Streamflow Response to Storm Dynamics, Journal of Hydrology
%   https://doi.org/10.1016/j.jhydrol.2023.130422


%% Width Function Unit Hydrograph

function [Q_conv,P_rain_conv,wf_conv,wf_time_Start,wf_time_End,t_final,wf_x,wf_ch,d_vector,t_vector] = WIUH(Table_river,i_rain,Duration_rain,Vc,Vh,D,t0,D_t)

%% Inputs:
% Table_river: Table describing the below HRUs characteristics. 
    % LinkID :      Unique indentifier
    % A_up :        Upstream area [km2]
    % x_c :         x-coordinate centroid of the link [decimal degrees]
    % y_c :         y-coordinate centroid of the link [decimal degrees]
    % A_h :         Hillslope area [km2]
    % L_out :       Length to the basin outlet [km]. This includes the local link length
    % L_h :         Link length [km]
% i_rain:        Effective Rainfall intensity [mm/hr]
% Duration_rain: Rainfall duration [s]
% Vc:            Channel celerity [m/s] 
% Vh:            Hillslope celerity [m/s] 
% D:             Dispersion coefficient [m2/s]
% t0:            Initial time of storm event [s]
% D_t:           Delta time to compute convolution [s].

%% Outputs: 
% Q_conv: Time series of streamflow at the basin outlet [m3/s]
% P_rain_conv: Time series of effective rainfall intensity [m/s]
% wf_conv: Transfer function based on the width function and a Deltat
% wf_time_Start: Vector of the initial travel time to get to the outlet for each link [s]
% wf_end_Start: Vector of the final travel time to get to the outlet for each link [s]
% t_final: Longest travel time in the system [s]
% wf_x:  Width function based on a D_x (Vc*D_t)
% d_vector: Distance related to each element in wf_x
% t_vector: Time related to each element in wf_conv

%% Set rainfall vector for the convolution
P_rain_conv = i_rain.*ones(size((t0:D_t:Duration_rain)))./(1000*60*60); % Hyetogram (Rainfall Depth) from [mm/hr] to [m/s]

% Set the TTD at the Hillslope scale. To this end, 
% an "equivalent hillslope distance" is defined for each hillslope.
% We assume that this distance can be approaximated as 2*{Hillslope Area}/{Link Length}
x_h = 2.*Table_river.A_h./Table_river.L_h; % [km]

% Rescaled distance in the hillslope [km]
if Vh == 0 
    x_h_rescaled = 0; 
else
    x_h_rescaled = (Vc/Vh).*x_h; 
end

% Distance x_ch %[km]
x_ch = Table_river.L_out + x_h; 

% Rescaled Distance x_ch to account the hillslope lenght [km]
x_ch_rescaled = Table_river.L_out + x_h_rescaled; 

% Reference distance to get width function
wf_distance_Start_ref = round((x_ch-x_h)*1000); % [m]
wf_distance_End_ref = round(x_ch*1000); % [m]

% Reference distance and times for convolution
% Compute reference travel times (initial and end) for each link
wf_time_Start = round((x_ch_rescaled*1000./Vc)-(x_h_rescaled*1000./Vc)); % WF from distance to time [s]
wf_time_End = round((x_ch_rescaled*1000./Vc)); % WF from distance to time [s]
% Compute reference distances (initial and end) for each link
wf_distance_Start = wf_time_Start.*Vc; % WF in distance to time (m)
wf_distance_End = wf_time_End.*Vc; % WF in distance to time (m)

% Get the longest distance to set width function vector
d_final = max(wf_distance_End); % [m]
D_x = Vc*D_t; % Delta x [m]. This should less than the smallest length link
d_vector = 0:D_x:round(d_final);
% Get the longest Travel Time to set vector
t_final = max(wf_time_End); % [s]
t_vector = d_vector./Vc; % [s]

%% Compute the Width Function
wf_x = zeros(size(0:D_x:round(d_final),2),1);
for i=1:size(Table_river,1)   
    X = [wf_distance_Start_ref(i) wf_distance_End_ref(i)]; % [s]
    Cte_TTD = 1/(wf_distance_End_ref(i)-wf_distance_Start_ref(i)); % Rectangle Height [1/s]
    Y = [Cte_TTD Cte_TTD]; % [1/s]
    Yi = interp1q(X',Y',d_vector').*(Table_river.A_h(i).*1000^2); % Assign weigths to each hillslope based on hillslope area [m2]
    Yi(isnan(Yi)) = 0; 
    wf_x = Yi+wf_x; 
end

wf_x = wf_x./sum(Table_river.A_h.*1000^2); % Normalize []
nansum(wf_x.*D_x) % Must be equal to 1

% Compute the Rescaled Width Function
wf_ch = zeros(size(0:D_x:round(d_final),2),1);
for i=1:size(Table_river,1)   
    X = [wf_distance_Start(i) wf_distance_End(i)]; % [s]
    Cte_TTD = 1/(wf_distance_End(i)-wf_distance_Start(i)); % Rectangle Height [1/s]
    Y = [Cte_TTD Cte_TTD]; % [1/s]
    Yi = interp1q(X',Y',d_vector').*(Table_river.A_h(i).*1000^2); % Assign weigths to each hillslope based on hillslope area [m2]
    Yi(isnan(Yi)) = 0; 
    wf_ch = Yi+wf_ch; 
end

wf_ch = wf_ch./sum(Table_river.A_h.*1000^2); % Normalize []
nansum(wf_ch.*D_x) % Must be equal to 1

if D == 0 % Compute the transfer function. Kinematic case
    wf_conv = wf_ch./sum(wf_ch.*D_t); % Transform  to time and normalize
    nansum(wf_conv.*D_t) % Must be equal to 1
else % Dispersion
    Delta_Fix = 100; % [m] This add a delta distance in the vector to avoid innestability in the equation for wf_conv
    Distance_X = Delta_Fix:D_x:round(d_final)+Delta_Fix; % [m]
    Distance_X = Distance_X';
    wf_conv = zeros(size(0:D_t:round(t_final),2),1);
    Pos = 1;
    for ii = 0:D_t:t_final
        wf_conv(Pos,1) = nansum(D_x.*((Distance_X.*wf_ch)./sqrt(4*pi*D*(ii^3))).*(exp(((Distance_X-Vc.*ii).^2)./(-4*D*ii))));        
        Pos = 1+Pos;
    end
    nansum(wf_conv.*D_t) % Must be equal to 1 
end

%% Compute hydrograph
Q_conv = 1000^2*sum(Table_river.A_h).*conv(wf_conv,P_rain_conv).*D_t; % [m3/s]

%% Compute mass balance error
MassP = (i_rain/(1000*60*60))*Duration_rain*sum(Table_river.A_h)*1000^2; % Rainfall volume [m3]
MassQ = sum(Q_conv.*D_t);
Error_Mass = 100*(MassQ-MassP)/MassP;
disp(['Mass balance error [%] = ' num2str(Error_Mass)])

if Error_Mass > 1
    disp('Decrease D_t to reduce mass balance error')
end