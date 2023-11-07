%   Coded by Gabriel Perez
%   Repository : Directional-UH
%   Email:   perezmesagj@ornl.gov
%   Last update: 06/07/2023,   MATLAB   2019b  version
%   IF  YOU PUBLISH  WORK  BENEFITING  FROM  THIS  M-FILE,   PLEASE  CITE  IT AS:
%   Perez. G., et al., (2023)  The Directional Unit Hydrograph Model: 
%   Connecting Streamflow Response to Storm Dynamics, Journal of Hydrology
%   https://doi.org/10.1016/j.jhydrol.2023.130422


%% Example to run the Directional-UH

% This code implements the D-UH for the Turkey River basin, a watershed
% with drainage area of 4,385 km2. The gauge USGS 05412500
% is at the watershed outlet.

% The D-UH can be formulated from Hydrological Response Units
% (HRU) at the pixel scale. For simplicity and computational efficiency,
% this example uses HRU as Hillslope-Link units. To this end, a Travel Time
% Distribution (TTD) at the hillslope scale is derived based on the
% hillslope area and reach (link) length.

%% Input Data: 
% A shapefile of the river network with the below attributes
% are the only geomorphological inputs required to implement the D-UH. 
% This can be derived by any GIS package.

% link_id :     Unique indentifier
% up_area :     Upstream area [km2]
% xc_link :     x-coordinate centroid of the link [decimal degrees]
% yc_link :     y-coordinate centroid of the link [decimal degrees]
% Area_H_km2 :  Hillslope area [km2]
% L_out_km :    Length to the basin outlet [km]. This includes the local link length
% L_km : Link   length [km]

%% Load data
Path=pwd;
Basin_Name='Turkey_IFCscale'; 
S = shaperead([Path '\Input_Data\GIS_Data\River_Network_' Basin_Name '.shp']); % Read river network

%% Get the basin shape.
S_basin = shaperead([Path '\Input_Data\GIS_Data\Basin_' Basin_Name '.shp']); % Basin shape

%% Create table with all the data
Table_river = table(cell2mat({S(:).link_id})',cell2mat({S(:).up_area})',cell2mat({S(:).xc_link})', cell2mat({S(:).yc_link})',...
    cell2mat({S(:).Area_H_km2})', cell2mat({S(:).L_out_km})', cell2mat({S(:).L_km})', ...
    'VariableNames',{'LinkID','A_up','x_c','y_c','A_h','L_out','L_h'});

x_centroid = cell2mat({S_basin(:).xc});  % [degrees]
y_centroid = cell2mat({S_basin(:).yc});  % [degrees]

%% Set parameters
RC = 0.2;               % Runoff Coefficient []
i_rain = RC*9.06;       % Effective Rainfall intensity [mm/hr]
Duration_rain=6*60*60;  % Rainfall duration [s]
Vc = 0.5947*(i_rain^0.3717)*(Duration_rain/(60*60))^0.2737;               % Channel celerity [m/s]  
Vh = 0.02;              % Hillslope celerity [m/s]
Vs = 5;                 % Storm speed [m/s]
th_storm=170*pi/180;    % Storm direction [degrees]. Zero degrees is East. The angles are measured counter clockwise direction
D = 1;                  % Dispersion Coefficient [m2/s]
t0 = 0;                 % Initial time [s]
D_t = 60;               % Delta time to compute convolution [s]. Use a D_t less than the travel time of the smallest link
A_total = sum(Table_river.A_h).*(1000^2);     % Basin area [m2]
Q_base = 0.003871*(A_total/(1000^2))^1.04;    % Baseflow [m3/s]. Regression from observations Q_base=Q_mean=aA^b; For Iowa using 92 sites, with an average of 65 years of data
Q_ref = 283;                                  % Flood Stage [m3/s] [1914 to 2020] [https://pubs.usgs.gov/of/2006/1067/pdf/OFR_2006-1067.pdf]

%% Compute the Width Function Unit Hydrograph (spatially uniform rainfall), and extract Width Function
[Q_conv,P_rain_conv,wf_conv,wf_time_Start,wf_time_End,t_final,wf_x,wf_ch,d_vector,t_vector] = WIUH(Table_river,i_rain,Duration_rain,Vc,Vh,D,t0,D_t);
Q_conv = Q_conv + Q_base; % Add baseflow

%% Compute the Directional-WIUH
[Q_conv_Storm,wf_conv_storm,wf_chs,t_final_storm,Geo_Rain,d_vector_storm] = DWIUH(Table_river,Vc,Vs,D,th_storm,i_rain,Duration_rain,wf_time_Start,wf_time_End,t0,D_t,x_centroid,y_centroid);
Q_conv_Storm = Q_conv_Storm + Q_base; % Add base flow optionally

%% Do some nice figures

%% Plot the basin with the respective locations of the Directional Unit Hydrograph
% Get ID outlet
[~, pos_outlet] = min(Table_river.L_out);
% Trim network for plotting
S_plot=S(Table_river.A_up>10);
figure
set(gcf,'color','white')
subplot(3,2,[1 3 5])
mapshow(S_basin); hold on
mapshow(S_plot); 
scatter(Table_river.x_c(pos_outlet),Table_river.y_c(pos_outlet),40,'MarkerFaceColor','r','MarkerEdgeColor','k');
plot(Geo_Rain.xunit, Geo_Rain.yunit,'color','k'); set(gca,'visible','off'); hold on
fill_h=fill([Geo_Rain.Pt1_x Geo_Rain.Pt2_x Geo_Rain.Pt2_x_w Geo_Rain.Pt1_x_w Geo_Rain.Pt1_x],...
    [Geo_Rain.Pt1_y Geo_Rain.Pt2_y Geo_Rain.Pt2_y_w Geo_Rain.Pt1_y_w Geo_Rain.Pt1_y],'b'); set(fill_h,'facealpha',.5); 
subplot(3,2,2)
plot((0:D_t:round(t_final_storm))./(24*60*60),wf_conv_storm, 'linewidth',2); 
xlabel('Time (days)'); ylabel('g_{Directional-UH} [1/s]'); 
xlim([0 (round(t_final_storm)+Duration_rain)]./(24*60*60)); 
set(gca,'fontsize',12); 
subplot(3,2,4)
plot(([t0 t0:D_t:Duration_rain Duration_rain])./(24*60*60),[0 P_rain_conv.*(1000*60*60) 0], 'linewidth',2); 
xlabel('Time (days)'); ylabel('i_e [mm/hr]'); 
xlim([0 (round(t_final_storm)+Duration_rain)]./(24*60*60)); 
set(gca,'fontsize',12); 
subplot(3,2,6)
plot((0:D_t:(round(t_final_storm)+Duration_rain))'./(24*60*60),Q_conv_Storm', 'linewidth',2); 
xlabel('Time (days)'); ylabel('Q [m^3/s]'); 
xlim([0 (round(t_final_storm)+Duration_rain)]./(24*60*60)); 
set(gca,'fontsize',12); 

%% Width Functions
figure; set(gcf,'color','white')
subplot(3,1,1)
plot(d_vector,wf_x);
ylabel('W [1/m]'); xlabel('Distance to outlet [m]');
subplot(3,1,2)
plot(d_vector,wf_ch);
ylabel('W_{c,h} [1/m]'); xlabel('Rescaled distance to outlet [m]');
subplot(3,1,3)
plot(d_vector_storm,wf_chs);
ylabel('W_{c,h,s} [1/m]'); xlabel('Rescaled distance to outlet [m]');

%% Hydrograph from Width Function (spatially uniform rainfall) and the Directional-UH
figure; set(gcf,'color','white')
plot((0:D_t:(round(t_final)+Duration_rain))'./(24*60*60),Q_conv,'b'); hold on
plot((0:D_t:(round(t_final_storm)+Duration_rain))'./(24*60*60),Q_conv_Storm,'k')
legend(['Spatial Uniform Rainfall'; 'Rectangular Moving Storm']);
xlabel('Time [days]'); ylabel('Q [m^3/s]'); 
set(gca,'fontsize',12); 

%% Compute the Directional-UH for "all" Directions From 0 to 360 degrees
Max_Q_Dir=zeros(36,1);
Angles=[0:10:40 45 50:10:130 135 140:10:220 225 230:10:310 315 320:10:350];
tic
parpool(4);
parfor i=1:length(Angles)
    i
    th_storm=Angles(i)*pi/180; % Degrees
    [Q_conv_Storm_T,wf_conv_storm_T,wf_chs_T,t_final_storm_T,Geo_Rain_T,d_vector_storm_T] = DWIUH(Table_river,Vc,Vs,D,th_storm,i_rain,Duration_rain,wf_time_Start,wf_time_End,t0,D_t,x_centroid,y_centroid);
    Max_Q_Dir(i)=nanmax(Q_conv_Storm_T);
    Q_conv_Storm_All{i}=Q_conv_Storm_T;
    wf_conv_storm_All{i}=wf_conv_storm_T;
    t_final_storm_All{i}=t_final_storm_T;
    Geo_Rain_All{i}=Geo_Rain_T;
end
delete(gcp);
toc

%% Plot in polar coordinates the Max_Q for all the directions
figure
set(gcf,'color','white')
subplot(1,2,1)
Angles_to_plot=[6 16 26 36];
Color={'k','b', [84,39,143]./256, [35,139,69]./256};
plot((0:D_t:(round(t_final_storm_All{Angles_to_plot(1)})+Duration_rain))./(24*60*60),Q_conv_Storm_All{Angles_to_plot(1)}'+Q_base,Color{1},'LineWidth',2); hold on
plot((0:D_t:(round(t_final_storm_All{Angles_to_plot(2)})+Duration_rain))./(24*60*60),Q_conv_Storm_All{Angles_to_plot(2)}'+Q_base,Color{2},'LineWidth',2);
plot((0:D_t:(round(t_final_storm_All{Angles_to_plot(3)})+Duration_rain))./(24*60*60),Q_conv_Storm_All{Angles_to_plot(3)}'+Q_base,'color',Color{3},'LineWidth',2);
plot((0:D_t:(round(t_final_storm_All{Angles_to_plot(4)})+Duration_rain))./(24*60*60),Q_conv_Storm_All{Angles_to_plot(4)}'+Q_base,'color',Color{4},'LineWidth',2);
plot((0:D_t:(round(t_final)+Duration_rain))'./(24*60*60),Q_conv,'--k'); hold on
plot([0 max((0:D_t:(round(t_final_storm_All{Angles_to_plot(1)})+Duration_rain))./(24*60*60))],[Q_ref Q_ref],'--r')
legend(['$\theta_s=$' num2str(Angles(Angles_to_plot(1))) '$^o$'],...
    ['$\theta_s=$' num2str(Angles(Angles_to_plot(2))) '$^o$'],...
    ['$\theta_s=$' num2str(Angles(Angles_to_plot(3))) '$^o$'],...
    ['$\theta_s=$' num2str(Angles(Angles_to_plot(4))) '$^o$'],...
    'Spatial Uniform Rainfall', 'Flood Stage','Interpreter','latex')
xlabel('Time [days]'); ylabel('Q [m^3/s]'); set(gca,'fontsize',12);

subplot(1,2,2)
Angles_rad=deg2rad(Angles);
ax1=polarplot([Angles_rad Angles_rad(1)],[Max_Q_Dir; Max_Q_Dir(1)]+Q_base,'LineWidth',3,'Color','k'); title('Peak flow (radius) vs Storm direction (angle)'); hold on
[Maximum_Q, Pos_Maximum_Q] = max(Max_Q_Dir+Q_base); 
[Minimum_Q, Pos_Minimum_Q] = min(Max_Q_Dir+Q_base);
polarplot(Angles_rad(Pos_Maximum_Q),Maximum_Q,'o','color',[0.5 0.5 0.5],'MarkerSize',10,'MarkerFaceColor','r');
polarplot(Angles_rad(Pos_Minimum_Q),Minimum_Q,'o','color',[0.5 0.5 0.5],'MarkerSize',10,'MarkerFaceColor','b');
text(Angles_rad(Pos_Maximum_Q),Maximum_Q+10,['Q_{max}=' num2str(Maximum_Q,'% 10.0f') 'm^3/s']);
text(Angles_rad(Pos_Minimum_Q),Minimum_Q+10,['Q_{min}=' num2str(Minimum_Q,'% 10.0f') 'm^3/s']);
set(gca,'RAxisLocation',120); 
set(gca,'fontsize',12);  %rlim([0 200]);



