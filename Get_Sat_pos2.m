% Dual-antenna data processing

% clear all;
close all;
warning('off','all');

%% Add path

% [homedir, mfile, ext] = fileparts(mfilename('fullpath'));
homedir= pwd;
addpath(homedir);
data_label = 'rooftop4_v2';
rinex_dir = sprintf('%s/data_rooftop/%s',homedir, data_label);
addpath(rinex_dir);
source_dir = sprintf('%s/source',homedir);
addpath(genpath(source_dir));
custom_dir = sprintf('%s/custom',homedir);
addpath(genpath(custom_dir));
navutils_dir = sprintf('%s/navutils',homedir);
addpath(genpath(navutils_dir));
output_dir = sprintf('%s/output/%s',homedir);
addpath(genpath(output_dir));

data_name = fullfile(output_dir, data_label);

% rooftop4
Station.name = 'station';
Rover.name = 'rover';

% rooftop5
% Station.name = 'BMHR19440195V_25-05-2021_13-17-56';
% Rover.name = 'BMHR21030714D_25-05-2021_13-18-20';

ID_prn.GPS = [1 32];
ID_prn.GLO = [33 59];
ID_prn.GAL = [60 95];
ID_prn.BDS = [96 158];

%% Station Data

filename = Station.name;    % 1 - Station
obs_file = strcat(filename,'.21O');
gpsNav_file = strcat(filename,'.21N');
% galNav_file = strcat(filename,'.21L');
% gloNav_file = strcat(filename,'.21G');
% bdsNav_file = strcat(filename,'.21F');
% qzsNav_file = strcat(filename,'.21Q');

processing_interval = 1;
prn1st_id = first_prn_id();

[pr1, ph1, pr2, ph2, pr3, ph3, dop1, dop2, dop3, snr1, snr2, snr3, time_zero, time_GPS, time, week, date,...
    pos, interval, antoff, antmod, codeC1, marker] = ...
          load_RINEX_obs(obs_file, prn1st_id, processing_interval);
      
            
[Eph, iono] = load_RINEX_nav(gpsNav_file, prn1st_id,0,2);

nSatTot = prn1st_id.getNumSat();
XS_tot1 = NaN(size(pr1,2), nSatTot,3);
VS_tot1 = NaN(size(pr1,2), nSatTot,3);
lambda = goGNSS.getGNSSWavelengths(Eph, [], nSatTot);

nsat = size(pr1,1);
time_rx = time_GPS + time_zero;
dtR = 0;
err_tropo = zeros(nsat,1);
err_iono  = zeros(nsat,1);
p_rate = 1e-6;

for i = 1:size(pr1,2)
    avail_sat1 = find(pr1(:,i));
    [XS, dtS, XS_tx, VS_tx, time_tx, no_eph, eclipsed, sys_idx] = satellite_positions(time_rx(i), pr1(:,i), avail_sat1, Eph, [], [], err_tropo, err_iono, dtR, 1, 'NONE', lambda, p_rate);
    for j = 1:size(XS_tx,1)
        XS_tot1(i,avail_sat1(j),:) = XS_tx(j,:); % at transmission time
%         XS_tot1(i,avail_sat1(j),:) = XS(j,:); % at reception time
        VS_tot1(i,avail_sat1(j),:) = VS_tx(j,:);
    end
end

[pr1, ph1, pr2, ph2, pr3, ph3, dop1, dop2, dop3, snr1, snr2, snr3, XS_tot1, VS_tot1] = ...
          zero_nan(pr1, ph1, pr2, ph2, pr3, ph3, dop1, dop2, dop3, snr1, snr2, snr3, XS_tot1, VS_tot1);

Station.pr1 = pr1.'; Station.ph1 = ph1.'; Station.pr2 = pr2.'; Station.ph2 = ph2.'; 
Station.pr3 = pr3.'; Station.ph3 = ph3.';
Station.dop1 = dop1.'; Station.dop2 = dop2.'; Station.dop3 = dop3.';
Station.snr1 = snr1.'; Station.snr2 = snr2.'; Station.snr3 = snr3.';
Station.time_zero = time_zero; Station.time_GPS = time_GPS; Station.time = time; Station.week = week; Station.date = date;
Station.approx_pos = pos;
Station.interval = interval;
Station.antoff = antoff;
Station.antmod = antmod;
Station.codeC1 = codeC1;
Station.marker = marker;
Station.SVpos_x = squeeze(XS_tot1(:,:,1)); Station.SVpos_y = squeeze(XS_tot1(:,:,2)); Station.SVpos_z = squeeze(XS_tot1(:,:,3));
Station.SVvel_x = squeeze(VS_tot1(:,:,1)); Station.SVvel_y = squeeze(VS_tot1(:,:,2)); Station.SVvel_z = squeeze(VS_tot1(:,:,3));

clearvars -except Rover Station data_name ID_prn rinex_dir

%% Rover Data


filename1 = Rover.name;    % 1 - Rover
obs_file = strcat(filename1,'.21O');
gpsNav_file = strcat(filename1,'.21N');
galNav_file = strcat(filename1,'.21L');
gloNav_file = strcat(filename1,'.21G');
bdsNav_file = strcat(filename1,'.21F');
qzsNav_file = strcat(filename1,'.21Q');

processing_interval = 1;
prn1st_id = first_prn_id();

[pr1, ph1, pr2, ph2, pr3, ph3, dop1, dop2, dop3, snr1, snr2, snr3, time_zero, time_GPS, time, week, date,...
    pos, interval, antoff, antmod, codeC1, marker] = ...
          load_RINEX_obs(obs_file, prn1st_id, processing_interval);
      
[Eph, iono] = load_RINEX_nav(gpsNav_file, prn1st_id,0,2);

nSatTot = prn1st_id.getNumSat();
XS_tot1 = NaN(size(pr1,2), nSatTot,3);
VS_tot1 = NaN(size(pr1,2), nSatTot,3);
lambda = goGNSS.getGNSSWavelengths(Eph, [], nSatTot);

nsat = size(pr1,1);
time_rx = time_GPS + time_zero;
dtR = 0;
err_tropo = zeros(nsat,1);
err_iono  = zeros(nsat,1);
p_rate = 1e-6;

for i = 1:size(pr1,2)
    avail_sat1 = find(pr1(:,i));
    [XS, dtS, XS_tx, VS_tx, time_tx, no_eph, eclipsed, sys_idx] = satellite_positions(time_rx(i), pr1(:,i), avail_sat1, Eph, [], [], err_tropo, err_iono, dtR, 1, 'NONE', lambda, p_rate);
    for j = 1:size(XS_tx,1)
        XS_tot1(i,avail_sat1(j),:) = XS_tx(j,:); % at transmitted time
%         XS_tot1(i,avail_sat1(j),:) = XS(j,:); % at received time
        VS_tot1(i,avail_sat1(j),:) = VS_tx(j,:); % at transmitted time
    end
end

[pr1, ph1, pr2, ph2, pr3, ph3, dop1, dop2, dop3, snr1, snr2, snr3, XS_tot1, VS_tot1] = ...
          zero_nan(pr1, ph1, pr2, ph2, pr3, ph3, dop1, dop2, dop3, snr1, snr2, snr3, XS_tot1, VS_tot1);
     
Rover.pr1 = pr1.'; Rover.ph1 = ph1.'; Rover.pr2 = pr2.';
Rover.ph2 = ph2.'; Rover.pr3 = pr3.'; Rover.ph3 = ph3.';
Rover.dop1 = dop1.'; Rover.dop2 = dop2.'; Rover.dop3 = dop3.';
Rover.snr1 = snr1.'; Rover.snr2 = snr2.'; Rover.snr3 = snr3.';
Rover.time_zero = time_zero; Rover.time_GPS = time_GPS; Rover.time = time; Rover.week = week; Rover.date = date;
Rover.approx_pos = pos;
Rover.interval = interval;
Rover.antoff = antoff;
Rover.antmod = antmod;
Rover.codeC1 = codeC1;
Rover.marker = marker;
Rover.SVpos_x = squeeze(XS_tot1(:,:,1)); Rover.SVpos_y = squeeze(XS_tot1(:,:,2)); Rover.SVpos_z = squeeze(XS_tot1(:,:,3));
Rover.SVvel_x = squeeze(VS_tot1(:,:,1)); Rover.SVvel_y = squeeze(VS_tot1(:,:,2)); Rover.SVvel_z = squeeze(VS_tot1(:,:,3));

clearvars -except Rover Station Eph iono lambda ID_prn data_name


%% Time sync
[Station, Rover] = timesync(Station, Rover); % With triple data, do this for all cyclic combination

save(data_name);