%% Start up ncorr
handles_ncorr = ncorr;

%% Load translation data
load translation_data;

%% Load rotation data
load rotation_data;

%% Set Parameters

% Set ref
handles_ncorr.set_ref(ref);
% Set cur
handles_ncorr.set_cur(cur);
% Set ROI
handles_ncorr.set_roi_ref(roi);