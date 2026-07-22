%% ECE 535 Meteor-M2-4 pass-prediction example
% This script predicts example passes over Rio Rancho, New Mexico.

clear;
clc;

% TLE downloaded for the report example. Replace this file with a newer TLE
% when predicting an actual future receiving pass.
tleFile = "meteor_m2_4.tle";

% Fixed interval used for the numerical example in the report.
startUTC = datetime(2026, 7, 21, 0, 0, 0, "TimeZone", "UTC");
durationHours = 48;

% AOS and LOS are defined when the satellite crosses 5 degrees elevation.
minElevationDeg = 5;

% Meteor LRPT downlink used for the example Doppler calculation.
downlinkHz = 137.900e6;

[passes, tleElements] = predict_tle_passes_rio_rancho( ...
    tleFile, startUTC, durationHours, minElevationDeg, downlinkHz);

disp("Parsed TLE elements and derived values:");
disp(tleElements);

disp("Predicted Rio Rancho passes:");
disp(passes);
