function [passTable, tleTable] = predict_tle_passes_rio_rancho( ...
    tleFile, startUTC, durationHours, minElevationDeg, downlinkHz)
%PREDICT_TLE_PASSES_RIO_RANCHO Predict TLE passes over Rio Rancho, NM.
%
%   [passTable,tleTable] = predict_tle_passes_rio_rancho(tleFile)
%   reads one or more satellites from a standard 2-line or 3-line TLE file,
%   propagates the orbit with the propagator selected by MATLAB for TLE data
%   (SGP4 for near-Earth objects or SDP4 for deep-space objects), and reports
%   every pass over a representative Rio Rancho ground-station location.
%
%   Optional inputs:
%       startUTC         Prediction start time as datetime. Default: now UTC.
%       durationHours    Prediction duration in hours. Default: 48.
%       minElevationDeg  AOS/LOS elevation mask in degrees. Default: 5.
%       downlinkHz       Carrier frequency for Doppler calculations. Set NaN
%                        to omit Doppler. Default: NaN.
%
%   Outputs:
%       passTable  AOS, TCA, LOS, azimuth, elevation, range, and Doppler data.
%       tleTable   Parsed TLE elements and derived two-body quantities.
%
%   Files written to the current folder:
%       Rio_Rancho_Passes.csv
%       TLE_Element_Summary.csv
%       rio_rancho_pass_table.tex
%       tle_element_table.tex
%       tle_calculation_block.tex
%       best_pass_elevation.png     (when at least one pass exists)
%       best_pass_doppler.png       (when downlinkHz is supplied)
%
%   Required toolbox:
%       Aerospace Toolbox or Satellite Communications Toolbox.
%
%   Example:
%       [passes,elements] = predict_tle_passes_rio_rancho( ...
%           "my_satellite.tle", datetime("now","TimeZone","UTC"), ...
%           48, 5, 400.502e6);
%
%   Notes:
%   1. The coordinates represent the approximate center of Rio Rancho.
%      Replace rrLat, rrLon, and rrAlt_m with the actual antenna location for
%      the most accurate pass times and azimuths.
%   2. The model does not include terrain, buildings, trees, or local RF
%      blockage. The elevation mask is a smooth-horizon approximation.

    arguments
        tleFile (1,1) string
        startUTC (1,1) datetime = datetime("now", "TimeZone", "UTC")
        durationHours (1,1) double {mustBePositive} = 48
        minElevationDeg (1,1) double = 5
        downlinkHz (1,1) double = NaN
    end

    validateattributes(minElevationDeg, {"numeric"}, ...
        {"real","finite",">=",-90,"<=",90}, mfilename, ...
        "minElevationDeg");
    if ~(isnan(downlinkHz) || (isfinite(downlinkHz) && downlinkHz > 0))
        error("downlinkHz must be a positive frequency in hertz or NaN.");
    end
    if ~isfile(tleFile)
        error("TLE file not found: %s", tleFile);
    end

    % Representative Rio Rancho city-center coordinates. Replace these with
    % the exact antenna coordinates when available.
    rrLat_deg = 35.2334;
    rrLon_deg = -106.6645;
    rrAlt_m   = 1610;

    % A one-second sample provides AOS/LOS and peak-elevation values suitable
    % for a student pass-prediction report. Increase for faster coarse scans.
    sampleTime_s = 1;

    % Convert the prediction interval to UTC without changing the instant.
    if isempty(startUTC.TimeZone)
        startUTC.TimeZone = "UTC";
    else
        startUTC.TimeZone = "UTC";
    end
    stopUTC = startUTC + hours(durationHours);

    % Parse the TLE separately so the report can show the element fields and
    % the derived mean-motion, period, and radius calculations.
    tleRecords = readTLEFile(tleFile);
    tleTable = struct2table(tleRecords);
    writetable(tleTable, "TLE_Element_Summary.csv");
    exportTLETableLatex(tleTable, "tle_element_table.tex");

    % Create the scenario and import the TLE. MATLAB selects the appropriate
    % general-perturbations propagator for the TLE orbit.
    sc = satelliteScenario(startUTC, stopUTC, sampleTime_s);
    sat = satellite(sc, tleFile);

    % R2025a and later prefer a visibility mask. The fallback supports older
    % releases that use MinElevationAngle.
    try
        gs = groundStation(sc, rrLat_deg, rrLon_deg, ...
            "Altitude", rrAlt_m, ...
            "Name", "Rio Rancho, NM", ...
            "MaskAzimuthEdges", [0 360], ...
            "MaskElevationAngle", minElevationDeg);
    catch
        gs = groundStation(sc, rrLat_deg, rrLon_deg, ...
            "Altitude", rrAlt_m, ...
            "Name", "Rio Rancho, NM", ...
            "MinElevationAngle", minElevationDeg);
    end

    totalSeconds = round(seconds(stopUTC - startUTC));
    timesUTC = (startUTC + seconds(0:sampleTime_s:totalSeconds)).';

    c_mps = 299792458;
    passRows = struct([]);
    row = 0;

    for k = 1:numel(sat)
        % Azimuth, elevation, and range of the satellite as observed from the
        % ground station. The NED frame makes azimuth clockwise from north.
        [az_deg, el_deg, range_m] = aer(gs, sat(k), timesUTC);
        az_deg = az_deg(:);
        el_deg = el_deg(:);
        range_m = range_m(:);

        % Positive range rate means the satellite is receding. Therefore,
        % the one-way received Doppler shift is negative for positive rate.
        rangeRate_mps = gradient(range_m, sampleTime_s);
        if isnan(downlinkHz)
            doppler_Hz = nan(size(rangeRate_mps));
        else
            doppler_Hz = -(rangeRate_mps ./ c_mps) .* downlinkHz;
        end

        ac = access(sat(k), gs);
        intervals = accessIntervals(ac);

        for p = 1:height(intervals)
            aosUTC = intervals.StartTime(p);
            losUTC = intervals.EndTime(p);

            inPass = find(timesUTC >= aosUTC & timesUTC <= losUTC);
            if isempty(inPass)
                continue;
            end

            [maxEl_deg, localPeak] = max(el_deg(inPass));
            peakIndex = inPass(localPeak);
            tcaUTC = timesUTC(peakIndex);

            aosIndex = nearestTimeIndex(timesUTC, aosUTC);
            losIndex = nearestTimeIndex(timesUTC, losUTC);

            aosLocal = aosUTC;
            tcaLocal = tcaUTC;
            losLocal = losUTC;
            aosLocal.TimeZone = "America/Denver";
            tcaLocal.TimeZone = "America/Denver";
            losLocal.TimeZone = "America/Denver";

            row = row + 1;
            passRows(row).SatelliteIndex = k; %#ok<AGROW>
            passRows(row).Satellite = string(sat(k).Name);
            passRows(row).PassNumber = p;
            passRows(row).AOS_UTC = aosUTC;
            passRows(row).TCA_UTC = tcaUTC;
            passRows(row).LOS_UTC = losUTC;
            passRows(row).AOS_Local = aosLocal;
            passRows(row).TCA_Local = tcaLocal;
            passRows(row).LOS_Local = losLocal;
            passRows(row).Duration_s = seconds(losUTC - aosUTC);
            passRows(row).MaxElevation_deg = maxEl_deg;
            passRows(row).AOSAzimuth_deg = az_deg(aosIndex);
            passRows(row).TCAAzimuth_deg = az_deg(peakIndex);
            passRows(row).LOSAzimuth_deg = az_deg(losIndex);
            passRows(row).AOSRange_km = range_m(aosIndex) / 1000;
            passRows(row).RangeAtTCA_km = range_m(peakIndex) / 1000;
            passRows(row).LOSRange_km = range_m(losIndex) / 1000;
            passRows(row).AOSRangeRate_mps = rangeRate_mps(aosIndex);
            passRows(row).TCARangeRate_mps = rangeRate_mps(peakIndex);
            passRows(row).LOSRangeRate_mps = rangeRate_mps(losIndex);
            passRows(row).AOSDoppler_Hz = doppler_Hz(aosIndex);
            passRows(row).TCADoppler_Hz = doppler_Hz(peakIndex);
            passRows(row).LOSDoppler_Hz = doppler_Hz(losIndex);
        end
    end

    if isempty(passRows)
        passTable = table;
        warning("No passes met the %.1f-degree elevation mask in the requested interval.", ...
            minElevationDeg);
        exportPassTableLatex(passTable, "rio_rancho_pass_table.tex", ...
            minElevationDeg, downlinkHz);
        exportCalculationLatex(tleTable, passTable, ...
            "tle_calculation_block.tex", minElevationDeg, downlinkHz);
        return;
    end

    passTable = struct2table(passRows);
    passTable = sortrows(passTable, "AOS_UTC");

    % Set readable formats for display and CSV output.
    utcFormat = "yyyy-MM-dd HH:mm:ss 'UTC'";
    localFormat = "yyyy-MM-dd HH:mm:ss z";
    passTable.AOS_UTC.Format = utcFormat;
    passTable.TCA_UTC.Format = utcFormat;
    passTable.LOS_UTC.Format = utcFormat;
    passTable.AOS_Local.Format = localFormat;
    passTable.TCA_Local.Format = localFormat;
    passTable.LOS_Local.Format = localFormat;

    writetable(passTable, "Rio_Rancho_Passes.csv");
    exportPassTableLatex(passTable, "rio_rancho_pass_table.tex", ...
        minElevationDeg, downlinkHz);
    exportCalculationLatex(tleTable, passTable, ...
        "tle_calculation_block.tex", minElevationDeg, downlinkHz);

    makeBestPassPlots(passTable, sat, gs, sampleTime_s, downlinkHz, c_mps);

    fprintf("\nGround station: %.4f deg N, %.4f deg E, %.0f m\n", ...
        rrLat_deg, rrLon_deg, rrAlt_m);
    fprintf("Prediction interval: %s through %s\n", ...
        formatDateTime(startUTC, utcFormat), formatDateTime(stopUTC, utcFormat));
    fprintf("Elevation mask: %.1f deg\n", minElevationDeg);
    fprintf("Passes found: %d\n\n", height(passTable));
    disp(passTable(:, ["Satellite","AOS_Local","TCA_Local","LOS_Local", ...
        "MaxElevation_deg","AOSAzimuth_deg","TCAAzimuth_deg", ...
        "LOSAzimuth_deg","RangeAtTCA_km"]));
end

function idx = nearestTimeIndex(times, target)
    [~, idx] = min(abs(seconds(times - target)));
end

function records = readTLEFile(tleFile)
%READTLEFILE Parse standard two-line and three-line element records.

    raw = strip(readlines(tleFile));
    raw(raw == "") = [];

    mu_km3_s2 = 398600.4418;
    earthRadius_km = 6378.137;

    records = struct([]);
    recordIndex = 0;
    pendingName = "";
    i = 1;

    while i <= numel(raw)
        line = raw(i);

        if startsWith(line, "1 ")
            if i == numel(raw) || ~startsWith(raw(i + 1), "2 ")
                error("TLE line 1 at file line %d is not followed by line 2.", i);
            end

            line1 = char(line);
            line2 = char(raw(i + 1));
            if numel(line1) < 61 || numel(line2) < 63
                error("TLE record near file line %d is shorter than expected.", i);
            end

            noradID = str2double(strtrim(line1(3:7)));
            epochYY = str2double(line1(19:20));
            epochDay = str2double(line1(21:32));
            if epochYY < 57
                epochYear = 2000 + epochYY;
            else
                epochYear = 1900 + epochYY;
            end
            epochUTC = datetime(epochYear, 1, 1, 0, 0, 0, ...
                "TimeZone", "UTC") + days(epochDay - 1);

            inclination_deg = str2double(line2(9:16));
            raan_deg = str2double(line2(18:25));
            eccentricity = str2double("0." + string(strtrim(line2(27:33))));
            argPerigee_deg = str2double(line2(35:42));
            meanAnomaly_deg = str2double(line2(44:51));
            meanMotion_rev_day = str2double(line2(53:63));
            bstar = parseTLEExponent(line1(54:61));

            meanMotion_rad_s = meanMotion_rev_day * 2*pi / 86400;
            period_min = (86400 / meanMotion_rev_day) / 60;
            semiMajorAxis_km = (mu_km3_s2 / meanMotion_rad_s^2)^(1/3);
            perigeeRadius_km = semiMajorAxis_km * (1 - eccentricity);
            apogeeRadius_km = semiMajorAxis_km * (1 + eccentricity);
            perigeeAltitude_km = perigeeRadius_km - earthRadius_km;
            apogeeAltitude_km = apogeeRadius_km - earthRadius_km;

            eccentricAnomaly_rad = solveKepler( ...
                deg2rad(meanAnomaly_deg), eccentricity);
            trueAnomaly_rad = 2 * atan2( ...
                sqrt(1 + eccentricity) * sin(eccentricAnomaly_rad/2), ...
                sqrt(1 - eccentricity) * cos(eccentricAnomaly_rad/2));
            trueAnomaly_deg = mod(rad2deg(trueAnomaly_rad), 360);

            if strlength(pendingName) == 0
                satelliteName = "NORAD " + string(noradID);
            else
                satelliteName = regexprep(pendingName, "^0\\s+", "");
            end

            recordIndex = recordIndex + 1;
            records(recordIndex).Satellite = satelliteName; %#ok<AGROW>
            records(recordIndex).NORAD_ID = noradID;
            records(recordIndex).Epoch_UTC = epochUTC;
            records(recordIndex).Inclination_deg = inclination_deg;
            records(recordIndex).RAAN_deg = raan_deg;
            records(recordIndex).Eccentricity = eccentricity;
            records(recordIndex).ArgumentOfPerigee_deg = argPerigee_deg;
            records(recordIndex).MeanAnomaly_deg = meanAnomaly_deg;
            records(recordIndex).EccentricAnomaly_deg = rad2deg(eccentricAnomaly_rad);
            records(recordIndex).TrueAnomalyAtEpoch_deg = trueAnomaly_deg;
            records(recordIndex).MeanMotion_rev_per_day = meanMotion_rev_day;
            records(recordIndex).MeanMotion_rad_per_s = meanMotion_rad_s;
            records(recordIndex).Period_min = period_min;
            records(recordIndex).SemiMajorAxis_km = semiMajorAxis_km;
            records(recordIndex).PerigeeRadius_km = perigeeRadius_km;
            records(recordIndex).ApogeeRadius_km = apogeeRadius_km;
            records(recordIndex).PerigeeAltitude_km = perigeeAltitude_km;
            records(recordIndex).ApogeeAltitude_km = apogeeAltitude_km;
            records(recordIndex).BSTAR = bstar;

            pendingName = "";
            i = i + 2;
        elseif startsWith(line, "2 ")
            error("Orphan TLE line 2 at file line %d.", i);
        else
            pendingName = line;
            i = i + 1;
        end
    end

    if isempty(records)
        error("No valid TLE records were found in %s.", tleFile);
    end
end

function value = parseTLEExponent(field)
%PARSETLEEXPONENT Parse implied-decimal TLE notation such as 29621-4.

    s = char(strtrim(string(field)));
    if isempty(s)
        value = NaN;
        return;
    end

    token = regexp(s, '^([+-]?)(\d{5})([+-])(\d+)$', 'tokens', 'once');
    if isempty(token)
        value = NaN;
        return;
    end

    mantissaSign = 1;
    if strcmp(token{1}, '-')
        mantissaSign = -1;
    end
    exponentSign = 1;
    if strcmp(token{3}, '-')
        exponentSign = -1;
    end

    mantissa = str2double(token{2}) * 1e-5;
    exponent = exponentSign * str2double(token{4});
    value = mantissaSign * mantissa * 10^exponent;
end

function E = solveKepler(M, e)
%SOLVEKEPLER Solve M = E - e sin(E) with Newton iteration.

    M = mod(M, 2*pi);
    if e < 0.8
        E = M;
    else
        E = pi;
    end

    for iteration = 1:50
        correction = (E - e*sin(E) - M) / (1 - e*cos(E));
        E = E - correction;
        if abs(correction) < 1e-13
            return;
        end
    end

    warning("Kepler iteration reached its iteration limit.");
end

function exportTLETableLatex(tleTable, fileName)
    fid = fopen(fileName, "w");
    assert(fid >= 0, "Could not create %s.", fileName);
    cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid, "%% Generated by predict_tle_passes_rio_rancho.m\n");
    fprintf(fid, "\\begin{table}[H]\n");
    fprintf(fid, "\\centering\n");
    fprintf(fid, "\\caption{TLE elements and derived orbital quantities.}\n");
    fprintf(fid, "\\label{tab:tle-elements}\n");
    fprintf(fid, "\\resizebox{\\textwidth}{!}{%%\n");
    fprintf(fid, "\\begin{tabular}{lrrrrrrrr}\n");
    fprintf(fid, "\\toprule\n");
    fprintf(fid, "Satellite & $i$ (deg) & $e$ & $\\Omega$ (deg) & $\\omega$ (deg) & $n$ (rev/day) & $T$ (min) & $a$ (km) & $h_p/h_a$ (km) \\\\\n");
    fprintf(fid, "\\midrule\n");

    for i = 1:height(tleTable)
        fprintf(fid, "%s & %.4f & %.7f & %.4f & %.4f & %.8f & %.3f & %.3f & %.1f/%.1f \\\\\n", ...
            latexEscape(tleTable.Satellite(i)), ...
            tleTable.Inclination_deg(i), ...
            tleTable.Eccentricity(i), ...
            tleTable.RAAN_deg(i), ...
            tleTable.ArgumentOfPerigee_deg(i), ...
            tleTable.MeanMotion_rev_per_day(i), ...
            tleTable.Period_min(i), ...
            tleTable.SemiMajorAxis_km(i), ...
            tleTable.PerigeeAltitude_km(i), ...
            tleTable.ApogeeAltitude_km(i));
    end

    fprintf(fid, "\\bottomrule\n");
    fprintf(fid, "\\end{tabular}}\n");
    fprintf(fid, "\\end{table}\n");
end

function exportPassTableLatex(passTable, fileName, minEl, downlinkHz)
    fid = fopen(fileName, "w");
    assert(fid >= 0, "Could not create %s.", fileName);
    cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid, "%% Generated by predict_tle_passes_rio_rancho.m\n");
    if isempty(passTable)
        fprintf(fid, "No passes exceeded the %.1f$^{\\circ}$ elevation mask during the prediction interval.\n", minEl);
        return;
    end

    fprintf(fid, "\\begin{landscape}\n");
    fprintf(fid, "\\begin{longtable}{llrrlrrr}\n");
    fprintf(fid, "\\caption{Predicted passes over Rio Rancho using a %.1f$^{\\circ}$ elevation mask.}\\label{tab:passes}\\\\\n", minEl);
    fprintf(fid, "\\toprule\n");
    fprintf(fid, "Satellite & AOS (local) & AOS az. & Max el. & TCA (local) & TCA az. & LOS az. & Range at TCA (km) \\\\\n");
    fprintf(fid, "\\midrule\n");
    fprintf(fid, "\\endfirsthead\n");
    fprintf(fid, "\\toprule\n");
    fprintf(fid, "Satellite & AOS (local) & AOS az. & Max el. & TCA (local) & TCA az. & LOS az. & Range at TCA (km) \\\\\n");
    fprintf(fid, "\\midrule\n");
    fprintf(fid, "\\endhead\n");

    for i = 1:height(passTable)
        fprintf(fid, "%s & %s & %.1f$^{\\circ}$ & %.1f$^{\\circ}$ & %s & %.1f$^{\\circ}$ & %.1f$^{\\circ}$ & %.1f \\\\\n", ...
            latexEscape(passTable.Satellite(i)), ...
            latexEscape(formatDateTime(passTable.AOS_Local(i), "yyyy-MM-dd HH:mm:ss z")), ...
            passTable.AOSAzimuth_deg(i), ...
            passTable.MaxElevation_deg(i), ...
            latexEscape(formatDateTime(passTable.TCA_Local(i), "yyyy-MM-dd HH:mm:ss z")), ...
            passTable.TCAAzimuth_deg(i), ...
            passTable.LOSAzimuth_deg(i), ...
            passTable.RangeAtTCA_km(i));
    end

    fprintf(fid, "\\bottomrule\n");
    fprintf(fid, "\\end{longtable}\n");
    if ~isnan(downlinkHz)
        fprintf(fid, "The Doppler columns are retained in the CSV output for a nominal carrier of %.6f MHz.\n", downlinkHz/1e6);
    end
    fprintf(fid, "\\end{landscape}\n");
end

function exportCalculationLatex(tleTable, passTable, fileName, minEl, downlinkHz)
    fid = fopen(fileName, "w");
    assert(fid >= 0, "Could not create %s.", fileName);
    cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

    mu = 398600.4418;
    earthRadius = 6378.137;

    fprintf(fid, "%% Generated by predict_tle_passes_rio_rancho.m\n");
    fprintf(fid, "\\subsection{TLE Element Calculations}\n");

    for i = 1:height(tleTable)
        satName = latexEscape(tleTable.Satellite(i));
        epochText = latexEscape(formatDateTime(tleTable.Epoch_UTC(i), ...
            "yyyy-MM-dd HH:mm:ss.SSS 'UTC'"));

        fprintf(fid, "\\subsubsection{%s}\n", satName);
        fprintf(fid, "The TLE epoch is %s. The principal values read from line 2 are\n", epochText);
        fprintf(fid, "\\[i=%.4f^{\\circ},\\quad \\Omega=%.4f^{\\circ},\\quad e=%.7f,\\quad \\omega=%.4f^{\\circ},\\quad M=%.4f^{\\circ}.\\]\n", ...
            tleTable.Inclination_deg(i), tleTable.RAAN_deg(i), ...
            tleTable.Eccentricity(i), tleTable.ArgumentOfPerigee_deg(i), ...
            tleTable.MeanAnomaly_deg(i));

        fprintf(fid, "\\begin{align}\n");
        fprintf(fid, "n &= %.8f\\left(\\frac{2\\pi}{86400}\\right) = %.10e\\ \\mathrm{rad/s},\\\\\n", ...
            tleTable.MeanMotion_rev_per_day(i), tleTable.MeanMotion_rad_per_s(i));
        fprintf(fid, "T &= \\frac{86400}{%.8f} = %.3f\\ \\mathrm{min},\\\\\n", ...
            tleTable.MeanMotion_rev_per_day(i), tleTable.Period_min(i));
        fprintf(fid, "a &= \\left(\\frac{%.4f}{(%.10e)^2}\\right)^{1/3} = %.3f\\ \\mathrm{km},\\\\\n", ...
            mu, tleTable.MeanMotion_rad_per_s(i), tleTable.SemiMajorAxis_km(i));
        fprintf(fid, "r_p &= a(1-e) = %.3f\\ \\mathrm{km},\\qquad h_p=r_p-%.3f=%.3f\\ \\mathrm{km},\\\\\n", ...
            tleTable.PerigeeRadius_km(i), earthRadius, tleTable.PerigeeAltitude_km(i));
        fprintf(fid, "r_a &= a(1+e) = %.3f\\ \\mathrm{km},\\qquad h_a=r_a-%.3f=%.3f\\ \\mathrm{km}.\n", ...
            tleTable.ApogeeRadius_km(i), earthRadius, tleTable.ApogeeAltitude_km(i));
        fprintf(fid, "\\end{align}\n");

        fprintf(fid, "Solving Kepler's equation, $M=E-e\\sin E$, gives $E=%.4f^{\\circ}$ and\n", ...
            tleTable.EccentricAnomaly_deg(i));
        fprintf(fid, "\\[\\nu=2\\tan^{-1}\\!\\left(\\sqrt{\\frac{1+e}{1-e}}\\tan\\frac{E}{2}\\right)=%.4f^{\\circ}.\\]\n", ...
            tleTable.TrueAnomalyAtEpoch_deg(i));
    end

    fprintf(fid, "\\subsection{Pass Results}\n");
    if isempty(passTable)
        fprintf(fid, "No access interval exceeded the %.1f$^{\\circ}$ elevation mask.\n", minEl);
        return;
    end

    [~, i] = max(passTable.MaxElevation_deg);
    fprintf(fid, "The highest-elevation pass in the prediction interval was selected for the detailed calculation.\n");
    fprintf(fid, "\\subsubsection{%s, Pass %d}\n", ...
        latexEscape(passTable.Satellite(i)), passTable.PassNumber(i));
    fprintf(fid, "\\begin{align}\n");
    fprintf(fid, "t_{\\mathrm{AOS}} &= \\text{%s}, & A_{\\mathrm{AOS}} &= %.2f^{\\circ},\\\\\n", ...
        latexEscape(formatDateTime(passTable.AOS_Local(i), "yyyy-MM-dd HH:mm:ss z")), ...
        passTable.AOSAzimuth_deg(i));
    fprintf(fid, "t_{\\max} &= \\text{%s}, & \\epsilon_{\\max} &= %.2f^{\\circ},\\quad A_{\\max}=%.2f^{\\circ},\\\\\n", ...
        latexEscape(formatDateTime(passTable.TCA_Local(i), "yyyy-MM-dd HH:mm:ss z")), ...
        passTable.MaxElevation_deg(i), passTable.TCAAzimuth_deg(i));
    fprintf(fid, "t_{\\mathrm{LOS}} &= \\text{%s}, & A_{\\mathrm{LOS}} &= %.2f^{\\circ},\\\\\n", ...
        latexEscape(formatDateTime(passTable.LOS_Local(i), "yyyy-MM-dd HH:mm:ss z")), ...
        passTable.LOSAzimuth_deg(i));
    fprintf(fid, "\\Delta t &= %.0f\\ \\mathrm{s}, & \\rho(t_{\\max}) &= %.2f\\ \\mathrm{km}.\n", ...
        passTable.Duration_s(i), passTable.RangeAtTCA_km(i));
    fprintf(fid, "\\end{align}\n");

    if ~isnan(downlinkHz)
        fprintf(fid, "Using $\\Delta f=-(\\dot{\\rho}/c)f_0$ with $f_0=%.6f$ MHz,\n", downlinkHz/1e6);
        fprintf(fid, "\\[\\Delta f_{\\mathrm{AOS}}=%.1f\\ \\mathrm{Hz},\\qquad \\Delta f_{\\max}=%.1f\\ \\mathrm{Hz},\\qquad \\Delta f_{\\mathrm{LOS}}=%.1f\\ \\mathrm{Hz}.\\]\n", ...
            passTable.AOSDoppler_Hz(i), passTable.TCADoppler_Hz(i), ...
            passTable.LOSDoppler_Hz(i));
    end
end

function makeBestPassPlots(passTable, sat, gs, sampleTime_s, downlinkHz, c_mps)
    [~, bestRow] = max(passTable.MaxElevation_deg);
    satIndex = passTable.SatelliteIndex(bestRow);
    aosUTC = passTable.AOS_UTC(bestRow);
    losUTC = passTable.LOS_UTC(bestRow);
    duration_s = round(seconds(losUTC - aosUTC));
    tUTC = (aosUTC + seconds(0:sampleTime_s:duration_s)).';

    [~, el_deg, range_m] = aer(gs, sat(satIndex), tUTC);
    el_deg = el_deg(:);
    range_m = range_m(:);
    tLocal = tUTC;
    tLocal.TimeZone = "America/Denver";

    f = figure("Visible", "off");
    plot(tLocal, el_deg, "LineWidth", 1.2);
    grid on;
    xlabel("Local Time");
    ylabel("Elevation (deg)");
    title("Best Predicted Pass: " + passTable.Satellite(bestRow));
    exportgraphics(f, "best_pass_elevation.png", "Resolution", 300);
    close(f);

    if ~isnan(downlinkHz)
        rangeRate_mps = gradient(range_m, sampleTime_s);
        doppler_Hz = -(rangeRate_mps ./ c_mps) .* downlinkHz;
        f = figure("Visible", "off");
        plot(tLocal, doppler_Hz, "LineWidth", 1.2);
        grid on;
        xlabel("Local Time");
        ylabel("One-Way Doppler Shift (Hz)");
        title("Predicted Doppler: " + passTable.Satellite(bestRow));
        exportgraphics(f, "best_pass_doppler.png", "Resolution", 300);
        close(f);
    end
end

function text = formatDateTime(dt, formatString)
    dt.Format = formatString;
    text = char(string(dt));
end

function escaped = latexEscape(text)
    escaped = string(text);
    bs = string(char(92));
    escaped = replace(escaped, bs, "<<<BACKSLASH>>>");
    escaped = replace(escaped, "{", bs + "{");
    escaped = replace(escaped, "}", bs + "}");
    escaped = replace(escaped, "_", bs + "_");
    escaped = replace(escaped, "%", bs + "%");
    escaped = replace(escaped, "&", bs + "&");
    escaped = replace(escaped, "#", bs + "#");
    escaped = replace(escaped, "$", bs + "$");
    escaped = replace(escaped, "~", bs + "textasciitilde{}");
    escaped = replace(escaped, "^", bs + "textasciicircum{}");
    escaped = replace(escaped, "<<<BACKSLASH>>>", bs + "textbackslash{}");
    escaped = char(escaped);
end
