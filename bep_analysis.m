function bep_analysis(dateStrings1, dateStrings2)
tic
% input example: dateStrings1 = ['20230511', '20230518', etc], dateStrings2 = ['2023-05-11', '2023-05-18', etc]
hit_eventPhase_cell = {};
miss_eventPhase_cell = {};
falsealarm_eventPhase_cell = {};
correctrejection_eventPhase_cell = {};
eventPhase_cell = {};
not_eventPhase_cell = {};
all_hit_eventMeanPhases = [];
all_miss_eventMeanPhases = [];
all_correctrejection_eventMeanPhases = [];
all_falsealarm_eventMeanPhases = [];
all_eventMeanPhases = [];
all_not_eventMeanPhases = [];
all_hit_eventMeanPhases_deg = [];
all_miss_eventMeanPhases_deg = [];
all_correctrejection_eventMeanPhases_deg = [];
all_falsealarm_eventMeanPhases_deg = [];
all_eventMeanPhases_deg = [];
all_not_eventMeanPhases_deg = [];
% loading files
dataDir = 'Y:\Projects\Pulv_bodysignal\Magnus_SDT\';
for num = 1:numel(dateStrings1)
    filePattern = [dataDir 'Magcombined' dateStrings2{num} '*.mat'];
    fileList = dir(filePattern);
    blockId = arrayfun(@(x) str2double(x.name(32:33)), fileList);
    Magcombined = cellfun(@(x) load([dataDir filesep x]), {fileList.name}, 'UniformOutput', false);
    
    ecgFilename = [dataDir dateStrings1{num} '_ecg.mat'];
    ecg = load(ecgFilename);
    ecg.out = ecg.out(blockId);
    ecg.Tab_outlier = ecg.Tab_outlier(blockId);
    
    % loop through blocks, compute event times throughout the whole session
    blockEndTime = 0;
%     B = [];
    for blockNum = 2:length(Magcombined)
        
        blockEndTime = blockEndTime + Magcombined{blockNum - 1}.trial(end).TDT_state_onsets_aligned_to_1st_INI(end);
        TDT_state_onsets_aligned_to_1st_INI = {Magcombined{blockNum}.trial.TDT_state_onsets_aligned_to_1st_INI};
        TDT_state_onsets_aligned_to_1st_INI = cellfun(@(x) x + blockEndTime, TDT_state_onsets_aligned_to_1st_INI, 'UniformOutput', false);
        [Magcombined{blockNum}.trial(:).TDT_state_onsets_aligned_to_1st_INI] = TDT_state_onsets_aligned_to_1st_INI{:};
        
        ecg.out(blockNum).Rpeak_t = ecg.out(blockNum).Rpeak_t + blockEndTime;
        ecg.out(blockNum).R2R_t = ecg.out(blockNum).R2R_t + blockEndTime;
        
%         A = cellfun(@transpose, TDT_state_onsets_aligned_to_1st_INI, 'UniformOutput', false);
%         B = [B, A];
%         figure, plot([B{:}])
        
    end
    
    % storing trials of each type in cells and converting to structures and arrays
    for blockNum = 1:length(Magcombined)
        
        state4times = arrayfun(@(x) x.TDT_state_onsets_aligned_to_1st_INI(x.TDT_states == 2), Magcombined{blockNum}.trial, 'UniformOutput', false);
        [Magcombined{blockNum}.trial(:).state4times] = state4times{:};
        
    end
    
    idx_hit_trials = cellfun(@(x) [x.trial.SDT_trial_type] == 1, Magcombined, 'UniformOutput', false);
    idx_miss_trials = cellfun(@(x) [x.trial.SDT_trial_type] == 2, Magcombined, 'UniformOutput', false);
    idx_falsealarm_trials = cellfun(@(x) [x.trial.SDT_trial_type] == 3, Magcombined, 'UniformOutput', false);
    idx_correctrejection_trials = cellfun(@(x) [x.trial.SDT_trial_type] == 4, Magcombined, 'UniformOutput', false);
    
    idx_completed_trials = cellfun(@(x) logical([x.trial.completed]), Magcombined, 'UniformOutput', false);
    idx_noncompleted_trials = cellfun(@not, idx_completed_trials, 'UniformOutput', false);
    
    hit_trials_cell = cellfun(@(x,y) x.trial(y), Magcombined, idx_hit_trials, 'UniformOutput', false);
    miss_trials_cell = cellfun(@(x,y) x.trial(y), Magcombined, idx_miss_trials, 'UniformOutput', false);
    falsealarm_trials_cell = cellfun(@(x,y) x.trial(y), Magcombined, idx_falsealarm_trials, 'UniformOutput', false);
    correctrejection_trials_cell = cellfun(@(x,y) x.trial(y), Magcombined, idx_correctrejection_trials, 'UniformOutput', false);
    
    completed_trials = cellfun(@(x,y) x.trial(y), Magcombined, idx_completed_trials, 'UniformOutput', false);
    noncompleted_trials = cellfun(@(x,y) x.trial(y), Magcombined, idx_noncompleted_trials, 'UniformOutput', false);
    
    completed_trials_array = [completed_trials{:}];
    noncompleted_trials_array = [noncompleted_trials{:}];
    hit_trials = [hit_trials_cell{:}];
    miss_trials = [miss_trials_cell{:}];
    falsealarm_trials = [falsealarm_trials_cell{:}];
    correctrejection_trials = [correctrejection_trials_cell{:}];
    
    completed_event_times = [completed_trials_array.state4times];
    not_event_times = [noncompleted_trials_array.state4times];
    
    hit_event_times = [hit_trials.state4times];
    miss_event_times = [miss_trials.state4times];
    correctrejection_event_times = [correctrejection_trials.state4times];
    falsealarm_event_times = [falsealarm_trials.state4times];
    
	R2R_t = [ecg.out.R2R_t];
    R2R_t = R2R_t([ecg.out.idx_valid_R2R_consec]); % take only consecutive RR-intervals
	R2R_durations = [ecg.out.R2R_valid];
    R2R_durations = R2R_durations([ecg.out.idx_valid_R2R_consec]); % take only consecutive RR-intervals
    cycleStart = R2R_t - R2R_durations;
    cycleEnd = R2R_t;
    cycleDuration = cycleEnd - cycleStart;
    eventTimesCycle = arrayfun(@(x,y) completed_event_times((completed_event_times >= x) & (completed_event_times < y)), cycleStart, cycleEnd, 'UniformOutput', false);
    eventTimesNorm = arrayfun(@(x,y,z) (x{1} - y) / z, eventTimesCycle, cycleStart, cycleDuration, 'UniformOutput', false);
    not_eventTimesCycle = arrayfun(@(x,y) not_event_times((not_event_times >= x) & (not_event_times < y)), cycleStart, cycleEnd, 'UniformOutput', false);
    not_eventTimesNorm = arrayfun(@(x,y,z) (x{1} - y) / z, not_eventTimesCycle, cycleStart, cycleDuration, 'UniformOutput', false);
    hit_eventTimesCycle = arrayfun(@(x,y) hit_event_times((hit_event_times >= x) & (hit_event_times < y)), cycleStart, cycleEnd, 'UniformOutput', false);
    hit_eventTimesNorm = arrayfun(@(x,y,z) (x{1} - y) / z, hit_eventTimesCycle, cycleStart, cycleDuration, 'UniformOutput', false);
    miss_eventTimesCycle = arrayfun(@(x,y) miss_event_times((miss_event_times >= x) & (miss_event_times < y)), cycleStart, cycleEnd, 'UniformOutput', false);
    miss_eventTimesNorm = arrayfun(@(x,y,z) (x{1} - y) / z, miss_eventTimesCycle, cycleStart, cycleDuration, 'UniformOutput', false);
    falsealarm_eventTimesCycle = ...
        arrayfun(@(x,y) falsealarm_event_times((falsealarm_event_times >= x) & (falsealarm_event_times < y)), cycleStart, cycleEnd, 'UniformOutput', false);
    falsealarm_eventTimesNorm = ...
        arrayfun(@(x,y,z) (x{1} - y) / z, falsealarm_eventTimesCycle, cycleStart, cycleDuration, 'UniformOutput', false);
    correctrejection_eventTimesCycle = ...
        arrayfun(@(x,y) correctrejection_event_times((correctrejection_event_times >= x) & (correctrejection_event_times < y)), cycleStart, cycleEnd, 'UniformOutput', false);
    correctrejection_eventTimesNorm = ...
        arrayfun(@(x,y,z) (x{1} - y) / z, correctrejection_eventTimesCycle, cycleStart, cycleDuration, 'UniformOutput', false);
    
    hit_eventPhase = 2*pi*[hit_eventTimesNorm{:}];
    miss_eventPhase = 2*pi*[miss_eventTimesNorm{:}];
    falsealarm_eventPhase = 2*pi*[falsealarm_eventTimesNorm{:}];
    correctrejection_eventPhase = 2*pi*[correctrejection_eventTimesNorm{:}];
    eventPhase = 2*pi*[eventTimesNorm{:}];
    not_eventPhase = 2*pi*[not_eventTimesNorm{:}];
    
    hit_eventPhase = mod(hit_eventPhase, 2*pi);
    miss_eventPhase = mod(miss_eventPhase, 2*pi);
    falsealarm_eventPhase = mod(falsealarm_eventPhase, 2*pi);
    correctrejection_eventPhase = mod(correctrejection_eventPhase, 2*pi);
    eventPhase = mod(eventPhase, 2*pi);
    not_eventPhase = mod(not_eventPhase, 2*pi);
    
    hit_eventPhase_cell{end+1} = hit_eventPhase;
    miss_eventPhase_cell{end+1} = miss_eventPhase;
    falsealarm_eventPhase_cell{end+1} = falsealarm_eventPhase;
    correctrejection_eventPhase_cell{end+1} = correctrejection_eventPhase;
    eventPhase_cell{end+1} = eventPhase;
    not_eventPhase_cell{end+1} = not_eventPhase;
    
    hit_eventMeanPhase_session = circ_mean(hit_eventPhase');
    miss_eventMeanPhase_session = circ_mean(miss_eventPhase');
    falsealarm_eventMeanPhase_session = circ_mean(falsealarm_eventPhase');
    correctrejection_eventMeanPhase_session = circ_mean(correctrejection_eventPhase');
    eventMeanPhase_session = circ_mean(eventPhase');
    not_eventMeanPhase_session = circ_mean(not_eventPhase');
    
    hit_eventMeanPhase_session(hit_eventMeanPhase_session < 0) = hit_eventMeanPhase_session(hit_eventMeanPhase_session < 0) + 2*pi;
    miss_eventMeanPhase_session(miss_eventMeanPhase_session < 0) = miss_eventMeanPhase_session(miss_eventMeanPhase_session < 0) + 2*pi;
    falsealarm_eventMeanPhase_session(falsealarm_eventMeanPhase_session < 0) = falsealarm_eventMeanPhase_session(falsealarm_eventMeanPhase_session < 0) + 2*pi;
    correctrejection_eventMeanPhase_session(correctrejection_eventMeanPhase_session < 0) = correctrejection_eventMeanPhase_session(correctrejection_eventMeanPhase_session < 0) + 2*pi;
    eventMeanPhase_session(eventMeanPhase_session < 0) = eventMeanPhase_session(eventMeanPhase_session < 0) + 2*pi;
    not_eventMeanPhase_session(not_eventMeanPhase_session < 0) = not_eventMeanPhase_session(not_eventMeanPhase_session < 0) + 2*pi;
    
    hit_eventMeanPhase_deg_session = hit_eventMeanPhase_session * (180/pi);
    miss_eventMeanPhase_deg_session = miss_eventMeanPhase_session * (180/pi);
    falsealarm_eventMeanPhase_deg_session = falsealarm_eventMeanPhase_session * (180/pi);
    correctrejection_eventMeanPhase_deg_session = correctrejection_eventMeanPhase_session * (180/pi);
    eventMeanPhase_deg_session = eventMeanPhase_session * (180/pi);
    not_eventMeanPhase_deg_session = not_eventMeanPhase_session * (180/pi);
    
    all_hit_eventMeanPhases = [all_hit_eventMeanPhases; hit_eventMeanPhase_session];
    all_miss_eventMeanPhases = [all_miss_eventMeanPhases; miss_eventMeanPhase_session];
    all_correctrejection_eventMeanPhases = [all_correctrejection_eventMeanPhases; correctrejection_eventMeanPhase_session];
    all_falsealarm_eventMeanPhases = [all_falsealarm_eventMeanPhases; falsealarm_eventMeanPhase_session];
    all_eventMeanPhases = [all_eventMeanPhases; eventMeanPhase_session];
    all_not_eventMeanPhases = [all_not_eventMeanPhases; not_eventMeanPhase_session];
    
    all_hit_eventMeanPhases_deg = [all_hit_eventMeanPhases_deg; hit_eventMeanPhase_deg_session];
    all_miss_eventMeanPhases_deg = [all_miss_eventMeanPhases_deg; miss_eventMeanPhase_deg_session];
    all_correctrejection_eventMeanPhases_deg = [all_correctrejection_eventMeanPhases_deg; correctrejection_eventMeanPhase_deg_session];
    all_falsealarm_eventMeanPhases_deg = [all_falsealarm_eventMeanPhases_deg; falsealarm_eventMeanPhase_deg_session];
    all_eventMeanPhases_deg = [all_eventMeanPhases_deg; eventMeanPhase_deg_session];
    all_not_eventMeanPhases_deg = [all_not_eventMeanPhases_deg; not_eventMeanPhase_deg_session];
end

circ_otest1_sessions = circ_otest(all_hit_eventMeanPhases,1);
circ_otest2_sessions = circ_otest(all_miss_eventMeanPhases,1);
circ_otest3_sessions = circ_otest(all_correctrejection_eventMeanPhases,1);
circ_otest4_sessions = circ_otest(all_falsealarm_eventMeanPhases,1);
circ_otest5_sessions = circ_otest(all_eventMeanPhases,1);
circ_otest6_sessions = circ_otest(all_not_eventMeanPhases,1);

figure;
subplot(2, 3, 1);
if ~isempty(all_hit_eventMeanPhases)
    polar(all_hit_eventMeanPhases, ones(size(all_hit_eventMeanPhases)), 'ro');
    title('Hit');
    hold on
    axis equal;
    hold off
end
subplot(2, 3, 2);
if ~isempty(all_miss_eventMeanPhases)
    polar(all_miss_eventMeanPhases, ones(size(all_miss_eventMeanPhases)), 'ro');
    title('Miss');
    hold on
    axis equal;
    hold off
end

subplot(2, 3, 3);
if ~isempty(all_correctrejection_eventMeanPhases)
    polar(all_correctrejection_eventMeanPhases, ones(size(all_correctrejection_eventMeanPhases)), 'ro');
    title('Correct rejection');
    hold on
    axis equal;
    hold off
end

subplot(2, 3, 4);
if ~isempty(all_falsealarm_eventMeanPhases)
    polar(all_falsealarm_eventMeanPhases, ones(size(all_falsealarm_eventMeanPhases)), 'ro');
    title('False alarms');
    hold on
    axis equal;
    hold off
end

subplot(2, 3, 5);
if ~isempty(all_eventMeanPhases)
    polar(all_eventMeanPhases, ones(size(all_eventMeanPhases)), 'ro');
    title('Completed');
    hold on
    axis equal;
    hold off
end

subplot(2, 3, 6);
if ~isempty(all_not_eventMeanPhases)
    polar(all_not_eventMeanPhases, ones(size(all_not_eventMeanPhases)), 'ro');
    title('Not completed');
    hold on
    axis equal;
    hold off
end

total = sum(all_hit_eventMeanPhases_deg);
mean1 = total / length(all_hit_eventMeanPhases_deg);
disp('Hit Event Mean Phases for all sessions:');
disp(mean1);

total = sum(all_miss_eventMeanPhases_deg);
mean2 = total / length(all_miss_eventMeanPhases_deg);
disp('Miss Event Mean Phases for all sessions:');
disp(mean2);

total = sum(all_correctrejection_eventMeanPhases_deg);
mean3 = total / length(all_correctrejection_eventMeanPhases_deg);
disp('Correct rejections Event Mean Phases for all sessions:');
disp(mean3);

total = sum(all_falsealarm_eventMeanPhases_deg);
mean4 = total / length(all_falsealarm_eventMeanPhases_deg);
disp('False alarm Event Mean Phases for all sessions:');
disp(mean4);

total = sum(all_eventMeanPhases_deg);
mean6 = total / length(all_eventMeanPhases_deg);
disp('All Completed Event Mean Phases for sessions:');
disp(mean6);

total = sum(all_not_eventMeanPhases_deg);
mean7 = total / length(all_not_eventMeanPhases_deg);
disp('All Not Completed Event Mean Phases for sessions:');
disp(mean7);

disp('Hit O-tests for sessions:');
disp(circ_otest1_sessions);

disp('All Miss O-tests for sessions:');
disp(circ_otest2_sessions);

disp('All Correct Rejection O-tests for sessions:');
disp(circ_otest3_sessions);

disp('All False alarms O-tests for sessions:');
disp(circ_otest4_sessions);

disp('All Completed O-tests for sessions:');
disp(circ_otest5_sessions);

disp('All Not Completed O-tests for sessions:');
disp(circ_otest6_sessions);

eventPhase_all = horzcat(eventPhase_cell{:});
not_eventPhase_all = horzcat(not_eventPhase_cell{:});
hit_eventPhase_all = horzcat(hit_eventPhase_cell{:});
miss_eventPhase_all = horzcat(miss_eventPhase_cell{:});
falsealarm_eventPhase_all = horzcat(falsealarm_eventPhase_cell{:});
correctrejection_eventPhase_all = horzcat(correctrejection_eventPhase_cell{:});

hit_eventMeanPhase = circ_mean(hit_eventPhase_all');
hit_mean_length = circ_r(hit_eventPhase_all')*100;
hit_x_end = hit_mean_length * cos(hit_eventMeanPhase);
hit_y_end = hit_mean_length * sin(hit_eventMeanPhase);

miss_eventMeanPhase = circ_mean(miss_eventPhase_all');
miss_mean_length = circ_r(miss_eventPhase_all')*100;
miss_x_end = miss_mean_length * cos(miss_eventMeanPhase);
miss_y_end = miss_mean_length * sin(miss_eventMeanPhase);

falsealarm_eventMeanPhase = circ_mean(falsealarm_eventPhase_all');
falsealarm_mean_length = circ_r(falsealarm_eventPhase_all')*100;
falsealarm_x_end = falsealarm_mean_length * cos(falsealarm_eventMeanPhase);
falsealarm_y_end = falsealarm_mean_length * sin(falsealarm_eventMeanPhase);

correctrejection_eventMeanPhase = circ_mean(correctrejection_eventPhase_all');
correctrejection_mean_length = circ_r(correctrejection_eventPhase_all')*100;
correctrejection_x_end = correctrejection_mean_length * cos(correctrejection_eventMeanPhase);
correctrejection_y_end = correctrejection_mean_length * sin(correctrejection_eventMeanPhase);

eventMeanPhase = circ_mean(eventPhase_all');
mean_length = circ_r(eventPhase_all')*100;
x_end = mean_length * cos(eventMeanPhase);
y_end = mean_length * sin(eventMeanPhase);

not_eventMeanPhase = circ_mean(not_eventPhase_all');
not_mean_length = circ_r(not_eventPhase_all')*100;
not_x_end = not_mean_length * cos(not_eventMeanPhase);
not_y_end = not_mean_length * sin(not_eventMeanPhase);

hit_eventMeanPhase(hit_eventMeanPhase < 0) = hit_eventMeanPhase(hit_eventMeanPhase < 0) + 2*pi;
miss_eventMeanPhase(miss_eventMeanPhase < 0) = miss_eventMeanPhase(miss_eventMeanPhase < 0) + 2*pi;
falsealarm_eventMeanPhase(falsealarm_eventMeanPhase < 0) = falsealarm_eventMeanPhase(falsealarm_eventMeanPhase < 0) + 2*pi;
correctrejection_eventMeanPhase(correctrejection_eventMeanPhase < 0) = correctrejection_eventMeanPhase(correctrejection_eventMeanPhase < 0) + 2*pi;
eventMeanPhase(eventMeanPhase < 0) = eventMeanPhase(eventMeanPhase < 0) + 2*pi;
not_eventMeanPhase(not_eventMeanPhase < 0) = not_eventMeanPhase(not_eventMeanPhase < 0) + 2*pi;

hit_total_count = numel(hit_eventPhase_all);
miss_total_count = numel(miss_eventPhase_all);
falsealarm_total_count = numel(falsealarm_eventPhase_all);
correctrejection_total_count = numel(correctrejection_eventPhase_all);
total_count = numel(eventPhase_all);
not_total_count = numel(not_eventPhase_all);

figure;
subplot(2, 2, 1);
if ~isempty(hit_eventPhase_all)
    ig_rose(hit_eventPhase_all, 20, true);
    title('Hit');
    hold on
    text(-1.5, 1.5, ['Total: ', num2str(hit_total_count)], 'Color', 'black', 'FontWeight', 'bold');
    plot([0, hit_x_end], [0, hit_y_end], 'r', 'LineWidth', 2);
    axis equal;
    hold off;
end
subplot(2, 2, 2);
if ~isempty(miss_eventPhase_all)
    ig_rose(miss_eventPhase_all, 20, true);
    title('Miss');
    hold on;
    plot([0, miss_x_end], [0, miss_y_end], 'r', 'LineWidth', 2);
    text(-1.5, 1.5, ['Total: ', num2str(miss_total_count)], 'Color', 'black', 'FontWeight', 'bold');
    axis equal;
    hold off;
end

subplot(2, 2, 3);
if ~isempty(falsealarm_eventPhase_all)
    ig_rose(falsealarm_eventPhase_all, 20, true);
    title('False Alarm');
    hold on;
    plot([0, falsealarm_x_end], [0, falsealarm_y_end], 'r', 'LineWidth', 2);
    text(-1.5, 1.5, ['Total: ', num2str(falsealarm_total_count)], 'Color', 'black', 'FontWeight', 'bold');
    axis equal;
    hold off;
end

subplot(2, 2, 4);
if ~isempty(correctrejection_eventPhase_all)
    ig_rose(correctrejection_eventPhase_all, 20, true);
    title('Correct Rejection');
    hold on;
    plot([0, correctrejection_x_end], [0, correctrejection_y_end], 'r', 'LineWidth', 2);
    text(-1.5, 1.5, ['Total: ', num2str(correctrejection_total_count)], 'Color', 'black', 'FontWeight', 'bold');
    axis equal;
    hold off;
end
spacing = 1;

figure;
subplot(1, 2, 1);
if ~isempty(eventPhase_all)
    ig_rose(eventPhase_all, 20, true);
    title('completed');
    hold on
    text(-1.5, 1.5, ['Total: ', num2str(total_count)], 'Color', 'black', 'FontWeight', 'bold');
    plot([0, x_end], [0, y_end], 'r', 'LineWidth', 2);
    axis equal;
    hold off;
end
subplot(1, 2, 2);
if ~isempty(not_eventPhase_all)
    ig_rose(not_eventPhase_all, 20, true);
    title('Not completed');
    hold on;
    plot([0, not_x_end], [0, not_y_end], 'r', 'LineWidth', 2);
    text(-1.5, 1.5, ['Total: ', num2str(not_total_count)], 'Color', 'black', 'FontWeight', 'bold');
    axis equal;
    hold off;
end

hit_eventMeanPhase_deg = hit_eventMeanPhase * (180/pi);
miss_eventMeanPhase_deg = miss_eventMeanPhase * (180/pi);
falsealarm_eventMeanPhase_deg = falsealarm_eventMeanPhase * (180/pi);
correctrejection_eventMeanPhase_deg = correctrejection_eventMeanPhase * (180/pi);
eventMeanPhase_deg = eventMeanPhase * (180/pi);
not_eventMeanPhase_deg = not_eventMeanPhase * (180/pi);

disp('Circular Mean Phase Angles (degrees):');
disp(['Completed: ', num2str(eventMeanPhase_deg)]);
disp(['Not completed: ', num2str(not_eventMeanPhase_deg)]);
disp(['Hits: ', num2str(hit_eventMeanPhase_deg)]);
disp(['Misses: ', num2str(miss_eventMeanPhase_deg)]);
disp(['False Alarms: ', num2str(falsealarm_eventMeanPhase_deg)]);
disp(['Correct Rejections: ', num2str(correctrejection_eventMeanPhase_deg)]);
circ_otest1 = circ_otest(hit_eventPhase_all,1);
circ_otest2 = circ_otest(miss_eventPhase_all,1);
circ_otest3 = circ_otest(correctrejection_eventPhase_all,1);
circ_otest4 = circ_otest(falsealarm_eventPhase_all,1);
circ_otest5 = circ_otest(eventPhase_all,1);
circ_otest6 = circ_otest(not_eventPhase_all,1);
disp('O-test p-values:');
disp(['Completed: ', num2str(circ_otest5)]);
disp(['Not completed: ', num2str(circ_otest6)]);
disp(['Hits: ', num2str(circ_otest1)]);
disp(['Misses: ', num2str(circ_otest2)]);
disp(['False Alarms: ', num2str(circ_otest3)]);
disp(['Correct Rejections: ', num2str(circ_otest4)]);

disp('Rayleigh Test p-value:');
[h, p] = circ_rtest(hit_eventPhase_all);
fprintf('Hits: %.4f\n', p);
[h1, p1] = circ_rtest(eventPhase_all);
fprintf('Completed: %.4f\n', h1);
[h2, p2] = circ_rtest(miss_eventPhase_all);
fprintf('Misses: %.4f\n', h2);
[h3, p3] = circ_rtest(not_eventPhase_all);
fprintf('Not completed: %.4f\n', h3);
[h4, p4] = circ_rtest(correctrejection_eventPhase_all);
fprintf('Correct rejections: %.4f\n', h4);
[h5, p5] = circ_rtest(falsealarm_eventPhase_all);
fprintf('False alarms: %.4f\n', h5);
total_trials = numel(eventPhase_all);
fprintf('Number of completed trials: %d\n', total_trials);
pHit = numel(hit_eventPhase_all) / total_trials;
pFA = numel(falsealarm_eventPhase_all) / total_trials;
testsim_dprime(pHit, pFA);
toc


