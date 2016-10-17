function [corrected_h, corrected_p, critical_t] = multcomp_blaire_karniski_permtest(cond1_data, cond2_data, varargin)

%__________________________________________________________________________
% Multiple comparisons correction function written by Daniel Feuerriegel 21/04/2016 
% to complement DDTBOX scripts written by Stefan Bode 01/03/2013.
%
% The toolbox was written with contributions from:
% Daniel Bennett, Daniel Feuerriegel, Phillip Alday
%
% The author (Stefan Bode) further acknowledges helpful conceptual input/work from: 
% Jutta Stahl, Simon Lilburn, Philip L. Smith, Elaine Corbett, Carsten Murawski, 
% Carsten Bogler, John-Dylan Haynes
%__________________________________________________________________________
%
% This script receives the original data and outputs corrected p-values and
% hypothesis test results based on a maximum statistic permutation test
% (Blaire & Karniski, 1993). The permutation test in this script is based
% on the t-statistic, but could be adapted to use with other statistics
% such as the trimmed mean.
%
% Blair, R. C., & Karniski, W. (1993). An alternative method for 
% significance testing of waveform difference potentials. 
% Psychophysiology, 30, 518-524. DOI: 10.1111/j.1469-8986.1993.tb02075.x
%
% requires:
% - cond1_data (data from condition 1, a subjects x time windows matrix)
% - cond2_data (data from condition 2, a subjects x time windows matrix)
% 
% optional inputs:
% - alpha (uncorrected alpha level for statistical significance, default is 0.05)
% - iterations (number of permutation samples to draw. Default is 5000.
% At least 1000 is recommended for the p = 0.05 alpha level, and at least 5000 is
% recommended for the p = 0.01 alpha level. This is due to extreme events
% at the tails being very rare, needing many random permutations to find
% enough of them).
%
%
% outputs:
% - corrected_h (vector of hypothesis tests in which statistical significance
% is defined by values above a threshold of the (alpha_level * 100)th percentile
% of the maximum statistic distribution.
% 1 = statistically significant, 0 = not statistically significant)
%
% - corrected_p (vector of p-values derived from assessing the t-value of
% each test relative to the distribution of maximum t-values across
% iterations in the permutation test. For example, if above the 99th
% percentile then p < .01.
%
% - critical_t (absolute critical t-value. t-values higher than this are
% counted as statistically significant).
%__________________________________________________________________________
%
% Variable naming convention: STRUCTURE_NAME.example_variable

%% Handling variadic inputs
% Define defaults at the beginning
options = struct(...
    'alpha', 0.05,...
    'iterations', 5000);

% Read the acceptable names
option_names = fieldnames(options);

% Count arguments
n_args = length(varargin);
if round(n_args/2) ~= n_args/2
   error([mfilename ' needs property name/property value pairs'])
end

for pair = reshape(varargin,2,[]) % pair is {propName;propValue}
   inp_name = lower(pair{1}); % make case insensitive

   % Overwrite default options
   if any(strcmp(inp_name, option_names))
      options.(inp_name) = pair{2};
   else
      error('%s is not a recognized parameter name', inp_name)
   end
end
clear pair
clear inp_name

% Renaming variables for use below:
alpha_level = options.alpha;
n_iterations = options.iterations;
clear options;


%% Permutation test

% Checking whether the number of steps of the first and second datasets are equal
if size(cond1_data, 2) ~= size(cond2_data, 2)
   error('Condition 1 and 2 datasets do not contain the same number of comparisons!');
end
if size(cond1_data, 1) ~= size(cond2_data, 1)
   error('Condition 1 and 2 datasets do not contain the same number of subjects!');
end

% Generate difference scores between conditions
diff_scores = cond1_data - cond2_data;

n_subjects = size(diff_scores, 1); % Calculate number of subjects
n_total_comparisons = size(diff_scores, 2); % Calculating the number of comparisons

% Perform t-tests at each step    
[~, ~, ~, extra_stats] = ttest(diff_scores, 0, 'Alpha', alpha_level);
uncorrected_t = extra_stats.tstat; % Vector of t statistics from each test

% Seed the random number generator based on the clock time
rng('shuffle');

% Generate t(max) distribution from the randomly-permuted data
t_max = zeros(1, n_iterations); % Preallocate
t_stat = zeros(n_total_comparisons, n_iterations); % Preallocate
        
    
for iteration = 1:n_iterations
    clear temp; % Clearing out temp variable
    temp_signs = zeros(n_subjects, n_total_comparisons);
    temp = zeros(n_subjects, n_total_comparisons);

    % Draw a random permutation sample for each test
    
    % Randomly switch the sign of difference scores (equivalent to
    % switching labels of conditions)
    for step = 1:n_total_comparisons
        temp_signs(1:n_subjects, step) = (rand(n_subjects, 1) > .5) * 2 - 1; % Switches signs of difference scores
        temp(1:n_subjects, step) = temp_signs(1:n_subjects, step) .* diff_scores(1:n_subjects, step);
    end % of for step
    [~, ~, ~, temp_stats] = ttest(temp, 0, 'Alpha', alpha_level);
    t_stat(:, iteration) = abs(temp_stats.tstat);   

    % Get the maximum t-value within the family of tests and store in a
    % vector. This is to create a null hypothesis distribution.
    t_max(iteration) = max(t_stat(:, iteration));  
end % of for iteration loop

% Calculating the 95th percentile of t_max values (two-tailed, used as decision
% critieria for statistical significance)
critical_t = prctile(t_max(1:n_iterations), ((1 - alpha_level) * 100));

corrected_h = zeros(1,n_total_comparisons); % Preallocate
corrected_p = zeros(1,n_total_comparisons); % Preallocate

% Compare each result with the t-value threshold
corrected_h(abs(uncorrected_t) > critical_t) = 1;

% Calculating a p-value for each step
for step = 1:n_total_comparisons
    corrected_p(step) = mean(t_max(:) >= abs(uncorrected_t(step)));
end % of for step loop