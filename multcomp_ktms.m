function [corrected_h] = multcomp_ktms(cond1_data, cond2_data, varargin)

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
% hypothesis test results based on control of the generalised family-wise error rate
% (Korn, 2004, method A in appendices). The permutation test in this script is based
% on the t-statistic, but could be adapted to use with other statistics
% such as the trimmed mean.
%
% Korn, E. L., Troendle, J. F., McShane, L. M., & Simon, R. (2004). Controlling
% the number of false discoveries: Application to high-dimensional genomic data.
% Journal of Statistical Planning and Inference, 124, 379-398. 
% doi 10.1016/S0378-3758(03)00211-8
%
%
% requires:
% - cond1_data (data from condition 1, a subjects x time windows matrix)
% - cond2_data (data from condition 2, a subjects x time windows matrix)
%
% optional:
% - alpha (uncorrected alpha level for statistical significance, default 0.05)
% - iterations (number of permutation samples to draw, default 5000. 
% At least 1000 is recommended for the p = 0.05 alpha level, and at least 5000 is
% recommended for the p = 0.01 alpha level. This is due to extreme events
% at the tails being very rare, needing many random permutations to find
% enough of them).
% - ktms_u (the u parameter of the procedure, or the number of hypotheses
% to automatically reject. Allowing for more false discoveries improves the
% sensitivity of the method. Default is 1).
%
%
% outputs:
% corrected_h (vector of hypothesis tests in which statistical significance
% is defined by values above a threshold of the (alpha_level * 100)th percentile
% of the maximum statistic distribution.
% 1 = statistically significant, 0 = not statistically significant)
%__________________________________________________________________________
%
% Variable naming convention: STRUCTURE_NAME.example_variable

%% Handling variadic inputs
% Define defaults at the beginning
options = struct(...
    'alpha', 0.05,...
    'iterations', 5000,...
    'ktms_u', 1);

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
ktms_u = options.ktms_u;
clear options;



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

p_values = zeros(1, n_total_comparisons); % Preallocate
uncorrected_t = zeros(1, n_total_comparisons); % Preallocate

% Perform t-tests at each step
[~, p_values, ~, extra_stats] = ttest(diff_scores, 0, 'Alpha', alpha_level);
uncorrected_t = extra_stats.tstat; % Vector of t statistics from each test


% Make a vector to denote statistically significant steps
ktms_sig_effect_locations = zeros(1, n_total_comparisons);

sorted_p = sort(p_values); % Sort p-values from smallest to largest

% Automatically reject the u smallest hypotheses (u is set by user as ktms_u variable).
if ktms_u > 0
    ktms_auto_reject_threshold = sorted_p(ktms_u);
elseif ktms_u == 0 % if ktms_u is set to zero
    ktms_auto_reject_threshold = 0;
end

ktms_sig_effect_locations(p_values <= ktms_auto_reject_threshold) = 1; % Mark tests with u smallest p-values as statistically significant.

% Run strong FWER control permutation test but use u + 1th most extreme
% test statistic.
ktms_t_max = zeros(1, n_iterations);
t_stat = zeros(n_total_comparisons, n_iterations);
temp_signs = zeros(n_subjects, n_total_comparisons);

for iteration = 1:n_iterations

    % Draw a random sample for each test
        % Randomly switch the sign of difference scores (equivalent to
        % switching labels of conditions)
        temp = zeros(n_subjects, n_total_comparisons);
    for step = 1:n_total_comparisons
        temp_signs(1:n_subjects, step) = (rand(n_subjects, 1) > .5) * 2 - 1; % Switches signs of labels
        temp(:, step) = temp_signs(1:n_subjects, step) .* diff_scores(1:n_subjects, step);
    end % of for step
        [~, ~, ~, temp_stats] = ttest(temp, 0, 'Alpha', alpha_level);
        t_stat(:, iteration) = abs(temp_stats.tstat);

    % Get the maximum t-value within the family of tests and store in a
    % vector. This is to create a null hypothesis distribution.
    t_sorted = sort(t_stat(:, iteration), 'descend');
    ktms_t_max(iteration) = t_sorted(ktms_u + 1);
end % of for iteration

% Calculating the 95th percentile of t_max values (used as decision
% critieria for statistical significance)
ktms_Null_Cutoff = prctile(ktms_t_max, ((1 - alpha_level) * 100));

% Checking whether each test statistic is above the specified threshold:
for step = 1:n_total_comparisons
    if abs(uncorrected_t(step)) > ktms_Null_Cutoff;
        ktms_sig_effect_locations(step) = 1;
    end
end % of for step

% Marking statistically significant tests
corrected_h = ktms_sig_effect_locations;    