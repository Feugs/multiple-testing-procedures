function [corrected_h] = multcomp_cluster_permtest(cond1_data, cond2_data, varargin)

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
% hypothesis test results based on a maximum cluster statistic permutation test.
% The permutation test in this script is based on the t-statistic, 
% but could be adapted to use with other statistics such as the trimmed mean.
%
% 
% Bullmore, E. T., Suckling, J., Overmeyer, S., Rabe-Hesketh, S., 
% Taylor, E., & Brammer, M. J. (1999). Global, voxel, and cluster tests, 
% by theory and permutation, for a difference between two groups of 
% structural MR images of the brain. IEEE Transactions on Medical Imaging,
% 18, 32-42. doi 10.1109/42.750253
%
% requires:
% - cond1_data (data from condition 1, a subjects x time windows matrix)
% - cond2_data (data from condition 2, a subjects x time windows matrix)
%
% optional:
% - alpha (uncorrected alpha level for statistical significance, default 0.05)
% - iterations (number of permutation samples to draw. At least 1000 is
% recommended for the p = 0.05 alpha level, and at least 5000 is
% recommended for the p = 0.01 alpha level. This is due to extreme events
% at the tails being very rare, needing many random permutations to find
% enough of them).
% - clusteringalpha (the significance threshold used to define
% individual points within a cluster. Setting this to larger values (e.g.
% 0.05) will detect broadly distributed clusters, whereas setting it to
% 0.01 will help detect smaller clusters that exhibit strong effects.
%
% outputs:
% corrected_h (vector of hypothesis tests in which statistical significance
% is defined by values above a threshold of the (alpha_level * 100)th percentile
% of the maximum statistic distribution.
% 1 = statistically significant, 0 = not statistically significant)
%__________________________________________________________________________
%
% Variable naming convention: STRUCTURE_NAME.example_variable

% alpha_level, n_iterations, clustering_alpha

%% Handling variadic inputs
% Define defaults at the beginning
options = struct(...
    'alpha', 0.05,...
    'iterations', 5000,...
    'clusteringalpha', 0.05);

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
clustering_alpha = options.clusteringalpha;
clear options;



%% Cluster-based permutation test

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
uncorrected_h = zeros(1, n_total_comparisons); % preallocate
uncorrected_t = zeros(1, n_total_comparisons); % preallocate

[uncorrected_h, ~, ~, extra_stats] = ttest(diff_scores, 0, 'Alpha', clustering_alpha);
uncorrected_t = extra_stats.tstat; % Vector of t statistics from each test

% Seed the random number generator based on the clock time
rng('shuffle');

% Generate t(max) distribution from the randomly-permuted data

t_stat = zeros(n_total_comparisons, n_iterations); % Preallocate
max_cluster_mass = zeros(1, n_iterations); % Preallocate
cluster_perm_test_h = zeros(n_total_comparisons, n_iterations); % Preallocate
t_sign = zeros(n_total_comparisons, n_iterations); % Preallocate

for iteration = 1:n_iterations
    
    % Draw a random permutation sample for each test
    temp = zeros(n_subjects, n_total_comparisons); % Preallocate
    for step = 1:n_total_comparisons 
        % Randomly switch the sign of difference scores (equivalent to
        % switching labels of conditions)
        temp_signs = (rand(n_subjects, 1) > .5) * 2 - 1; % Switches signs of labels
        temp(:,step) = temp_signs .* diff_scores(1:n_subjects, step);
    end % of for step
    % Run t tests
    [cluster_perm_test_h(:, iteration), ~, ~, temp_stats] = ttest(temp, 0, 'Alpha', clustering_alpha);
    t_stat(:, iteration) = temp_stats.tstat; % Get t statistics
        
        % Marking the sign of each t statistic to avoid clustering pos
        % and neg direction significant results
        for step = 1:n_total_comparisons
            if t_stat(step, iteration) < 0;
                t_sign(step, iteration) = -1; 
            else
                t_sign(step, iteration) = 1; 
            end % of if t_stat
        end % of for step
    % Identify clusters and generate a maximum cluster statistic
    cluster_mass_vector = [0]; % Resets vector of cluster masses
    cluster_counter = 0;

    for step = 1:n_total_comparisons    
        if cluster_perm_test_h(step, iteration) == 1
            if step == 1 % If the first test in the set
                cluster_counter = cluster_counter + 1;
                cluster_mass_vector(cluster_counter) = abs(t_stat(step, iteration));
            else
                % Add to the cluster if there are consecutive
                % statistically significant tests with the same sign.
                % Otherwise, make a new cluster.
                if cluster_perm_test_h(step - 1, iteration) == 1 && t_sign(step - 1, iteration) == t_sign(step, iteration)
                    cluster_mass_vector(cluster_counter) = cluster_mass_vector(cluster_counter) + abs(t_stat(step, iteration));
                else
                    cluster_counter = cluster_counter + 1;
                    cluster_mass_vector(cluster_counter) = abs(t_stat(step, iteration));
                end 
            end % of if test == 1
        end % of if clusterPermTest
    end % of for steps = 1:n_total_steps

    % Find the maximum cluster mass
    max_cluster_mass(iteration) = max(cluster_mass_vector);
end % of iterations loop


% Calculating the 95th percentile of maximum cluster mass values (used as decision
% critieria for statistical significance)
cluster_mass_null_cutoff = prctile(max_cluster_mass, ((1 - alpha_level) * 100));


% Calculate cluster masses in the actual (non-permutation) tests
cluster_mass_vector = [0]; % Resets vector of cluster masses
cluster_counter = 0;
cluster_locations = zeros(1, n_total_comparisons);
cluster_corrected_sig_steps = zeros(1, n_total_comparisons);
clear t_sign;

for step = 1:n_total_comparisons   
    if uncorrected_h(step) == 1
        if step == 1 % If the first test in the set
            cluster_counter = cluster_counter + 1;
            cluster_mass_vector(cluster_counter) = abs(uncorrected_t(step));
            cluster_locations(step) = cluster_counter;
            % Tagging as positive or negative sign effect
            if uncorrected_t < 0
                t_sign(step) = -1;
            else
                t_sign(step) = 1;
            end
        elseif step > 1
            % Tagging as positive or negative sign effect
            if uncorrected_t < 0
                t_sign(step) = -1;
            else
                t_sign(step) = 1;
            end

            % Add to the same cluster only if the previous test was sig.
            % and of the same sign (direction).
            if uncorrected_h(step - 1) == 1 && t_sign(step - 1) == t_sign(step)
                cluster_mass_vector(cluster_counter) = cluster_mass_vector(cluster_counter) + abs(uncorrected_t(step));
                cluster_locations(step) = cluster_counter;
            else
                cluster_counter = cluster_counter + 1;
                cluster_mass_vector(cluster_counter) = abs(uncorrected_t(step));
                cluster_locations(step) = cluster_counter;
            end 
        end % of if step == 1
    end % of if ANALYSIS.RES.h_ttest_uncorrected(na,step) == 1  
end % of for step = 1:n_total_steps

for cluster_no = 1:length(cluster_mass_vector);
    if cluster_mass_vector(cluster_no) > cluster_mass_null_cutoff
        cluster_corrected_sig_steps(cluster_locations == cluster_no) = 1;
    end
end

% Update analysis structure with cluster-corrected significant time
% windows
corrected_h = cluster_corrected_sig_steps;