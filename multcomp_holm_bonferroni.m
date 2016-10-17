function [holm_corrected_h, holm_corrected_alpha] = multcomp_holm_bonferroni(p_values, varargin)

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
% This script receives a vector of p-values and outputs
% Holm-Bonferroni corrected null hypothesis test results. The number of tests is
% determined by the length of the vector of p-values.
%
% Holm, S. (1979). A simple sequentially rejective multiple test procedure. 
% Scandinavian Journal of Statistics 6 (2): 65?70.
%
% requires:
% - p_values (vector of p-values from the hypothesis tests of interest)
%
% optional:
% - alpha (uncorrected alpha level for statistical significance, default 0.05)
%
%
% outputs:
% - holm_corrected_h (vector of Holm-Bonferroni corrected hypothesis tests 
% derived from comparing p-values to Holm-Bonferroni adjusted critical alpha level. 
% 1 = statistically significant, 0 = not statistically significant)
%
% - holm_corrected_alpha (the adjusted alpha threshold)
%__________________________________________________________________________
%
% Variable naming convention: STRUCTURE_NAME.example_variable

%% Handling variadic inputs
% Define defaults at the beginning
options = struct(...
    'alpha', 0.05);

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
clear options;


%% Holm-Bonferroni Correction
n_total_comparisons = length(p_values); % Get the number of comparisons
holm_corrected_h = zeros(1, length(p_values)); % preallocate

sorted_p = sort(p_values); % Sort p-values from smallest to largest
found_crit_alpha = 0; % Reset to signify that we have not found the Holm-Bonferroni corrected critical alpha level

for holm_step = 1:n_total_comparisons
   % Iteratively look for the critical alpha level
   if sorted_p(holm_step) > alpha_level / (n_total_comparisons + 1 - holm_step) && found_crit_alpha == 0
       holm_corrected_alpha = sorted_p(holm_step);
       found_crit_alpha = 1;
   end  
end

if ~exist('holm_corrected_alpha', 'var') % If all null hypotheses are rejected
    holm_corrected_alpha = alpha_level;
end

holm_corrected_h(p_values < holm_corrected_alpha) = 1; % Compare each p-value to the corrected threshold.