function [bonferroni_corrected_h, bonferroni_corrected_alpha] = multcomp_bonferroni(p_values, varargin)

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
% Bonferroni-corrected null hypothesis test results. The number of tests is
% determined by the length of the vector of p-values.
%
%
% Dunn, O. J. (1959). Estimation of the medians for dependent variables. 
% Annals of Mathematical Statistics, 30(1), 192-197. doi 10.1214/aoms/1177706374
%
% Dunn, O.J. (1961). Multiple comparisons among means. Journal of the 
% American Statistical Association, 56(293), 52-64. doi 10.1080/01621459.1961.10482090
%
%
% requires:
% - p_values (vector of p-values from the hypothesis tests of interest)
% 
% optional:
% - alpha (uncorrected alpha level for statistical significance, default 0.05)
%
%
% outputs:
% - bonferroni_corrected_h (vector of Bonferroni-corrected hypothesis tests 
% derived from comparing p-values to Bonferroni adjusted critical alpha level. 
% 1 = statistically significant, 0 = not statistically significant)
%
% - bonferroni_corrected_alpha (the corrected critical alpha level)
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



%% Bonferroni correction
n_total_comparisons = length(p_values); % Get the number of comparisons
bonferroni_corrected_alpha = alpha_level / n_total_comparisons; % Calculate bonferroni-corrected alpha

bonferroni_corrected_h = zeros(1, length(p_values)); % preallocate

bonferroni_corrected_h(p_values < bonferroni_corrected_alpha) = 1; % Compare each p-value to the corrected threshold.
