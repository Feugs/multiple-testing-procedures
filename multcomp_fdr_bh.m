function [fdr_corrected_h, benhoch_critical_alpha] = multcomp_fdr_bh(p_values, varargin)

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
% false discovery rate corrected null hypothesis test results (Benjamin-Hochberg procedure).
% The number of tests is determined by the length of the vector of p-values.
%
% Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: 
% A practical and powerful approach to multiple testing. Journal of the 
% Royal Statistical Society. Series B (Methodological), 57, 289-300. 
% Stable link:http://www.jstor.org/stable/2346101 
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
% - fdr_corrected_h (vector of false discovery rate corrected hypothesis tests 
% derived from comparing p-values to false discovery rate adjusted critical alpha level. 
% 1 = statistically significant, 0 = not statistically significant)
%
% - benhoch_critical_alpha (the adjusted critical alpha for the false
% discovery rate procedure. p-values smaller or equal to this value are
% declared statistically significant. This value is 0 if no tests reached 
% statistical significance).
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



%% Benjamini-Hochberg False Discovery Rate Correction
n_total_comparisons = length(p_values); % Get the number of comparisons
fdr_corrected_h = zeros(1, length(p_values)); % preallocate

sorted_p = sort(p_values); % Sort p-values from smallest to largest

% Find critical k value
for benhoch_step = 1:n_total_comparisons
    if sorted_p(benhoch_step) <= (benhoch_step / n_total_comparisons) * alpha_level
        benhoch_critical_alpha = sorted_p(benhoch_step);
    end
end

% If no steps are significant set critical alpha to zero
if ~exist('benhoch_critical_alpha', 'var')
    benhoch_critical_alpha = 0;
end

% Declare tests significant if they are smaller than or equal to the adjusted critical alpha
fdr_corrected_h(p_values <= benhoch_critical_alpha) = 1;