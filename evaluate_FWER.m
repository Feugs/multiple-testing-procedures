% evaluate_FWER.m
%
% Estimates the FWER of several multiple comparisons correction functions, to
% assess whether these functions and methods are controlling the FWER in
% the expected way.
%
% All tests are evaluating paired-samples effects, as implemented in DDTBox.
% However, several MC correction methods also apply to between-groups
% comparisons.
%
% Note: False Discovery Rate (FDR) corresponds to the expected(average) proportion of all
% discoveries (rejected null hypotheses) that are false positives. For
% example, with a FDR of 0.05 this means that 5% of all rejected null
% hypotheses are tolerated to be false positives. (FDR / nFalseDiscoveries / nTotalDiscoveries)
%
% Written by DF 10/16


% Houskeeping
clear all;
close all;

% Seed random number generator based on computer clock
rng('shuffle');

% Simulation Settings
Settings.sampleSizesToUse = [50]; % Sample size per test
Settings.nTestsToUse = [11:20]; % Number of tests
Settings.nIterations = 1000; % Number of iterations
Settings.alphaLevel = 0.05; % Nominal alpha level

% Resampling method settings
Settings.blaireKarniskiIterations = 1000; % Number of permutation samples to draw for the Blaire-Karniski correction
Settings.clusterIterations = 1000; % Number of permutation samples to draw for the maximum cluster mass null distribution
Settings.clusteringAlphaLevel = 0.05; % Alpha level for detecting individual tests to include within a cluster
Settings.ktmsIterations = 1000; % Number of iterations for ktms GFWER control procedure
Settings.ktms_u = 0; % u parameter for ktms GFWER control procedure

% Preallocate FWER false positive count matrices
FWER.FalsePositives.uncorrectedAllTests = zeros(Settings.nIterations, length(Settings.sampleSizesToUse), length(Settings.nTestsToUse));
FWER.FalsePositives.bonferroniAllTests = zeros(Settings.nIterations, length(Settings.sampleSizesToUse), length(Settings.nTestsToUse));
FWER.FalsePositives.holmAllTests = zeros(Settings.nIterations, length(Settings.sampleSizesToUse), length(Settings.nTestsToUse));
FWER.FalsePositives.bhAllTests = zeros(Settings.nIterations, length(Settings.sampleSizesToUse), length(Settings.nTestsToUse));
FWER.FalsePositives.bkyAllTests = zeros(Settings.nIterations, length(Settings.sampleSizesToUse), length(Settings.nTestsToUse));
FWER.FalsePositives.byAllTests = zeros(Settings.nIterations, length(Settings.sampleSizesToUse), length(Settings.nTestsToUse));
FWER.FalsePositives.blaireKarniskiAllTests = zeros(Settings.nIterations, length(Settings.sampleSizesToUse), length(Settings.nTestsToUse));
FWER.FalsePositives.clusterAllTests = zeros(Settings.nIterations, length(Settings.sampleSizesToUse), length(Settings.nTestsToUse));
FWER.FalsePositives.ktmsAllTests = zeros(Settings.nIterations, length(Settings.sampleSizesToUse), length(Settings.nTestsToUse));

% Generate random samples for multiple hypothesis tests and record FWER for
% each method:
for nTests = 1:length(Settings.nTestsToUse)    
    for sampleSize = 1:length(Settings.sampleSizesToUse)
        for i = 1:Settings.nIterations
            fprintf('Running for %i tests with sample size %i iteration %i \n', Settings.nTestsToUse(nTests), Settings.sampleSizesToUse(sampleSize), i);
            
            % Generate a random samples from a normal distribution and perform
            % paired-samples t tests
            clear tempSample1; clear tempSample2; % Clear out temporary samples from previous iteration
            for j = 1:Settings.nTestsToUse(nTests)
                tempSample1(:,j) = randn(Settings.sampleSizesToUse(sampleSize), 1);
                tempSample2(:,j) = randn(Settings.sampleSizesToUse(sampleSize), 1);
            end

            [temp_h, temp_p] = ttest(tempSample1, tempSample2, 'Alpha', Settings.alphaLevel); 

            
            % Blaire-Karniski Maximum Statistic Permutation-Based
            % Correction
            [blaireKarniski_corrected_h, bkp_corrected_p, bkp_critical_t(i, sampleSize, nTests)] = multcomp_blaire_karniski_permtest(tempSample1, tempSample2, 'alpha', Settings.alphaLevel, 'iterations', Settings.blaireKarniskiIterations);
            
            % Cluster-based correction
            [cluster_corrected_h] = multcomp_cluster_permtest(tempSample1, tempSample2, 'alpha', Settings.alphaLevel, 'iterations', Settings.clusterIterations, 'clusteringalpha', Settings.clusteringAlphaLevel);

            % Generalised FWER control procedure (KTMS)
            [ktms_h] = multcomp_ktms(tempSample1, tempSample2, 'alpha', Settings.alphaLevel, 'iterations', Settings.ktmsIterations, 'ktms_u', Settings.ktms_u);
            
            % Bonferroni correction
            [bonferroni_corrected_h, bonferroni_corrected_alpha(i, sampleSize, nTests)] = multcomp_bonferroni(temp_p, 'alpha', Settings.alphaLevel);

            % Holm-Bonferroni correction
            [holm_corrected_h, holm_corrected_alpha(i, sampleSize, nTests)] = multcomp_holm_bonferroni(temp_p, 'alpha', Settings.alphaLevel);

            % Benjamini-Hochberg FDR control procedure
            [fdr_bh_corrected_h, benhoch_critical_alpha(i, sampleSize, nTests)] = multcomp_fdr_bh(temp_p, 'alpha', Settings.alphaLevel);

            % Benjamini-Krieger-Yekutieli FDR control procedure
            [fdr_bky_corrected_h, bky_stage2_critical_alpha(i, sampleSize, nTests)] = multcomp_fdr_bky(temp_p, 'alpha', Settings.alphaLevel);

            % Benjamini-Yekutieli FDR control procedure
            [fdr_by_corrected_h, benyek_critical_alpha(i, sampleSize, nTests)] = multcomp_fdr_by(temp_p, 'alpha', Settings.alphaLevel);
            

            % Calculate the proportion of false positives using each method
            FWER.FalsePositives.uncorrected(i, sampleSize, nTests) = sum(temp_h) / Settings.nTestsToUse(nTests);
            FWER.FalsePositives.bonferroni(i, sampleSize, nTests) = sum(bonferroni_corrected_h) / Settings.nTestsToUse(nTests);
            FWER.FalsePositives.holm(i, sampleSize, nTests) = sum(holm_corrected_h) / Settings.nTestsToUse(nTests);
            FWER.FalsePositives.bh(i, sampleSize, nTests) = sum(fdr_bh_corrected_h) / Settings.nTestsToUse(nTests);
            FWER.FalsePositives.bky(i, sampleSize, nTests) = sum(fdr_bky_corrected_h) / Settings.nTestsToUse(nTests);
            FWER.FalsePositives.by(i, sampleSize, nTests) = sum(fdr_by_corrected_h) / Settings.nTestsToUse(nTests);
            FWER.FalsePositives.blaireKarniski(i, sampleSize, nTests) = sum(blaireKarniski_corrected_h) / Settings.nTestsToUse(nTests);
            FWER.FalsePositives.cluster(i, sampleSize, nTests) = sum(cluster_corrected_h) / Settings.nTestsToUse(nTests);
            FWER.FalsePositives.ktms(i, sampleSize, nTests) = sum(ktms_h) / Settings.nTestsToUse(nTests);


        end % of for i = 1:Settings.nIterations loop

        FWER.FalsePositives.uncorrectedAllTests(FWER.FalsePositives.uncorrected > 0) = 1;
        FWER.FalsePositives.bonferroniAllTests(FWER.FalsePositives.bonferroni > 0) = 1;
        FWER.FalsePositives.holmAllTests(FWER.FalsePositives.holm > 0) = 1;
        FWER.FalsePositives.bhAllTests(FWER.FalsePositives.bh > 0) = 1;
        FWER.FalsePositives.bkyAllTests(FWER.FalsePositives.bky > 0) = 1;
        FWER.FalsePositives.byAllTests(FWER.FalsePositives.by > 0) = 1;
        FWER.FalsePositives.blaireKarniskiAllTests(FWER.FalsePositives.blaireKarniski > 0) = 1;
        FWER.FalsePositives.clusterAllTests(FWER.FalsePositives.cluster > 0) = 1;
        FWER.FalsePositives.ktmsAllTests(FWER.FalsePositives.ktms > 0) = 1; % No. of false positives needs to be above ktms_u parameter for GFWER

        % Calculate the FWER
        FWER.uncorrected(sampleSize, nTests) = sum(FWER.FalsePositives.uncorrectedAllTests(:, sampleSize, nTests)) / Settings.nIterations;
        FWER.bonferroni(sampleSize, nTests) = sum(FWER.FalsePositives.bonferroniAllTests(:, sampleSize, nTests)) / Settings.nIterations;
        FWER.holm(sampleSize, nTests) = sum(FWER.FalsePositives.holmAllTests(:, sampleSize, nTests)) / Settings.nIterations;
        FWER.bh(sampleSize, nTests) = sum(FWER.FalsePositives.bhAllTests(:, sampleSize, nTests)) / Settings.nIterations;
        FWER.bky(sampleSize, nTests) = sum(FWER.FalsePositives.bkyAllTests(:, sampleSize, nTests)) / Settings.nIterations;
        FWER.by(sampleSize, nTests) = sum(FWER.FalsePositives.byAllTests(:, sampleSize, nTests)) / Settings.nIterations;
        FWER.blaireKarniski(sampleSize, nTests) = sum(FWER.FalsePositives.blaireKarniskiAllTests(:, sampleSize, nTests)) / Settings.nIterations;
        FWER.cluster(sampleSize, nTests) = sum(FWER.FalsePositives.clusterAllTests(:, sampleSize, nTests)) / Settings.nIterations;
        FWER.GFWER_ktms(sampleSize, nTests) = sum(FWER.FalsePositives.ktmsAllTests(:, sampleSize, nTests)) / Settings.nIterations;

        % Estimate the expected (average) number of false positives within
        % a given family of tests
        FalsePosRate.uncorrected(sampleSize, nTests) = nanmean(FWER.FalsePositives.uncorrected(:, sampleSize, nTests));
        FalsePosRate.bonferroni(sampleSize, nTests) = nanmean(FWER.FalsePositives.bonferroni(:, sampleSize, nTests));
        FalsePosRate.holm(sampleSize, nTests) = nanmean(FWER.FalsePositives.holm(:, sampleSize, nTests));
        FalsePosRate.bh(sampleSize, nTests) = nanmean(FWER.FalsePositives.bh(:, sampleSize, nTests));
        FalsePosRate.bky(sampleSize, nTests) = nanmean(FWER.FalsePositives.bky(:, sampleSize, nTests));
        FalsePosRate.by(sampleSize, nTests) = nanmean(FWER.FalsePositives.by(:, sampleSize, nTests));
        FalsePosRate.blaireKarniski(sampleSize, nTests) = nanmean(FWER.FalsePositives.blaireKarniski(:, sampleSize, nTests));
        FalsePosRate.cluster(sampleSize, nTests) = nanmean(FWER.FalsePositives.cluster(:, sampleSize, nTests));
        FalsePosRate.ktms(sampleSize, nTests) = nanmean(FWER.FalsePositives.ktms(:, sampleSize, nTests));

    end % of for sampleSize

end % of for nTests

% Save all results and settings to the workspace
save(['Workspace Saves/workspace_save' datestr(now, 30)]);


% Plot FWER of each method
figure('Name', 'FWER for each method');
plot(FWER.uncorrected);
hold on;
plot(FWER.bonferroni);
plot(FWER.holm);
plot(FWER.bh);
plot(FWER.bky);
plot(FWER.by);
plot(FWER.blaireKarniski);
plot(FWER.cluster);
plot(FWER.GFWER_ktms);


% Plot False Positive Rate of each method
figure('Name', 'False positive rate for each method');
plot(FalsePosRate.uncorrected);
hold on;
plot(FalsePosRate.bonferroni);
plot(FalsePosRate.holm);
plot(FalsePosRate.bh);
plot(FalsePosRate.bky);
plot(FalsePosRate.by);
plot(FalsePosRate.blaireKarniski);
plot(FalsePosRate.cluster);