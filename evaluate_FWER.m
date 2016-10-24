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
% Conversely, the False Null Rate (FNR) corresponds to the expeted
% proportion of all accepted null hypotheses that are false negatives. 
%
% Written by DF 10/16


% Houskeeping
clear all;
close all;

% Seed random number generator based on computer clock
rng('shuffle');

% Simulation Settings
Settings.sampleSizesToUse = [30]; % Sample size per test
Settings.nTestsToUse = [10 50 100]; % Number of tests
Settings.trueEffectProportion = 0.2; % Proportion of hypotheses that are real effects (between 0 and 1)
Settings.meanEffect = 1; % Mean magnitude of the effect (note: SD is on average 1)
Settings.nIterations = 3; % Number of iterations
Settings.alphaLevel = 0.05; % Nominal alpha level

% Resampling method settings
Settings.blaireKarniskiIterations = 1000; % Number of permutation samples to draw for the Blaire-Karniski correction
Settings.clusterIterations = 1000; % Number of permutation samples to draw for the maximum cluster mass null distribution
Settings.clusteringAlphaLevel = 0.05; % Alpha level for detecting individual tests to include within a cluster
Settings.ktmsIterations = 1000; % Number of iterations for ktms GFWER control procedure
Settings.ktms_u = 0; % u parameter for ktms GFWER control procedure

% Preallocate FWER false positive count matrices
FWER.FalsePositives.AllTests.uncorrected = zeros(Settings.nIterations, length(Settings.sampleSizesToUse), length(Settings.nTestsToUse));
FWER.FalsePositives.AllTests.bonferroni = zeros(Settings.nIterations, length(Settings.sampleSizesToUse), length(Settings.nTestsToUse));
FWER.FalsePositives.AllTests.holm = zeros(Settings.nIterations, length(Settings.sampleSizesToUse), length(Settings.nTestsToUse));
FWER.FalsePositives.AllTests.bh = zeros(Settings.nIterations, length(Settings.sampleSizesToUse), length(Settings.nTestsToUse));
FWER.FalsePositives.AllTests.bky = zeros(Settings.nIterations, length(Settings.sampleSizesToUse), length(Settings.nTestsToUse));
FWER.FalsePositives.AllTests.by = zeros(Settings.nIterations, length(Settings.sampleSizesToUse), length(Settings.nTestsToUse));
FWER.FalsePositives.AllTests.blaireKarniski = zeros(Settings.nIterations, length(Settings.sampleSizesToUse), length(Settings.nTestsToUse));
FWER.FalsePositives.AllTests.cluster = zeros(Settings.nIterations, length(Settings.sampleSizesToUse), length(Settings.nTestsToUse));
FWER.FalsePositives.AllTests.ktms = zeros(Settings.nIterations, length(Settings.sampleSizesToUse), length(Settings.nTestsToUse));

% Calculate number of true effects and true negatives within all tests
Settings.nTrueEffects = round(Settings.nTestsToUse * Settings.trueEffectProportion);
Settings.nTrueNegatives = Settings.nTestsToUse - Settings.nTrueEffects;

% Generate random samples for multiple hypothesis tests and record FWER for
% each method:
for nTests = 1:length(Settings.nTestsToUse)    
    for sampleSize = 1:length(Settings.sampleSizesToUse)
        for i = 1:Settings.nIterations
            fprintf('Running for %i tests with sample size %i iteration %i \n', Settings.nTestsToUse(nTests), Settings.sampleSizesToUse(sampleSize), i);
            
            % Generate a random samples from a normal distribution (SD = 1)
            clear tempSample1; clear tempSample2; % Clear out temporary samples from previous iteration
            for j = 1:Settings.nTestsToUse(nTests)
                tempSample1(:,j) = randn(Settings.sampleSizesToUse(sampleSize), 1);
                tempSample2(:,j) = randn(Settings.sampleSizesToUse(sampleSize), 1);
            end

            % Make vector of the ground truths (null or alternative hypothesis). 
            % 0 = true null, 1 = true alternative
            trueNullOrAlt = zeros(1, Settings.nTestsToUse(nTests));
            
            % Randomly allocate true effects
            trueEffectLocations = randi(Settings.nTestsToUse(nTests), 1, Settings.nTrueEffects(nTests));
            tempSample1(:, trueEffectLocations) = tempSample1(:, trueEffectLocations) + Settings.meanEffect;
            trueNullOrAlt(1, trueEffectLocations) = 1;
            
            % Perform paired-samples t test
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
            
            % Calculate the number of "hits" (true positives) using each method
            FWER.TruePositives.uncorrected(i, sampleSize, nTests) = sum(temp_h(trueNullOrAlt == 1));
            FWER.TruePositives.bonferroni(i, sampleSize, nTests) = sum(bonferroni_corrected_h(trueNullOrAlt == 1));
            FWER.TruePositives.holm(i, sampleSize, nTests) = sum(holm_corrected_h(trueNullOrAlt == 1));
            FWER.TruePositives.bh(i, sampleSize, nTests) = sum(fdr_bh_corrected_h(trueNullOrAlt == 1));
            FWER.TruePositives.bky(i, sampleSize, nTests) = sum(fdr_bky_corrected_h(trueNullOrAlt == 1));
            FWER.TruePositives.by(i, sampleSize, nTests) = sum(fdr_by_corrected_h(trueNullOrAlt == 0));
            FWER.TruePositives.blaireKarniski(i, sampleSize, nTests) = sum(blaireKarniski_corrected_h(trueNullOrAlt == 1));
            FWER.TruePositives.cluster(i, sampleSize, nTests) = sum(cluster_corrected_h(trueNullOrAlt == 1));
            FWER.TruePositives.ktms(i, sampleSize, nTests) = sum(ktms_h(trueNullOrAlt == 1));
            
            % Calculate the number of "correct rejections" (true negatives) using each method
            FWER.TrueNegatives.uncorrected(i, sampleSize, nTests) = Settings.nTrueNegatives(nTests) - sum(temp_h(trueNullOrAlt == 0));
            FWER.TrueNegatives.bonferroni(i, sampleSize, nTests) = Settings.nTrueNegatives(nTests) - sum(bonferroni_corrected_h(trueNullOrAlt == 0));
            FWER.TrueNegatives.holm(i, sampleSize, nTests) = Settings.nTrueNegatives(nTests) - sum(holm_corrected_h(trueNullOrAlt == 0));
            FWER.TrueNegatives.bh(i, sampleSize, nTests) = Settings.nTrueNegatives(nTests) - sum(fdr_bh_corrected_h(trueNullOrAlt == 0));
            FWER.TrueNegatives.bky(i, sampleSize, nTests) = Settings.nTrueNegatives(nTests) - sum(fdr_bky_corrected_h(trueNullOrAlt == 0));
            FWER.TrueNegatives.by(i, sampleSize, nTests) = Settings.nTrueNegatives(nTests) - sum(fdr_by_corrected_h(trueNullOrAlt == 0));
            FWER.TrueNegatives.blaireKarniski(i, sampleSize, nTests) = Settings.nTrueNegatives(nTests) - sum(blaireKarniski_corrected_h(trueNullOrAlt == 0));
            FWER.TrueNegatives.cluster(i, sampleSize, nTests) = Settings.nTrueNegatives(nTests) - sum(cluster_corrected_h(trueNullOrAlt == 0));
            FWER.TrueNegatives.ktms(i, sampleSize, nTests) = Settings.nTrueNegatives(nTests) - sum(ktms_h(trueNullOrAlt == 0));
            
            % Calculate the number of false positives using each method
            FWER.FalsePositives.uncorrected(i, sampleSize, nTests) = sum(temp_h(trueNullOrAlt == 0));
            FWER.FalsePositives.bonferroni(i, sampleSize, nTests) = sum(bonferroni_corrected_h(trueNullOrAlt == 0));
            FWER.FalsePositives.holm(i, sampleSize, nTests) = sum(holm_corrected_h(trueNullOrAlt == 0));
            FWER.FalsePositives.bh(i, sampleSize, nTests) = sum(fdr_bh_corrected_h(trueNullOrAlt == 0));
            FWER.FalsePositives.bky(i, sampleSize, nTests) = sum(fdr_bky_corrected_h(trueNullOrAlt == 0));
            FWER.FalsePositives.by(i, sampleSize, nTests) = sum(fdr_by_corrected_h(trueNullOrAlt == 0));
            FWER.FalsePositives.blaireKarniski(i, sampleSize, nTests) = sum(blaireKarniski_corrected_h(trueNullOrAlt == 0));
            FWER.FalsePositives.cluster(i, sampleSize, nTests) = sum(cluster_corrected_h(trueNullOrAlt == 0));
            FWER.FalsePositives.ktms(i, sampleSize, nTests) = sum(ktms_h(trueNullOrAlt == 0));

            % Calculate the number of false negatives using each method
            FWER.FalseNegatives.uncorrected(i, sampleSize, nTests) = Settings.nTrueEffects(nTests) - sum(temp_h(trueNullOrAlt == 1));
            FWER.FalseNegatives.bonferroni(i, sampleSize, nTests) = Settings.nTrueEffects(nTests) - sum(bonferroni_corrected_h(trueNullOrAlt == 1));
            FWER.FalseNegatives.holm(i, sampleSize, nTests) = Settings.nTrueEffects(nTests) - sum(holm_corrected_h(trueNullOrAlt == 1));
            FWER.FalseNegatives.bh(i, sampleSize, nTests) = Settings.nTrueEffects(nTests) - sum(fdr_bh_corrected_h(trueNullOrAlt == 1));
            FWER.FalseNegatives.bky(i, sampleSize, nTests) = Settings.nTrueEffects(nTests) - sum(fdr_bky_corrected_h(trueNullOrAlt == 1));
            FWER.FalseNegatives.by(i, sampleSize, nTests) = Settings.nTrueEffects(nTests) - sum(fdr_by_corrected_h(trueNullOrAlt == 1));
            FWER.FalseNegatives.blaireKarniski(i, sampleSize, nTests) = Settings.nTrueEffects(nTests) - sum(blaireKarniski_corrected_h(trueNullOrAlt == 1));
            FWER.FalseNegatives.cluster(i, sampleSize, nTests) = Settings.nTrueEffects(nTests) - sum(cluster_corrected_h(trueNullOrAlt == 1));
            FWER.FalseNegatives.ktms(i, sampleSize, nTests) = Settings.nTrueEffects(nTests) - sum(ktms_h(trueNullOrAlt == 1));
            
            % Calculate the False Discovery Proportion (FDP) for each
            % method. This is the proportion of false discoveries compared to
            % the number of total discoveries (rejections of the null
            % hypothesis).
            FDP.uncorrected(i, sampleSize, nTests) = FWER.FalsePositives.uncorrected(i, sampleSize, nTests) / (FWER.FalsePositives.uncorrected(i, sampleSize, nTests) + FWER.TruePositives.uncorrected(i, sampleSize, nTests));
            FDP.bonferroni(i, sampleSize, nTests) = FWER.FalsePositives.bonferroni(i, sampleSize, nTests) / (FWER.FalsePositives.bonferroni(i, sampleSize, nTests) + FWER.TruePositives.bonferroni(i, sampleSize, nTests));
            FDP.holm(i, sampleSize, nTests) = FWER.FalsePositives.holm(i, sampleSize, nTests) / (FWER.FalsePositives.holm(i, sampleSize, nTests) + FWER.TruePositives.holm(i, sampleSize, nTests));
            FDP.bh(i, sampleSize, nTests) = FWER.FalsePositives.bh(i, sampleSize, nTests) / (FWER.FalsePositives.bh(i, sampleSize, nTests) + FWER.TruePositives.bh(i, sampleSize, nTests));
            FDP.bky(i, sampleSize, nTests) = FWER.FalsePositives.bky(i, sampleSize, nTests) / (FWER.FalsePositives.bky(i, sampleSize, nTests) + FWER.TruePositives.bky(i, sampleSize, nTests));
            FDP.by(i, sampleSize, nTests) = FWER.FalsePositives.by(i, sampleSize, nTests) / (FWER.FalsePositives.by(i, sampleSize, nTests) + FWER.TruePositives.by(i, sampleSize, nTests));
            FDP.blaireKarniski(i, sampleSize, nTests) = FWER.FalsePositives.blaireKarniski(i, sampleSize, nTests) / (FWER.FalsePositives.blaireKarniski(i, sampleSize, nTests) + FWER.TruePositives.blaireKarniski(i, sampleSize, nTests));
            FDP.cluster(i, sampleSize, nTests) = FWER.FalsePositives.cluster(i, sampleSize, nTests) / (FWER.FalsePositives.cluster(i, sampleSize, nTests) + FWER.TruePositives.cluster(i, sampleSize, nTests));
            FDP.ktms(i, sampleSize, nTests) = FWER.FalsePositives.ktms(i, sampleSize, nTests) / (FWER.FalsePositives.ktms(i, sampleSize, nTests) + FWER.TruePositives.ktms(i, sampleSize, nTests));

            
        end % of for i = 1:Settings.nIterations loop

        FWER.FalsePositives.AllTests.uncorrected(FWER.FalsePositives.uncorrected > 0) = 1;
        FWER.FalsePositives.AllTests.bonferroni(FWER.FalsePositives.bonferroni > 0) = 1;
        FWER.FalsePositives.AllTests.holm(FWER.FalsePositives.holm > 0) = 1;
        FWER.FalsePositives.AllTests.bh(FWER.FalsePositives.bh > 0) = 1;
        FWER.FalsePositives.AllTests.bky(FWER.FalsePositives.bky > 0) = 1;
        FWER.FalsePositives.AllTests.by(FWER.FalsePositives.by > 0) = 1;
        FWER.FalsePositives.AllTests.blaireKarniski(FWER.FalsePositives.blaireKarniski > 0) = 1;
        FWER.FalsePositives.AllTests.cluster(FWER.FalsePositives.cluster > 0) = 1;
        FWER.FalsePositives.AllTests.ktms(FWER.FalsePositives.ktms > 0) = 1; % No. of false positives needs to be above ktms_u parameter for GFWER

        % Calculate the FWER (type 1 error rate)
        FWER.uncorrected(sampleSize, nTests) = sum(FWER.FalsePositives.AllTests.uncorrected(:, sampleSize, nTests)) / Settings.nIterations;
        FWER.bonferroni(sampleSize, nTests) = sum(FWER.FalsePositives.AllTests.bonferroni(:, sampleSize, nTests)) / Settings.nIterations;
        FWER.holm(sampleSize, nTests) = sum(FWER.FalsePositives.AllTests.holm(:, sampleSize, nTests)) / Settings.nIterations;
        FWER.bh(sampleSize, nTests) = sum(FWER.FalsePositives.AllTests.bh(:, sampleSize, nTests)) / Settings.nIterations;
        FWER.bky(sampleSize, nTests) = sum(FWER.FalsePositives.AllTests.bky(:, sampleSize, nTests)) / Settings.nIterations;
        FWER.by(sampleSize, nTests) = sum(FWER.FalsePositives.AllTests.by(:, sampleSize, nTests)) / Settings.nIterations;
        FWER.blaireKarniski(sampleSize, nTests) = sum(FWER.FalsePositives.AllTests.blaireKarniski(:, sampleSize, nTests)) / Settings.nIterations;
        FWER.cluster(sampleSize, nTests) = sum(FWER.FalsePositives.AllTests.cluster(:, sampleSize, nTests)) / Settings.nIterations;
        FWER.GFWER_ktms(sampleSize, nTests) = sum(FWER.FalsePositives.AllTests.ktms(:, sampleSize, nTests)) / Settings.nIterations;

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

        % Estimate the expected (average) number of false negatives within
        % a given family of tests
        FalseNegRate.uncorrected(sampleSize, nTests) = nanmean(FWER.FalseNegatives.uncorrected(:, sampleSize, nTests));
        FalseNegRate.bonferroni(sampleSize, nTests) = nanmean(FWER.FalseNegatives.bonferroni(:, sampleSize, nTests));
        FalseNegRate.holm(sampleSize, nTests) = nanmean(FWER.FalseNegatives.holm(:, sampleSize, nTests));
        FalseNegRate.bh(sampleSize, nTests) = nanmean(FWER.FalseNegatives.bh(:, sampleSize, nTests));
        FalseNegRate.bky(sampleSize, nTests) = nanmean(FWER.FalseNegatives.bky(:, sampleSize, nTests));
        FalseNegRate.by(sampleSize, nTests) = nanmean(FWER.FalseNegatives.by(:, sampleSize, nTests));
        FalseNegRate.blaireKarniski(sampleSize, nTests) = nanmean(FWER.FalseNegatives.blaireKarniski(:, sampleSize, nTests));
        FalseNegRate.cluster(sampleSize, nTests) = nanmean(FWER.FalseNegatives.cluster(:, sampleSize, nTests));
        FalseNegRate.ktms(sampleSize, nTests) = nanmean(FWER.FalseNegatives.ktms(:, sampleSize, nTests));
        
        
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


% Plot average number of false positives for each method
figure('Name', 'Average number of false positives for each method');
plot(FalsePosRate.uncorrected);
hold on;
plot(FalsePosRate.bonferroni);
plot(FalsePosRate.holm);
plot(FalsePosRate.bh);
plot(FalsePosRate.bky);
plot(FalsePosRate.by);
plot(FalsePosRate.blaireKarniski);
plot(FalsePosRate.cluster);

% Plot average number of false negatives for each method
figure('Name', 'Average number of false negatives for each method');
plot(FalseNegRate.uncorrected);
hold on;
plot(FalseNegRate.bonferroni);
plot(FalseNegRate.holm);
plot(FalseNegRate.bh);
plot(FalseNegRate.bky);
plot(FalseNegRate.by);
plot(FalseNegRate.blaireKarniski);
plot(FalseNegRate.cluster);