function results = run_identification(trainingData, validationData, config)
%RUN_IDENTIFICATION Identify two-axis antenna dynamics from numeric arrays.
% trainingData and validationData follow the column schema in README.md.
% Pass [] for validationData when no separate validation run is available.
% Physical consistency estimation requires YALMIP and SDPT3.

arguments
    trainingData (:,:) double
    validationData (:,:) double = []
    config.sampleRate (1,1) double {mustBePositive} = 200
    config.cutoffFrequency (1,1) double {mustBePositive} = 0.5
    config.initialAngle (2,1) double = [0; 50]
    config.gearRatio (2,1) double {mustBePositive} = [560; 700]
    config.trainingRange (1,2) double {mustBeInteger, mustBePositive} = [1, 1]
    config.validationRange (1,2) double {mustBeInteger, mustBePositive} = [1, 1]
    config.method (1,1) string {mustBeMember(config.method,["physical","ols"])} = "physical"
    config.randomSeed (1,1) double {mustBeInteger} = 0
end

if size(trainingData,2) < 7
    error('Training input must have at least seven columns.');
end
if ~isempty(validationData) && size(validationData,2) < 7
    error('Validation input must have at least seven columns.');
end
if config.cutoffFrequency >= config.sampleRate/2
    error('Cutoff frequency must be below the Nyquist frequency.');
end
trainingRange = config.trainingRange;
if isequal(trainingRange, [1, 1])
    trainingRange = [1, size(trainingData,1)];
end
validationRange = config.validationRange;
if isequal(validationRange, [1, 1])
    validationRange = [1, size(validationData,1)];
end
if trainingRange(2) > size(trainingData,1) || trainingRange(1) >= trainingRange(2)
    error('trainingRange is outside trainingData or too short.');
end
if ~isempty(validationData) && (validationRange(2) > size(validationData,1) || validationRange(1) >= validationRange(2))
    error('validationRange is outside validationData or too short.');
end

train = Filtered_data(trainingData, trainingRange(1), trainingRange(2), ...
    config.cutoffFrequency, config.sampleRate, config.gearRatio, config.initialAngle);
[~, baseQR] = Base_paras_QR(config.randomSeed);
[measuredTorque, observationMatrix] = Build_Observation_Matrices(train, baseQR, config.gearRatio);

solution = struct;
switch config.method
    case "physical"
        [solution.pi_b, solution.pi_fr, solution.pi_full, solverInfo] = ...
            PhysicallyConsistentEstimation(measuredTorque, observationMatrix, baseQR);
        if solverInfo.problem ~= 0
            error('Physical consistency optimization failed: %s', solverInfo.info);
        end
    case "ols"
        [solution.pi_b, solution.pi_fr] = OrdinaryLeastSquareEstimation( ...
            measuredTorque, observationMatrix, baseQR.numberOfBaseParameters);
        solverInfo = [];
end

[~, trainError1, trainError2] = Display_result(solution, train, observationMatrix, config.gearRatio);
results = struct;
results.solution = solution;
results.baseQR = baseQR;
results.solverInfo = solverInfo;
results.trainingRMSE = [sqrt(mean(trainError1.^2)); sqrt(mean(trainError2.^2))];

if ~isempty(validationData)
    validation = Filtered_data(validationData, validationRange(1), validationRange(2), ...
        config.cutoffFrequency, config.sampleRate, config.gearRatio, config.initialAngle);
    [~, validationMatrix] = Build_Observation_Matrices(validation, baseQR, config.gearRatio);
    [~, validationError1, validationError2] = Display_result( ...
        solution, validation, validationMatrix, config.gearRatio);
    results.validationRMSE = [sqrt(mean(validationError1.^2)); sqrt(mean(validationError2.^2))];
end
end
