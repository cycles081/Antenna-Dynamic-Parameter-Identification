# Two-Axis Antenna Dynamic Parameter Identification

MATLAB code for identifying the dynamic parameters of a two-axis (azimuth and elevation) antenna. The workflow preprocesses measured joint signals, builds a Newton–Euler dynamics regressor, selects a base parameter set by pivoted QR decomposition, and estimates parameters using ordinary least squares or a physically consistent semidefinite formulation. A separate recording can be used for torque prediction validation.

This repository contains **code only**. Experimental recordings, Excel workbooks, and saved MATLAB results are not included. No numerical accuracy claim is made without those recordings.

## Pipeline

1. `Filtered_data.m` derives velocity and acceleration from recorded position, applies a Butterworth low-pass filter, and converts joint angles to radians.
2. `RegressorWithMotor.m` and `Base_paras_QR.m` construct the inertial regressor and select independent base parameters.
3. `RegressorWithoutMotor.m`, `FrictionRegressor.m`, and `Build_Observation_Matrices.m` assemble the sampled observation matrix and measured torque vector.
4. `OrdinaryLeastSquareEstimation.m` provides an unconstrained baseline. `PhysicallyConsistentEstimation.m` adds mass, inertia, and friction constraints with YALMIP and SDPT3.
5. `Display_result.m` plots measured torque, predicted torque, and residuals. `run_identification.m` connects these steps and reports per-joint RMSE.

## Requirements

- MATLAB with support for `arguments` blocks and `tiledlayout`.
- Signal Processing Toolbox for `butter` and `filtfilt`.
- For `method="physical"`: [YALMIP](https://yalmip.github.io/) and the [SDPT3](https://github.com/sqlp/sdpt3) solver, available on the MATLAB path.
- `method="ols"` does not use YALMIP or SDPT3.

The regressor uses a fixed link offset of **0.132 m** and gravitational acceleration of **9.8 m/s²** in `RegressorWithMotor.m` and `RegressorWithoutMotor.m`. Adapt these values to another mechanism before interpreting its parameters.

Despite its legacy filename, `RegressorWithMotor.m` currently returns a 20-column rigid-body regressor and does **not** contain separate motor inertia columns.

## Input format

Pass a numeric matrix with at least seven columns. The code reads the following columns from the original acquisition layout:

| Column | Meaning in the current code |
| --- | --- |
| 1 | Sample index; time is calculated as column 1 divided by `sampleRate` |
| 2 | Elevation position in degrees; converted as `85 - value` |
| 4 | Elevation motor torque signal; converted as `-value * 0.147 / 1000` |
| 5 | Azimuth position in degrees |
| 7 | Azimuth motor torque signal; converted as `value * 0.147 / 1000` |

These conversions are specific to the original measurement setup. Confirm calibration, signs, and gear ratios before using other recordings. Columns 3 and 6 are not used.

## Run

Place these `.m` files on the MATLAB path. Load your own recordings into numeric arrays using an appropriate importer, then call:

```matlab
train = readmatrix("path/to/training_recording.xlsx");
valid = readmatrix("path/to/validation_recording.xlsx");

result = run_identification(train, valid, ...
    sampleRate=200, cutoffFrequency=0.5, ...
    initialAngle=[0;50], gearRatio=[560;700], ...
    trainingRange=[4000,16000], validationRange=[1,12000], ...
    method="physical", randomSeed=0);

disp(result.trainingRMSE)
disp(result.validationRMSE)
```

The ranges above illustrate the original experiment configuration; adjust them to your recordings. Omit the range options to use all rows. Pass `[]` as the second argument if no separate validation recording is available. For the unconstrained baseline, use `method="ols"`.

```matlab
result = run_identification(train, [], method="ols");
```

The entry point returns the parameter estimates, QR decomposition information, solver status (when used), and torque RMSE for each joint. It also opens diagnostic plots.

## Scope and limitations

- This is a research implementation for one antenna geometry and acquisition format, not a general robot identification package.
- The original recording files are not distributed, so the published code cannot reproduce the original numerical results on its own.
- The physical consistency constraints are formulated in YALMIP and require a successful SDPT3 solve. Check `result.solverInfo` and the returned parameters before drawing physical conclusions.
- The separate excitation trajectory optimization experiments were left out of this repository because their archived versions have inconsistent function signatures and are not part of the verified identification path.

## Repository contents

| File | Role |
| --- | --- |
| `run_identification.m` | Data-independent entry point for the identification and validation workflow |
| `Filtered_data.m` | Signal preparation and plotting |
| `RegressorWithMotor.m`, `RegressorWithoutMotor.m` | Two-axis Newton–Euler dynamics regressors |
| `Base_paras_QR.m` | Base parameter selection by pivoted QR |
| `FrictionRegressor.m`, `Build_Observation_Matrices.m` | Friction and stacked observation matrix |
| `OrdinaryLeastSquareEstimation.m` | Least-squares baseline |
| `PhysicallyConsistentEstimation.m` | Constrained parameter estimation |
| `Display_result.m` | Torque comparison and residual plotting |

## Background

Developed as part of a two-axis antenna dynamics identification project. If you use this code in research, cite the related publication or thesis once its bibliographic details are available.
