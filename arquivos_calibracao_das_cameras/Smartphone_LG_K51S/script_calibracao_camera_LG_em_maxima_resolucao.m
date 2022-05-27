% Auto-generated by cameraCalibrator app on 11-May-2022
%-------------------------------------------------------


% Define images to process
imageFileNames = {'C:\Users\ruben\OneDrive\Documentos\PDI_Doutorado\mosaicagem_2\calibracao\Smartphone_LG_K51S\20220511_204526.jpg',...
    'C:\Users\ruben\OneDrive\Documentos\PDI_Doutorado\mosaicagem_2\calibracao\Smartphone_LG_K51S\20220511_204557.jpg',...
    'C:\Users\ruben\OneDrive\Documentos\PDI_Doutorado\mosaicagem_2\calibracao\Smartphone_LG_K51S\20220511_204604.jpg',...
    'C:\Users\ruben\OneDrive\Documentos\PDI_Doutorado\mosaicagem_2\calibracao\Smartphone_LG_K51S\20220511_204608.jpg',...
    'C:\Users\ruben\OneDrive\Documentos\PDI_Doutorado\mosaicagem_2\calibracao\Smartphone_LG_K51S\20220511_204617.jpg',...
    'C:\Users\ruben\OneDrive\Documentos\PDI_Doutorado\mosaicagem_2\calibracao\Smartphone_LG_K51S\20220511_204624.jpg',...
    'C:\Users\ruben\OneDrive\Documentos\PDI_Doutorado\mosaicagem_2\calibracao\Smartphone_LG_K51S\20220511_204626.jpg',...
    'C:\Users\ruben\OneDrive\Documentos\PDI_Doutorado\mosaicagem_2\calibracao\Smartphone_LG_K51S\20220511_204643.jpg',...
    'C:\Users\ruben\OneDrive\Documentos\PDI_Doutorado\mosaicagem_2\calibracao\Smartphone_LG_K51S\20220511_204653.jpg',...
    'C:\Users\ruben\OneDrive\Documentos\PDI_Doutorado\mosaicagem_2\calibracao\Smartphone_LG_K51S\20220511_204722.jpg',...
    'C:\Users\ruben\OneDrive\Documentos\PDI_Doutorado\mosaicagem_2\calibracao\Smartphone_LG_K51S\20220511_204728.jpg',...
    'C:\Users\ruben\OneDrive\Documentos\PDI_Doutorado\mosaicagem_2\calibracao\Smartphone_LG_K51S\20220511_204738.jpg',...
    'C:\Users\ruben\OneDrive\Documentos\PDI_Doutorado\mosaicagem_2\calibracao\Smartphone_LG_K51S\20220511_204742.jpg',...
    'C:\Users\ruben\OneDrive\Documentos\PDI_Doutorado\mosaicagem_2\calibracao\Smartphone_LG_K51S\20220511_204754.jpg',...
    'C:\Users\ruben\OneDrive\Documentos\PDI_Doutorado\mosaicagem_2\calibracao\Smartphone_LG_K51S\20220511_204826.jpg',...
    };
% Detect calibration pattern in images
detector = vision.calibration.monocular.CheckerboardDetector();
[imagePoints, imagesUsed] = detectPatternPoints(detector, imageFileNames, 'HighDistortion', true);
imageFileNames = imageFileNames(imagesUsed);

% Read the first image to obtain image size
originalImage = imread(imageFileNames{1});
[mrows, ncols, ~] = size(originalImage);

% Generate world coordinates for the planar pattern keypoints
squareSize = 15;  % in units of 'centimeters'
worldPoints = generateWorldPoints(detector, 'SquareSize', squareSize);

% Calibrate the camera
[cameraParams, imagesUsed, estimationErrors] = estimateCameraParameters(imagePoints, worldPoints, ...
    'EstimateSkew', true, 'EstimateTangentialDistortion', true, ...
    'NumRadialDistortionCoefficients', 2, 'WorldUnits', 'centimeters', ...
    'InitialIntrinsicMatrix', [], 'InitialRadialDistortion', [], ...
    'ImageSize', [mrows, ncols]);

% View reprojection errors
h1=figure; showReprojectionErrors(cameraParams);

% Visualize pattern locations
h2=figure; showExtrinsics(cameraParams, 'CameraCentric');

% Display parameter estimation errors
displayErrors(estimationErrors, cameraParams);

% For example, you can use the calibration data to remove effects of lens distortion.
undistortedImage = undistortImage(originalImage, cameraParams);

% See additional examples of how to use the calibration data.  At the prompt type:
% showdemo('MeasuringPlanarObjectsExample')
% showdemo('StructureFromMotionExample')
