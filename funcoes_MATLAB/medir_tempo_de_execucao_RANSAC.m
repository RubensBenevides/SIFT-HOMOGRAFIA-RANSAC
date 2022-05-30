%% IMPORTAR IMAGENS
clc
clear
close all
format long g

% Importar par imagens A e B
imArgb = im2single(imread('dataset_Nikon_P600/e.jpg')) ;
imBrgb = im2single(imread('dataset_Nikon_P600/f.jpg')) ;
imA = rgb2gray(imArgb) ;
imB = rgb2gray(imBrgb) ;

%% DISTORCER E AMOSTRAR IMAGENS

% Carregar parametros de calibracao da camera.
data = load(fullfile('calibracao\Nikon_P900\cameraParams.mat'));
cameraParams = data.cameraParams;
intrinsics = cameraParams.Intrinsics ;

% Distorcer o par de imagens atuais.
imA = undistortImage(imA, intrinsics) ; 
imB = undistortImage(imB, intrinsics) ; 

% Amostrar imagens para obter desempenho
% Imagens maiores que 3000x3000 exigem mais de 8 GB de memoria RAM
r_max = 2000 ;              
if any(size(imA) > r_max) || any(size(imB) > r_max)
    imA = imresize(imA, 1/mean(size(imA)/r_max), "bilinear") ;
    imB = imresize(imB, 1/mean(size(imB)/r_max), "bilinear") ;
    imArgb = imresize(imArgb, 1/mean(size(imArgb,1,2)/r_max), "bilinear") ;
    imBrgb = imresize(imBrgb, 1/mean(size(imBrgb,1,2)/r_max), "bilinear") ;
end

%% Calcular tempo de execucao de SIFT + HOMOGRAFIA + RANSAC

% f = funcao que mede o tempo de exacucao da funcao
f = @() SIFT_HOMOGRAFIA_RANSAC(imA,imB); 
timeit(f)


