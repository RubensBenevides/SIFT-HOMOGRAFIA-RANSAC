clc
clear
close all
% Ler par de imagens e colocar no intervalo [0,1]
im1 = im2single(imread('fotos_diversas/s.jpg')) ;
im2 = im2single(imread('fotos_diversas/t.jpg')) ;

% Converter para niveis de cinza
if size(im1,3) > 1, im1g = rgb2gray(im1) ; end
if size(im2,3) > 1, im2g = rgb2gray(im2) ; end

% Carregar parametros de calibracao da camera
data = load(fullfile('C:\Users\ruben\OneDrive\Documentos\PDI_Doutorado\mosaicagem_2\calibracao\cameraParams.mat'));
cameraParams = data.cameraParams;

% Distorcer o par de imagens atuais.
intrinsics = cameraParams.Intrinsics ;
Iprev = undistortImage(im1g, intrinsics) ; 
Inext = undistortImage(im2g, intrinsics) ; 

% Detectar pontos SIFT e fazer o match com RANSAC e calcular H_otima.
[prevPoints, nextPoints, matchedPoints1, matchedPoints2, ...
    ~, ~, indexPairs, ok, ~] = sift_ransac(Iprev,Inext) ;

% Mostre as corresp. nas duas primeiras imagens com e sem filtragem
show_matches(prevPoints, nextPoints, indexPairs, ok, im1, im2) ;

% Filtrar correspondencias dos indices obtidos na saida RANSAC
idx_bons = indexPairs.*double(ok) ;
idx_bons = nonzeros(idx_bons) ;
idx_bons = reshape(idx_bons, 2, size(idx_bons,1)/2) ;

% Calcular Homografias antes e depois do RANSAC
H_ruim  = calcular_H2(prevPoints, nextPoints, indexPairs) ;
H_boa  = calcular_H2(prevPoints, nextPoints, idx_bons) ;

% Mostre o mosaico feito com e sem filtragem
show_mosaico(im1, im2, H_boa) ;
%show_mosaico(im1, im2, H_ruim) ;
