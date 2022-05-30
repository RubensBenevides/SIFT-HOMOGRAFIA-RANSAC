%% IMPORTAR IMAGENS & DISTORCER
clc
clear
close all
format long g

% Importar par imagens A e B
imArgb = im2single(imread('dataset_LG_K51s_GEENG/e.jpg')) ;
imBrgb = im2single(imread('dataset_LG_K51s_GEENG/f.jpg')) ;
imA = rgb2gray(imArgb) ;
imB = rgb2gray(imBrgb) ;

%% DISTORCER IMAGENS
% Carregar parametros de calibracao da camera.
data = load(fullfile('calibracao\Smartphone_LG_K51S\cameraParams.mat'));
cameraParams = data.cameraParams_LG_K51s;
intrinsics = cameraParams.Intrinsics ;

% Distorcer o par de imagens atuais.
imA = undistortImage(imA, intrinsics) ; 
imB = undistortImage(imB, intrinsics) ; 

% Amostrar imagens para obter desempenho
r_max = 2000 ;               % Se passar de 3000 exige mais de 8 GB de RAM
if any(size(imA) > r_max) || any(size(imB) > r_max)
    imA = imresize(imA, 1/mean(size(imA)/r_max), "bilinear") ;
    imB = imresize(imB, 1/mean(size(imB)/r_max), "bilinear") ;
    imArgb = imresize(imArgb, 1/mean(size(imArgb,1,2)/r_max), "bilinear") ;
    imBrgb = imresize(imBrgb, 1/mean(size(imBrgb,1,2)/r_max), "bilinear") ;
end

%% OBTER PONTOS HOMÓLOGOS POR MEIO DO DESCRITOR SIFT
% Descrever pontos com SIFT
pts_A_sift = detectSIFTFeatures(imA);
pts_B_sift = detectSIFTFeatures(imB);

% Extrair apenas pts. relevantes com orientacao.
[descritores_A, pts_A_sift_ori] = extractFeatures(imA, pts_A_sift);
[descritores_B, pts_B_sift_ori] = extractFeatures(imB, pts_B_sift);

% Match para obter indices dos pares.
ind_pares = matchFeatures(descritores_A,descritores_B);
numMatches = length(ind_pares);

% Selecao das coordenadas na imA que fazem match com imB e vice-versa.
pts_corr_A = pts_A_sift_ori(ind_pares(:,1), :);
pts_corr_B = pts_B_sift_ori(ind_pares(:,2), :);

%% NORMALIZAR PONTOS E TRANSFORMAR EM COORDENADAS HOMOGENEAS
% Transformar objeto SIFTPoint para matriz de coordenadas.
pts_A_xy = double(pts_corr_A.Location) ;
pts_B_xy = double(pts_corr_B.Location) ;

% Obter centroide das coordenadas.
medias_A = mean(pts_A_xy) ; txA=medias_A(1) ; tyA=medias_A(2) ;  
medias_B = mean(pts_B_xy) ; txB=medias_B(1) ; tyB=medias_B(2) ;

% Centralizar.
pts_A_centro = pts_A_xy - medias_A ;
pts_B_centro = pts_B_xy - medias_B ;

% Obter escalas.
aux_A = sum(sqrt(pts_A_centro(:,1).^2 + pts_A_centro(:,2).^2))/numMatches ;
aux_B = sum(sqrt(pts_B_centro(:,1).^2 + pts_B_centro(:,2).^2))/numMatches ;
sA = sqrt(2)/aux_A ;
sB = sqrt(2)/aux_B ;

% Montar transformacoes que normalizam pts. em A e pts. em B.
TN_A = [sA   0 -sA*txA;   
        0   sA -sA*tyA; 
        0    0      1];
TN_B = [sB   0 -sB*txB; 
        0   sB -sB*tyB; 
        0    0      1];

% Transformar em coordenadas homogeneas.
pts_A = [pts_A_xy, ones(numMatches,1)]' ;
pts_B = [pts_B_xy, ones(numMatches,1)]' ;

% Aplicar Normalizacao.
pts_A_Norm = TN_A*pts_A ;
pts_B_Norm = TN_B*pts_B ;

%% CALCULAR HOMOGRAFIA COM RANSAC ADAPTATIVO
% Iniciar RANSAC com nº de iteracoes elevado (k = 1005 significa 
% capacidade de lidar com ate 75 % de outliers nos pares de correspondenc.)

erros = 1:30 ;
matriz_desvios_dos_inliers = zeros(100, length(erros)) ;
% Loop que testa 10 limiares de erro(e) no RANSAC
for j=1:1:length(erros)
    j
    desvios_dos_inliers = zeros(100,1) ;
    % Loop que roda o RANSAC 100x
    for i=1:1:100
        % Inicializar celulas que receberam Homografias e escores do RANSAC
        clear H score ok
        H = cell(1,1005) ; 
        ok = cell(1,1005) ; 
        escore = zeros(1,1005) ;
        % Inicializar parametros do RANSAC
        k = 1005 ;  t = 0 ; iteracoes_RANSAC = 0 ;
        % Loop do RANSAC adaptativo nas iteracoes k
        while (t <= k)
            t = t + 1 ;
            % Seleciona 4 pares de pts. aleatorios (min. para a Homografia).
            random_ind_pairs = datasample(1:numMatches, 4) ; 
        
            % Monta a matriz de HOMOGRAFIA por produto tensorial de Kronecker
            % Utiliza o vetor em forma de matriz anti-simétrica (skew-simetric)
            idx1 = random_ind_pairs(1);
            idx2 = random_ind_pairs(2);
            idx3 = random_ind_pairs(3);
            idx4 = random_ind_pairs(4);
            A1 = kron(pts_A_Norm(:,idx1)', vec2skew(pts_B_Norm(:,idx1)) ) ;
            A2 = kron(pts_A_Norm(:,idx2)', vec2skew(pts_B_Norm(:,idx2)) ) ;
            A3 = kron(pts_A_Norm(:,idx3)', vec2skew(pts_B_Norm(:,idx3)) ) ;
            A4 = kron(pts_A_Norm(:,idx4)', vec2skew(pts_B_Norm(:,idx4)) ) ;
            A = [A1;A2;A3;A4] ;
            [~,~,V] = svd(A) ;
            
            % A ultima coluna de V corresponde ao autovetor de menor autovalor
            % associado. Neste vetor estao os parametros da Homografia normalizados.
            H{t} = reshape(V(:,9),3,3);
            
            % Testar Homografia estimada transformando coordenadas de A para B.
            pts_B_Norm_ = H{t}*pts_A_Norm ;
            
            % Normaliza pelo fator de escala imaterial z e subtrai das 
            % coordenadas originais as transformadas de A para B.
            e_x = pts_B_Norm_(1,:)./pts_B_Norm_(3,:) - pts_B_Norm(1,:) ; % erro em x: (x'/z' - x)
            e_y = pts_B_Norm_(2,:)./pts_B_Norm_(3,:) - pts_B_Norm(2,:) ; % erro em y: (y'/z' - y)
            
            % Vetor logico para informar quantos pontos obedecem a restricao do 
            % limiar em pixels na escala. Rodar 100x para fazer o fine-tunning.
            % Testar 10 valores de erro. e = limiar de erro em pixels
            e = erros(j) ;

            % Distancia euclidiana
            ok{t} =  (sqrt(e_x.*e_x + e_y.*e_y)) < e*min(sA,sB) ;
            
            % Dist. de Manhattan, metrica mais rapida e menos precisa.
            %ok{t} = (abs(e_x) + abs(e_y)) < e*min(sA,sB) ;  
        
            % Escore de H = cardinalidade dos inliers. Total de 1's em ok
            escore(t) = sum(ok{t}) ;
        
            % Estimar taxa omega (inliers/pares).
            omega = escore(t)/numMatches ;
          
            % Novo nº de iteracoes do RANSAC.
            iteracoes = round(log(1-0.99)/log(1-omega^4)) ;
        
            % Se o novo nº de iteracoes for menor que o restante, atualize k.
            if (iteracoes < k-t) && (k-t > 0)
                k = iteracoes ; 
                t = 1;
            end
            % Salvar quantas vezes o RANSAC iterou ate a convergencia.
            iteracoes_RANSAC = iteracoes_RANSAC + 1 ;
        end
        % best = indice da matriz com maior score dentro da celula H.
        [escore, best] = max(escore) ;
        H = H{best} ;
        ok = ok{best} ;  % escolhe o vetor logico de maior cardinalidade
        omega = sum(ok)/numMatches ;

        % Mostrar metricas do resultado do RANSAC
        %fprintf("\nTaxa inliers: %.2f porcento\n", omega) ;
        %fprintf("Iteracoes do RANSAC: %d\n", iteracoes_RANSAC) ;
        %fprintf("Qualidade da Homografia: %d inliers de %d matches\n", escore, numMatches) ;

        % Calcular o erro com i-esima rodada do RANSAC
        pts_B_Norm_ = H*pts_A_Norm ;
        e_x = pts_B_Norm_(1,:)./pts_B_Norm_(3,:) - pts_B_Norm(1,:) ; % erro em x: (x'/z' - x)
        e_y = pts_B_Norm_(2,:)./pts_B_Norm_(3,:) - pts_B_Norm(2,:) ; % erro em y: (y'/z' - y)

        % Calcular Coeficiente de Variacao do desvio padrao da media.
        erros_inliers = nonzeros(sqrt(e_x.*e_x + e_y.*e_y).*ok) ;
        CV = std(erros_inliers)/(sqrt(length(erros_inliers))*mean(erros_inliers)) ;
        desvios_dos_inliers(i) = CV ;
    end
    % Salvar os 100 resultados do escore para aquele erro (e)
    matriz_desvios_dos_inliers(:,j) = desvios_dos_inliers ;
end

%% Análise estatística do erro dos inliers segundo limiar de erro (ε)

%matriz_desvios_dos_inliers = log(matriz_desvios_dos_inliers) ;
% Plot de todos os dados (waterfall)
subplot(2,2,1)
waterfall(matriz_desvios_dos_inliers')
title('CV dos inliers em cada limiar ε')
xlabel('100 rodadas do RANSAC')
ylabel('limiar de erro ε')
zlabel('CV dos inliers)')

% Boxplot das distribuicoes de omega(ω) segundo o limiar de erro(ε)
subplot(2,2,2)
boxplot(matriz_desvios_dos_inliers)
title('Boxplot do CV em cada limiar ε')
xlabel('limiar de erro ε')
ylabel('CV dos inliers)')

% Plot da variacao da media de -log(C.V) e do C.V. de -log(C.V) vs ε
subplot(2,2,3)
errorbar(mean(matriz_desvios_dos_inliers), std(matriz_desvios_dos_inliers)) ;
title('Média do CV dos inliers em cada limiar ε. Barra de erro: desvio padrão do CV dos inliers')
xlabel('limiar de erro ε')
ylabel('CV dos inliers)')

% Plot da variacao da mediana de -log(C.V) e da MAD de -log(C.V) segundo ε
subplot(2,2,4)
errorbar(median(matriz_desvios_dos_inliers), mad(matriz_desvios_dos_inliers,1)) ;
title('Mediana do CV dos inliers em cada limiar ε. Barra de erro: Median Absolute Error (MAD) do CV dos inliers')
xlabel('limiar de erro ε')
ylabel('CV dos inliers')

%% PLOT DAS CORRESPONDENCIAS FILTRADAS NO RANSAC

% Mostra as correspondencias filtradas entre o par de imagens.
desenhar_correspondencias(pts_A_sift_ori, pts_B_sift_ori, ...
    ind_pares, ok, imArgb, imBrgb) ;

%% MOSAICO COM REGISTRO DE A -> B ANTES/DEPOIS DA FILTRAGEM COM RANSAC

% Filtrar pares de indices obtidos na saida RANSAC com o vetor logico ok
ind_bons = ind_pares(ok,:) ;

% Calcular H_ruim com todas as correspondencias 
H_ruim = calcular_Homografia(pts_A_sift_ori, pts_B_sift_ori, ind_pares') ;

% Calcular H_boa com correspondencias filtradas
H_boa  = calcular_Homografia(pts_A_sift_ori, pts_B_sift_ori, ind_bons') ;

% Construir mosaico RGB (resultado do registro) antes e apos RANSAC 
mosaico_bom = desenhar_mosaico_RGB(imArgb, imBrgb, H_boa) ;
mosaico_ruim = desenhar_mosaico_RGB(imArgb, imBrgb, H_ruim) ;

% Mostrar mosaicos
subplot(1,2,1), imshow(mosaico_bom) % RANSAC ON
title('Registro de B em A - RANSAC ON')

subplot(1,2,2), imshow(mosaico_ruim) % RANSAC OFF
title('Registro de B em A - RANSAC OFF')