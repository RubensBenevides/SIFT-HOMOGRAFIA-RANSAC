% ENTRADA: duas imagens imA e imB
% SAIDA: homografia H e pares de correspondencias inliers: (H, ok)

function [H, ok] = SIFT_HOMOGRAFIA_RANSAC(imA, imB)
    % Descrever pontos com SIFT
    pts_A_sift = detectSIFTFeatures(imA);
    pts_B_sift = detectSIFTFeatures(imB);
    
    % Extrair apenas pts. relevantes com orientacao.
    [descritores_A, pts_A_sift_ori] = extractFeatures(imA, pts_A_sift);
    [descritores_B, pts_B_sift_ori] = extractFeatures(imB, pts_B_sift);
    
    % Match para obter indices dos pares.
    ind_pares = matchFeatures(descritores_A,descritores_B);
    total_pares = length(ind_pares);
    
    % Selecao das coordenadas na imA que fazem match com imB e vice-versa.
    pts_corr_A = pts_A_sift_ori(ind_pares(:,1), :);
    pts_corr_B = pts_B_sift_ori(ind_pares(:,2), :);
    
    % NORMALIZAR PONTOS E TRANSFORMAR EM COORDENADAS HOMOGENEAS
    
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
    aux_A = sum(sqrt(pts_A_centro(:,1).^2 + pts_A_centro(:,2).^2))/total_pares ;
    aux_B = sum(sqrt(pts_B_centro(:,1).^2 + pts_B_centro(:,2).^2))/total_pares ;
    sA = sqrt(2)/aux_A ;
    sB = sqrt(2)/aux_B ;
    
    % Montar transformacoes que normalizam pts. em A e B.
    TN_A = [sA   0 -sA*txA;   
            0   sA -sA*tyA; 
            0    0      1];
    TN_B = [sB   0 -sB*txB; 
            0   sB -sB*tyB; 
            0    0      1];
    
    % Transformar em coordenadas homogeneas.
    pts_A = [pts_A_xy, ones(total_pares,1)]' ;
    pts_B = [pts_B_xy, ones(total_pares,1)]' ;
    
    % Aplicar Normalizacao.
    pts_A_Norm = TN_A*pts_A ;
    pts_B_Norm = TN_B*pts_B ;
    
    % CALCULAR HOMOGRAFIA COM RANSAC ADAPTATIVO
    % Iniciar RANSAC com nº de iteracoes elevado (k = 1005 significa 
    % capacidade de lidar com ate 74 % de outliers nos pares)
    
    % Inicializar celulas que receberao Homografias, vetores logicos e escores
    % do RANSAC
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
        random_ind_pairs = datasample(1:total_pares, 4) ; 
    
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
        % associado. Neste vetor estao os parametros da Homografia normalizados
        H{t} = reshape(V(:,9),3,3);
        
        % Aplicar Homografia estimada transformando coordenadas de A para B.
        pts_B_Norm_ = H{t}*pts_A_Norm ;
        
        % Normalizar pelo fator de escala da homografia e subtrair das 
        % coordenadas originais as transformadas de A para B.
        e_x = pts_B_Norm_(1,:)./pts_B_Norm_(3,:) - pts_B_Norm(1,:) ; % erro em x: (xB'/zB' - xB)
        e_y = pts_B_Norm_(2,:)./pts_B_Norm_(3,:) - pts_B_Norm(2,:) ; % erro em y: (yB'/zB' - yB)
        
        % Limiar em pixels*escala
        erro_max = 5*mean([sA,sB]) ;
    
        % Vetor logico para informar quantos pontos obedecem a restricao
        ok{t} = (sqrt(e_x.*e_x + e_y.*e_y)) < erro_max ; % Dist. euclidiana
        
        % Dist. de Manhattan, mais rapida, mas menos precisa.
        %ok{t} = (abs(e_x) + abs(e_y)) < erro_max ;  
    
        % Escore de H = cardinalidade dos inliers. Total de 1's em ok
        escore(t) = sum(ok{t}) ;
    
        % Estimar taxa omega (inliers/total_de_pares).
        omega = escore(t)/total_pares ;
      
        % Estimar novo nº de iteracoes do RANSAC.
        iteracoes = round(log(1-0.99)/log(1-omega^4)) ;
    
        % Se o novo nº de iteracoes for menor que o restante, atualize k.
        if (iteracoes < k-t) && (k-t > 0)
            k = iteracoes ; 
            t = 1;
        end
        % Salvar quantas vezes o RANSAC iterou ate a convergencia.
        iteracoes_RANSAC = iteracoes_RANSAC + 1 ;
    end
    % Salvar matriz e indice com maior score dentro da celula H
    [~, best] = max(escore) ;
    H = H{best} ;
   
    % Salvar vetor logico de maior cardinalidade
    ok = ok{best} ;
end