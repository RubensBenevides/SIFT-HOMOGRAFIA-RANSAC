% Funcao para mostrar pontos correspondentes entre duas imagens.
%     ENTRADA: 
%         Pares de pontos correspondentes na forma de SIFTPoints.
%         prevPoints = pontos na img. A; nextPoints = pontos na img. B).
%         indexPairs = indices dos pontos correspondentes.
%         ok = indices dos pontos correspondentes filtrados por RANSAC.
%         im1 = imagem A; im2 = imagem B.
%     SAIDA:
%         Vazio, apenas desenha o par de imagens com as correspondencias


function [] = desenhar_correspondencias(prevPoints, nextPoints, indexPairs, ok, imA, imB)
    % Transformar siftPoint em matriz de pontos/colunas
    pts_A = prevPoints.Location' ;
    pts_B = nextPoints.Location' ;
    

    % Transpor vetor de indices
    matches = indexPairs' ;
    numMatches = size(matches,2) ; % quantidade de matches

    % Verificar dimensoes das imagens
    dh1 = max(size(imB,1)-size(imA,1), 0) ; % max na direcao horizontal
    dh2 = max(size(imA,1)-size(imB,1), 0) ;
    
    % Desenhar uma figura dividida em duas (uma ao lado da outra)
    figure(1) ; clf ;
    subplot(2,1,1) ;

    % Adicionar um 'padding' (espacamento) nas imagens
    x = padarray(imA,dh1,'post') ; y = padarray(imB,dh2,'post') ;
    imagesc([x y]) ;
    
    % Criar linhas entre correspondencias nao filtradas
    o = size(imA,2) ;  % tamanho da imagem A em X
    line([pts_A(1,matches(1,:));pts_B(1,matches(2,:))+o], ...
         [pts_A(2,matches(1,:));pts_B(2,matches(2,:))]) ;

    % Mostrar titulo e quantidade
    title(sprintf('%d Correspondências', numMatches)) ; axis image off ;
    
    % Segunda parte do subplot (2,1,2) uma imagem com padding  em cima
    % outra imagem com padding em baixo
    subplot(2,1,2) ;
    imagesc([padarray(imA,dh1,'post') padarray(imB,dh2,'post')]) ;

    % Criar linhas entre correspondencias filtradas
    o = size(imA,2) ; % tamanho da imagem A em X
    line([pts_A(1,matches(1,ok));pts_B(1,matches(2,ok))+o], ...
         [pts_A(2,matches(1,ok));pts_B(2,matches(2,ok))]) ;
    
    % Calcular porcentagem omega=inliers/dados
    taxa = 100*sum(ok)/numMatches ;
    
    % Incluir titulo com % de inliers/dados e desenhar imagem
    title(sprintf('%d (%.2f%%) pares-inlier de %d', sum(ok), taxa, numMatches)) ;
    axis image off ;
    drawnow ;
end