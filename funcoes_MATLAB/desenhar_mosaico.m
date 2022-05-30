% Funcao que desenha o mosaico das imagens sobrepostas utilizando
% interpoalacao bicúbica dos valores da imagem. 
%     ENTRADA: 
%         Par de imagens [A,B] e transformacao H para aplicar na imagem B. 
%         Funcao feita para aplicar o resultado da homografia calculada
%         por RANSAC utilizando pares de correspondencias entre SIFT (A,B).
%         A homografia H registra A(source) em B(target). Nesta funcao ela
%         eh invertida para ser aplicada em B. A saída mantém a imagem A 
%         inalterada, a imagem transformada eh B.
%     SAIDA:
%         Vazio, apenas desenha o par de imagens (A,B) sobrepostas
%         projetadas num plano com origem no sistema da imagem A

function [mosaico] = desenhar_mosaico(imA, imB, H)
    
    % Caixa da imagem 2 em coordenadas homogeneas
    boxB = [1  size(imB,2) size(imB,2)  1 ;
            1  1           size(imB,1)  size(imB,1) ;
            1  1           1            1 ] ;
    
    % Aplica homografia inversa na caixa da imagem 2.
    boxB_ = inv(H)*boxB ;
    
    % Normaliza (obter coordenadas da caixa com o fator de
    % escala da homografia)
    boxB_(1,:) = boxB_(1,:) ./ boxB_(3,:) ; 
    boxB_(2,:) = boxB_(2,:) ./ boxB_(3,:) ;
    
    % Cria um grid com a caixa transformada utilizando os limites da imA
    ur = min([1 boxB_(1,:)]):max([size(imA,2) boxB_(1,:)]) ;
    vr = min([1 boxB_(2,:)]):max([size(imA,1) boxB_(2,:)]) ;
    [u,v] = meshgrid(ur,vr) ;

    % Aplica interpolacao bicubica na im1
    imA_ = interp2(im2double(imA), u, v, 'linear') ;

    % Aplica homografia nas coordenadas do grid
    z_ = H(3,1) * u + H(3,2) * v + H(3,3) ;
    u_ = (H(1,1) * u + H(1,2) * v + H(1,3)) ./ z_ ; % normaliza x/z
    v_ = (H(2,1) * u + H(2,2) * v + H(2,3)) ./ z_ ; % normaliza y/z

    % Aplica interpolacao bicubica na im2
    imB_ = interp2(im2double(imB), u_, v_, 'cubic') ;

    % Define como 0 (preto) onde for NaN no mosaico
    imA_(isnan(imA_)) = 0 ;
    imB_(isnan(imB_)) = 0 ;
    mosaico = (imA_ + imB_);
end