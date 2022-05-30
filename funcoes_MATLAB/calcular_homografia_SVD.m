% Funcao que calcula a homografia com as correspondencias encontradas no 
% SIFT do MATLAB. Utiliza o produto tensorial na construcao da matriz. Essa
% funcao utiliza uma das saidas do RANSAC, onde ha os indices entre todas
% as correspondencias inliers. Todas elas sao utilizadas nesta funcao para
% construir a matriz A de (2*n_linhas por 12 colunas) n = numero de pares
% de correspondencias. Calcula a matriz de A(source) -> B (target).
% Ou seja, H eh a matriz que registra a imagem A em B.

% ENTRADAS:
%           indexPairs = vetor de indices correspondentes com pares na
%           forma de colunas: [pa1,pb1; pa2,pb2; pa3,pb3; ...].
% SAIDAS:   H = matriz de Homografia de norma-2 = 1.

function H = calcular_Homografia(prevPoints, nextPoints, indexPairs)
    
    % Extrair pontos do objeto siftPoints  
    f1 = double(prevPoints.Location') ;
    f2 = double(nextPoints.Location') ;

    % Selecao das coordenadas em f1 que fazem match com f2 e vice-versa
    XaYa = f1(1:2,indexPairs(1,:)) ; XaYa(3,:) = 1 ;
    XbYb = f2(1:2,indexPairs(2,:)) ; XbYb(3,:) = 1 ;
    
    % Montar matriz A
    L = length(XaYa);
    A = [] ;
    for i = 1:5:L
        A = cat(1, A, kron(XaYa(:,i)', vec2skew(XbYb(:,i)))) ;
    end
 
    % Soluciona A por SVD
    [~,~,V] = svd(A) ;
    H = reshape(V(:,9),3,3) ; % Monta matriz H de homografia
end