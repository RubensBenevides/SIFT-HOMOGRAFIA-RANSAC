% ENTRADAS: 
%           X1 -> pontos obtidos pelo descritor na imagem 1
%           X2 -> pontos obtidos pelo descritor na imagem 2
% SAIDAS: 
%           H -> Matriz 3x3 de homografia
%           ok -> Vetor logico dos pontos inliers (0 = ourlier, 1 = inlier)
%           score = [0,1] % porcentagem de inliers/pares na imagem, quanto 
%                   menor, mais vai demorarar, pois significa que os dados
%                   estao muito corrompidos por outliers: score =~ 0.

function [ok, H, score] = RANSAC_Adaptativo(X1,X2)
    % Pre alocar arrays e celulas para melhorar desempenho. A quantidade
    % inicial considera o pior caso, apenas 26 % de inliers -> 1005 itera
    H = cell(1,1005) ;
    ok = cell(1,1005) ;
    score = zeros(1,1005) ;

    % RANSAC ADAPTATIVO loop
    k = 1005 ;  t = 0 ; iteracoes_RANSAC = 0 ; 
    while (t <= k)
        t = t + 1 ;
        % Estimara transformacao de homografia selecionando 4 pares de pts.
        subset = datasample(1:numMatches, 4) ;
        A = [] ;
        
        % MODELO DE HOMOGRAFIA com produto tensorial
        for i = subset
            A = cat(1, A, kron(X1(:,i)', vec2skew(X2(:,i)))) ;
        end
        [~,~,V] = svd(A) ;
        
        % A ultima coluna de autovetores eh salva como a matriz homografia
        % 3x3.
        H{t} = reshape(V(:,9),3,3) ;
        
        % Para saber se a estimacao utilizando os 4 pontos aleatorios foi boa
        % calcula-se o score da transformacao H aplicando ela nos pontos A,
        % o erro eh calculado com as distancias entre pontos na imagem B.
        X2_ = H{t} * X1 ;
        du = X2_(1,:)./X2_(3,:) - X2(1,:) ; % erro em x: (x'/z' - x)
        dv = X2_(2,:)./X2_(3,:) - X2(2,:) ; % erro em y: (y'/z' - y)
        
        % Vetor logico para informar pts. inliers
        ok{t} = (abs(du) + abs(dv)) < 6  ;

        % Score de H sera dado pela soma dos erros dos inliers ponderado
        % pela quantidade de pts
        score(t) = sum(ok{t}) ;

        % Estimar taxa omega (inliers/dados)
        omega = score(t)/numMatches ;
        
        % Novo nº de iteracoes do RANSAC
        iteracoes = round(log(1-0.99)/log(1-omega^4)) ;

        % Se o novo nº de iteracoes for menor que o restante, atualize k
        if (iteracoes < k-t) && (k-t > 0)
            k = iteracoes ; 
            t = 1;
        end
        iteracoes_RANSAC = iteracoes_RANSAC + 1 ;
    end
    
    % Varias H foram calculadas, cada uma com seu score, a melhor sera 
    % escolhida, best = indice da melhor matriz H na celula H.
    [score, best] = max(score) ;
    H = H{best} ;
    ok = ok{best} ;

    % Utiliza o melhor vetor logico ok para recuperar a posicao dos pontos
    % em indexPairs que obedecem a restricao
    prevPoints = [f1(1,indexPairs(1,ok)); f1(2,indexPairs(1,ok))] ;
    nextPoints = [f2(1,indexPairs(2,ok)); f2(2,indexPairs(2,ok))] ;
    
    % Transformar pontos em objeto SIFTPoints
    prevPoints = SIFTPoints(prevPoints');
    nextPoints = SIFTPoints(nextPoints');
    fprintf("Iteracoes feitas no RANSAC: %d", iteracoes_RANSAC)
end


