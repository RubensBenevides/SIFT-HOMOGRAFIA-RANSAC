% ENTRADAS: 
%           v -> vetor 1x3 
% SAIDAS: 
%           S -> Matriz skew simetrica (anti-simetrica) 3x3 do vetor
function S = vec2skew(v)
    S = [0    -v(3)  v(2);
         v(3)    0  -v(1);
        -v(2)  v(1)    0];