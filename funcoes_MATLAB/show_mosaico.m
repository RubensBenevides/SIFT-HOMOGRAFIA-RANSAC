function [mosaico] = show_mosaico(im1, im2, H)
    
    % Caixa da imagem 2 em coordenadas homogeneas
    box2 = [1  size(im2,2) size(im2,2)  1 ;
            1  1           size(im2,1)  size(im2,1) ;
            1  1           1            1 ] ;
    
    % Aplica homografia inversa na caixa da imagem 2
    box2_ = inv(H) * box2 ;
    
    % Normaliza (obter coordenadas da caixa com o fator de
    % escala da homografia)
    box2_(1,:) = box2_(1,:) ./ box2_(3,:) ; 
    box2_(2,:) = box2_(2,:) ./ box2_(3,:) ;
    
    % Cria um grid com a caixa transformada (im2 -> im1)
    ur = min([1 box2_(1,:)]):max([size(im1,2) box2_(1,:)]) ;
    vr = min([1 box2_(2,:)]):max([size(im1,1) box2_(2,:)]) ;
    [u,v] = meshgrid(ur,vr) ;

    % Aplica interpolacao bilinear na imagem 1
    im1_ = vl_imwbackward(im2double(im1),u,v) ;

    % Aplica homografia nas coordenadas do grid
    z_ = H(3,1) * u + H(3,2) * v + H(3,3) ;
    u_ = (H(1,1) * u + H(1,2) * v + H(1,3)) ./ z_ ; % normaliza x/z
    v_ = (H(2,1) * u + H(2,2) * v + H(2,3)) ./ z_ ; % normaliza y/z

    % Aplica interpolacao bilinear
    im2_ = vl_imwbackward(im2double(im2),u_,v_) ;
    
    % Seta pra 0 (preto) onde for NaN no mosaico
    mass = ~isnan(im1_) + ~isnan(im2_) ;
    im1_(isnan(im1_)) = 0 ;
    im2_(isnan(im2_)) = 0 ;
    mosaico = (im1_ + im2_) ./ mass ;
end