function [t_k,x_k,img_update] = fista_update(t_k,x_k,img_update)

    t_kplus1        = (1 + sqrt(1 + 4 * t_k^2)) / 2;
    coef_kneg1      = -(t_k - 1) / t_kplus1;
    coef_k          = (t_kplus1 + t_k - 1) / t_kplus1;
    y_kplus1        = coef_k * img_update + coef_kneg1 * x_k;

    t_k             = t_kplus1;
    x_k             = img_update;
    img_update      = y_kplus1;
                    
end

