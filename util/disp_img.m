function disp_img(img, imax)

[nx,ny,nz] = size(img);

img_tmp = zeros(nx,ny*nz);

for t = 1:nz
    ind = (t-1) * ny + 1 : t * ny ;
    img_tmp(:,ind) = img(:,:,t);
end

figure, imagesc( abs(img_tmp) , [0, imax] ), colormap gray, axis off image;

end

