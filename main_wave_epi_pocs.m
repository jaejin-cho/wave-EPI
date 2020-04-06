% addpath / load data

addpath util;
load('data.mat');

[nx,ny,nz,nc] = size(sens);

%% size parameter

params      = [];
params.nx   = nx;
params.ny   = ny;
params.nc   = nc;
params.ns   = ns;
params.rz   = rz;

% deblur matrix
sms_blur            = ones([nx,ny,rz,nc,ns,2]);
ry_eff              = ry / ns;

for ss = 1:rz
    for irz = 1:FOV_SHIFT*ry_eff 
        sms_blur(:,irz:FOV_SHIFT*ry_eff:end,ss,:,:,:,:,:,:) = ...
                 exp( (ss-1).*sqrt(-1).*(irz-1).* 2.*pi./ FOV_SHIFT ./ry_eff);
    end
end

% low-rank constraints
winSize     = [7,7];
lambda      = 0.75;
keep        = 1:floor(lambda*prod(winSize));


%% image reconstruction - epi

im_epi = zeros(nx,ny,nz,1,ns);

tic;
for slc = 1:nz/rz    
    ind_z           = slc:nz/rz:nz;
    sns_lsqr        = repmat(sens(:,:,ind_z,:),[1,1,1,1,ns,2]);    
    kdata_epi_lsqr 	= kdata_epi(:,:,slc,:,:,:,:);    
    psf_epi_lsqr    = psf_epi(:,:,ind_z,:,:,:);    
    im_epi_lsqr 	= im_epi(:,:,ind_z,:,:);
    
    params.mask     = repmat(abs(sum(sum(sum(kdata_epi_lsqr,8),7),3)) ~= 0,[1,1,rz,1,1,1]);
  
    t_epi_k      	= 1;
    x_epi_k         = im_epi_lsqr;
            
    while (1)        
        im_epi_prev     = im_epi_lsqr;

        % SENSE - EPI
        im_epi_coil  	= repmat(  im_epi_lsqr , [1,1,1,nc,1,2]) .* sns_lsqr; % repmat(sns_lsqr,[1,1,1,1,1,2]);
        kd_epi_coil     = mfft(mfft(im_epi_coil,1) .* psf_epi_lsqr,2) .* sms_blur ;        
        im_epi_coil     = kconsistency(mifft2(kd_epi_coil), kdata_epi_lsqr ,params);
        im_epi_coil     = mifft(mifft(mfft2(im_epi_coil) .* conj(sms_blur),2) .* conj(psf_epi_lsqr) ,1);
        im_epi_lsqr     = mean( sum(im_epi_coil .* conj(sns_lsqr),4)  ./ (eps + sum(abs(sns_lsqr).^2, 4))  ,6);

        % low-rank constraints
        for ss = 1:nz            
            Am_epi  = im2row(mfft2(im_epi_lsqr(:,:,ss,:,:,:,:)), winSize);            
            [U_e,S_e,V_e] = svd(Am_epi,'econ');            
            Am_epi  = U_e(:,keep) * S_e(keep,keep) * V_e(:,keep)';            
            k_epi  = Row2im(Am_epi,[nx,ny,1,1,ns,1],winSize);            
            im_epi_lsqr(:,:,ss,:,:,:,:)   = mifft2(k_epi);
        end        
        % fista update
        [t_epi_k,x_epi_k,im_epi_lsqr] = fista_update(t_epi_k,x_epi_k,im_epi_lsqr);
        
        % check break condition
        rmse_epi = rmse(im_epi_prev,im_epi_lsqr);
        
        if rmse_epi < 1
            break;
        end
    end    
    im_epi(:,:,ind_z,:,:) = im_epi_lsqr;    
end
toc;
disp_img(rot90(mean(im_epi(0.25*end+1:0.75*end,:,:,:,:),5)),1e-3);


%% image reconstruction - wav

im_wav = zeros(nx,ny,nz,1,ns);

tic;

for slc = 1:nz/rz    
    ind_z           = slc:nz/rz:nz;    
    sns_lsqr        = repmat(sens(:,:,ind_z,:),[1,1,1,1,ns,2]);    
    kdata_wav_lsqr 	= kdata_wav(:,:,slc,:,:,:,:);
    psf_wav_lsqr    = psf_wav(:,:,ind_z,:,:,:);
    im_wav_lsqr 	= im_wav(:,:,ind_z,:,:);
    
    params.mask     = repmat(abs(sum(sum(sum(kdata_wav_lsqr,8),7),3)) ~= 0,[1,1,rz,1,1,1]);
  
    t_wav_k      	= 1;
    x_wav_k         = im_wav_lsqr;
            
    while (1)        
        im_wav_prev     = im_wav_lsqr;

        % SENSE - Wave-EPI
        im_wav_coil  	= repmat(  im_wav_lsqr , [1,1,1,nc,1,2]) .* sns_lsqr; % repmat(sns_lsqr,[1,1,1,1,1,2]);
        kd_wav_coil     = mfft(mfft(im_wav_coil,1) .* psf_wav_lsqr,2) .* sms_blur ;        
        im_wav_coil     = kconsistency(mifft2(kd_wav_coil), kdata_wav_lsqr ,params);
        im_wav_coil     = mifft(mifft(mfft2(im_wav_coil) .* conj(sms_blur),2) .* conj(psf_wav_lsqr) ,1);
        im_wav_lsqr     = mean( sum(im_wav_coil .* conj(sns_lsqr),4)  ./ (eps + sum(abs(sns_lsqr).^2, 4))  ,6);
        
        % low-rank constraints
        for ss = 1:nz            
            Am_wav  = im2row(mfft2(im_wav_lsqr(:,:,ss,:,:,:,:)), winSize);
            
            [U_w,S_w,V_w] = svd(Am_wav,'econ');
            
            Am_wav  = U_w(:,keep) * S_w(keep,keep) * V_w(:,keep)';
            k_wav  = Row2im(Am_wav,[nx,ny,1,1,ns,1],winSize);            
            im_wav_lsqr(:,:,ss,:,:,:,:)   = mifft2(k_wav);
        end
        
        % fista update
        [t_wav_k,x_wav_k,im_wav_lsqr] = fista_update(t_wav_k,x_wav_k,im_wav_lsqr);
        
        % check break condition
        rmse_wav = rmse(im_wav_prev,im_wav_lsqr);
        
        if rmse_wav < 1
            break;
        end
    end    
    im_wav(:,:,ind_z,:,:) = im_wav_lsqr;    
end
toc;

disp_img(rot90(mean(im_wav(0.25*end+1:0.75*end,:,:,:,:),5)),1e-3);
