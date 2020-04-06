function img_coil = kconsistency(img_coil,kdata,params)

    img_concat	= zeros([params.nx*params.rz,params.ny,1,params.nc,params.ns,2]); 
    
    % reshape to params.rz*params.nx cases
    for ss = 1:params.rz
        ind = (ss-1)*params.nx + 1 : ss*params.nx;
        img_concat(ind,:,:,:,:,:,:,:) = img_coil(:,:,ss,:,:,:,:,:,:);
    end        
    if mod(params.rz,2) == 0
        img_concat = circshift(img_concat,[-params.nx/2,0,0,0,0,0,0,0,0,0]);
    end        
    kdata_concat = mfft2(img_concat);
    % coparams.nsisteparams.ncy
    kdata_concat(1:params.rz:end,:,:,:,:,:)	= (1-params.mask(:,:,1,:,:,:,:)) .* kdata_concat(1:params.rz:end,:,:,:,:,:) + params.mask(:,:,1,:,:,:,:) .* kdata;

    img_concat  = mifft2(kdata_concat);

    if mod(params.rz,2) == 0
        img_concat = circshift(img_concat,[+params.nx/2,0,0,0,0,0,0,0,0,0]);
    end
    for ss = 1:params.rz
        ind = (ss-1)*params.nx + 1 : ss*params.nx;
        img_coil(:,:,ss,:,:,:,:) = img_concat(ind,:,:,:,:,:,:);
    end

end

