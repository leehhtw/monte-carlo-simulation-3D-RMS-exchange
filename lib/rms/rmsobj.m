classdef rmsobj < handle
    
properties (GetAccess = public, SetAccess = public)
    sig0;
    sig;
    TD;
    dtime; 
    Tstep; 
    NPar; 
    res;
    dx2t; 
    dx4t;
    NParbin;
    bval; 
    bvec; 
    Nbvec;
    iPoints;
    Din; 
    Dex; 
    kappa;
    T2in;
    T2ex;
    iFlag;
    pgse;
    gval;
    DEL;
    del;
    TE;
end % public prop

properties (GetAccess = private, SetAccess = private)
    D_cnt = [1 2 2 1 2 1];
    W_cnt = [1 4 4 6 12 6 4 12 12 4 1 4 6 4 1];
    D_ind = [1 1; 1 2; 1 3; 2 2; 2 3; 3 3];
    W_ind = [1 1 1 1; 1 1 1 2; 1 1 1 3; 1 1 2 2; 1 1 2 3;
            1 1 3 3; 1 2 2 2; 1 2 2 3; 1 2 3 3; 1 3 3 3;
            2 2 2 2; 2 2 2 3; 2 2 3 3; 2 3 3 3; 3 3 3 3];
end % private prop

methods (Access = public)
    function this = rmsobj(root,fbtab)
        if(nargin > 0 && ~isempty(root) && ~isempty(fbtab))
            this.readJob(root,fbtab);
        else
        this.sig0 = 0; this.sig = 0;
        this.TD = 0; this.dtime = 0; this.Tstep = 0; this.NPar = 0; this.res = 0;
        this.dx2t = 0; this.dx4t = 0; this.bval = 0; this.bvec = 0; this.Nbvec = 0;
        this.iPoints = 0; this.Din = 0; this.Dex = 0; this.kappa = 0; this.iFlag = 0;
        end
    end % constructor

    function saveSubstrate(this,filename,fiber,vs)
        % filename: path, name and extension of the file to save
        % fiber   : 3D binary (0/1) volume contaning the fibers/axons (1), empty space (0) 
        % vs      : voxel size specified in nanometers (integer)

        % It has to be uint8 (cpp code reads data as char)
        fb         = uint8(fiber);
        [nx,ny,nz] = size(fiber); 

        % open file .bin
        fileID = fopen(filename,'w');
        % save substrate properties 
        fwrite(fileID,nx,'uint32',0,'n');
        fwrite(fileID,ny,'uint32',0,'n');
        fwrite(fileID,nz,'uint32',0,'n');
        fwrite(fileID,vs,'uint32',0,'n');
        % save substrate
        fwrite(fileID,fb(:),'uint8',0,'n');
        % save file .bin
        fclose(fileID);
    end % save 

    function [fiber,vs] = readSubstrate(this,filename)
        % load parameters of simulation and data
        fID    = fopen(filename,'rb');
        params = fread(fID,4,'uint32');
        nx     = params(1);
        ny     = params(2);
        nz     = params(3);
        vs     = params(4);
        % load substrate
        fbraw  = fread(fID,nx*ny*nz,'uint8');  
        fiber  = reshape(fbraw,[nx ny nz]);
        fclose(fID);   
    end % read

    function this = readJob(this,pathout,pathbtab)
        % Load btable
        btab = load(pathbtab);
        this.bval = btab(:,4);
        this.bvec = btab(:,1:3);            

        % load parameters of simulation and data
        fID       = fopen(pathout,'rb');
        nparams   = fread(fID,1,'double'); 
        params    = fread(fID,(nparams-1),'double');
        nmat      = params(1);
        ncolumns  = params(2);
        ntp       = params(3); % n saved time points
        Rpoints   = fread(fID,(nmat*3),'double');
        nvals     = nmat*ncolumns*ntp;
        Rdata     = fread(fID,nvals,'double');
        fclose(fID);

        % parameters of the simulation
        this.dtime = params(4);
        this.Tstep = params(5);
        this.NPar  = params(6);
        this.Nbvec = params(7);
        this.Din   = params(8);
        this.Dex   = params(9);
        this.kappa = params(10);
        this.iFlag = params(11); 
        this.res   = params(12);
        this.T2in  = params(13);
        this.T2ex  = params(14);

        % Initialization points (in case multiple initializations were used)
        this.iPoints = reshape(Rpoints,[3 nmat])';

        % Simulation Results as 3D Volume
        Results = zeros(ntp,ncolumns,nmat);
        for i = 1:nmat
        iini = (i-1)*ntp*ncolumns + 1;   
        iend = (i-1)*ntp*ncolumns + ntp*ncolumns;
        Results(:,:,i) = reshape(Rdata(iini:iend),[ntp ncolumns]);
        end % i

        % Obtaining diffusion time, cumulants and signal
        this.TD    = zeros(ntp,nmat);
        this.dx2t  = zeros(ntp,6,nmat);
        this.dx4t  = zeros(ntp,15,nmat);
        this.sig0  = zeros(ntp,nmat);
        this.sig   = zeros(ntp,this.Nbvec,nmat);
        for k = 1:nmat
        this.TD(:,k)      = Results(:,1,k);
        this.sig0(:,k)    = Results(:,2,k);
        this.dx2t(:,:,k)  = Results(:,3:8,k).*this.D_cnt./this.sig0(:,k);
        this.dx4t(:,:,k)  = Results(:,9:23,k).*this.W_cnt./this.sig0(:,k);
        this.sig(:,:,k)   = abs(Results(:,24:23+this.Nbvec,k))./this.sig0(:,k);
        end % k
        
        if this.iFlag==4
            Nbin = ncolumns-(1+6+15+1+this.Nbvec);
            this.NParbin = zeros(ntp,Nbin,nmat);
            for k = 1:nmat
                this.NParbin(:,:,k) = Results(:,24+this.Nbvec:end,k);
            end
        end

    end

    function this = readPGSE(this,pathout,pathbtab,pathpgse)
        % Load btable
        btab = load(pathbtab);
        this.bval = btab(:,4);
        this.bvec = btab(:,1:3);
        this.DEL  = btab(:,5);
        this.del  = btab(:,6);
        this.TE   = btab(:,7);
        this.gval = sqrt(this.bval ./ this.del.^2 ./ (this.DEL-this.del/3));

        % load parameters of simulation and data
        fID       = fopen(pathout,'rb');
        nparams   = fread(fID,1,'double'); 
        params    = fread(fID,(nparams-1),'double');
        nmat      = params(1);
        ncolumns  = params(2);
        ntp       = params(3); % n saved time points
        Rpoints   = fread(fID,(nmat*3),'double');
        nvals     = nmat*ncolumns*ntp;
        Rdata     = fread(fID,nvals,'double');
        fclose(fID);

        % parameters of the simulation
        this.dtime = params(4);
        this.Tstep = params(5);
        this.NPar  = params(6);
        this.Nbvec = params(7);
        this.Din   = params(8);
        this.Dex   = params(9);
        this.kappa = params(10);
        this.iFlag = params(11); 
        this.res   = params(12);

        % Initialization points (in case multiple initializations were used)
        this.iPoints = reshape(Rpoints,[3 nmat])';

        % Simulation Results as 3D Volume
        Results = zeros(ntp,ncolumns,nmat);
        for i = 1:nmat
            iini = (i-1)*ntp*ncolumns + 1;   
            iend = (i-1)*ntp*ncolumns + ntp*ncolumns;
            Results(:,:,i) = reshape(Rdata(iini:iend),[ntp ncolumns]);
        end % i

        % Obtaining diffusion time, cumulants and signal
        this.TD    = zeros(ntp,nmat);
        this.dx2t  = zeros(ntp,6,nmat);
        this.dx4t  = zeros(ntp,15,nmat);
        this.sig0  = zeros(ntp,nmat);
        this.sig   = zeros(ntp,this.Nbvec,nmat);
        for k = 1:nmat
            this.TD(:,k)      = Results(:,1,k);
            this.sig0(:,k)    = Results(:,2,k);
            this.dx2t(:,:,k)  = Results(:,3:8,k).*this.D_cnt./this.sig0(:,k);
            this.dx4t(:,:,k)  = Results(:,9:23,k).*this.W_cnt./this.sig0(:,k);
            this.sig(:,:,k)   = abs(Results(:,24:end,k))./this.sig0;
        end % k

        fID       = fopen(pathpgse,'rb');
        this.pgse = fread(fID,numel(this.bval),'double');
        fclose(fID);
    end

%% Analysis
        
    function dt = dki_shell(this,sig,bvec,bval)
        Nt = size(sig,1);
        bvecu = unique(bvec,'row');
        n2 = this.ndir2(bvecu);
        n4 = this.ndir4(bvecu);

        nbvec = size(bvecu,1);
        nbval = numel(unique(bval));
        dt = zeros(Nt,21);
        for i = 1:Nt
            sigi = sig(i,:);
            Di = zeros(nbvec,1);
            Ki = zeros(nbvec,1);
            for j = 1:nbvec
                Ij = ismember(bvec,bvecu(j,:),'rows');
                sigj = sigi(Ij); sigj = sigj(:);
                bvalj = bval(Ij); bvalj = bvalj(:);
                A = [-bvalj 1/6*bvalj.^(2:nbval)];
                X = A\log(sigj + eps);
                Di(j) = X(1);
                Ki(j) = X(2)/X(1)^2;
            end
            dx2g = Di*2.*this.TD(i);
            dx4g = (Ki+3).*dx2g.^2;
            dt(i,1:6) = ( n2\dx2g ).';
            dt(i,7:21) = ( n4\dx4g ).';
        end
    end

    function dt = dti_lls(this,sig,bvec,bval)
        % Prepare b table
        grad = [bvec bval];
        normgrad = sqrt(sum(grad(:, 1:3).^2, 2)); normgrad(normgrad == 0) = 1;
        grad(:, 1:3) = grad(:, 1:3)./repmat(normgrad, [1 3]);

        b = -(grad(:,ones(1,6)*4).*grad(:,this.D_ind(1:6,1)).*grad(:,this.D_ind(1:6,2)))*diag(this.D_cnt);

        % Linear least square
        dt = b\log(sig(:) + eps);

    end

    function dt = dki_lls(this,sig,bvec,bval)
        % Prepare b table
        grad = [bvec bval];
        normgrad = sqrt(sum(grad(:, 1:3).^2, 2)); normgrad(normgrad == 0) = 1;
        grad(:, 1:3) = grad(:, 1:3)./repmat(normgrad, [1 3]);

        b = -(grad(:,ones(1,6)*4).*grad(:,this.D_ind(1:6,1)).*grad(:,this.D_ind(1:6,2)))*diag(this.D_cnt);

        bsqd6 = grad(:, 4).^2 * ones(1,15)/6;
        b = [b, (bsqd6 .* prod(reshape(grad(:,this.W_ind),[],15,4),3))*diag(this.W_cnt)];

        % Linear least square
        dtTmp = b\log(sig + eps).';
        dt = this.W2K(dtTmp);
    end

    function dt = W2K(this, dt)
        % dt = DKI.W2K(dt)
        %
        % Conversion from W to K, W being MD^2*K
        % W is the estimated tensor, K is the actual kurtosis tensor
        D_apprSq = 1./(sum(dt([1 4 6],:),1)/3).^2;
        dt(7:21,:) = dt(7:21,:) .* D_apprSq(ones(15,1),:);
    end

    function [K,D] = akc_mom(this,n)
        n2 = this.ndir2(n);
        n4 = this.ndir4(n);
        x2 = this.dx2t*n2.';
        x4 = this.dx4t*n4.';
        D = x2/2./this.TD;
        K = x4./x2.^2-3;
    end

    function [K,D] = akc_shell(this,dt,n)
        n2 = this.ndir2(n);
        n4 = this.ndir4(n);
        x2 = dt(:,1:6)*n2.';
        x4 = dt(:,7:21)*n4.';
        D = x2/2./this.TD;
        K = x4./x2.^2-3;
    end

    function [K,D] = akc_lls(this,dt,n)
        [K,D] = DKI.akc(dt,n);
        K = K.'; D = D.';
    end

    function n4i = ndir4(this,ni)
        n4i = prod(reshape(ni(:,this.W_ind),[],15,4),3);
    end

    function n2i = ndir2(this,ni)
        n2i = prod(reshape(ni(:,this.D_ind),[],6,2),3);
    end

    function [fb,sk] = bwaxon2cyl(this,BW,un,cv)
        [nx,ny,nz] = size(BW);
        sk = zeros(nz,2);
        Ai = zeros(nz,1);
        for k = 1:nz
            [I,J] = ind2sub([nx,ny],find(BW(:,:,k)));
            sk(k,:) = [mean(I(:)) mean(J(:))] - 0.5;
            Ai(k,1) = nnz(BW(:,:,k));
        end
        sk = (sk-mean(sk,1))*un + mean(sk,1);
        Am = nnz(BW)/nz;
        Ai = (Ai-Am)*cv + Am;
        ri = sqrt(Ai/pi);

        fb = zeros(nx,ny,nz);
        [y,x] = meshgrid(1:ny,1:nx);
        for k = 1:nz
            fb(:,:,k) = (x-sk(k,1)).^2 + (y-sk(k,2)).^2 <= ri(k)^2;
        end
    end

    function [h,TR] = plotaxon(this,BW)
        [f,v] = isosurface(BW,0);
        TR = triangulation(f,v);
        h = trisurf(TR);
    end

    function cm = skeleton(this,BW,vox,span)
        cm = this.centermass(BW,vox);
        cm(:,1) = imgaussfilt(cm(:,1),span);
        cm(:,2) = imgaussfilt(cm(:,2),span);
    end

    function cm = centermass(this,BW,vox)
        if numel(vox)==1
            vox = vox*[1 1 1];
        end
        [nx,ny,~] = size(BW);
        zlist = find(sum(sum(BW,1),2));
        cm = zeros(numel(zlist),3);
        for i = 1:numel(zlist)
            fiberi = BW(:,:,zlist(i));
            [I,J] = ind2sub([nx,ny],find(fiberi));
            cm(i,:) = [mean(I-0.5)*vox(1),mean(J-0.5)*vox(2),(zlist(i)-0.5)*vox(3)];
        end
    end

    function [w0,lambda,ru,lambda2] = undulation(this,cm)
        wx = cm(1:end-1,1)/2+cm(2:end,1)/2;
        wy = cm(1:end-1,2)/2+cm(2:end,2)/2;
        wx = wx-mean(wx);
        wy = wy-mean(wy);
        wkx = 2*ifft([wx; wx(end:-1:1)]);
        wky = 2*ifft([wy; wy(end:-1:1)]);
        wkx = real(wkx);
        wky = real(wky);
        wk = sqrt(wkx.^2+wky.^2);
        wk = wk(2:10+1);
        w0 = sqrt(sum(wk.^2));

        dx = diff(cm,1,1);
        dl = sqrt(sum(dx.^2,2));
        sin2 = sum(dx(:,1:2).^2,2)./sum(dx.^2,2);
        sin2 = sum(dl/sum(dl).*sin2);
        lambda = pi*w0*sqrt(2/sin2);
        lambda2 = pi*w0*2/sqrt(1-sqrt(1-4*sin2));

        ru = 0.5429*sqrt(w0*lambda);
    end
                
end
        
end % class