function [brain_area,lipid_ring] = createBrainArea(MRSIStructspectral,origin_x,origin_y,options)
    arguments
        MRSIStructspectral (1,1) struct %For mask generation, PRE-LIPID REMOVAL!
        origin_x (1,:) double
        origin_y (1,:) double
        options.isPlotMask (1,1) logical = false;
        options.isPGAi1 (1,1) logical = false;
        options.MRSIStructspatial (1,1) struct %For mask generation, PRE-LIPID REMOVAL!
    end

    %Mask generation is from the spectral data, and is true for all time
    %points
    format long
    lipidmap = op_CSIintegrate(MRSIStructspectral,0.9,1.8,'mag');
    
    %figure
    %imagesc(lipidmap)
    mask = zeros(size(lipidmap));
    matrix_size = size(mask,1);
    max_length = sqrt((matrix_size/2)^2 + (matrix_size/2)^2);
    distance_to_ring = size(lipidmap,1)/2;
    total_slices= ceil((2*pi) / atan(1/distance_to_ring));
    

    angle = 0;
    for slice = 1:total_slices
        vector_1_x = (max_length) * cos(-angle) + origin_x ; %rcos(theta)+x_shift
        vector_1_y = (max_length) * sin(-angle) + origin_y ;
        angle = angle + (2*pi) / total_slices;

        [cx,cy,c] = improfile(abs(lipidmap),[origin_x,vector_1_x],[origin_y,vector_1_y]);
        cx = cx(~isnan(c));
        cy = cy(~isnan(c));
        c = c(~isnan(c));
        Zscores = zscore(c);
        if any(Zscores > 3)
            Zscores(Zscores < 3) = 0;
            Zscores(Zscores > 3 ) = 1;
        else
            indx = find(Zscores == max(Zscores));
            Zscores = zeros(size(Zscores));
            Zscores(indx) = 1;
        end


        mask_profile = Zscores.*c;
        mask_profile(mask_profile > 0) =1;


        %plot3(cx,cy,mask_profile)

        for x_here = 1:size(cx,1)
            mask(round(cy(x_here)),round(cx(x_here))) = mask_profile(x_here);
        end
    end

    
    [border_y,border_x] = find(mask==1);
    k = convhull(border_x,border_y);


    brain_mask_y = repmat([1:size(mask,1)],size(mask,1),1);
    brain_mask_x = permute(brain_mask_y,[2,1]);

    [in,on] = inpolygon(brain_mask_x,brain_mask_y,border_x(k),border_y(k));
    lipid_mask_x_on = brain_mask_x(on);
    lipid_mask_y_on = brain_mask_y(on);
    lipid_ring = mask;


    brain_mask_x_in = brain_mask_x(in);
    brain_mask_y_in = brain_mask_y(in);
    brain_area = zeros(size(mask));


    for i=1:size(brain_mask_x_in,1)
        brain_area(brain_mask_y_in(i),brain_mask_x_in(i)) =1;
    end

    for i=1:size(lipid_mask_x_on,1)
        brain_area(lipid_mask_y_on(i),lipid_mask_x_on(i)) =0;
    end
    
    brain_area = imfill(brain_area,'holes');
    brain_area_flip = brain_area;
    brain_area_flip(brain_area==1) = 0;
    brain_area_flip(brain_area==0) = 1;
    lipid_ring = lipid_ring .* brain_area_flip;


    [B,~] = bwboundaries(brain_area,'noholes');
    for k = 1:length(B)
        lipid_k = B{k};
    end
    lipid_ring_x = lipid_k(:,1);
    lipid_ring_y = lipid_k(:,2);
    for i=1:size(lipid_ring_x)
        lipid_ring(lipid_ring_x(i),lipid_ring_y(i))=1;
    end

    for x=1:size(mask,1)
        for y=1:size(mask,2)
            if brain_area(x,y)==1 & lipid_ring(x,y)==1
                brain_area(x,y)=0;
            end
        end
    end


    if options.isPlotMask == true
        figure
        imagesc(abs(squeeze(options.MRSIStructspatial.data(1,:,:))))
        hold on
        plot(border_x(k),border_y(k))
        hold on
        scatter(brain_mask_x(in),brain_mask_y(in));
        hold off
        title("Image Data with Border and Brain Voxels")
    end
end
% function [brain_mask, lipid_mask] = createBrainArea(MRSIStructSpectral, origin_x, origin_y, varargin)
% % CREATEBRAINAREA  Build logical brain & lipid masks from a CSI data set.
% %
% %   [brain_mask, lipid_mask] = createBrainArea(specStruct, x0, y0)
% %   [...] = createBrainArea(specStruct, x0, y0, 'Name',Value, ...)
% %
% % INPUTS
% %   specStruct   : FID‑A CSI structure *already* Fourier‑transformed in space & freq.
% %   x0 , y0      : Floating‑point coordinates (in voxel units) of the brain centre.
% %
% % NAME–VALUE OPTIONS
% %   'Plot'       : true  → show QC figure (default = false)
% %   'SpatialStruct': additional CSI struct in the spatial domain (for nicer QC plots)
% %   'ZThresh'    : Z‑score threshold for lipid detection (default = 3)
% %
% % OUTPUTS
% %   brain_mask   : logical matrix (Ny×Nx)  — voxels inside the convex hull.
% %   lipid_mask   : logical matrix (Ny×Nx)  — voxels along the outer lipid ring.
% %
% % NOTES / IMPROVEMENTS OVER ORIGINAL VERSION
% %   • Accepts Name–Value syntax instead of an options struct.
% %   • Works whether or not an extras/averages dimension exists.
% %   • Uses built‑in *bwconvhull* to obtain the convex hull (faster & safer).
% %   • Ensures lipid ring is exactly a 1‑pixel perimeter via *bwperim*.
% %   • All intermediate indices are bounds‑checked.
% % -------------------------------------------------------------------------
% 
% %% ------------------ 0. Parse inputs -----------------------------------------
%     p = inputParser;
%     p.addRequired ( 'MRSIStructSpectral', @(s)isstruct(s)&&isfield(s,'data') );
%     p.addRequired ( 'origin_x', @isnumeric );
%     p.addRequired ( 'origin_y', @isnumeric );
%     p.addParameter( 'Plot',          false, @islogical );
%     p.addParameter( 'SpatialStruct', struct(), @(s)isstruct(s) );
%     p.addParameter( 'ZThresh',       3,     @(x)isnumeric(x)&&x>0 );
%     p.parse(MRSIStructSpectral, origin_x, origin_y, varargin{:});
%     args = p.Results;
% 
% %% ------------------ 1. 2‑D lipid integral map -------------------------------
%     lipid3D = op_CSIintegrate(args.MRSIStructSpectral, 0.9, 1.8, 'mag');  % y×x×extras
%     lipid2D = squeeze(lipid3D(:,:,1));                                    % force 2‑D
%     [Ny, Nx] = size(lipid2D);
% 
% %% ------------------ 2. Set up ray geometry ----------------------------------
%     maxR     = hypot(Nx/2, Ny/2);               % radius to an image corner
%     nSlices  = ceil( 2*pi / atan(1/(max(Nx,Ny)/2)) );
%     angles   = linspace(0, 2*pi, nSlices+1); angles(end) = [];
% 
%     rayMask  = false(Ny, Nx);
% 
% %% ------------------ 3. Ray‑wise lipid voxel detection ------------------------
%     for ang = angles
%         % tip point in 2‑D image coordinates
%         x_tip =  origin_x + maxR*cos(-ang);
%         y_tip =  origin_y + maxR*sin(-ang);
% 
%         [cx, cy, c] = improfile(abs(lipid2D), [origin_x, x_tip], [origin_y, y_tip]);
%         if isempty(c), continue, end
%         valid = isfinite(c(:));  c=c(valid);  cx=cx(valid);  cy=cy(valid);
%         if isempty(c), continue, end
% 
%         zc = zscore(double(c));
%         candidate = zc > args.ZThresh;
%         if ~any(candidate)          % fallback: pick global max
%             [~, idx] = max(c);
%             candidate = false(size(c)); candidate(idx)=true;
%         end
%         % write into mask (bounds‑checked)
%         for k = find(candidate).'
%             xi = min(max(round(cx(k)),1), Nx);
%             yi = min(max(round(cy(k)),1), Ny);
%             rayMask(yi,xi) = true;
%         end
%     end
% 
% %% ------------------ 4. Build lipid & brain masks ----------------------------
%     if ~any(rayMask(:))
%         brain_mask = false(Ny,Nx);
%         lipid_mask = false(Ny,Nx);
%         return
%     end
% 
%     lipid_mask = bwperim(bwconvhull(rayMask));   % 1‑pixel lipid perimeter
%     brain_mask = imfill(~lipid_mask,'holes');    % everything inside hull
%     brain_mask(lipid_mask) = false;             % ensure disjoint sets
% 
% %% ------------------ 5. Optional QC plot -------------------------------------
%     if args.Plot
%         figure('Color','k'); hold on
%         if ~isempty(fieldnames(args.SpatialStruct))
%             imagesc(abs(squeeze(args.SpatialStruct.data(1,:,:))))
%         else
%             imagesc(abs(lipid2D))
%         end
%         colormap(gray); axis image ij off
%         scatter(find(brain_mask), brain_mask(brain_mask), 15, 'g', 'filled');
%         scatter(find(lipid_mask), lipid_mask(lipid_mask), 15, 'r', 'filled');
%         title('createBrainArea – QC','Color','w')
%         legend({'Brain','Lipid'},'TextColor','w');
%     end
% end
