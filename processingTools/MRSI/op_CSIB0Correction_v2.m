% op_CSIB0Correction_v2.m
% Emile Kadalie, Sunnybrook 2025.
%
% USAGE:
% [out, out_w, freqMap, R2Map]=op_CSIB0Correction_v2(in, in_w);
%
% DESCRIPTION:
% Corrects for slight B0 shifts by subtracting water unsuppressed phase data 
% from water suppressed phase data. Also outputs B0 shift map with the slope values 
% of a linear fit of the water unsuppressed phase data. 
%
% INPUTS:
% in   = MRSI struct
% in_w   = MRSI struct water unsuppressed (time domain)
%
% OUTPUTS:
% out          = MRSI water suppr. struct with B0 corrected data
% out_w        = MRSI water unsuppr. struct with B0 corrected data
% freqMap      = 2D array of frequency shift intensity at each specific voxel
% R2Map      = 2D array of fitting rsquared

function [MRSIStruct, MRSIStruct_w, freqMap, R2Map] = op_CSIB0Correction_v2(MRSIStruct, MRSIStruct_w)
    arguments
        MRSIStruct (1, 1) struct
        MRSIStruct_w (1, 1) struct
    end
    %only correct if coils have be combined
    if(getFlags(MRSIStruct, 'addedrcvrs') == 0)
        error("MRSIStruct Error. Add coils first in MRSIStruct");
    end
    if(getFlags(MRSIStruct, 'spatialft') == 0)
        error("Please fourier transform along the spatial dimension");
    end
    if(getFlags(MRSIStruct, 'spectralft') == 0)
        error("Please fourier transform along the spectral dimension");
    end
    if(getFlags(MRSIStruct_w, 'spatialft') == 0)
        error("Please fourier transform water unsuppr. along the spatial dimension");
    end
    if isfield(MRSIStruct,'mask')
        mask = MRSIStruct.mask.brainmasks(:,:,1);
    else
        mask = ones(MRSIStruct.sz(2:3));
    end

    data_w=reshape(MRSIStruct_w.data,MRSIStruct_w.sz(1),[]);
    %figure;
    idx=find(mask>0);

    % Cap the phase-fit window at the REAL (non-zero) FID length. Zero-filling
    % the time domain (e.g. 576 -> 4096) pads with zeros whose phase is
    % meaningless; fitting over them adds nothing and makes this O(Nt^2) double
    % loop ~50x slower. Detect where the FID actually ends from a masked voxel.
    if ~isempty(idx)
        refFid   = abs(fft(ifftshift(data_w(:, idx(1)))));
        lastReal = find(refFid > 1e-6*max(refFid), 1, 'last');
    else
        lastReal = MRSIStruct_w.sz(1);
    end
    if isempty(lastReal), lastReal = MRSIStruct_w.sz(1); end
    max_tpt_range = 10:min(MRSIStruct_w.sz(1), lastReal);

    freqMap=zeros([prod(MRSIStruct_w.sz(2:3)) length(max_tpt_range)]);
    R2Map=zeros([prod(MRSIStruct_w.sz(2:3)) length(max_tpt_range)]);

    % Precompute each masked voxel's FID + unwrapped phase ONCE. The FID does
    % not depend on max_tpt, so recomputing fft/unwrap inside the max_tpt loop
    % (as before) repeated the same transform length(max_tpt_range) times per
    % voxel for nothing.
    nMax   = length(max_tpt_range);
    phiAll = zeros(size(data_w,1), length(idx));   % [Nspec, nVox] unwrapped phase
    for i = 1:length(idx)
        fid_i        = fft(ifftshift(data_w(:,idx(i))));
        phiAll(:,i)  = unwrap(angle(fid_i));
    end

    for m=1:nMax
        max_tpt = max_tpt_range(m);
        for i=1:length(idx)
            try
                [P,S]=polyfit(MRSIStruct_w.spectralTime(1:max_tpt),phiAll(1:max_tpt,i),1);
                freqMap(idx(i),m)=P(1)/(2*pi);
                R2Map(idx(i),m)=S.rsquared;
                % uncomment next part if you want to see voxels that are not
                % fitting well:
                % if S.rsquared<.5
                %     subplot(121); imagesc(squeeze(abs(ccav_w.data(1,:,:)))); colormap gray; hold on;
                %     [x,y]=ind2sub([40 40],i);
                %     rectangle('Position',[y-0.5 x-0.5 1 1],'LineWidth',2,'FaceColor','b','EdgeColor','b')
                %     hold off
                %     subplot(122); plot(ccav_w.spectralTime(1:max_tpt),unwrap(angle(data(1:max_tpt,i))))
                %     pause(1);
                % end
            catch
                continue
            end
        end
    end
    freqMap=reshape(freqMap,[size(MRSIStruct_w.data,2:3) length(max_tpt_range)]);
    R2Map=reshape(R2Map,[size(MRSIStruct_w.data,2:3) length(max_tpt_range)]);

    % Take best B0map w.r.t R2maps
    meanR2=squeeze(sum(R2Map,1:2)/length(idx));
    [R2val,R2idx]=max(meanR2);
    freqMap=freqMap(:,:,R2idx);
    R2Map=R2Map(:,:,R2idx);

    % B0 correction
    % First option:
    data=reshape(MRSIStruct.data,MRSIStruct.sz(1),[]);
    data_w=reshape(MRSIStruct_w.data,MRSIStruct_w.sz(1),[]);
    idx=find(mask>0);

    for i=1:length(idx)
        data_fid_w = fft(ifftshift(data_w(:,idx(i))));
        phase_offset=unwrap(angle(data_fid_w));

        data_fid=fft(ifftshift(data(:,idx(i))));
        data_fid=data_fid.*exp(1i*-phase_offset);
        data(:,idx(i))=fftshift(ifft(data_fid));

        data_fid_w=data_fid_w.*exp(1i*-phase_offset);
        data_w(:,idx(i))=fftshift(ifft(data_fid_w));
    end

    data=reshape(data,MRSIStruct.sz);
    MRSIStruct = setData(MRSIStruct, data);

    data_w=reshape(data_w,MRSIStruct_w.sz);
    MRSIStruct_w = setData(MRSIStruct_w, data_w);

    % Other option:
    %[MRSIStruct, MRSIStruct_w] = op_CSIecc_klose(MRSIStruct, MRSIStruct_w);
  
end

