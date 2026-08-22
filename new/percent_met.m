  %% compute_map_accuracy.m  — phantom metabolite-map accuracy vs nominal
  S = load('map.mat');  map  = S.map;     % .Lac/.Cho/.Cr/.Act  (40x40, ~mM)
  C = load('crlb.mat'); crlb = C.crlb;    % same fields, CRLB in %  (999 = failed)

  mets = {'Lac','Cho','Cr','Act'};
  nom.Lac=[60 0]; nom.Cho=[54 18]; nom.Cr=[60 60]; nom.Act=[18 90];  % [inner outer] mM
  TH = 20;                                            % CRLB cutoff (%)

  % --- ROIs (data-driven, no orientation assumptions) ---
  phantom = crlb.Cr < TH & map.Cr > 0;               % Cr present in both -> phantom
  inner   = phantom & crlb.Lac < TH;                 % lactate only inner -> inner bottle
  outer   = phantom & ~inner;
  fprintf('phantom=%d  inner=%d  outer=%d\n\n', nnz(phantom), nnz(inner), nnz(outer));

  % --- (1) ABSOLUTE accuracy ---
  errs=[]; meas_all=[]; nom_all=[];
  roi={inner,outer}; rn={'inner','outer'};
  fprintf('%-5s%-7s%8s%10s%8s\n','Metab','ROI','nominal','measured','%err');
  for m=1:numel(mets), M=mets{m};
    for r=1:2
      sel = roi{r} & crlb.(M)<TH;
      if ~nnz(sel), continue; end
      meas=mean(map.(M)(sel)); nv=nom.(M)(r);
      if nv>0, e=abs(meas-nv)/nv*100; errs(end+1)=e; meas_all(end+1)=meas; nom_all(end+1)=nv;
        fprintf('%-5s%-7s%8.0f%10.1f%7.1f%%\n',M,rn{r},nv,meas,e);
      else, fprintf('%-5s%-7s%8.0f%10.2f%8s\n',M,rn{r},nv,meas,'(skip)'); end
    end
  end
  XX=mean(errs);
  fprintf('\n(1) Absolute: XX = %.1f%% mean abs error (median %.1f%%, max %.1f%%)\n',XX,median(errs),max(errs));

  % --- (2) Scaling-independent: linear fit measured vs nominal ---
  p=polyfit(nom_all,meas_all,1); yf=polyval(p,nom_all);
  R2=1-sum((meas_all-yf).^2)/sum((meas_all-mean(meas_all)).^2);
  fprintf('(2) Linearity: slope=%.2f, R^2=%.3f  (slope~1 & high R^2 = good relative accuracy)\n',p(1),R2);

  % --- (3) Inner:outer CONTRAST recovery (cancels global scale) ---
  fprintf('(3) Inner:outer ratios (measured vs nominal):\n');
  for m=1:numel(mets), M=mets{m};
    if nom.(M)(2)>0
      ri=mean(map.(M)(inner & crlb.(M)<TH)); ro=mean(map.(M)(outer & crlb.(M)<TH));
      fprintf('    %-3s %.2f vs %.2f\n',M,ri/ro,nom.(M)(1)/nom.(M)(2));
    end
  end

  % --- (4) Cr-referenced ratios in the inner bottle (cancels water scaling) ---
  fprintf('(4) Inner-bottle [X]/Cr (measured vs nominal):\n');
  crI = mean(map.Cr(inner & crlb.Cr<TH));
  for m=1:numel(mets), M=mets{m};
    mi=mean(map.(M)(inner & crlb.(M)<TH));
    fprintf('    %-3s/Cr %.2f vs %.2f\n',M,mi/crI,nom.(M)(1)/60);
  end
    fprintf('\nMedian CRLB (%%)  [failed fits CRLB>99%% excluded]\n');
  fprintf('  %-4s %8s %8s %9s\n','Met','inner','outer','phantom');
  for m = 1:numel(mets)
      M = mets{m};
      mi = median(crlb.(M)(inner   & crlb.(M) < 99));
      mo = median(crlb.(M)(outer   & crlb.(M) < 99));
      mp = median(crlb.(M)(phantom & crlb.(M) < 99));
      fprintf('  %-4s %7.1f%% %7.1f%% %8.1f%%\n', M, mi, mo, mp);
  end