% function spec = gft1dRealInterpolate(g)
% Produce interpolated spectrogram from 
% real-partitioned GFT

N = length(g);
pars = gft.gft1dRealPartitions(N);
Npar = length(pars);

pardat = [];

fstart = 0;
for i=1:Npar
    fend = pars(i);
    fwidth = fend - fstart;
    pardat(i).Z = g( (fstart:(fstart+fwidth-1)) + 1);
    if fstart < N/2
        fcentre = fstart + floor(fwidth/2);
    else
        % !!! different to windowsFromPars
        fcentre = abs(-N+fstart-1) - floor(fwidth/2);
    end
    pardat(i).fcentre = fcentre;
    pardat(i).time = linspace(1,N,length(pardat(i).Z));
    if i==Npar
        break
    end
    fstart = fend;
end


F = cell(1,Npar);
T = cell(1,Npar);
Z = cell(1,Npar);
for i=1:Npar
    F{i} = pardat(i).fcentre * ones(1,length(pardat(i).Z));
    T{i} = pardat(i).time;
    Z{i} = pardat(i).Z;
end

F = cell2mat(F);
T = cell2mat(T);
Z = cell2mat(Z);

%%
% interpolate time
Z = cell(1,Npar);
for i=1:Npar
    if length(pardat(i).Z)==1
        Z{i} = pardat(i).Z.*ones(1,N);
    else
        Z{i} = interp1(pardat(i).time, pardat(i).Z, 1:N, 'nearest');
    end
end
freqs = [pardat.fcentre];

Z = cell2mat(Z');
% interpolate freqs
spec = zeros(N/2, N);
for i=1:N
    spec(:,i) = interp1(freqs, Z(:,i), 1:(N/2), 'nearest');
end

figure
imagesc(abs(spec))
set(gca,'YDir','normal')
colorbar

%%
% TSI = TriScatteredInterp(F',T',Z','nearest');
% [X Y] = meshgrid(1:(N/2), 1:N);
% spec = TSI(X,Y);
% 
% figure
% imagesc(abs(spec)')
% set(gca,'YDir','normal')
% colorbar

%%
st = gft.gft1d(y,'gaussian');
im = gft.gft1dInterpolate(st);
figure
imagesc(abs(im(:,1:N/2)'));
set(gca,'YDir','normal')
colorbar