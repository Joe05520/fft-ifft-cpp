root = fileparts(fileparts(mfilename('fullpath')));
resultsDir = fullfile(root, 'results');
outputDir = fullfile(root, 'figures_matlab');
showFigures = usejava('desktop');  % false: save figures without showing figure windows

sizes = [3, 32, 64, 127, 128];
numberOfSizes = numel(sizes);
patternNames = {'P1: Constant 1+j', 'P2: Alternating signs', ...
    'P3: Constellation traversal', 'P4: Seeded 16-QAM'};
blue = [0.141, 0.388, 0.651];
orange = [0.839, 0.455, 0.133];
olive = [0.416, 0.475, 0.231];
colors = [0.55,0.32,0.66; blue; orange; 0.12,0.52,0.48; olive];
required = {'index', 'input_re', 'input_im', 'fft_re', 'fft_im', ...
    'ifft_re', 'ifft_im', 'rx_re', 'rx_im'};

% Check every input before creating output files.
D = cell(numberOfSizes, 4);
for s = 1:numberOfSizes
    for p = 1:4
        path = fullfile(resultsDir, sprintf('N%d_pattern%d.csv', sizes(s), p));
        assert(isfile(path), 'Missing input CSV: %s. Run ./bin/ofdm first.', path);
        T = readtable(path);
        assert(all(ismember(required, T.Properties.VariableNames)), ...
            'Missing required columns: %s', path);
        assert(height(T) == sizes(s) && isequal(T.index, (0:sizes(s)-1)'), ...
            'Unexpected row count or indices: %s', path);
        assert(all(isfinite(T{:, required}), 'all'), 'Non-finite values: %s', path);
        assert(all(ismember(T.input_re, [-3,-1,1,3])) && ...
            all(ismember(T.input_im, [-3,-1,1,3])), 'Invalid 16-QAM input: %s', path);
        D{s,p} = T;
    end
end

% All reference computations use the original C++ input values.
referenceFFT = cell(numberOfSizes,4);
fftError = cell(numberOfSizes,4);
maxFFT = zeros(numberOfSizes,4);
maxIFFT = zeros(numberOfSizes,4);
maxOFDM = zeros(numberOfSizes,4);
symbolErrors = zeros(numberOfSizes,4);
levels = [-3,-1,1,3];
for s = 1:numberOfSizes
    for p = 1:4
        T = D{s,p};
        x = complex(T.input_re, T.input_im);
        actual = complex(T.fft_re, T.fft_im);
        tx = complex(T.ifft_re, T.ifft_im);
        rx = complex(T.rx_re, T.rx_im);
        referenceFFT{s,p} = fft(x);
        fftError{s,p} = abs(actual - referenceFFT{s,p});
        maxFFT(s,p) = max(fftError{s,p});
        maxIFFT(s,p) = max(abs(tx - ifft(x)));
        maxOFDM(s,p) = max(abs(rx - x));
        [~, reIndex] = min(abs(real(rx) - levels), [], 2);
        [~, imIndex] = min(abs(imag(rx) - levels), [], 2);
        detected = complex(reshape(levels(reIndex), [], 1), ...
            reshape(levels(imIndex), [], 1));
        symbolErrors(s,p) = nnz(detected ~= x);
    end
end
assert(all(maxFFT(:) < 1e-10) && all(maxIFFT(:) < 1e-10) && ...
    all(maxOFDM(:) < 1e-10), 'An absolute error exceeds the 1e-10 criterion.');
assert(all(symbolErrors(:) == 0), 'OFDM symbol errors were detected.');

pngDir = fullfile(outputDir, 'png');
figDir = fullfile(outputDir, 'fig');
if ~isfolder(pngDir), mkdir(pngDir); end
if ~isfolder(figDir), mkdir(figDir); end
detailedSizes = [64, 3, 127];
for detailIndex = 1:numel(detailedSizes)
    n = detailedSizes(detailIndex);
    s = find(sizes == n, 1);
    k = (0:n-1)';
    ticks = unique(round(linspace(0, n-1, min(n,5))));
    padding = max(0.1, n/64);
    bounds = [-padding, n-1+padding];
    if n == 64
        fileNames = {'Figure_01_FFT_MATLAB_comparison', ...
            'Figure_02_OFDM_reconstruction', 'Figure_03_FFT_per_bin_error'};
    elseif n == 3
        fileNames = {'Figure_05_N3_FFT_MATLAB_comparison', ...
            'Figure_06_N3_OFDM_reconstruction', 'Figure_07_N3_FFT_per_bin_error'};
    else
        fileNames = {'Figure_08_N127_FFT_MATLAB_comparison', ...
            'Figure_09_N127_OFDM_reconstruction', 'Figure_10_N127_FFT_per_bin_error'};
    end

% FFT comparison for the selected size.
f = newFigure(fileNames{1}, [100,60,1150,1050], showFigures);
t = tiledlayout(f, 4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(t, sprintf('%d-point FFT: C++ implementation vs. MATLAB', n), 'FontSize', 17);
for p = 1:4
    T = D{s,p};
    actual = complex(T.fft_re, T.fft_im);
    reference = referenceFFT{s,p};
    parts = {real(actual), imag(actual)};
    refParts = {real(reference), imag(reference)};
    partNames = {'Real part', 'Imaginary part'};
    for component = 1:2
        ax = nexttile(t);
        h1 = plot(ax, k, parts{component}, '-', 'Color', blue, 'LineWidth', 1.2);
        hold(ax, 'on');
        h2 = plot(ax, k, refParts{component}, 'o', 'Color', orange, ...
            'MarkerSize', 4, 'LineWidth', 0.9);
        title(ax, [patternNames{p}, ' | ', partNames{component}]);
        ylabel(ax, 'Amplitude'); xlim(ax, bounds); xticks(ax, ticks);
        lowerAmplitude = min([0; parts{component}; refParts{component}]);
        upperAmplitude = max([0; parts{component}; refParts{component}]);
        amplitudeSpan = max(upperAmplitude - lowerAmplitude, 1);
        ylim(ax, [lowerAmplitude - 0.1*amplitudeSpan, upperAmplitude + 0.1*amplitudeSpan]);
        if p == 4, xlabel(ax, 'Frequency bin k'); end
        if p == 1 && component == 1
            lg = legend(ax, [h1,h2], {'Custom C++ FFT','MATLAB fft()'}, ...
                'Orientation', 'horizontal');
            lg.Layout.Tile = 'north';
        end
    end
end
saveFigure(f, pngDir, figDir, showFigures);

% OFDM reconstruction for the selected size.
T = D{s,4};
x = complex(T.input_re, T.input_im);
rx = complex(T.rx_re, T.rx_im);
roundTripError = abs(rx-x);
f = newFigure(fileNames{2}, [120,80,1150,800], showFigures);
t = tiledlayout(f, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(t, sprintf('OFDM reconstruction: %d-point seeded 16-QAM', n), 'FontSize', 17);
ax = nexttile(t);
h1 = plot(ax, real(x), imag(x), 'o', 'Color', blue, 'MarkerSize', 10, 'LineWidth', 1.5);
hold(ax, 'on');
h2 = plot(ax, real(rx), imag(rx), 'x', 'Color', orange, 'MarkerSize', 7, 'LineWidth', 1.4);
axis(ax, 'equal'); xlim(ax, [-4,4]); ylim(ax, [-4,4]);
xticks(ax, levels); yticks(ax, levels);
xlabel(ax, 'In-phase component'); ylabel(ax, 'Quadrature component');
title(ax, '(a) 16-QAM constellation');
lg = legend(ax, [h1,h2], {'Original symbols','Recovered symbols'}, 'Orientation', 'horizontal');
lg.Layout.Tile = 'north';
ax = nexttile(t);
plot(ax, k, roundTripError / 1e-15, 'o-', 'Color', blue, 'MarkerSize', 3, 'LineWidth', 1);
xlabel(ax, 'Subcarrier index k'); ylabel(ax, 'Absolute error / 10^{-15}');
title(ax, '(b) Complex reconstruction error'); xlim(ax, bounds); xticks(ax, ticks);
ylim(ax, [0,max(1, max(roundTripError / 1e-15)*1.2)]);
parts = {real(x), imag(x)};
rxParts = {real(rx), imag(rx)};
partNames = {'(c) Real components before / after OFDM', ...
    '(d) Imaginary components before / after OFDM'};
for component = 1:2
    ax = nexttile(t);
    plot(ax, k, parts{component}, '-', 'Color', blue, 'LineWidth', 1);
    hold(ax, 'on');
    plot(ax, k, rxParts{component}, 'x', 'Color', orange, 'MarkerSize', 4);
    title(ax, partNames{component}); xlabel(ax, 'Subcarrier index k');
    ylabel(ax, 'Amplitude'); xlim(ax, bounds); xticks(ax, ticks); ylim(ax, [-3.6,3.6]); yticks(ax, levels);
end
saveFigure(f, pngDir, figDir, showFigures);

% Per-bin FFT error for the selected size.
f = newFigure(fileNames{3}, [140,100,1150,700], showFigures);
t = tiledlayout(f, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(t, sprintf('%d-point FFT: absolute complex error at every bin', n), 'FontSize', 17);
upper = max(1, max(maxFFT(s,:)) / 1e-14 * 1.3);
for p = 1:4
    ax = nexttile(t);
    plot(ax, k, fftError{s,p} / 1e-14, 'o-', 'Color', blue, 'MarkerSize', 3, 'LineWidth', 1);
    title(ax, patternNames{p}); xlabel(ax, 'Frequency bin k');
    ylabel(ax, 'Absolute error / 10^{-14}');
    xlim(ax, bounds); ylim(ax, [0,upper]); xticks(ax, ticks);
    text(ax, 0.03, 0.9, sprintf('Maximum = %.3e', maxFFT(s,p)), 'Units', 'normalized');
end
saveFigure(f, pngDir, figDir, showFigures);

end

% Figure 4: comparison of all required and additional sizes.
f = newFigure('Figure_04_Bonus_error_summary', [160,120,1100,950], showFigures);
t = tiledlayout(f, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(t, 'FFT / IFFT / OFDM errors for N = 3, 32, 64, 127, 128', 'FontSize', 16);
metrics = {maxFFT, maxIFFT, maxOFDM};
scales = [1e-14, 1e-16, 1e-15];
ylabels = {'Maximum error / 10^{-14}', 'Maximum error / 10^{-16}', 'Maximum error / 10^{-15}'};
titles = {'(a) FFT vs. MATLAB', '(b) IFFT vs. MATLAB', '(c) OFDM reconstruction'};
for metric = 1:3
    ax = nexttile(t);
    values = metrics{metric}.' / scales(metric);
    b = bar(ax, 1:4, values, 'grouped');
    hold(ax, 'on');
    for s = 1:numberOfSizes
        b(s).FaceColor = colors(s,:);
        b(s).EdgeColor = 'white';
        zero = values(:,s) == 0;
        if any(zero)
            plot(ax, b(s).XEndPoints(zero), zeros(1,nnz(zero)), '_', ...
                'Color', colors(s,:), 'MarkerSize', 10, 'LineWidth', 2, 'Clipping', 'off');
        end
    end
    top = max(1, max(values(:))*1.25);
    ylim(ax, [0,top]); xlim(ax, [0.5,4.5]); xticks(ax, 1:4);
    xticklabels(ax, {'P1 Constant','P2 Alternating','P3 Traversal','P4 Seeded'});
    ylabel(ax, ylabels{metric}); title(ax, titles{metric});
    for p = 1:4
        if all(values(p,:) == 0)
            text(ax, p, top*0.07, 'all sizes: 0', 'HorizontalAlignment', 'center', 'FontSize', 9);
        end
    end
    if metric == 1
        lg = legend(ax, b, arrayfun(@(n) sprintf('N = %d', n), sizes, 'UniformOutput', false), 'Orientation', 'horizontal');
        lg.Layout.Tile = 'north';
    end
end
xlabel(t, 'Input pattern');
saveFigure(f, pngDir, figDir, showFigures);

fprintf('\nSaved ten MATLAB figures (PNG at 300 dpi and editable FIG).\n');
fprintf('Output folder: %s\n', outputDir);
fprintf('Max FFT error: %.6e\n', max(maxFFT(:)));
fprintf('Max IFFT error: %.6e\n', max(maxIFFT(:)));
fprintf('Max OFDM error: %.6e\n', max(maxOFDM(:)));

function f = newFigure(name, position, showFigures)
    visibility = 'off';
    if showFigures, visibility = 'on'; end
    f = figure('Name', name, 'NumberTitle', 'off', 'Color', 'white', ...
        'Position', position, 'Visible', visibility);
end

function saveFigure(f, pngDir, figDir, showFigures)
    axesList = findall(f, 'Type', 'axes');
    for k = 1:numel(axesList)
        set(axesList(k), 'FontName', 'Arial', 'FontSize', 10, 'Box', 'off', ...
            'YGrid', 'on', 'GridAlpha', 0.15, 'Layer', 'top');
    end
    drawnow;
    exportgraphics(f, fullfile(pngDir, [f.Name,'.png']), 'Resolution', 300);
    savefig(f, fullfile(figDir, [f.Name,'.fig']));
    fprintf('Saved: %s\n', f.Name);
    if ~showFigures, close(f); end
end
