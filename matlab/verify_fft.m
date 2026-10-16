% Independent MATLAB FFT/IFFT reference verification.
% fft()/ifft() are used ONLY as independent reference implementations here.
root = fileparts(fileparts(mfilename('fullpath')));
out = fullfile(root, 'results');
files = dir(fullfile(out, 'N*_pattern*.csv'));
rows = zeros(numel(files), 7);
for q = 1:numel(files)
    token = regexp(files(q).name, '^N(\d+)_pattern(\d+)\.csv$', 'tokens', 'once');
    N = str2double(token{1}); p = str2double(token{2});
    T = readtable(fullfile(out, files(q).name));
    x = complex(T.input_re,T.input_im);
    actual = complex(T.fft_re,T.fft_im);
    tx = complex(T.ifft_re,T.ifft_im);
    rx = complex(T.rx_re,T.rx_im);
    reference = fft(x); reference_ifft = ifft(x);
    fft_error = abs(actual-reference);
    ifft_error = abs(tx-reference_ifft);
    rows(q,:) = [N,p,max(fft_error),norm(actual-reference)/norm(reference), ...
                 max(ifft_error),max(abs(rx-x)),sqrt(sum(abs(rx-x).^2)/sum(abs(x).^2))];
    assert(height(T)==N && all(T.index==(0:N-1)'));
    assert(all(ismember(real(x),[-3,-1,1,3])) && all(ismember(imag(x),[-3,-1,1,3])));
    assert(max(fft_error)<1e-10 && max(ifft_error)<1e-10 && max(abs(rx-x))<1e-10);
    B = table(T.index, real(reference), imag(reference), fft_error, ...
        real(reference_ifft), imag(reference_ifft), ifft_error, ...
        'VariableNames', {'index','matlab_fft_re','matlab_fft_im','fft_abs_error', ...
        'matlab_ifft_re','matlab_ifft_im','ifft_abs_error'});
    writetable(B,fullfile(out,sprintf('matlab_N%d_pattern%d.csv',N,p)));
end
rows = sortrows(rows,[1,2]);
S = array2table(rows,'VariableNames',{'N','pattern','fft_max_error','fft_relative_l2', ...
    'ifft_max_error','ofdm_max_error','evm_rms'});
writetable(S,fullfile(out,'matlab_summary.csv'));
disp(S);
f = fopen(fullfile(out,'matlab_version.txt'),'w');
fprintf(f,'%s\n%s\n',version,computer); fclose(f);
fprintf('PASS: %d cases; MATLAB fft/ifft comparison threshold 1e-10.\n',height(S));
