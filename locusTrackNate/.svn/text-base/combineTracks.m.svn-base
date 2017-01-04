function [track_com, err_com] = combineTracks(track_i , track_o, err_i, err_o)

ii = find(~isnan(track_i),1,'last');

ss_i = numel(track_i);
ss_o = numel(track_o);

ss_err = numel(err_o);

track_com(1:ii) = track_i(1:ii);

track_com(ii+1:(ii+ss_o)) = track_o+track_i(ii);

err_com(1:ii) = err_i(1:ii);

err_com(ii+1:(ii+ss_err)) = err_o+err_i(ii);

end