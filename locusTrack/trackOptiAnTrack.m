function tracks_array = trackOptiAnTrack(dirname,spots_flag)
%batch function of running eSpot on several directories

NUM_CELLS = 1000;
tracks_array = cell(1, NUM_CELLS);



dirseperator = filesep;
if(nargin<1 || isempty(dirname))
    
    dirname = '.';
    dirname = [dirname,dirseperator];
else
    if dirname(length(dirname))~=dirseperator
        dirname = [dirname,dirseperator];
    end
end

contents=dir([dirname '*ell*.mat']);


num_im = numel(contents);


h = waitbar( 0, 'Loading Date.');

for i = 1:min(num_im,NUM_CELLS);
    
    waitbar(i/num_im,h,['Loading Date--File: ',num2str(i),'/',num2str(num_im)]);
    
    data = load( [dirname, contents(i).name] );
    
    data = getLocusTracks( data );
    
    tracks_array{i} = data.track;
    
end

close(h);

histStepSize( tracks_array );

end

function data = loaderInternal( filename );
data = load( filename );
end


function data = histStepSize( tracks_array )


steps_x1 = zeros(1,1000000);
steps_x2 = zeros(1,1000000);
steps_y1 = zeros(1,1000000);
steps_y2 = zeros(1,1000000);

count_x1 = 0;
count_x2 = 0;
count_y1 = 0;
count_y2 = 0;

num_files = numel( tracks_array );

for j = 1:num_files
try    
    dx1 = tracks_array{j}.x1(2:end,:)-tracks_array{j}.x1(1:end-1,:);
    dx2 = tracks_array{j}.x2(2:end,:)-tracks_array{j}.x2(1:end-1,:);
    dy1 = tracks_array{j}.y1(2:end,:)-tracks_array{j}.y1(1:end-1,:);
    dy2 = tracks_array{j}.y2(2:end,:)-tracks_array{j}.y2(1:end-1,:);
    
    dx1 = dx1(~isnan(dx1));
    dx2 = dx2(~isnan(dx2));
    dy1 = dy1(~isnan(dy1));
    dy2 = dy2(~isnan(dy2));
    
    nx1 = numel(dx1);
    nx2 = numel(dy2);
    ny1 = numel(dx1);
    ny2 = numel(dy2);
    
    steps_x1(count_x1+(1:nx1)) = dx1;
    steps_x2(count_x2+(1:nx2)) = dx2;
    steps_y1(count_y1+(1:ny1)) = dy1;
    steps_y2(count_y2+(1:ny2)) = dy2;
    
    count_x1 = count_x1+nx1;
    count_x2 = count_x2+nx2;
    count_y1 = count_y1+ny1;
    count_y2 = count_y2+ny2;
catch
    
end
end


steps_x1 = steps_x1(1:count_x1);
steps_x2 = steps_x2(1:count_x2);
steps_y1 = steps_y1(1:count_y1);
steps_y2 = steps_y2(1:count_y2);

data = [];
data.steps_x1 = steps_x1;
data.steps_x2 = steps_x2;
data.steps_y1 = steps_y1;
data.steps_y2 = steps_y2;

xx = -10:.3333:10;

hist_x1 = hist( steps_x1, xx);
hist_x2 = hist( steps_x2, xx);
hist_y1 = hist( steps_y1, xx);
hist_y2 = hist( steps_y2, xx);


x1m = mean(steps_x1);
x2m = mean(steps_x2);
y1m = mean(steps_y1);
y2m = mean(steps_y2);

x1d = sqrt(mean((steps_x1-x1m).^2));
x2d = sqrt(mean((steps_x2-x2m).^2));
y1d = sqrt(mean((steps_y1-y1m).^2));
y2d = sqrt(mean((steps_y2-y2m).^2));

disp(['Green long Axis mean: ', num2str(x1m,3), ' std: ', num2str(x1d,3)]);
disp(['Green short Axis mean: ', num2str(y1m,3), ' std: ', num2str(y1d,3)]);
disp(['Red long Axis mean: ', num2str(x2m,3), ' std: ', num2str(x2d,3)]);
disp(['Red short Axis mean: ', num2str(y2m,3), ' std: ', num2str(y2d,3)]);

clf;
semilogy(xx,hist_x1/sum(hist_x1)/(xx(2)-xx(1)),'g');
hold on;
semilogy(xx,hist_y1/sum(hist_y1)/(xx(2)-xx(1)),'g--');
semilogy(xx,hist_x2/sum(hist_x2)/(xx(2)-xx(1)),'r');
semilogy(xx,hist_y2/sum(hist_y2)/(xx(2)-xx(1)),'r--');

ylabel( 'Probability Density');
xlabel( 'Set size (pixels)');

end