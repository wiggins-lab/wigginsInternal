%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function for making super-resolution construction
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = makeCountim( data_mask, data, CONST )
%%

ss = size(data_mask.phase);
jetter = 0.5*jet(256);
jetter = 0.5*ones(size(jetter)) + jetter;


%% Show sum image

% this image has the gaussian amplitudes
figure(1);
clf;
set(gcf, 'Position', [634 402 1000 600] );
    

subplot( 2,3,1 );
tmp = 0.8*ag(~bwmorph(data_mask.mask_cell,'dilate'));
color_im = cat(3, ag(data.sum_im)+tmp, ag(data.sum_im)+tmp, tmp);
color_im1 = cat(3, ag(data.im1)+tmp, ag(data.im1)+tmp, tmp);

imshow( [color_im,color_im1] );
hold on;
title( data.header, 'Interpreter', 'none');

subplot( 2,3,2 );
imshow( data.kymo, [] );
hold on;
colormap( jetter );
title( data.header, 'Interpreter', 'none');

subplot( 2,3,3 );
hold on;
set( gca, 'Box', 'on' );
title( data.header, 'Interpreter', 'none');
%%


% this image has summed intensity
figure(7);
clf;
set(gcf, 'Position', [634 402 1000 600] );


subplot( 2,3,1 );
tmp = 0.8*ag(~bwmorph(data_mask.mask_cell,'dilate'));
color_im  = cat(3, ag(data.sum_im)+tmp, ag(data.sum_im)+tmp, tmp);

imshow( [color_im,color_im1] );
hold on;
title( data.header, 'Interpreter', 'none');

subplot( 2,3,2 );
imshow( data.kymo, [] );
hold on;
colormap( jetter );
title( data.header, 'Interpreter', 'none');

subplot( 2,3,3 );
hold on;
set( gca, 'Box', 'on' );
title( data.header, 'Interpreter', 'none');

%%


n = [data.regs(1).trace(:).n];

[junk,ord] = sort( n, 'descend' );

cc = {'r','y','g','c','b','m','w'};


nc = numel(cc);

% this is the gaussian intensity
I     = [];

% this is the estimated summed intensity from gaussian fit
Is1   = [];

% this is summed intensity using entire cell background
Is2   = [];

% this is summed intensity using local background
Is3   = [];

% filtered I and Is2
If    = [];
Is2f  = [];

% initial values for each trace 
Iinit     = [];
Is1init   = [];
Is2init   = [];
Is3init   = [];

% initial values for each filtered trace 
Ifinit    = [];
Is2finit  = [];

% initial values for each trace 
Iinit0     = [];
Is1init0   = [];
Is2init0   = [];
Is3init0   = [];

% initial values for each filtered trace 
Ifinit0    = [];
Is2finit0  = [];



kk = 0

I0 = 250;

% loop through all traces
for ii = 1:1000
    
    
    % only consider traces with more than five frames
    if data.regs(1).trace(ii).n > 5
        
        % Show intensities

        kk = kk+1;
        ic = mod( kk, nc ) + 1;
        
        
        % these are the Intensity traces for trace ii
        Iii = [data.regs(1).trace(ii).I,0,0,0,0,0,0,0,0];
        
        Is1ii = [data.regs(1).trace(ii).Isum1,0,0,0,0,0,0,0,0];
        Is2ii = [data.regs(1).trace(ii).Isum2,0,0,0,0,0,0,0,0];
        Is3ii = [data.regs(1).trace(ii).Isum3,0,0,0,0,0,0,0,0];

        Ifii   = sed2( Iii, 5, 10);
        Is2fii = sed2( Is2ii, 5, 10);

        % time in frames
        nnii= data.regs(1).trace(ii).nn;
        nnii = [nnii, nnii(end)+(1:8)];
        
        figure(1);
        
        % plot traces
        subplot( 2,3,3 );
        plot( nnii, Iii, ['.',cc{ic}] );
        plot( nnii, Ifii, ['-',cc{ic}] );
        
        figure(7);
        % plot traces
        subplot( 2,3,3 );
        plot( nnii, Is2ii, ['.',cc{ic}] );
        plot( nnii, Is2fii, ['-',cc{ic}] );
        
        % make cum lists for intensities
        I   = [I,  nan(1,5), Iii  ];        
        If  = [If, nan(1,5), Ifii ];
        
        
        Is1  = [ Is1,  nan(1,5), Is1ii];
        Is2  = [ Is2,  nan(1,5), Is2ii];
        Is2f = [ Is2f, nan(1,5), Is2fii];
        Is3  = [ Is3,  nan(1,5), Is3ii];
        
        
        Iinit   = [Iinit,   Iii(1)  ];        
        Ifinit  = [Ifinit,  Ifii(1) ];
        
        
        Is1init  = [ Is1init,   Is1ii(1)];
        Is2init  = [ Is2init,   Is2ii(1)];
        Is2finit = [ Is2finit,  Is2fii(1)];
        Is3init  = [ Is3init,   Is3ii(1)];
        
        if nnii(1) == 0
        Iinit0   = [Iinit0,   Iii(1)  ];        
        Ifinit0  = [Ifinit0,  Ifii(1) ];
        
        
        Is1init0  = [ Is1init0,   Is1ii(1)];
        Is2init0  = [ Is2init0,   Is2ii(1)];
        Is2finit0 = [ Is2finit0,  Is2fii(1)];
        Is3init0  = [ Is3init0,   Is3ii(1)];
        
        
        end
        
        
        x = data.regs(1).trace(ii).x;
        y = data.regs(1).trace(ii).y;

        % convert into kymo coords
        theta = -data_mask.regs.props(1).Orientation*pi/180;
        RR = [cos(theta), sin(theta);-sin(theta), cos(theta)];
        r =  [x-ss(2)/2;y-ss(2)/2];
        rp = RR*r+[ss(2)/2*ones(1,numel(x));ss(2)/2*ones(1,numel(x))];
        
        figure(1);
        % show kymo track
        subplot( 2,3,2 );
        plot( data.regs(1).trace(ii).nn, rp(1,:), ['.-',cc{ic}] );
        
        % plot position of fluors
        subplot( 2,3,1 );
        plot( mean(x), mean(y), ['.r']);
        text( 3 + mean(x), mean(y), num2str(kk) );
        
        figure(7);
        % show kymo track
        subplot( 2,3,2 );
        plot( data.regs(1).trace(ii).nn, rp(1,:), ['.-',cc{ic}] );
        
        % plot position of fluors
        subplot( 2,3,1 );
        plot( mean(x), mean(y), ['.r']);
        text( 3 + mean(x), mean(y), num2str(kk), 'Color', 'r' );
        
        % plot position of fluors
        if nnii(1) == 0
        subplot( 2,3,1 );
        plot( mean(x)+ss(2), mean(y), ['.r']);
        text( 3 + mean(x)+ss(2), mean(y), num2str(kk), 'Color', 'r' );
        end
    end
end

% label axes
figure(7);
subplot( 2,3,3 );
xlabel( 'Time (Frames)' );
ylabel('Filter Intensity (AU)' );
%set( gca, 'YGrid', 'on' );


%% Make Histograms

nI = numel(I);



dI = I'*ones(1,nI);
dI = dI - dI';
dI = dI(:)';
dI = dI(dI~=0);

dIf = If'*ones(1,nI);
dIf = dIf - dIf';
dIf = dIf(:)';
dIf = dIf(dIf~=0);

dIs2 = Is2'*ones(1,nI);
dIs2 = dIs2 - dIs2';
dIs2 = dIs2(:)';
dIs2 = dIs2(dIs2~=0);

dIs2f = Is2f'*ones(1,nI);
dIs2f = dIs2f - dIs2f';
dIs2f = dIs2f(:)';
dIs2f = dIs2f(dIs2f~=0);



I_     = I(    I ~=0   );
If_    = If(   If~=0   );
Is2_   = Is2(  Is2 ~=0 );
Is2f_  = Is2f( Is2f ~=0);


Ix =  -2000:50:2000;
Ix2 = -30000:500:30000;

[y]   = hist( I_(:),  Ix );
[dy ] = hist( dI(:),  Ix );

[yf]  = hist( If_(:), Ix );
[dyf] = hist( dIf(:), Ix );

[y2]   = hist( Is2_(:),    Ix2 );
[dy2]  = hist( dIs2(:),    Ix2 );

[y2f]   = hist( Is2f_(:),  Ix2 );
[dy2f]  = hist( dIs2f(:),  Ix2 );

figure(1)
subplot( 2,3,4 );
semilogy( Ix,y/sum(y),   '.-r' );
hold on;
title( data.header, 'Interpreter', 'none');

%semilogy( Ix,yf/sum(yf), '.-r' );
semilogy( Ix,dy/sum(dy), '.-b' );
%semilogy( Ix,dyf/sum(dyf), '.-b' );

legend( {'I','dI'} );


xlabel( 'Intensity (AU)' );
ylabel( 'Probability' );

%%
figure(7);
subplot( 2,3,4 );

semilogy( Ix2,y2/sum(y2f),   '.-r' );
hold on;
title( data.header, 'Interpreter', 'none');

%semilogy( Ix,yf/sum(yf), '.-r' );
semilogy( Ix2,dy2/sum(dy2f), '.-b' );
%semilogy( Ix,dyf/sum(dyf), '.-b' );

legend( {'I','dI'} );

%% Fit the histograms
figure(7);
subplot( 2,3,5 );
Ps2 = countFitHist( dy2f, Ix2 , 2200 );
title( data.header, 'Interpreter', 'none');

I0s2 = Ps2(1)

%% Fit Histogram
figure(1);
subplot( 2,3,5 );
P = countFitHist( dy, Ix ,300 );
title( data.header, 'Interpreter', 'none');

I0 = P(1)

%% Plot traces end to end

figure(1)
subplot( 2,3,6 );
plot( I/I0, '.b' );
hold on;
title( data.header, 'Interpreter', 'none');

plot( If/I0, 'r-' );


xlabel( 'Time (Frames)' );
ylabel('Intensity (AU)' );
set( gca, 'YGrid', 'on' );

maxtick =  ceil(max( [I,If]/I0 ));
mintick =  floor(min( [I,If]/I0 ));
mintick = -1;
set( gca, 'YTick', mintick:maxtick )

%% Plot traces end to end
I0s__ = data.I0s__;

if isempty( I0s__ )
    I0s__ = I0s2;
end

figure(7)
subplot( 2,3,6 );
plot( Is2/I0s2, '.b' );
hold on;
title( data.header, 'Interpreter', 'none');

plot( Is2f/I0s2, 'r-' );


plot( Is2(1)/I0s__, '.b', 'MarkerSize', 30 );

ind = and( ~isnan(Is2),[1,isnan( Is2(1:end-1))]);
nn = 1:numel(Is2);
plot( nn(ind), Is2(ind)/I0s__, '.g', 'MarkerSize', 18 );


xlabel( 'Time (Frames)' );
ylabel('Intensity (AU)' );
set( gca, 'YGrid', 'on' );

maxtick =  ceil(max( [Is2]/I0s2 ));
mintick =  floor(min( [Is2]/I0s2 ));
mintick = -1;
set( gca, 'YTick', mintick:maxtick )


%%

figure(1);
clf;

plot( Is2/I0s2, '.b' );
hold on;
title( data.header, 'Interpreter', 'none');

plot( Is2f/I0s2, 'r-' );


plot( Is2(1)/I0s__, '.b', 'MarkerSize', 30 );

ind = and( ~isnan(Is2),[1,isnan( Is2(1:end-1))]);
nn = 1:numel(Is2);
plot( nn(ind), Is2(ind)/I0s__, '.g', 'MarkerSize', 18 );


xlabel( 'Time (Frames)' );
ylabel('Intensity (AU)' );
set( gca, 'YGrid', 'on' );

maxtick =  ceil(max( [Is2]/I0s2 ));
mintick =  floor(min( [Is2]/I0s2 ));
mintick = -1;
set( gca, 'YTick', mintick:maxtick )


%% Load all the fit params into data

% gaussian amplitudes for all traces
data.I  = I;

% all differences of I
data.dI = dI;

% smoothed I
data.If = If;

% fit I per fluor
data.I0 = I0;

% Max initial I value for all traces
data.Imax = max(Iinit);

% Max initial If value for all traces
data.Ifmax = max(Ifinit);

% Estimated number of fluors for brightest trace
data.N = data.Imax/data.I0;

% same as above but for summed intensities
data.Is  = Is2;
data.dIs = dIs2;
data.Isf = Is2f;
data.I0s = I0s2;
data.Isinit = Is2init;
data.Isinit0 = Is2init0;

data.Ismax = max(Is2init);
data.Isfmax = max(Is2finit);
data.Ns = data.Ismax/data.I0s;
data.Nsf = data.Isfmax/data.I0s

%%
%figure(7)
%set( gcf, 'Position', [66 1 1855 1121] );

%pause;

%%
% 
% figure(6);
% clf;
% hold on;% 
% Itmp = I/mean(I(~isnan(I)));
% plot( Itmp, 'r.-' );
% 
% %Itmp = If/mean(If(~isnan(If)));
% %plot( Itmp, 'r-' );
% 
% Itmp = Is1/mean(Is1(~isnan(Is1)));
% plot( Itmp, 'b.-' );
% 
% Itmp = Is2/mean(Is2(~isnan(Is2)));
% plot( Itmp, 'g.-' );
% 
% Itmp = Is3/mean(Is3(~isnan(Is3)));
% plot( Itmp, 'y.-' );




 

end