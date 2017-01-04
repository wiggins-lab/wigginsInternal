%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function for making super-resolution construction
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = makeSRim( data_mask, data, CONST )
%%
'hi'

%% Plot the positions of the super-res focus positions
figure(11);
clf;

% Show mean fluorescent image
tmp = data.sum_im;
tmp(tmp<mean(tmp(:))) = mean(tmp(:));
phase  = 0.1*ag(data_mask.phase);

masktmp = logical(data_mask.mask_cell);
outer   = imdilate(masktmp,strel('square',3));
outer(masktmp) = false;

tmp = data.sum_im;
tmpm = mean(tmp(:)) + std(tmp(:));

tmp( tmp < tmpm ) = tmpm;
tmp = tmp - tmpm;

imshow( cat(3,       phase, ...
    0.50*ag(tmp)   + phase, ...
    0.25*ag(outer) + phase));
hold on;


figure(13);
clf;
hold on;

I_MIN        = CONST.findFocusSR.I_MIN;

% Draw circle around focus position 
theta = (0:8)/8*2*pi;
xt = cos( theta );
yt = sin( theta ); 

% parfor ii = 1:numel(data.regs)
%     for jj = 1:data.regs(ii).numTrace
%         ind_ = (data.regs(ii).trace(jj).I>I_MIN);
%         nind = sum(double(ind_));
%         % foci contribute only in the lnegth of the trace is 
%         % longer than CONST.findFocusSR.MIN_TRACE_LEN
%         if nind > 0%CONST.findFocusSR.MIN_TRACE_LEN
%             % Compute mean locus position
%             xm   = mean( data.regs(ii).trace(jj).x(ind_));
%             ym   = mean( data.regs(ii).trace(jj).y(ind_));
% 
%             plot( xm,ym,'g.');
% 
%             % Compute uncertanty in locus position.
%             if nind > 3
%                 xs   = std( data.regs(ii).trace(jj).x(ind_))/sqrt(numel(ind_));
%                 ys   = std( data.regs(ii).trace(jj).y(ind_))/sqrt(numel(ind_));
%                         
%             %plot( xm,ym,'g.','MarkerSize',1);
%                 plot( xm + xs*xt, ym + ys*yt, 'g');
%             end
%         end
%     end
% end

%% Reconstruct super-resolution image

mag   = CONST.findFocusSR.mag;

sum_im = ag(imresize( tmp, mag));
phase  = 0.0*ag(imresize( data_mask.phase, mag));
outer = imresize( outer, mag );

pixelsize = 16/200/mag*1000;
ssim = size(phase);

ss    = size( data_mask.phase ) * mag;
SRim  = zeros( ss );
c2    = mag*CONST.findFocusSR.crop*10;

n_list = zeros(1, 1000);
n_count = 0;

h = waitbar(0,  'working on reconstruction' );

n_nn = 0;
Ifmax = 0;



for ii = 1:numel(data.regs)

    
    
    DIsum0 = data.regs(ii).DIsum0;

    
    waitbar((ii-1)/numel(data.regs), h );

    
    IsMax = 1.6*DIsum0;
    IsMin = 0.4*DIsum0;

    for jj = 1:data.regs(ii).numTrace
        
        
        Isf = data.regs(ii).trace(jj).Isum_f;
        
        ind_ = and( Isf > IsMin,...
                    Isf < IsMax );
                
        ind_ =  and(  ind_(1:end), ...
                and( [ind_(2:end),0],...
                     [0,ind_(1:end-1)] ));
        
        % foci contribute only in the lnegth of the trace is 
        % longer than CONST.findFocusSR.MIN_TRACE_LEN
        
        nind = sum(double(ind_));
        singular_num = 0;
        
        if  nind > 3
            %CONST.findFocusSR.MIN_TRACE_LEN
            
            n_count = n_count + 1;
            n_list(n_count) = nind;
            
            % Compute mean locus position
            xm   = mean( data.regs(ii).trace(jj).x(ind_)-0.5)*mag+0.5;
            ym   = mean( data.regs(ii).trace(jj).y(ind_)-0.5)*mag+0.5;

            
            figure(11);
            plot( data.regs(ii).trace(jj).x(ind_),...
                  data.regs(ii).trace(jj).y(ind_), 'w.-' );
            
              
              figure(13);
              num_jj = numel(data.regs(ii).trace(jj).Isum_f);
              nn = (1:num_jj) + n_nn;
              n_nn = n_nn + num_jj;
              
              plot( nn, Isf/DIsum0, 'r-' );
              plot( nn(ind_), Isf(ind_)/DIsum0, '.w' );
              
              
              Ifmax = max( [ Ifmax, data.regs(ii).trace(jj).Isum_f] );
              
            % Compute uncertanty in locus position. 
            %xs   = std( data.regs(ii).trace(jj).x(ind_))/sqrt(nind)*mag;
            %ys   = std( data.regs(ii).trace(jj).y(ind_))/sqrt(nind)*mag;
            
            C = zeros(2,2);
            C(1,1) = sum(( xm/mag - data.regs(ii).trace(jj).x(ind_) ).^2)/nind^2*mag^2;
            C(1,2) = sum(( xm/mag - data.regs(ii).trace(jj).x(ind_) ).* ...
                         ( ym/mag - data.regs(ii).trace(jj).y(ind_) ))/nind^2*mag^2;
            C(2,1) = C(1,2);
            C(2,2) = sum(( ym/mag - data.regs(ii).trace(jj).y(ind_) ).^2)/nind^2*mag^2;
            
            cc = det(C);
            K = inv(C);
            
            % Set focus intensity by summing Im. Make the integrated
            % intensity Im.
            Im   = sum( data.regs(ii).trace(jj).I(ind_));
            I0   = Im/(2*pi*sqrt(det(C)));
            
            % crop out region corresponding to the focus position
            x0   = round(xm);
            y0   = round(ym);
            
            xmin = max([x0-c2,1]);
            xmax = min([x0+c2,ss(2)]);
            
            ymin = max([y0-c2,1]);
            ymax = min([y0+c2,ss(1)]);
            
            xind = xmin:xmax;
            yind = ymin:ymax;
            
            [X,Y] = meshgrid( xind,yind );
            
            % add locus to summed image. Gaussian width is determined by
            % position uncertainty 
            
            if cc > 0.01
            SRim(yind,xind) = SRim(yind,xind) + log(I0)*exp( - (  K(1,1)*(X-xm).*(X-xm) + ...
                                                           2*K(1,2)*(X-xm).*(Y-ym) + ...
                                                             K(2,2)*(Y-ym).*(Y-ym))/2 );          
            else
              singular_num = singular_num + 1;  
            end
            
            if any(isnan(SRim))
                'hi'
            end
        end
        
        
    end
    
    SR = ag( (SRim) );
comp = cat(3, 0.5*sum_im + phase, SR + phase, phase + ag(outer)*.25 );
figure(2);
clf;
imagesc( pixelsize*(1:ssim(2)), pixelsize*(1:ssim(1)), comp );
axis equal;
hold on;

xlabel('nm');
ylabel('nm');

drawnow;
    
end

close(h);

disp( ['There were ', num2str(singular_num), ' singular entries.'] );


%% Make histogram
figure(12);
[y,x] = hist( n_list, 10 );

semilogy( x, y, 'r.-' );

%% Make composite image and show it

figure(13);
set(gca, 'YGrid', 'on' );
set(gca, 'YTick', [-1:1:ceil(Ifmax/DIsum0)]  );



data.comp = comp;





SR = ag( log(log(500) + SRim) );

figure(14);
clf;
comp = cat(3, 0.5*sum_im + phase, SR + phase, phase + ag(outer)*.25 );
imagesc( pixelsize*(1:ssim(2)), pixelsize*(1:ssim(1)), comp );
axis equal;
hold on;

xlabel('nm');
ylabel('nm');


end



