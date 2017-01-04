%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function for making super-resolution construction
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = makeCountImSteppy( data_mask, data, CONST )
%%

DIsum0_guess = data.DIsum0_guess;

% SED2 filter parameters
SED_WINDOW = CONST.findFocusSR.SED_WINDOW;
SED_P      = CONST.findFocusSR.SED_P;


ss = size(data_mask.phase);

% colormap for kymo
jetter = jet(256);


numRegs = data_mask.regs.num_regs;

hh = waitbar( 0, 'Step fitting.' );

for jj = 1:numRegs
    
    waitbar( jj/numRegs, hh );
    
    data.regs(jj).include = false;
    
    if  data.regs(jj).numTrace
        
        data.regs(jj).include = true;
        
        ndisk    = max(data.regs(jj).ndisk);
        
        dI_noise = mean(data.regs(jj).I0_std);
        DI0 = dI_noise*sqrt(ndisk);
        
        %% Show sum image

        mask_ii = (data_mask.regs.regs_label==jj);
        outline = 0.5*ag(and( ~mask_ii, bwmorph( mask_ii, 'dilate' ) ));
        
        if isfield( data, 'sum_im_proc' )
            back  = ag(data.sum_im_proc);
        else
            back  = ag(data.sum_im);
        end
        back1 = ag(data.im1);
        
        color_im  = cat(3, back,  back,  back  + outline );
        color_im1 = cat(3, back1, back1, back1 + outline );
        
        
        % this image has summed intensity
        figure(7);
        clf;
        set(gcf, 'Position', [634 402 1000 600] );
        
        
        % show summed image
        subplot( 2,3,1 );
        imshow( [color_im] );
        hold on;
        title( data.header, 'Interpreter', 'none');
        
        % Show first image
        subplot( 2,3,6 );
        
        imshow( [color_im1] );
        hold on;
        title( data.header, 'Interpreter', 'none');
        
        % Show kymo
        subplot( 2,3,2 );
        imagesc( data.regs(jj).kymo );
        hold on;
        colormap( jetter );
        title( data.header, 'Interpreter', 'none');
        ylabel( 'Position (Pixels)' );
        xlabel( 'Time (Frames)' );
        
        subplot( 2,3,3 );
        hold on;
        set( gca, 'Box', 'on' );
        title( data.header, 'Interpreter', 'none');
        
        %%
        
        n = [data.regs(jj).trace(:).n];
        
        [junk,ord] = sort( n, 'descend' );
        
        cc = {'r','y','g','c','b','m','w'};
        
        
        nc = numel(cc);
        
        
        %% Fit intensity of all traces to a exponential
        
        %
        %         % loop through all traces to compile the summed intensity
        %         ttmax = 0;
        %         Isum_All = zeros( [1,1e6] );
        %         for ii = 1:data.regs(jj).numTrace
        %             % only consider traces with more than five frames
        %             if data.regs(jj).trace(ii).n > 5
        %                 IsumG = [data.regs(jj).trace(ii).IsumG];
        %
        %                 n_ii = numel(IsumG);
        %
        %                 ttmax = max( [ttmax,n_ii] );
        %
        %                 Isum_All(1:n_ii) = Isum_All(1:n_ii) + IsumG;
        %             end
        %         end
        %
        %         tt = 1:ttmax;
        %         Isum_All = Isum_All(tt);
        %
        %
        %         % fit the intensity until it is reduced by a factor of 5
        %         indInt = Isum_All(tt)>Isum_All(1)/5;
        %         indInt = find( ~indInt, 1, 'first' )-1;
        %
        %         tt0 = 1:indInt;
        %
        %         % fit to an exponential
        %         figure(8);
        %         clf;
        %         PP = countFitInt( Isum_All(tt0), tt0, tt );
        %         tau = PP(2);
        %
        %         data.regs(jj).tau = tau;
        %         %disp( ['Decay time ', num2str( tau), ' frames.'] );
        %
        %
        %
        
        
        
        
        
        %% Loop through all traces and fit trace by trace
        for ii = 1:data.regs(jj).numTrace
            
            kk = 0;
            
            
            
            
            data.regs(jj).trace(ii).DIsum0_guess = DIsum0_guess;
            data.regs(jj).trace(ii).include      = true;
            
            % Show intensities
            
            kk = kk+1;
            ic = mod( kk, nc ) + 1;
            
            
            
            
            
            %Isum_ii = [data.regs(jj).trace(ii).IsumG,0,0,0,0,0,0,0,0];
            Isum  = [data.regs(jj).trace(ii).Isum];
            IsumL = [data.regs(jj).trace(ii).IsumL];
            IsumG = [data.regs(jj).trace(ii).IsumG];
            
            
            %                 I0   = [data.regs(jj).trace(ii).I0]*data.regs(jj).ndisk(ii);
            %                 I0_smooth  = smooth( I0, 10 );
            %                 Isum_smooth  = smooth( Isum, 10 );
            %
            
            
            
            
            %% Make Models
            lmin = 1e-1;
            lmax = 1e4;
            
            llguess = DIsum0_guess;
            x = 1.2*(max(Isum)-min(Isum))*(-1:.005:1);
            opt =  optimset('MaxIter',25,'Display','off', 'TolX', 1/10);

            
            %% Make unprocess data model
            data.regs(jj).trace(ii).IsumG_model = ...
                intMakeModel( data.regs(jj).trace(ii).IsumG,...
                lmin, lmax, x );
            
            %% Model Global subtraction IsumG
            %model = steppi3( IsumG', {[],[],0,0}, false, 200 );
            CONST.count.steppi_rad = 3;
            CONST.count.steppi_kmax = [];
            model = steppi_muv( IsumG', {[],[]}, false, [], ...
                CONST.count.steppi_rad,CONST.count.steppi_kmax );
            
            %model = steppi4( IsumG', {[],(makeElectVar(Isum).^-1),0,0}, true );
            model = fftmodel( model, lmin, lmax, false, llguess, inf, x );
            model = intMakeImModel( model, ...
                data.regs(jj).trace(ii).micro_im,...
                data.regs(jj).trace(ii).mask,...
                data.regs(jj).trace(ii).I0,...
                data.regs(jj).trace(ii).dI0,...
                CONST,...
                opt );
            
            % Discrete level model
            IsumGF = model.Xideal;
            data.regs(jj).trace(ii).IsumGF       = IsumGF;
            data.regs(jj).trace(ii).IsumGF_model = model;
            
            if ~isempty( model.ll0 )
                data.regs(jj).trace(ii).DIsum0 = model.ll0;
            else
                data.regs(jj).trace(ii).DIsum0 = [];
            end
            
            data.regs(jj).trace(ii).IsumGF_model_20 = ...
                fftmodel( model, lmin, lmax, false, llguess, 20, x );
            
            data.regs(jj).trace(ii).IsumGF_model_100 = ...
                fftmodel( model, lmin, lmax, false, llguess, 100, x );
            
            data.regs(jj).trace(ii).IsumGF_model_500 = ...
                fftmodel( model, lmin, lmax, false, llguess, 500, x );
            
            % sed model
            data.regs(jj).trace(ii).IsumGs = sed2( ...
                IsumG, 10, 10 );
            data.regs(jj).trace(ii).IsumGs_model = ...
                intMakeModel( data.regs(jj).trace(ii).IsumGs,...
                lmin, lmax, x );
            
            
            %% Model Local subtraction IsumL
            model = steppi3( IsumL', {[],[],0,0}, false );
            model = fftmodel( model, lmin, lmax, false, llguess, inf, x );
            
            % Discrete level model
            IsumLF = model.Xideal;
            data.regs(jj).trace(ii).IsumLF = IsumLF;
            
            % sed model
            data.regs(jj).trace(ii).IsumLs = sed2( ...
                IsumL, 10, 10 );
            data.regs(jj).trace(ii).IsumLs_model = ...
                intMakeModel( data.regs(jj).trace(ii).IsumLs,...
                lmin, lmax, x );
            
            
            %% Model raw Isum
            data.regs(jj).trace(ii).Isums = sed2( ...
                Isum, 10, 10 );
            data.regs(jj).trace(ii).Isums_model = ...
                intMakeModel( data.regs(jj).trace(ii).Isums,...
                lmin, lmax, x );
            
            
            
            
            %% Show Model FFTs
            figure(7);
            subplot(2,3,4 );
            cla;
            
            loglog( data.regs(jj).trace(ii).IsumGF_model_20.ll, ...
                data.regs(jj).trace(ii).IsumGF_model_20.rhop, 'm' );
            hold on;
            loglog( data.regs(jj).trace(ii).IsumGF_model_100.ll, ...
                data.regs(jj).trace(ii).IsumGF_model_100.rhop, 'y' );
            loglog( data.regs(jj).trace(ii).IsumGF_model_500.ll, ...
                data.regs(jj).trace(ii).IsumGF_model_500.rhop, 'r' );
            loglog( data.regs(jj).trace(ii).IsumGF_model.ll, ...
                data.regs(jj).trace(ii).IsumGF_model.rhop, 'g' );
            loglog( data.regs(jj).trace(ii).IsumGs_model.ll, ...
                data.regs(jj).trace(ii).IsumGs_model.rhop, 'b' );
            loglog( data.regs(jj).trace(ii).IsumG_model.ll, ...
                data.regs(jj).trace(ii).IsumG_model.rhop, 'b:' );
            
            yy = [data.regs(jj).trace(ii).IsumGF_model_20.rhop,...
                data.regs(jj).trace(ii).IsumGF_model_100.rhop,...
                data.regs(jj).trace(ii).IsumGF_model_500.rhop,...
                data.regs(jj).trace(ii).IsumGF_model.rhop, ...
                data.regs(jj).trace(ii).IsumGs_model.rhop, ...
                data.regs(jj).trace(ii).IsumG_model.rhop ];
            
            yylim = max(yy)*[1e-3,2];
            axis( [1e2,1e5,yylim] );
             
            % take care of step size here.
            DIsum0 = data.regs(jj).trace(ii).DIsum0;
            if ~isempty( DIsum0 )
                loglog(  DIsum0 + [0,0], yylim, ':w' );
            end
            title( 'Transform Distribution' );
            ylabel( 'Power (AU)' );
            xlabel( '\Delta I (AU)' );

            %% Show stepsize distribution
            figure(7);
            subplot(2,3,5 );
            cla;
            
            semilogy( data.regs(jj).trace(ii).IsumGF_model_20.x,...
                data.regs(jj).trace(ii).IsumGF_model_20.rho,  'm' );
            hold on;
            semilogy( data.regs(jj).trace(ii).IsumGF_model_100.x,...
                data.regs(jj).trace(ii).IsumGF_model_100.rho, 'y' );
            semilogy( data.regs(jj).trace(ii).IsumGF_model_500.x,...
                data.regs(jj).trace(ii).IsumGF_model_500.rho, 'r' );
            semilogy( data.regs(jj).trace(ii).IsumGF_model.x,...
                data.regs(jj).trace(ii).IsumGF_model.rho,     'g' );
            semilogy( data.regs(jj).trace(ii).IsumGs_model.x,...
                data.regs(jj).trace(ii).IsumGs_model.rho,     'b' );
            semilogy( data.regs(jj).trace(ii).IsumG_model.x,...
                data.regs(jj).trace(ii).IsumG_model.rho,      'b:' );
            
            xx = data.regs(jj).trace(ii).IsumG_model.x;
            yy = [data.regs(jj).trace(ii).IsumGF_model_20.rho,...
                data.regs(jj).trace(ii).IsumGF_model_100.rho,...
                data.regs(jj).trace(ii).IsumGF_model_500.rho,...
                data.regs(jj).trace(ii).IsumGF_model.rho, ...
                data.regs(jj).trace(ii).IsumGs_model.rho, ...
                data.regs(jj).trace(ii).IsumG_model.rho ];
            
            axis( [min(xx), max(xx) ,max(yy)*[1e-3,2]] );
            title( 'Stepsize Distribution' );
            ylabel( 'Probability Desnity (AU)' );
            xlabel( '\Delta I (AU)' );
            
            %% Show models for bleaching curve
            figure(7);
            subplot( 2,3,3 );
            cla;
            
            plot( IsumG,  ['.',cc{ic}] );
            hold on;
            plot( IsumGF, ['-',cc{ic}] );
            ylabel( 'Sum Intensity' );
            xlabel( 'Time (Frames)' );
            
            
            
%             %% Fit the exponential decay (Method 1)
%             indInt = IsumG>(IsumG(1)/5);
%             indInt = find( ~indInt, 1, 'first' )-1;
%             
%             tt0 = 1:indInt;
%             PP = countFitInt2( IsumG(tt0), tt0, tau, false );
            
            
            
            %% Fit the exponential decay (Method 2)
            indd = (IsumG>0);
            tt0 = 1:numel(IsumG);
            PP2 = countFitInt( IsumG(indd), tt0(indd), tt0(indd), false );
            
            
            data.regs(jj).trace(ii).IsumExpFit = PP2;
            data.regs(jj).trace(ii).IsumMax = max( IsumG(:) );

                       
            
            
            %% Show kymo plot
            x = data.regs(jj).trace(ii).x;
            y = data.regs(jj).trace(ii).y;
            
            % convert into kymo coords
            theta = -data_mask.regs.props(jj).Orientation*pi/180;
            RR = [cos(theta), sin(theta);-sin(theta), cos(theta)];
            xpp = x - data.regs(jj).xx(1)+1;
            ypp = y - data.regs(jj).yy(1)+1;
            ss_ = [numel(data.regs(jj).yy), numel(data.regs(jj).xx)];
            
            r =  [xpp-ss_(2)/2;ypp-ss_(1)/2];
            
            ss_ = data.regs(jj).kymo_ss;
            rp = RR*r+[ss_(2)/2*ones(1,numel(x));ss_(1)/2*ones(1,numel(x))];
            
            % show kymo track
            subplot( 2,3,2 );
            plot( data.regs(jj).trace(ii).nn/data.skip, rp(1,:), ['.-',cc{ic}] );
            
            % plot position of fluors
            subplot( 2,3,1 );
            plot( mean(x), mean(y), ['.',cc{ic}]);
            text( 3 + mean(x), mean(y), num2str(kk), 'Color', cc{ic} );
            
            
            
        end
        
        %% House Keeping
        figure(7);
        subplot( 2,3,3 );
        xlabel( 'Time (Frames)' );
        ylabel('Filter Intensity (AU)' );
        %set( gca, 'YGrid', 'on' );
        
        
        %%
        
        
        
        
        
    end
    
    %pause;
end

close(hh);


end

%%
function model = intMakeModel( I, lmin, lmax, x )

dI = pairwise( I );
ndI = numel( dI );

yy = hist( dI, x );



[rhop,ll] = fftmod( I, lmin, lmax, false, [] );


model.rho = yy/(x(2)-x(1))/ndI;
model.x = x;


model.rhop = rhop;
model.ll   = ll;
model.lmin = lmin;
model.lmax = lmax;

model.tau  = [];


end

%%
function model = intMakeImModel( model, im, mask_ii, I0, dI0, CONST, opt )
disp_flag = false;
MAX_FOCI  = 5;

nind = numel( model.ind )-1;

for ii = 1:nind
%     if ii==1
%         T1 = model.ind(ii)-1;
%     else
%         T1 = model.ind(ii);
%     end
    
     T1 = model.ind(ii);
    
    T2 = model.ind(ii+1)-1;
        
    mean_im = mean( reshape( cell2mat( im(T1:T2) ),...
        [size(im{1}),T2-T1+1] ), 3 ) - ...
        mean( I0(T1:T2));
    
    model.micro_im{ii} = mean_im;
    
    model.t_mid(ii) = (T1+T2)/2;
        
    I0ii  = mean( I0(T1:T2));
    dI0ii = mean( dI0(T1:T2))/sqrt(T2-T1);
    
    focus(ii) = intCountFitFocus( mean_im, mask_ii, 0, dI0ii, CONST, disp_flag, opt );
end

model.focus = focus;


end
