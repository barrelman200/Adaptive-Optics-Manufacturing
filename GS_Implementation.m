%%  Parameters 
lambda       = 800e-9;     % wavelength: 800 nm
f_obj        = 8.0e-5;     % effective focal length of the objective
NA           = 1.49;       % numerical aperture
magnification= 100;        % magnification of objective lens / total optical train between slm and focal plane

% SLM native specs
Nx           = 4160;       
Ny           = 2464;       
pixelPitch   = 3.74e-6;    

% Target (real-space) pattern in the final plane
numSpotsX    = 2;          
numSpotsY    = 2;
spotSpacing  = 117e-6;     
numIter      = 20;         % Gerchberg-Saxton iterations
% In the interest of computation time, keep the iterations <20; convergence
% typically occurs in <10 iterations anyways.

%% Define single-axis tilt about vertical axis 
alphaDeg     = 0;          % tilt angle in degrees (anticlockwise from normal)
alphaRad     = alphaDeg * pi/180;
tilt_x_phys  = (4*pi / lambda) * sin(alphaRad);  
tilt_x_comp = -tilt_x_phys;
fprintf('Applying physical tilt=%.1f deg => tilt_x=%.2g rad/m\n', alphaDeg, tilt_x_phys);

%% ---- define a circular aperture radius in the SLM plane  ----
% This is just essentially applying a circular mask on the slm to ensure a
% (truncated) gaussian beam output
radius_circ  = 0.0045;  % 4.5 mm
fprintf('Circular aperture radius = %.2f mm\n', radius_circ*1e3);

%% Effective Pixel Pitch & Aperture 
pixelPitch_eff = pixelPitch / magnification;  
Lx_eff = Nx * pixelPitch_eff;   
Ly_eff = Ny * pixelPitch_eff;

fprintf('Effective pixel pitch: %.2f nm\n', pixelPitch_eff*1e9);
fprintf('Effective aperture: %.2f mm x %.2f mm\n', Lx_eff*1e3, Ly_eff*1e3);

% Check if desired offset fits in the NA:
r_max = NA * f_obj;  
fprintf('Max radial extent by NA=%.2f at f=%.2f mm => %.1f um.\n',...
    NA, f_obj*1e3, r_max*1e6);

%% Define Real-Space Targets 
mx = linspace(-(numSpotsX-1)/2, (numSpotsX-1)/2, numSpotsX);
my = linspace(-(numSpotsY-1)/2, (numSpotsY-1)/2, numSpotsY);
[Mgrid, Ngrid] = meshgrid(mx, my);

X_spots = Mgrid * spotSpacing;  
Y_spots = Ngrid * spotSpacing;  

% Check if they exceed lens NA
spotRadius = sqrt(X_spots.^2 + Y_spots.^2);
if any(spotRadius(:) > r_max)
    warning('Some spots lie outside NA=%.2f. They may be clipped', NA);
end

%% Convert Spots to Spatial Frequencies 
% fraunhofer relation
u_target = X_spots / (lambda * f_obj);
v_target = Y_spots / (lambda * f_obj);

%% Frequency Grid (with Effective Pitch) 
fx = ((-Nx/2):(Nx/2-1)) / Lx_eff;
fy = ((-Ny/2):(Ny/2-1)) / Ly_eff;
[FX, FY] = meshgrid(fx, fy);

%% Build Target Field (Fourier Domain) 
targetField = zeros(Ny, Nx);  
spotSize_pix = 4;  
for i = 1:numel(u_target)
    [~, ix] = min(abs(fx - u_target(i)));
    [~, iy] = min(abs(fy - v_target(i)));
    for dy = -spotSize_pix:spotSize_pix
        for dx_ = -spotSize_pix:spotSize_pix
            xx = ix + dx_;
            yy = iy + dy;
            if xx>=1 && xx<=Nx && yy>=1 && yy<=Ny
                r2 = (dx_^2 + dy^2)/(spotSize_pix^2);
                if r2<=1
                    targetField(yy, xx) = exp(-4*r2);
                end
            end
        end
    end
end

%% SLM plane coords & circular mask 
x_slm = ((-Nx/2):(Nx/2-1)) * pixelPitch;  
y_slm = ((-Ny/2):(Ny/2-1)) * pixelPitch;  
[X_slm, Y_slm] = meshgrid(x_slm, y_slm);

R_slm = sqrt(X_slm.^2 + Y_slm.^2);
circMask = (R_slm <= radius_circ);  

%% Gerchberg-Saxton Iterations 
% We'll keep two fields:
%   1) "algoField" => the field the algorithm sees (no tilt)
%   2) "physField" => the actual physically displayed field (with tilt)
% Start them the same.
algoField = exp(1i * 2*pi * rand(Ny, Nx));
err = zeros(1, numIter);

% For easy DC removal:
cx = floor(Nx/2)+1;
cy = floor(Ny/2)+1;

for iter = 1:numIter
    
    % (1) Build the actual physical field = algoField with tilt:
    %     amplitude from algoField, but add tilt ramp in phase
    physPhase = angle(algoField) + tilt_x_comp.*X_slm;  
    physAmp   = abs(algoField);  % or just 1 inside mask, zero outside
    % We'll enforce the circular mask here
    outField = physAmp .* exp(1i*physPhase);
    outField(~circMask) = 0;
    
    % (2) Forward transform => freq domain
    F_field = fftshift( fft2( ifftshift(outField) ) );
    
    % Zero out the DC pixel to remove the bright center 
    F_field(cy, cx) = 0;
    
    % (3) Enforce amplitude = targetField in freq domain
    F_phase  = angle(F_field);
    F_field2 = targetField .* exp(1i * F_phase);
    
    % (4) Inverse transform => back to SLM plane
    backField = fftshift( ifft2( ifftshift(F_field2) ));
    
    % (5) Phase-only => amplitude=1 inside circle, 0 outside
    newPhase = angle(backField);
    newAmp   = ones(size(newPhase));  
    newAmp(~circMask) = 0;
    algoField = newAmp .* exp(1i*newPhase);
    
    % (6) Compute error: compare "targetField" vs. current freq amplitude
    checkF = fftshift( fft2( ifftshift( outField ) ) );
    %   outField is the physically displayed field from step (A)
    err(iter) = image_error( targetField.^2, abs(checkF).^2 );
    
    if mod(iter, 5)==0
        fprintf('Iteration %2d, Error: %e\n', iter, err(iter));
    end
end

finalPhase   = angle(algoField) + tilt_x_comp.*X_slm ; %adding linear phase ramp to the g-s produced phase mask
finalHologram = exp(1i * finalPhase);
finalHologram(~circMask) = 0;
hologram_phase = angle(finalHologram);

%% Reconstruct & Show in the Final (Sample) Plane 
dxp_eff = lambda * f_obj / Lx_eff;
dyp_eff = lambda * f_obj / Ly_eff;

x_recon = ((-Nx/2):(Nx/2-1))*dxp_eff;
y_recon = ((-Ny/2):(Ny/2-1))*dyp_eff;

% Simulate the reconstruction from finalHologram
reconField = fftshift( fft2( ifftshift( finalHologram ) ) );
reconI = abs(reconField).^2;
reconI = reconI / max(reconI(:));

%% ========================= Plotting Sections =========================
figure('Position',[100,100,1500,450]);

% 1) Hologram Phase
subplot(1,3,1);
x_slm_mm = x_slm*1e3; 
y_slm_mm = y_slm*1e3;
imagesc(x_slm_mm,y_slm_mm, mod(hologram_phase,2*pi));
axis image; axis xy;
colormap(parula); colorbar;
xlabel('x (mm)'); ylabel('y (mm)');
title('Final SLM Phase');

% 2) Reconstructed Intensity
subplot(1,3,2);
imagesc(x_recon*1e6, y_recon*1e6, reconI);
axis image; axis xy;
colormap(parula); colorbar;
xlabel('x (\mum)'); ylabel('y (\mum)');
title('Reconstruction');

% 3) Target amplitude in frequency domain
fx_mm = fx*1e-3; 
fy_mm = fy*1e-3;
subplot(1,3,3);
imagesc(fx_mm, fy_mm, targetField);
axis image; axis xy;
colormap(parula); colorbar;
xlabel('f_x (cyc/mm)'); ylabel('f_y (cyc/mm)');
title('Target Amplitude (Fourier Domain)');

%% figure  with linear scale + a zoom around a chosen spot
figure('Position',[100,600,1500,400]);

% (1) Full reconstruction in linear scale
subplot(1,3,1);
imagesc(x_recon*1e6, y_recon*1e6, reconI);
axis image; axis xy;
colormap(parula); colorbar;
xlabel('x (\mum)'); ylabel('y (\mum)');
title('Reconstructed Intensity (linear scale)');

% (2) Zoom on the top-right spot
spotCoord = +58.5;  
[~, iX] = min(abs(x_recon*1e6 - spotCoord));
[~, iY] = min(abs(y_recon*1e6 - spotCoord));

zoomRange = 50;
x_zoom = x_recon(iX-zoomRange : iX+zoomRange)*1e6;
y_zoom = y_recon(iY-zoomRange : iY+zoomRange)*1e6;
roi = reconI(iY-zoomRange : iY+zoomRange, iX-zoomRange : iX+zoomRange);

subplot(1,3,2);
imagesc(x_zoom, y_zoom, roi);
axis image; axis xy;
colormap(parula); colorbar;
xlabel('x (\mum)'); ylabel('y (\mum)');
title('Zoom on top-right spot (linear)');

% (3) A log-scale view
subplot(1,3,3);
imagesc(x_zoom, y_zoom, log10(roi + 1e-10));
axis image; axis xy;
colormap(parula); colorbar;
xlabel('x (\mum)'); ylabel('y (\mum)');
title('Zoom on top-right spot (log_{10} scale)');

%% Plotting error vs iteration for g-s algorithm
figure('Position',[800,100,500,300]);
plot(1:numIter, err, '-o','LineWidth',1.5);
grid on; xlabel('Iteration'); ylabel('Error');
title('G-S Residual Error vs. Iteration');

%% 3D spectral domain
figure;
% Convert to amplitude:
ax = axes('Parent',gcf);
amplitudeSpectrum = abs(reconI);

% Create a 3D surface plot:
surf(ax, fx_mm, fy_mm, amplitudeSpectrum);

shading interp;        
axis tight;            
view(3);               %3D perspective
% Turn on full 3D box style & customize axis colors
ax.BoxStyle = 'full';     
ax.Box = 'on';
ax.Color = 'k';           
ax.XColor = 'k';          
ax.YColor = 'k';          
ax.ZColor = 'k';          
ax.GridColor = 'k';       
colormap(parula);         
colorbar;              
xlabel('x (\mum)');
ylabel('y (\mum)');
zlabel('Amplitude');
title('3D Spectral Intensity of the Reconstructed Field');








