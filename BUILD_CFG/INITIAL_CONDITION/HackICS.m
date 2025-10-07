close all 
clear all 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                            USER PARAMETERS
% SET PATHS AND FILENMAMES  

% INPUTS    
file = '<Path_to_Original_ICS>/UK500_ICS2.nc';
vosaline = ncread(file, 'vosaline'); 
votemper = ncread(file,'votemper');

bathyfile = '<Path_to_UK500_Bathymetry>/bathy_meter.nc';
bathy1 = ncread(bathyfile,'Bathymetry');

dom500 = '<Path_to_UK500_Domain>/domain_cfg_UK500.nc';
lat500 = ncread(dom500,'nav_lat'); % full uk500 coords
lon500 = ncread(dom500,'nav_lon');
bathy2 = ncread(dom500,'bathy_meter');

% outputs
fileout = '<Path_for_Output_File>/UK500_ICS2_hacked.nc';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                      CODE to modify variables for UK500

% fins the 3d masked area from salinity
msk = find(vosaline==0); 
figure; pcolor(vosaline(:,:,1)); shading flat; colorbar

% Load a 2D horizontal slice for reference
slice = vosaline(:, :, 1);  % or any level

% Display the slice
figure;
imagesc(slice); axis image; set(gca, 'YDir', 'normal');
title('Draw polygon and double-click to confirm');
colormap jet;
colorbar;

% Step 1: Manually draw a polygon (returns logical mask)
h = drawpolygon('LineWidth',1.5);
mask = poly2mask(h.Position(:,1), h.Position(:,2), size(slice,1), size(slice,2));


% Step 2: Change values within the polygon for all levels
for k = 1:31
    %Modify salinity
    sal = vosaline(:, :, k);
    mean(sal(mask))
    mask_cond = mask & (sal < 34.5);
    sal(mask_cond) = 34.5;  % or assign any value you want, e.g., tmp(mask) = 35
    vosaline(:, :, k) = sal;
    
    %Modify temperature
    % tmp = votemper(:, :, k);
    % %mnr=nanmean(tmp(mask)); %mean of values over the masked region
    % tmp(mask) = mnr;  % same here
    % votemper(:, :, k) = tmp;
end

vosaline(msk)=0;
%votemper(msk)=0; 

figure; pcolor(vosaline(:,:,1)); shading flat; colorbar
%caxis([25 35]

%figure; pcolor(votemper(:,:,1)); shading flat; colorbar

%% copy new variable in updated netcdf 
copyfile(file, fileout);%makes a copy of the original ICS file
ncwrite(fileout, 'vosaline', vosaline); %copies the modified varibable into the copy

%% TEST plot
figure; pcolor(lon500,lat500,bathy1); shading flat; colorbar; hold on 
scatter(lon500(347,1109), lat500(347,1109),'rx') % crahses here with hacked ICS (after almost 2days)
scatter(lon500(374,1159), lat500(374,1159),'kx') % crahses here with hacked ICS (after almost 2days)
%scatter(lon500(1275,2498), lat500(1275,2498),'kx') % crashes here when using 'regular' un-hacked ICS (at time step 1)
caxis([0 300])
