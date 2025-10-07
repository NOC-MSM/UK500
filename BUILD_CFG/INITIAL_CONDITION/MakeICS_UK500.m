close all; clear all;
%% MAKE ICS
% Script to generate UK500 ICS for T&S from AMM15 files
% Because the AMM15 grid overlaps with UK500 grid horizontally,the scrip does:
% vertical interp first and then the horizontal interpolation. 

% NOTE: use rn_wd_ref_depth = 0, having and offset in the e3t leads to issues the vertical interpolation. 
% NOTE: you can select 2 method to mask the final outputs in the section
% 'mask out'. One uses the mesh_mask.nc one uses nemo_mask.nc.
% NOTE: you need the Inpaint_nans function 
% NOTE: The code takes about 4h30 to run. If you run it once, it will save some *.mat files after the vertical interpolation
% that you can re-load in line 203 to speed things up the second time, if you want to change stuff in the horizontal interp. 
t= tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       USER PARAMETERS
%% FILES
%child model
dom500 = '<Path_to_UK500_Domain_file>/domain.0m_cfg.nc';
mask500='<Path_to_UK500_MeshMask_File>/mesh_mask.nc';
% nemo mask
bdy_msk = ncread('<Path_To_NemoMask>/nemo_mask.nc','bdy_msk');
bdy_msk=double(bdy_msk);

%parent model
dom15='<Path_to_AMM15_Domain_File>/GEG_SF12.nc'; %file used by ryan for chamfer
mask15='<Path_to_AMM15_MeshMask>/mesh_mask.nc'; % file gen by R. for chamfer

file15T='<Path_to_AMM15input>/CHAMFER_1d_19940501_19940531_25hourm_grid_T.nc'; %'/projectsa/CHAMFER/MakeRestart/AMM15input/CHAMFER_1d_19930701_19930731_25hourm_grid_T.nc';
file15U='<Path_to_AMM15input>/CHAMFER_1d_19940501_19940531_25hourm_grid_U.nc'; %'/projectsa/CHAMFER/MakeRestart/AMM15input/CHAMFER_1d_19930701_19930731_25hourm_grid_U.nc';
file15V='<Path_to_AMM15input>/CHAMFER_1d_19940501_19940531_25hourm_grid_V.nc'; %'/projectsa/CHAMFER/MakeRestart/AMM15input/CHAMFER_1d_19930701_19930731_25hourm_grid_V.nc';

%% SELECT PERIOD to use for generating ICS 
day = 1; %pick first day of file

%% OUTPUT FILE NAME 
filename = '<Path_For_the_Output_File>/UK500_ICS_y1994m05d01.nc'; %name of the ICS file that will be generated

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       CODE - DON'T TOUCH BELOW
%% LOAD VARIABLES 
%CHILD MODEL - UK500
% domain and e3t 
lat500 = ncread(dom500,'nav_lat'); % full uk500 coords
lon500 = ncread(dom500,'nav_lon');
e3t500=ncread(dom500,'e3t_0');
e3u500=ncread(dom500,'e3u_0');
e3v500=ncread(dom500,'e3v_0');

%mask 
mask_t500=double(ncread(mask500,'tmask'));
mask_u500=double(ncread(mask500,'umask'));
mask_v500=double(ncread(mask500,'vmask'));

lon500ms = lon500; 
lon500ms(isnan(mask_t500(:,:,1))) = Inf;% set to Inf to ingnore land point in calulations
lat500ms = lat500; 
lat500ms(isnan(mask_t500(:,:,1))) = Inf;

%PARENT MODEL - AMM15 
% domain and e3t
lat15=ncread(dom15,'nav_lat');
lon15=ncread(dom15,'nav_lon');
e3t15=ncread(dom15,'e3t_0');
e3u15=ncread(dom15,'e3u_0');
e3v15=ncread(dom15,'e3v_0');

%mask
mask_t15=double(ncread(mask15,'tmask')); 
mask_u15=double(ncread(mask15,'umask'));
mask_v15=double(ncread(mask15,'vmask'));

lon15ms = lon15;
lon15ms(isnan(mask_t15(:,:,1))) = Inf;% set to Inf to ingnore land point in calulations ??NANS?? or zeros?
lat15ms = lat15;
lat15ms(isnan(mask_t15(:,:,1))) = Inf;

sal =  ncread(file15T,'vosaline');salt1=sal(:,:,:,day);
clear sal 
tmp = ncread(file15T,'votemper'); tmp1=tmp(:,:,:,day);
clear tmp 

%% DEPTHS 
% child - estimate depths from e3 level thickness
Depth500(:,:,1)=(e3t500(:,:,1)./2).*mask_t500(:,:,1); 
for zz=2:size(e3t500,3)
    Depth500(:,:,zz)=nansum((e3t500(:,:,1:zz-1).*mask_t500(:,:,1:zz-1)),3)+(e3t500(:,:,zz)./2).*mask_t500(:,:,zz);
end
Depth500masked = Depth500.*mask_t500;
Depth500naned = Depth500masked; Depth500naned(Depth500naned==0)=NaN;

% parent - estimate depths from e3 level thickness
Depth15(:,:,1)=(e3t15(:,:,1)./2).*mask_t15(:,:,1); % 
for zz=2:size(e3t15,3)
    Depth15(:,:,zz)=nansum((e3t15(:,:,1:zz-1).*mask_t15(:,:,1:zz-1)),3)+(e3t15(:,:,zz)./2).*mask_t15(:,:,zz);
end
Depth15masked = Depth15.*mask_t15;
Depth15naned = Depth15masked; Depth15naned(Depth15naned==0)=NaN;

%% Load functions 
addpath '/projectsa/ecowind/Restart_Make/Inpaint_nans'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%        fun stuff begins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 3D fields interpolation - vertical then inpaint.  
% VERTICAL INTERPOLATION - 1D at each point
sal_vinterp=zeros(size(Depth500));
tmp_vinterp=zeros(size(Depth500));
matched_ii = zeros(size(lon500));
matched_jj = zeros(size(lon500));
minidx= zeros(size(lon500));

 for i = 1:1668 %50:53 
     disp(['loop i =',num2str(i)])
     for j = 1:2934 %80:82
         %find closest amm15 point to uk500 (points should be overlapping)
         d = sqrt((lon15-lon500(i,j)).^2+(lat15-lat500(i,j)).^2);
         minim = find(d==min(min(d)));
         [ii,jj]=ind2sub(size(lon15),minim);

         %save idx of closest points
         matched_ii(i,j) = ii;
         matched_jj(i,j) = jj;
         minidx(i,j) = minim;

         %interpolate
         if all(squeeze(Depth500(i,j,:))==-25) % if it's a land point in uk500 just set to Nan. you lowered from 0 to -25, so land is now -25. 
             sal_vinterp(i,j,:)=NaN;
             tmp_vinterp(i,j,:)=NaN;
         else
            x=squeeze(Depth15masked(ii,jj,:));
            v=squeeze(salt1(ii,jj,:));
            v1 = squeeze(tmp1(ii,jj,:));
            xq=squeeze(Depth500(i,j,:));

            if  all(squeeze(Depth15masked(ii,jj,:))==0) % if it's a land point in AMM15
                 sal_vinterp(i,j,:)=NaN;
                 tmp_vinterp(i,j,:)=NaN;

            elseif x(end)== x(end-1) % check if there are regions below seafloor level in amm15
                last_nonzero_idx = find(x ~= 0, 1, 'last');
                x= x(1:last_nonzero_idx);
                v=v(1:last_nonzero_idx);%salinity
                vq = interp1(x(~isnan(v)),v(~isnan(v)),xq,'linear',NaN); %extrap for depth where uk500 is deeper than amm15?
                sal_vinterp(i,j,:)=vq;

                v1=v1(1:last_nonzero_idx);%temperature
                vq1 = interp1(x(~isnan(v1)),v1(~isnan(v)),xq,'linear',NaN); %extrap for depth where uk500 is deeper than amm15?
                tmp_vinterp(i,j,:)=vq1;

                % and if the surface is nans, copy the nearest surface value
                if isnan(vq(1))
                    first_nonnan_idx = find(~isnan(vq),1,'first');
                    vq(1:first_nonnan_idx)=vq(first_nonnan_idx);
                    sal_vinterp(i,j,:)=vq;
                end
                if isnan(vq1(1))
                    first_nonnan_idx = find(~isnan(vq1),1,'first');
                    vq1(1:first_nonnan_idx)=vq1(first_nonnan_idx);
                    tmp_vinterp(i,j,:)=vq1;
                end

            else
                vq = interp1(x(~isnan(v)),v(~isnan(v)),xq,'linear',NaN);
                sal_vinterp(i,j,:)=vq;

                vq1 = interp1(x(~isnan(v1)),v1(~isnan(v1)),xq,'linear',NaN);
                tmp_vinterp(i,j,:)=vq1;

                % and if the surface is nans, copy the nearest surface value
                if isnan(vq(1))
                    first_nonnan_idx = find(~isnan(vq),1,'first');
                    vq(1:first_nonnan_idx)=vq(first_nonnan_idx);
                    sal_vinterp(i,j,:)=vq;
                end
                if isnan(vq1(1))
                    first_nonnan_idx = find(~isnan(vq1),1,'first');
                    vq1(1:first_nonnan_idx)=vq1(first_nonnan_idx);
                    tmp_vinterp(i,j,:)=vq1;
                end

            end
         end
     end
 end

save('sal_vinterp.mat', 'sal_vinterp', '-v7.3')
save('tmp_vinterp.mat', 'tmp_vinterp', '-v7.3')
save('matched_ii.mat', 'matched_ii', '-v7.3')
save('matched_jj.mat', 'matched_jj', '-v7.3')
save('minidx.mat', 'minidx', '-v7.3')

%  % test plot
lev=1; 
figure
pcolor(lon500, lat500, sal_vinterp(:,:,8)); shading flat; colorbar; hold on
scatter(lon500(1298,2472), lat500(1298,2472),'rx')
figure
pcolor(lon500, lat500, tmp_vinterp(:,:,1)); shading flat; colorbar; hold on
scatter(lon500(1000,500), lat500(1000,500),'rx')

%% Load the pre-save matrix. if not, re-run above. 
% load('sal_vinterp.mat');
% load('tmp_vinterp.mat');
% load('matched_ii.mat');
% load('matched_jj.mat');
% load('minidx.mat');

% for some reason you get an odd line of zeros at what looks like the
% northen bdy of amm15. Just nan it. 

tmp_vinterp(sal_vinterp==0)=NaN; 
sal_vinterp(sal_vinterp==0)=NaN; 

% figure
% pcolor(nav_lon,nav_lat,tmp_vinterp(:,:,30));shading flat; colorbar;hold on
% scatter(nav_lon(4908917), nav_lat(4908917),'rx')

%% add edges of the domain to avoid inpaint giving odd value near boudaries
for kk=1:size(sal_vinterp,3)
Z = sal_vinterp(:, :, kk);

% Define edge mask (1-pixel edge)
edgeWidth = 1;
edgeMask = false(size(Z));
edgeMask(1:edgeWidth, :) = true;                       % top
edgeMask(end-edgeWidth+1:end, :) = true;               % bottom
edgeMask(:, 1:edgeWidth) = true;                       % left
edgeMask(:, end-edgeWidth+1:end) = true;               % right

% Extract coordinates and values of edge pixels
latEdge = lat500(edgeMask);
lonEdge = lon500(edgeMask);
zEdge = Z(edgeMask);
zind = find(edgeMask==1);

% Make sure Z is 2D, and zind contains valid linear indices
assert(isvector(zind) && all(zind > 0 & zind <= numel(Z)), 'Invalid zind indices');

% Create mask of valid (non-NaN) points
validMask = ~isnan(Z);

% Use bwdist to get index of nearest non-NaN for every point
[~, idxNearest] = bwdist(validMask);

% For each index in zind, find nearest non-NaN value
nearestIdx = idxNearest(zind);
nearestVals = Z(nearestIdx);
Z(zind)=nearestVals;
sal_vinterp(:,:,kk)=Z;
end

% same with temperature (would be more efficent to put them into one loop)
for kk=1:size(sal_vinterp,3)
Z = tmp_vinterp(:, :, kk);

% Define edge mask (1-pixel edge)
edgeWidth = 1;
edgeMask = false(size(Z));
edgeMask(1:edgeWidth, :) = true;                       % top
edgeMask(end-edgeWidth+1:end, :) = true;               % bottom
edgeMask(:, 1:edgeWidth) = true;                       % left
edgeMask(:, end-edgeWidth+1:end) = true;               % right

% Extract coordinates and values of edge pixels
latEdge = lat500(edgeMask);
lonEdge = lon500(edgeMask);
zEdge = Z(edgeMask);
zind = find(edgeMask==1);

% Make sure Z is 2D, and zind contains valid linear indices
assert(isvector(zind) && all(zind > 0 & zind <= numel(Z)), 'Invalid zind indices');

% Create mask of valid (non-NaN) points
validMask = ~isnan(Z);

% Use bwdist to get index of nearest non-NaN for every point
[~, idxNearest] = bwdist(validMask);

% For each index in zind, find nearest non-NaN value
nearestIdx = idxNearest(zind);
nearestVals = Z(nearestIdx);
Z(zind)=nearestVals;
tmp_vinterp(:,:,kk)=Z;
end

%% Inpaint the gaps. 
for kk = 1:size(sal_vinterp,3)

   %flood only for T and S not for velocities
   salinity_3d(:,:,kk)=inpaint_nans(sal_vinterp(:,:,kk),2);
   temperature_3d(:,:,kk)=inpaint_nans(tmp_vinterp(:,:,kk),2);

   %account for shallow domain that may have no values at last depths
   LL=salinity_3d(:,:,kk);
   if isempty(find(LL(:)~=0))
       salinity_3d(:,:,kk)=nan;
       temperature_3d(:,:,kk)=nan;
   end  
end

%% maks out 
% method 1 with mesh_mask (makes more sense, but nemo is outputting it
%weired)
% salinity_3d=salinity_3d.*mask_t500;
% temperature_3d=temperature_3d.*mask_t500;

% method 2 with nemo_mask.nc
salnew = salinity_3d;
tempnew = temperature_3d;
for lev = 1:31
salnew(:,:,lev) = salinity_3d(:,:,lev).*bdy_msk;
tempnew(:,:,lev) = temperature_3d(:,:,lev).*bdy_msk;
end
salinity_3d=salnew;
temperature_3d = tempnew;

figure;
pcolor(lon500, lat500, salinity_3d(:,:,30));shading flat; colorbar; hold on 
%scatter(lon500(1182,2934), lat500(1182,2934),'rx')

figure;
pcolor(lon500, lat500, temperature_3d(:,:,30));shading flat; colorbar;

% figure;
% pcolor(lon500, lat500, tempnew(:,:,30));shading flat; colorbar;

figure;
pcolor(lon500, lat500, mask_t500(:,:,30));shading flat; colorbar; hold on 
%scatter(lon500(1316,2443), lat500(1316,2443),'rx')
scatter(lon500(325,101), lat500(325,101),'rx')


%% TESTs
x=1000; y=500; %uk500 point 

%closest amm15 point
d = sqrt((lon15-lon500(x,y)).^2+(lat15-lat500(x,y)).^2);
minim = find(d==min(min(d)));
[x15,y15]=ind2sub(size(lon15),minim);

figure;
plot(squeeze(temperature_3d(x,y,:)),squeeze(Depth500(x,y,:)),'r.','MarkerSize', 12); hold on
plot(squeeze(tmp1(x15,y15,:)),squeeze(Depth15(x15,y15,:)),'k.','MarkerSize', 12); hold on
legend('depth 500 (rn wd ref depth = 0)','depth amm15')
title('temperature x depth at random point after interpolating')
%xlim([14.5 15])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate Netcdf restart file

% Read levels and set sizes
lev = ncread(dom500,'nav_lev');     % source name; output var will be 'gdep'
x = size(lon500,1);
y = size(lon500,2);
z = length(lev);

if ~isfile(filename)
    % generate all the 'basic' info for the netcdf file
    time = 1;
    ntime = 0;
    e3t_new = e3t500; % ncread('/projectsa/CHAMFER/domain_cfg_UK500_noWAD.nc','e3t_0');

    % --- coordinates (already with final names/dims) ---
    % 2D lon/lat on {'lon','lat'}
    nccreate(filename,'lat', 'Dimensions',{'lon',x,'lat',y})
    ncwrite(filename,'lat',lat500);

    nccreate(filename,'lon', 'Dimensions',{'lon',x,'lat',y})
    ncwrite(filename,'lon',lon500);

    % vertical coord as 'gdep'
    nccreate(filename,'gdep', 'Dimensions',{'gdep',z})
    ncwrite(filename,'gdep',lev);

    % time axis as 'time_counter' (unlimited)
    nccreate(filename,'time_counter', 'Dimensions',{'time_counter',Inf}, 'Datatype','double')
    ncwrite(filename,'time_counter',time);

    % keep 'ntime' on the time axis as before
    nccreate(filename,'ntime', 'Dimensions',{'time_counter',Inf}, 'Datatype','double')
    ncwrite(filename,'ntime',ntime);

    % e3t_* with the final dim names (lon,lat,gdep,time_counter)
    nccreate(filename,'e3t_b', 'Dimensions',{'lon',x,'lat',y,'gdep',z,'time_counter',Inf}, 'Datatype','double')
    ncwrite(filename,'e3t_b',e3t_new);

    nccreate(filename,'e3t_n', 'Dimensions',{'lon',x,'lat',y,'gdep',z,'time_counter',Inf}, 'Datatype','double')
    ncwrite(filename,'e3t_n',e3t_new);
end

% 3D fields written with final CF-like names and dims
nccreate(filename,'votemper', 'Dimensions',{'lon',x,'lat',y,'gdep',z,'time_counter',Inf}, 'Datatype','double')
ncwrite(filename,'votemper',temperature_3d);

nccreate(filename,'vosaline', 'Dimensions',{'lon',x,'lat',y,'gdep',z,'time_counter',Inf}, 'Datatype','double')
ncwrite(filename,'vosaline',salinity_3d);

% add global attribute
ncwriteatt(filename,'/','Note:','Initial conditions made from AMM15 - CHAMFER run.');

%% test plots 1

z500=1; z15=1;
figure
pcolor(salt1(:,:,z15)'); shading flat; colorbar
caxis([0 35])
xlim([500 1200])
ylim([220 1200])
title(sprintf('cropped AMM15 - salinity, (z = %.0f)', z15))

figure
pcolor(salinity_3d(:,:,z500)'); shading flat ; colorbar
caxis([0 35])
title(sprintf('new restart - salinity(z = %.0f)', z500))

figure
pcolor(tmp1(:,:,z15)'); shading flat; colorbar
caxis([0 20])
xlim([500 1200])
ylim([220 1200])
title(sprintf('cropped AMM15 - temperature, (z = %.0f )', z15))

figure
pcolor(temperature_3d(:,:,z500)'); shading flat ; colorbar
caxis([0 20])
title(sprintf('new restart - temperature(z = %.0f )', z500))

%% check run time 
elapsed = toc(t);
fprintf('Tempo trascorso: %.3f s\n', elapsed);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% salinity_3d(1297,2478,:)
Bathy = ncread('/projectsa/CHAMFER/jrule/RemakeUK500/FromArcher/bathy_meter.nc','Bathymetry');

 % figure;
 % pcolor(lon500, lat500, mask_t500(:,:,30));shading flat; colorbar; hold on 
 % %scatter(lon500(1316,2443), lat500(1316,2443),'rx')
 % scatter(lon500(1297,2478), lat500(1297,2478),'rx')
% 
%  figure;
%  pcolor(lon500, lat500, salinity_3d(:,:,1));shading flat; colorbar; hold on 
%  %scatter(lon500(1316,2443), lat500(1316,2443),'rx')
%  scatter(lon500(1297,2478), lat500(1297,2478),'rx')
%  caxis([25 35])
% 
 figure;
 pcolor(lon500, lat500, Bathy);shading flat; colorbar; hold on 
 %scatter(lon500(1316,2443), lat500(1316,2443),'rx')
 %scatter(lon500(1297,2478), lat500(1297,2478),'rx')
 scatter(lon500(8,2190), lat500(8,2190),'rx')
 caxis([0 300])
