clear;
close all;
clc;

%% Load input files

Case = 'hole';

if strcmp(Case,'crack')
    filedirectory = 'Crack_Images';
    ref = imread(strcat(filedirectory,'/ref.TIF'));
    ref = double(ref);
    maxref = max(max(ref));
    ref = ref/maxref;
    cur = imread(strcat(filedirectory,'/cur.TIF'));
    cur = double(cur);
    cur = cur/maxref;
    gridType = 'rectangular';
    candidateGridDimensionX = 20;
    candidateGridDimensionY = 20;
    width = 60;
    height = 60;
    curSubimageBuffer = 300;
    imageEdgeBuffer = 400;
    smoothingIterations
elseif strcmp(Case,'hole')
    filedirectory = 'Hole_Plate_Images';
    filename = 'ohtcfrp_';
    ref = read(Tiff(strcat(filedirectory,'/',filename,'01','.tif')));
    ref = double(ref)/256.;
    cur = read(Tiff(strcat(filedirectory,'/',filename,'11','.tif')));
    cur = double(cur)/256.;
    gridType = 'rectangular';
    candidateGridDimensionX = 20;
    candidateGridDimensionY = 40;
    width = 20;
    height = 20;
    curSubimageBuffer = 20;
    imageEdgeBuffer = 75;
elseif strcmp(Case,'holeRandomGrid')
    filedirectory = 'Hole_Plate_Images';
    filename = 'ohtcfrp_';
    ref = read(Tiff(strcat(filedirectory,'/',filename,'01','.tif')));
    ref = double(ref)/256.;
    cur = read(Tiff(strcat(filedirectory,'/',filename,'11','.tif')));
    cur = double(cur)/256.;
    gridType = 'random';
    candidateGridSize = 800;
    width = 20;
    height = 20;
    curSubimageBuffer = 20;
    imageEdgeBuffer = 75;
elseif strcmp(Case,'rotation')
    load('rotation_data.mat')
    candidateGridDimensionX = 10;
    candidateGridDimensionY = 10;
    gridType = 'rectangular';
    width = 20;
    height = 20;
    curSubimageBuffer = 20;
    imageEdgeBuffer = 75;
elseif strcmp(Case, 'arbitraryF')
    load('rotation_data.mat')
    Fxx = 1.01;
    Fxy = -0.05;
    Fyx = -.001;
    Fyy = .99;
    tform = affine2d([Fxx Fxy 0; Fyx Fyy 0; 0 0 1]);
    cur = imwarp(ref,tform);
    gridType = 'rectangular';
    candidateGridDimensionX = 10;
    candidateGridDimensionY = 10;
    width = 20;
    height = 20;
    curSubimageBuffer = 20;
    imageEdgeBuffer = 75;
elseif strcmp(Case, 'translation')
    load('translation_data.mat')
    candidateGridDimensionX = 10;
    candidateGridDimensionY = 10;
    gridType = 'rectangular';
    width = 20;
    height = 20;
    curSubimageBuffer = 20;
    imageEdgeBuffer = 75;
else
    error("Invalid case string")
end

%% Create unstructured grid
backgroundIntensityCutoff = .1;

if strcmp(gridType,'random') 
    % Generate random x and y coordinates on the image
    buffer = 41;
    candidateGridX = round((size(ref,2)-2*buffer)*rand(candidateGridSize,1)+buffer);
    candidateGridY = round((size(ref,1)-2*buffer)*rand(candidateGridSize,1)+buffer);
    % Discard points in the background
    grid = [];
    for i=1:candidateGridSize
        if ref(candidateGridY(i),candidateGridX(i))>backgroundIntensityCutoff
            grid = [grid; [candidateGridX(i) candidateGridY(i)]];
        end
    end
    % Discard duplicates
    grid = unique(grid,'rows');
end

if strcmp(gridType,'rectangular')
    % Create rectangular grid
    
    candidateGridX = round(linspace(imageEdgeBuffer, size(ref,2)-imageEdgeBuffer, candidateGridDimensionX));
    candidateGridY = round(linspace(imageEdgeBuffer, size(ref,1)-imageEdgeBuffer, candidateGridDimensionY));
    % Discard points in the background
    grid = [];
    for i=1:size(candidateGridY,2)
        for j=1:size(candidateGridX,2)
            if ref(candidateGridY(i),candidateGridX(j))>backgroundIntensityCutoff
                grid = [grid; [candidateGridX(j) candidateGridY(i)]];
            end
        end
    end
end

%% Get Delaunay triangulation of unstructured grid
gridX = grid(:,1);
gridY = grid(:,2);
DT = delaunay(grid(:,1),grid(:,2));

%Get rid of triangles whose centroids are background
trianglesToRemove = [];
for i = 1:size(DT,1)
    triangle = DT(i,:);
    centroid = round([(grid(triangle(1),1)+grid(triangle(2),1)+grid(triangle(3),1))/3,(grid(triangle(1),2)+grid(triangle(2),2)+grid(triangle(3),2))/3]);
    if ref(centroid(2),centroid(1))<backgroundIntensityCutoff
        trianglesToRemove = [trianglesToRemove; i];
    end
end
for i = 1:size(trianglesToRemove)
    DT(trianglesToRemove(end+1-i),:)=[];
end

%% Get displacements with normxcorr2
displacementsList = [];
for i=1:length(grid)
    subImageX = grid(i,1);
    subImageY = grid(i,2);

    refSubImageTopLeftY = round(subImageY-height/2);
    refSubImageTopLeftX = round(subImageX-width/2);
    refSubImageBottomRightY = round(subImageY+height/2);
    refSubImageBottomRightX = round(subImageX+width/2);
    refSubimage = ref(refSubImageTopLeftY:refSubImageBottomRightY, refSubImageTopLeftX:refSubImageBottomRightX,1);
    
    curSubImageTopLeftY = round(subImageY-height/2-curSubimageBuffer);
    curSubImageTopLeftX = round(subImageX-width/2-curSubimageBuffer);
    curSubImageBottomRightY = round(subImageY+height/2+curSubimageBuffer);
    curSubImageBottomRightX = round(subImageX+width/2+curSubimageBuffer);
    curSubimage = cur(curSubImageTopLeftY:curSubImageBottomRightY, curSubImageTopLeftX:curSubImageBottomRightX,1);
    
    c = normxcorr2(refSubimage, curSubimage);        
    [ypeak,xpeak] = find(c==max(c(:)));
    ytopleft = ypeak-size(refSubimage,1)+curSubImageTopLeftY;
    xtopleft = xpeak-size(refSubimage,2)+curSubImageTopLeftX;

    rectangleCur = [xtopleft,ytopleft,width,height];
    rectangleRef = [refSubImageTopLeftX,refSubImageTopLeftY,width,height];
    ref = insertShape(ref,'Rectangle',rectangleRef,'LineWidth',1,'Color',[1,0,0]);
    cur = insertShape(cur,'Rectangle',rectangleCur,'LineWidth',1,'Color',[1,0,0]);

    displacementsList = [displacementsList;[subImageX,subImageY,xtopleft+width/2-subImageX,ytopleft+height/2-subImageY]];
end

%% Tune displacements with cpcorr
movingPoints = [displacementsList(:,1)+displacementsList(:,3),displacementsList(:,2)+displacementsList(:,4)];
fixedPoints = [displacementsList(:,1),displacementsList(:,2)];
newPoints = cpcorr(movingPoints, fixedPoints, cur(:,:,1), ref(:,:,1));
displacementsList(:,3) = newPoints(:,1)-displacementsList(:,1);
displacementsList(:,4) = newPoints(:,2)-displacementsList(:,2);

%% Convert Displacements to useful units
lambda = 1.0;  % [length unit/pixel]
displacementsList = displacementsList*lambda;

%% Plot images with regions
figure();
tiledlayout(1,2);
nexttile
imshow(ref)
title("Subregions in reference image")
nexttile
imshow(cur)
title("Subregions in deformed image")

%% Scatterplot of displacements
% Multiply y by -1 since image starts in the upper left
figure();
scatter(displacementsList(:,1),-displacementsList(:,2)+size(ref,1), 40, 'red')
hold on;
scatter(displacementsList(:,1)+displacementsList(:,3), -displacementsList(:,2)-displacementsList(:,4)+size(ref,1), 40, 'green')
hold off;
xlim([0,size(ref,2)])
ylim([0,size(ref,1)])
title("Grid displacement")
xlabel("x (pixels)")
ylabel("y (pixels)")

%% Countour plots of displacements
figure();
colormap jet;
tiledlayout(1,2);
nexttile
sm = surfaceMesh([gridX, -gridY+size(ref,1), displacementsList(:,3)],DT);
smOut = smoothSurfaceMesh(sm,10,Method="Laplacian",ScaleFactor=.5);
[~,h] = tricontf(gridX, -gridY+size(ref,1), double(smOut.Faces), smOut.Vertices(:,3));
set(h,'edgecolor','none');
colorbar
xlim([0,size(ref,2)])
ylim([0,size(ref,1)])
title("x displacement")
xlabel("x (pixels)")
ylabel("y (pixels)")
nexttile
sm = surfaceMesh([gridX, -gridY+size(ref,1), displacementsList(:,4)],DT);
smOut = smoothSurfaceMesh(sm,10,Method="Laplacian",ScaleFactor=.5);
[~,h] = tricontf(gridX, -gridY+size(ref,1), double(smOut.Faces), smOut.Vertices(:,3));
set(h,'edgecolor','none');
colorbar
xlim([0,size(ref,2)])
ylim([0,size(ref,1)])
title("y displacement")
xlabel("x (pixels)")
ylabel("y (pixels)")

%% Calculate F and strain measures
[uxx,uxy] = trigradient(displacementsList(:,1),displacementsList(:,2),displacementsList(:,3),DT);
[uyx,uyy] = trigradient(displacementsList(:,1),displacementsList(:,2),displacementsList(:,4),DT);

I = [1 0;0,1];
for i=1:size(grid,1)
    gradientu(i,:,:) = [uxx(i),uxy(i);uyx(i),uyy(i)];
    F(i,:,:) = reshape(gradientu(i,:,:),[2,2])+I;
    F_local = reshape(F(i,:,:),[2,2]);
    [R_local,U_local,V_local] = poldecomp(F_local);
    R(i,:,:) = R_local;
    U(i,:,:) = U_local;
    V(i,:,:) = V_local;
    J(i) = det(F_local);
    E(i,:,:) = 1/2*(F_local'*F_local-I);
    C(i,:,:) = F_local'*F_local;
    B(i,:,:) = F_local*F_local';
    Estar(i,:,:) = 1/2*(I-inv(F_local)'*inv(F_local));
    epsilon(i,:,:) = 1/2*(reshape(gradientu(i,:,:),[2,2])+reshape(gradientu(i,:,:),[2,2])');
    omega(i,:,:) = 1/2*(reshape(gradientu(i,:,:),[2,2])'-reshape(gradientu(i,:,:),[2,2]));
end

%% Plot infinitesimal strain epsilon and infinitesimal rotation omega
figure()
colormap jet;
tiledlayout(2,2);
nexttile
sm = surfaceMesh([gridX, -gridY+size(ref,1), epsilon(:,1,1)],DT);
smOut = smoothSurfaceMesh(sm,10,Method="Laplacian",ScaleFactor=.5);
[~,h] = tricontf(gridX, -gridY+size(ref,1), double(smOut.Faces), smOut.Vertices(:,3));
set(h,'edgecolor','none');
colorbar
title("\epsilonxx")
xlim([0,size(ref,2)])
ylim([0,size(ref,1)])
nexttile
sm = surfaceMesh([gridX, -gridY+size(ref,1), epsilon(:,2,2)],DT);
smOut = smoothSurfaceMesh(sm,10,Method="Laplacian",ScaleFactor=.5);
[~,h] = tricontf(gridX, -gridY+size(ref,1), double(smOut.Faces), smOut.Vertices(:,3));
set(h,'edgecolor','none');
colorbar
title("\epsilonyy")
xlim([0,size(ref,2)])
ylim([0,size(ref,1)])
nexttile
sm = surfaceMesh([gridX, -gridY+size(ref,1), epsilon(:,1,2)],DT);
smOut = smoothSurfaceMesh(sm,10,Method="Laplacian",ScaleFactor=.5);
[~,h] = tricontf(gridX, -gridY+size(ref,1), double(smOut.Faces), smOut.Vertices(:,3));
set(h,'edgecolor','none');
colorbar
title("\epsilonxy")
xlim([0,size(ref,2)])
ylim([0,size(ref,1)])
nexttile
sm = surfaceMesh([gridX, -gridY+size(ref,1), omega(:,1,2)],DT);
smOut = smoothSurfaceMesh(sm,10,Method="Laplacian",ScaleFactor=.5);
[~,h] = tricontf(gridX, -gridY+size(ref,1), double(smOut.Faces), smOut.Vertices(:,3));
set(h,'edgecolor','none');
colorbar
title("\omegaxy")
xlim([0,size(ref,2)])
ylim([0,size(ref,1)])