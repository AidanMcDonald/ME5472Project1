clear;
close all;
clc;


%% Load input files
filedirectory = 'Hole_Plate_Images';
filename = 'ohtcfrp_';
ref = read(Tiff(strcat(filedirectory,'/',filename,'01','.tif')));
ref = double(ref)/256.;
cur = read(Tiff(strcat(filedirectory,'/',filename,'11','.tif')));
cur = double(cur)/256.;


%load('rotation_data.mat')

% cur = imrotate(cur, 3);
% cur = cur(round(size(cur)/2)-188:round(size(cur)/2)+187,round(size(cur)/2)-188:round(size(cur)/2)+187);
% 
% cur = imresize(cur,1.15);
% cur = cur(round(size(cur)/2)-188:round(size(cur)/2)+187,round(size(cur)/2)-188:round(size(cur)/2)+187);
% 
% Fxx = 1.01;
% Fxy = -0.05;
% Fyx = -.001;
% Fyy = .99;
% tform = affine2d([Fxx Fxy 0; Fyx Fyy 0; 0 0 1]);
% cur = imwarp(ref,tform);
% cur = cur(round(size(cur)/2)-188:round(size(cur)/2)+187,round(size(cur)/2)-188:round(size(cur)/2)+187);

%% Create unstructured grid
backgroundIntensityCutoff = .1;

% candidateGridSize = 1000;
% buffer = 30;
% 
% candidateGridX = round((size(ref,2)-2*buffer)*rand(candidateGridSize,1)+buffer);
% candidateGridY = round((size(ref,1)-2*buffer)*rand(candidateGridSize,1)+buffer);
% grid = [];
% for i=1:candidateGridSize
%     if ref(candidateGridY(i),candidateGridX(i))>backgroundIntensityCutoff
%         grid = [grid; [candidateGridX(i) candidateGridY(i)]];
%     end
% end


candidateGridDimensionX = 30;
candidateGridDimensionY = 50;
%candidateGridX = round(linspace(1, size(ref,2), candidateGridDimensionX));
%candidateGridY = round(linspace(1, size(ref,1), candidateGridDimensionY));
candidateGridX = round(linspace(75, size(ref,2)-75, candidateGridDimensionX));
candidateGridY = round(linspace(75, size(ref,1)-75, candidateGridDimensionY));
% candidateGridY = round(linspace(400, 640, candidateGridDimensionX));
% candidateGridX = round(linspace(80, 320, candidateGridDimensionY));

grid = [];
for i=1:size(candidateGridY,2)
    for j=1:size(candidateGridX,2)
        if ref(candidateGridY(i),candidateGridX(j))>backgroundIntensityCutoff
            grid = [grid; [candidateGridX(j) candidateGridY(i)]];
        end
    end
end


%% Get displacements with normxcorr2
displacementsList = [];
for i=1:length(grid)
    subImageX = grid(i,1);
    subImageY = grid(i,2);
    width = 20;
    height = 20;
    curSubimageBuffer = 20;

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

gridX = grid(:,1);
gridY = grid(:,2);

%DT = delaunayTriangulation(grid(:,1),grid(:,2));
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

% k = 1;
% for i=1:length(grid)
%    displacementsMatrix(j,i,:) = [displacementsList(k,3:4)];
%    k=k+1;
% end

%% Smooth displacement matrix
%displacementsMatrix(:,:,1) = smooth2a(displacementsMatrix(:,:,1), 2);
%displacementsMatrix(:,:,2) = smooth2a(displacementsMatrix(:,:,2), 2);

%% Convert Displacements to useful units
lambda = 1.0;  % [length unit/pixel]
% displacementsMatrix = displacementsMatrix*lambda;

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
tiledlayout(1,2);
nexttile
[~,h] = tricontf(gridX, -gridY+size(ref,1), DT, displacementsList(:,3));
set(h,'edgecolor','none');
colorbar
xlim([0,size(ref,2)])
ylim([0,size(ref,1)])
title("x displacement")
xlabel("x (pixels)")
ylabel("y (pixels)")
nexttile
[~,h] = tricontf(gridX, -gridY+size(ref,1), DT, -displacementsList(:,4));
set(h,'edgecolor','none');
colorbar
xlim([0,size(ref,2)])
ylim([0,size(ref,1)])
title("y displacement")
xlabel("x (pixels)")
ylabel("y (pixels)")

%% Calulate displacement gradient

% For each triangle, get the displacment gradient inside it by finite
% element interpolation

%% Calculate F and strain measures
[uxx,uxy] = trigradient(displacementsList(:,1),displacementsList(:,2),displacementsList(:,3),DT);
[uyx,uyy] = trigradient(displacementsList(:,1),displacementsList(:,2),displacementsList(:,4),DT);
% figure()
% quiver(gridX, -gridY+size(ref,2), Fx, Fy)

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

figure()
tiledlayout(2,2);
nexttile
[~,h] = tricontf(gridX, -gridY+size(ref,1), DT, epsilon(:,1,1));
set(h,'edgecolor','none');
colorbar
title("\epsilonxx")
xlim([0,size(ref,2)])
ylim([0,size(ref,1)])
nexttile
[~,h] = tricontf(gridX, -gridY+size(ref,1), DT, epsilon(:,2,2));
set(h,'edgecolor','none');
colorbar
title("\epsilonyy")
xlim([0,size(ref,2)])
ylim([0,size(ref,1)])
nexttile
[~,h] = tricontf(gridX, -gridY+size(ref,1), DT, epsilon(:,1,2));
set(h,'edgecolor','none');
colorbar
title("\epsilonxy")
xlim([0,size(ref,2)])
ylim([0,size(ref,1)])
nexttile
[~,h] = tricontf(gridX, -gridY+size(ref,1), DT, omega(:,1,2));
set(h,'edgecolor','none');
colorbar
title("\omegaxy")
xlim([0,size(ref,2)])
ylim([0,size(ref,1)])