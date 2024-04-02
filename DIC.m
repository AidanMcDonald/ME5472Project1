clear;
close all;
clc;

load('rotation_data.mat')

% cur = imrotate(cur, 3);
% cur = cur(round(size(cur)/2)-188:round(size(cur)/2)+187,round(size(cur)/2)-188:round(size(cur)/2)+187);
% 
% cur = imresize(cur,1.15);
% cur = cur(round(size(cur)/2)-188:round(size(cur)/2)+187,round(size(cur)/2)-188:round(size(cur)/2)+187);

Fxx = 1.04;
Fxy = -0.05;
Fyx = -.01;
Fyy = .96;
tform = affine2d([Fxx Fxy 0; Fyx Fyy 0; 0 0 1]);
cur = imwarp(ref,tform);
% cur = cur(round(size(cur)/2)-188:round(size(cur)/2)+187,round(size(cur)/2)-188:round(size(cur)/2)+187);

%% Get displacements with normxcorr2
nGridPoints = 15;
gridX = linspace(50, 325, nGridPoints);
gridY = linspace(50, 325, nGridPoints);
displacementsList = [];
for i=1:length(gridX)
    for j=1:length(gridY)
        gridCoords(i,j,:) = [gridX(i), gridY(j)];
        subImageX = gridX(i);
        subImageY = gridY(j);
        width = 10;
        height = 10;
        curSubimageBuffer = 39;

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
end

%% Tune displacements with cpcorr
movingPoints = [displacementsList(:,1)+displacementsList(:,3),displacementsList(:,2)+displacementsList(:,4)];
fixedPoints = [displacementsList(:,1),displacementsList(:,2)];
newPoints = cpcorr(movingPoints, fixedPoints, cur(:,:,1), ref(:,:,1));
displacementsList(:,3) = newPoints(:,1)-displacementsList(:,1);
displacementsList(:,4) = newPoints(:,2)-displacementsList(:,2);

k = 1;
for i=1:length(gridX)
    for j=1:length(gridY)
       displacementsMatrix(j,i,:) = [displacementsList(k,3:4)];
       k=k+1;
    end
end

%% Smooth displacement matrix
%displacementsMatrix(:,:,1) = smooth2a(displacementsMatrix(:,:,1), 2);
%displacementsMatrix(:,:,2) = smooth2a(displacementsMatrix(:,:,2), 2);

%% Convert Displacements to useful units
lambda = 1.0;  % [length unit/pixel]
displacementsMatrix = displacementsMatrix*lambda;

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
scatter(displacementsList(:,1),-displacementsList(:,2)+size(ref,2), 40, 'red')
hold on;
scatter(displacementsList(:,1)+displacementsList(:,3), -displacementsList(:,2)-displacementsList(:,4)+size(ref,2), 40, 'green')
hold off;
xlim([0,size(ref,1)])
ylim([0,size(ref,2)])
title("Grid displacement")
xlabel("x (pixels)")
ylabel("y (pixels)")

%% Countour plots of displacements
figure();
tiledlayout(1,2);
nexttile
contourf(gridX, -gridY+size(ref,2), displacementsMatrix(:,:,1), "ShowText",true)
xlim([0,size(ref,1)])
ylim([0,size(ref,2)])
title("x displacement")
xlabel("x (pixels)")
ylabel("y (pixels)")
nexttile
contourf(gridX, -gridY+size(ref,2), -displacementsMatrix(:,:,2), "ShowText",true)
xlim([0,size(ref,1)])
ylim([0,size(ref,2)])
title("y displacement")
xlabel("x (pixels)")
ylabel("y (pixels)")

%% Calculate F and strain measures
[uxx,uxy] = gradient(displacementsMatrix(:,:,1),gridX(2)-gridX(1),gridY(2)-gridY(1));
[uyx,uyy] = gradient(displacementsMatrix(:,:,2),gridX(2)-gridX(1),gridY(2)-gridY(1));
% figure()
% quiver(gridX, -gridY+size(ref,2), Fx, Fy)

I = [1 0;0,1];
for i=1:nGridPoints
    for j=1:nGridPoints
        gradientu(i,j,:,:) = [uxx(i,j),uxy(i,j);uyx(i,j),uyy(i,j)];
        F(i,j,:,:) = reshape(gradientu(i,j,:,:),[2,2])+I;
        F_local = reshape(F(i,j,:,:),[2,2]);
        [R_local,U_local,V_local] = poldecomp(F_local);
        R(i,j,:,:) = R_local;
        U(i,j,:,:) = U_local;
        V(i,j,:,:) = V_local;
        J(i,j) = det(F_local);
        E(i,j,:,:) = 1/2*(F_local'*F_local-I);
        C(i,j,:,:) = F_local'*F_local;
        B(i,j,:,:) = F_local*F_local';
        Estar(i,j,:,:) = 1/2*(I-inv(F_local)'*inv(F_local));
        epsilon(i,j,:,:) = 1/2*(reshape(gradientu(i,j,:,:),[2,2])+reshape(gradientu(i,j,:,:),[2,2])');
        omega(i,j,:,:) = 1/2*(reshape(gradientu(i,j,:,:),[2,2])'-reshape(gradientu(i,j,:,:),[2,2]));
    end
end

figure()
tiledlayout(2,2);
nexttile
contourf(gridX, -gridY+size(ref,2), smooth2a(epsilon(:,:,1,1),4), "ShowText",true)
title("\epsilonxx")
nexttile
contourf(gridX, -gridY+size(ref,2), smooth2a(epsilon(:,:,2,2),4), "ShowText",true)
title("\epsilonyy")
nexttile
contourf(gridX, -gridY+size(ref,2), smooth2a(epsilon(:,:,1,2),4), "ShowText",true)
title("\epsilonxy")
nexttile
contourf(gridX, -gridY+size(ref,2), smooth2a(omega(:,:,1,2),4), "ShowText",true)
title("\omegaxy")

%plotPrincipalComponents(epsilon, gridCoords)