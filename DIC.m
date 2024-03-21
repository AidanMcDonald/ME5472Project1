clear;
close all;
clc;

load('translation_data.mat')
%cur = ref;


%% Get displacements with normxcorr2
nGridPoints = 20;
gridX = linspace(50, 325, nGridPoints);
gridY = linspace(50, 325, nGridPoints);
displacementsList = [];
for i=1:length(gridX)
    for j=1:length(gridY)

        subImageX = gridX(i);
        subImageY = gridY(j);
        width = 10;
        height = 10;
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

        displacementsList = [displacementsList;[gridX(i),gridY(j),xtopleft+width/2-subImageX,ytopleft+height/2-subImageY]];
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
displacementsMatrix(:,:,1) = smooth2a(displacementsMatrix(:,:,1), 1);
displacementsMatrix(:,:,2) = smooth2a(displacementsMatrix(:,:,2), 1);

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
        F(i,j,:,:) = [uxx(i,j),uxy(i,j);uyx(i,j),uyy(i,j)]+I;
        J(i,j) = det(F(i,j));
        E(i,j,:,:) = 1/2*(F(i,j)'*F(i,j)-I);
        C(i,j,:,:) = F(i,j)'*F(i,j);
        B(i,j,:,:) = F(i,j)*F(i,j)';
        Estar(i,j,:,:) = 1/2*(I-inv(F(i,j))'*inv(F(i,j)));
    end
end