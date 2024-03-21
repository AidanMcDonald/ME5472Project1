clear;
close all;
clc;

load('rotation_data.mat')
%cur = ref;

gridX = linspace(50, 325, 10);
gridY = linspace(50, 325, 10);
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
        %I want the x coordinates to be along the horizontal dimension i.e.
        % columns
        displacementsMatrix(j,i,:) = [xtopleft+width/2-subImageX,ytopleft+height/2-subImageY];

    end
end

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
contourf(gridX, -gridY+size(ref,2), displacementsMatrix(:,:,2), "ShowText",true)
xlim([0,size(ref,1)])
ylim([0,size(ref,2)])
title("y displacement")
xlabel("x (pixels)")
ylabel("y (pixels)")