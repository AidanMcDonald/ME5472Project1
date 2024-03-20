clear;
close all;
clc;

load('rotation_data.mat')
%cur = ref;

gridX = linspace(50, 325, 10);
gridY = linspace(50, 325, 10);
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

        displacements(i,j,:) = [round(xtopleft+width/2)-round(subImageX),round(ytopleft+height/2)-round(subImageY)];

    end
end

montage({ref, cur})
