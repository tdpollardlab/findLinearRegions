function [intervalStarts,slopes,yInts] = findLinearRegions(data, intervalSize,showLineSegments,dataName)

%Parameters: intervalSize must be an odd number that is smaller than data.
timeStep = 0.1;
%data = concSmoothed;
%intervalSize = 35; %odd number

if nargin < 4
    dataName = 'Sample';
    if nargin == 2
        showLineSegments = false;
    end
end

%halfInterval = floor(intervalSize/2);
dataSize = length(data);
rSquaredArray = zeros(dataSize-intervalSize+1,1);
coeffsArray = zeros(dataSize-intervalSize+1,2);
xAxisSegment = flip(rot90(1:intervalSize));

for i = 1:dataSize-intervalSize+1
   [line,goodnessOfFit] = fit(xAxisSegment,data(i:i+intervalSize-1),'poly1');
   coeffsArray(i,:) = coeffvalues(line);
   rSquaredArray(i) = goodnessOfFit.rsquare;
end
%Plotting data, r-squared values
if showLineSegments
    % figure
    xAxis = 0:timeStep:(length(data)-1)*timeStep;
    xlim([0 xAxis(end)]);
    ylim([0 1.1*max(data)]);
    hold on
end

%%plotting r-squared values such that they align with center point of
%%interval
%rSquaredArrayShifted = NaN(dataSize,1);
%rSquaredArrayShifted(halfInterval+1:dataSize-halfInterval) = rSquaredArray;
%plot(rSquaredArrayShifted*max(data))

%finding maximal r-squared values and eliminating peaks that are too close
%to each other
[peaks,locs] = findpeaks(rSquaredArray);
%finalPeaks = peaks;
intervalStarts = locs;
for i = 1:length(locs) %eliminating linearity peaks where there's a higher point within intervalSize
    %check points to the left
    for j = 1:locs(i)-1
        if j > intervalSize
            break
        elseif rSquaredArray(locs(i)-j) >= peaks(i)
            %finalPeaks(i) = 0; %peaks that are too close are flagged with 0
            intervalStarts(i) = 0;
        end
    end
    %check points to the right
    for j = 1:dataSize-intervalSize+1-locs(i)
        if j > intervalSize
            break
        elseif rSquaredArray(locs(i)+j) > peaks(i)
            %finalPeaks(i) = 0;
            intervalStarts(i) = 0;
        end
    end
end
%finalPeaks(finalPeaks==0) = []; %removing flagged peaks
intervalStarts(intervalStarts==0) = [];
%%plotting interval starts so they align with center points of interval
%plot(intervalStarts+halfInterval+1,finalPeaks*max(data),'o')

%Plotting the best fit line segments at these points
numLocs = length(intervalStarts);
slopes = zeros(numLocs,1);
yInts = zeros(numLocs,1);
for i = 1:numLocs
    slopes(i) = coeffsArray(intervalStarts(i),1);
    yInts(i) = coeffsArray(intervalStarts(i),2);
end
if showLineSegments %Plot lines at these points
    predictedY = NaN(dataSize,numLocs);
    for i = 1:numLocs
        slope = slopes(i);
        yInt = yInts(i);
        xCurrent = intervalStarts(i);
        for x = xCurrent:xCurrent+intervalSize-1
            predictedY(x,i) = slope*(x-xCurrent+1)+yInt;
        end
        plot(xAxis,predictedY(:,i),'LineWidth',4); 
    end
    plot(xAxis,data,'color','black');
    % hold off
    title(['Linear regions of ' dataName])
    set(gca,'FontSize',20);
    xlabel('Time (s)')
    ylabel('Actin patch brightness (AU)')
    legend('Assembly phase','Disassembly phase')
end



%NEW STRATEGY: Find region of maximal linearity using a very large interval
%(i.e. maybe 2/3 the length of the assembly phase). Find the slope
%of the best fit line in that region. Then, find absolute sum of residuals
%between 

%Finding points of maximal linearity. For patch finding, there should only
%be two.
 
% %Test: first point
% fitArray = zeros(dataSize,1);
% i = 1;
% slope = finalCoeffs(i,1);
% yInt = finalCoeffs(i,2);
% fractions = zeros(dataSize,1);
% residuals = zeros(dataSize,1);
% predictedY = zeros(dataSize,1);
% %Finding residuals
% for j = 1:dataSize
%     Y = data(j);
%     deltaX = j-intervalStarts(i)+1;
%     predictedY(j) = slope*deltaX+yInt;
%     %difference(j) = 
%     fractions(j) = Y/predictedY(j);
%     residuals(j) = Y-predictedY(j);
% end
% %summing residuals/multiplying fractions
% sumResiduals = zeros(dataSize,1);
% productFractions = zeros(dataSize,1);
% for j = 1:dataSize
%     if j < intervalStarts(i) 
%         sumResiduals(j) = sum(residuals(j:intervalStarts(i)));
% %         productFractions(j) = prod(fractions(j:intervalStarts(i)));
% %     elseif j > intervalStarts(i)+intervalSize-1;
% %         sumResiduals(j) = sum(residuals(intervalStarts(i)+intervalSize-1:j));
% %         productFractions(j) = prod(fractions(intervalStarts(i)+intervalSize-1:j));
% %     else
% %         sumResiduals(j) = 0;
% %         productFractions(j) = 1;
%     else
%         sumResiduals(j) = sum(residuals(intervalStarts(i):j));
%         productFractions(j) = prod(fractions(intervalStarts(i):j));       
%     end
% end
% plot(predictedY);
% plot(varSquared*max(data))
% plot(sumResiduals+150);
% plot(productFractions*100+50);
% hold off;
% %ylim([-50 200])
% % %figure(3);
% % %plot(residua);
% % 
% % 
% % 
% % %varSquaredBinary = false(dataSize,1);
% % %varSquaredBinary(varSquared>threshold) = true;
% % %plot(varSquaredBinary*max(data));
% % %hold off
% % 
