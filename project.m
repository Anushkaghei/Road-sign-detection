function project(~) % function doesn't take an input argument
clc %clear command window

% Open a dialog box for the user to select an image
[filename, folder] = uigetfile({'*.jpg'; '*.png'; '*.bmp'; '*.tif'}, 'Select an image file');

% Check if the user clicked 'Cancel' (empty file selection)
if isequal(filename, 0)
    disp('No image selected. Program terminated.');
    return; % Exit the script if no image is selected
end

% Load the selected image
I1 = imread(fullfile(folder, filename));

% Continue with the rest of the code as before
% ...

% Create the figure with the image name as the title
figure('Name', filename, 'NumberTitle', 'off');

% Continue with the rest of the code as before
% ...




%dimensions of image I1
xdim = size(I1, 2); %number of rows
ydim = size(I1, 1); %number of columns
%New 3 channel(rgb) image R1 of the same size as the original with all
%pixels set to 0.
R1 = zeros(ydim, xdim, 3, 'uint8'); %data type- unsigned 8 bit integers

% Image with only red channel 
T1 = zeros(ydim, xdim, 1);

% store infoabout overlap data during polygon matching
%2x4 matrix with all 0s
matchData = zeros(2, 4); 

% Store information on overlap - converts from rgb to hsv
hsvI = rgb2hsv(I1);

N_sides = 4; %sets Number of sides in polygon
t = (1 /(N_sides * 2):1 / N_sides:1)' * 2 * pi; % t vector for polygon that 
% represents angles for polygon vertices. equal spaced angles 0 to 2pi
% radians. corresponding to number of sides in a polygon.

iteration = 1;% number of times images has been checked for polygon fit

%HSV criteria to filter images. pixels here will be included in R1 image
HueLower = 20 / 255;
HueUpper = 150 / 255;
SaturationMin = 50 / 255;
ValueMin =  50 / 255;

%nested loop to check if the pixel's hue is outside the specified range
%and if its saturation and value are greater than the specified thresholds 
for x = 1:1:xdim
    for y = 1:1:ydim
        if (hsvI(y, x, 1) > HueUpper || hsvI(y, x, 1) < HueLower) && (hsvI(y, x, 2) > SaturationMin) && (hsvI(y, x, 3) > ValueMin)
            %copies the RGB pixel values from the original to new image
            % preserving only the pixels that satisfy the HSV criteria. 
            % The rest of the pixels in R1 remain zero, as initialized.
            R1(y, x, 1) = I1(y, x, 1);
            R1(y, x, 2) = I1(y, x, 2);
            R1(y, x, 3) = I1(y, x, 3);
        end
    end
end
%binary edge map-edge detection on R1, by converting to grayscale then applying sobel
E1 = edge(rgb2gray(R1), 'sobel', 0.1);

%nested loop to iterate through binary edge map
% For each (x, y) it expands by setting the neighboring pixels to 1. 
% forms a filled polygon shape around the detected edges.
for x = 2:1:xdim
    for y = 2:1:ydim
        if E1(y, x) == 1
           T1(y + 1, x - 1) = 1;
           T1(y + 0, x - 1) = 1;
           T1(y - 1, x - 1) = 1;
           T1(y + 1, x + 0) = 1;
           T1(y - 0, x + 0) = 1;
           T1(y + 1, x + 1) = 1;
           T1(y + 0, x + 1) = 1;
           T1(y - 1, x + 1) = 1;
        end
    end
end

s = strel('square', 5);% structure element of 5x5 used in morphological operations
% helps fill gaps and smoothen the binary mask to create the largest
% connected object
C1 = imclose(T1, s); 
BW2 = imfill(C1,'holes');

% filter binary mask Get 3 largest objects from BW2 and T1
F1 = bwareafilt(logical(BW2), 3); 
N1 = bwareafilt(logical(T1), 3); 
% Image with blank pixels removed

% Get boundaries of edges, ~ used to discard unnecessary outputs
[~, ~, ~, A] = bwboundaries(E1, 'noholes'); 
% B = Boundaries, N = number of objects, 4 = connectivity

% Only grab slices of images that contain pixels to improve efficiency
columns = any(N1);%check if it contains atleast 1 non zero pixel(logical true)
Loc = 1;%current position while extracting columns from image

%extract columns with only non zero pixels
tic
for x = 1:1:xdim
        if columns(x) == 1
            Slice1(:, Loc) = N1(:, x);
            Slice1Image(:, Loc, 1) = I1(:, x, 1);
            Slice1Image(:, Loc, 2) = I1(:, x, 2);
            Slice1Image(:, Loc, 3) = I1(:, x, 3);
            Loc = Loc + 1;
        end
end
toc

rows = any(Slice1, 2);
Loc = 1;

tic
for y = 1:1:ydim
        if rows(y) == 1  
            Slice2(Loc, :) = Slice1(y, :);
            Slice2Image(Loc, :, 1) = Slice1Image(y, :, 1);
            Slice2Image(Loc, :, 2) = Slice1Image(y, :, 2);
            Slice2Image(Loc, :, 3) = Slice1Image(y, :, 3);
            Loc = Loc + 1;
        end
end
toc

% Reduce resolution of slice 2 to further improve efficiency
multiplier = 128 / size(Slice2, 1);
if multiplier < 1 %if size is less than 1 resolution of slice 2 needs to be reduced
    Slice2 = imresize(Slice2, multiplier);
end %done to avoid overprocessing the large images and speed up computation

xdim = size(Slice2, 2);
ydim = size(Slice2, 1);

% check different polygon positions and scales for the best-fit shape matching.
for scale = (xdim / 5):3:(xdim)
    for xOffset = abs(min(scale*sin(t))):3:xdim - abs(min(scale*sin(t)))
        for yOffset = abs(min(scale*cos(t)))+ 1:3:ydim - max(scale*cos(t))         
            % Get vertices
            polyCoord = [round(scale*sin(t) + xOffset) round(scale*cos(t) + yOffset)];
            polyCoord(N_sides + 1, 1) = polyCoord(1, 1);
            polyCoord(N_sides + 1, 2) = polyCoord(1, 2);
            
            % Create polygon using vertices
            j = 1;
            k = 1;
            
            for i = 1:1:(N_sides * 2) 
                if mod(i, 2) == 1 
                    poly(i) = polyCoord(j, 1);
                    j = j + 1;
                else
                    poly(i) = polyCoord(k, 2);
                    k = k + 1;
                end
            end
            
            % Draw Polygon on fresh image 
            S1 = zeros(ydim, xdim, 3, 'uint8');
            S1 = insertShape(S1, 'Polygon', poly, 'Color', 'white');

            % Convert image to logical type
            S2 =  imbinarize(rgb2gray(S1));

            % Get overlap
            Slice2 = imcrop(Slice2,[0 0 xdim ydim]);
            M1 = Slice2 .* S2;

            % Get number of overlapped pixels in images
            commonPixels = sum(M1(:) == 1);
            
            % Save data
            matchData(iteration, :) = [commonPixels scale xOffset yOffset];
            iteration = iteration + 1;     
        end
    end
end

% Calculate best match
[~, I] = max(matchData); % ~ discards the unnecessary first output
scale = matchData(I(1), 2);
xOffset = matchData(I(1), 3);
yOffset = matchData(I(1), 4);

% Create polygon for Best Fit Shape overlapped on Sliced/Processed Image
% Get vertices
polyCoord = [round(scale*sin(t) + xOffset) round(scale*cos(t) + yOffset)];
polyCoord(N_sides + 1, 1) = polyCoord(1, 1);
polyCoord(N_sides + 1, 2) = polyCoord(1, 2);

% Create polygon using vertices
j = 1;
k = 1;

for i = 1:1:(N_sides * 2) 
    if mod(i, 2) == 1 % odd
        poly(i) = polyCoord(j, 1);
        j = j + 1;
    else
        poly(i) = polyCoord(k, 2);
        k = k + 1;
    end
end

% Create another polygon for Best Fit Shape overlapped on Sliced Raw Image
% Only needed if multiplier < 1
if multiplier < 1
    polyCoord2 = [round((scale*sin(t) + xOffset)/multiplier) round((scale*cos(t) + yOffset)/multiplier)];
    polyCoord2(N_sides + 1, 1) = polyCoord2(1, 1);
    polyCoord2(N_sides + 1, 2) = polyCoord2(1, 2);
else
    polyCoord2 = polyCoord;
end

% Create polygon using vertices
j = 1;
k = 1;

for i = 1:1:(N_sides * 2) 
    if mod(i, 2) == 1 % odd
        poly(i) = polyCoord2(j, 1);
        j = j + 1;
    else
        poly(i) = polyCoord2(k, 2);
        k = k + 1;
    end
end

sum(A(:) == 1);

%figure('Name',sprintf('BIM472_Image%02d.jpg', 5),'NumberTitle','off');

subplot(2, 3, 1), imshow(I1);
title('Original');

subplot(2, 3, 2), imshow(R1);
title('Extract Red Pixels');


subplot(2, 3, 3), imshow(T1);
title('Fill in Edges');

subplot(2, 3, 4), imshow(F1);
title('Largest Object');

subplot(2, 3, 5)
hold on;
imshow(Slice2Image);
plot(polyCoord2(:,1), polyCoord2(:,2), 'g')
hold off;
title('Best Fit Shape');
