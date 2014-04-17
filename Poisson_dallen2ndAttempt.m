function [error] = Poisson_dallen2ndAttempt(tFileName, sFileName, mFileName, iFileName, nonC, iterations, dimReturnsThreshold, verbose)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Start timer
tic
fprintf('====\n');
fprintf('Poisson_dallen2ndAttempt :: Executing\n');
fprintf('====\n');

%% Handle default arguments
for i = 1 % For loop is for code collapsing only (so I don't have to look at these)
    if( nargin < 1 )
        tFileName = 'test_T.jpg';
    end
    
    if( nargin < 2 )
        sFileName = 'test_S.jpg';
    end
    
    if( nargin < 3 )
        mFileName = 'test_M.jpg';
    end
    
    if( nargin < 4 )
        iFileName = 'Poisson_dallen2ndAttempt_test_I'; % set for unique final composit output names
    end
    
    if( nargin < 5 )
        nonC = false; % set to true for non-conservative field compositing
    end
    
    if( nargin < 6 )
        iterations = 0.5E4;
    end
    
    if( nargin < 7 )
        dimReturnsThreshold = 1E-4;
    end
    
    if( nargin < 8 )
        %         verbose = true;
        verbose = false;
    end
end

%% Load images
% TODO: Generalize this into a GUI that allows you to select images
T0 = double(imread(tFileName)); % Target image that makes up largest area of composited image
S0 = double(imread(sFileName)); % Source image for desired elements to be composited into target image
m = double(imread(mFileName)); % Mask image defining area ? where composited image should "look like" S0 (source image)

% Save dimensions for quick reference
w = size(T0, 1);
h = size(T0, 2);
chan = size(T0, 3);

% Check dimensions of all images to make sure they match, as well as bit
% depth of S0 & T0 to make sure they match
if size(S0,1) ~= w | size(S0,2) ~= h | size(S0,3) ~= chan | size(m,1) ~= w | size(m,2) ~= h 
    % Throw error and return from function
    fprintf('====\n'); 
    fprintf('Error - Poisson_dallen2ndAttempt:\n'); 
    fprintf('--\n'); 
    fprintf('source, target, and mask image dimensions must match,\n'); 
    fprintf('and source, target must have same number of channels\n'); 
    fprintf(' - Target is %d x %d x %d\n', size(T0,1), size(T0,2), size(T0,3)); 
    fprintf(' - Source is %d x %d x %d\n', size(S0,1), size(S0,2), size(S0,3)); 
    fprintf(' - Mask is %d x %d x %d\n', size(m,1), size(m,2), size(m,3));
    fprintf('====\n'); 
    return;
end

% Convert Mask to 1-bit (based on first channel, if it has more than one)
m1bit = (m(:,:,1) > 127); 
% Expand out to match channels in T for matrix math
if chan > 1 
    M0 = zeros(w, h, chan); 
    for i = 1:(chan)
        M0(:,:,i) = m1bit(:,:);
    end
else
    M0 = m1bit;
end

% Clean up
clear tFileName sFileName mFileName m m1bit


%% Create Laplacian0
% Calculate Laplacians for T & S images
% K = [0,1,0; 1,-4,1; 0,1,0]; % Lapacian filter Kernel
K = -[0,1,0; 1,-4,1; 0,1,0]; % Lapacian filter Kernel (inverting it fixed channels that were inverted)
Lt = imfilter(T0,K,'replicate'); % Target Laplacian
Ls = imfilter(S0,K,'replicate'); % Source Laplacian

% Create a guidance vector field to keep gradients that are stronger in
% target image (T0) when nonC is true
LtAmp = (Lt(:,:,1) .* Lt(:,:,2) .* Lt(:,:,3)) .^(1/3); % Greyscale amplitudes of T gradient 
LsAmp = (Ls(:,:,1) .* Ls(:,:,2) .* Ls(:,:,3)) .^(1/3); % Greyscale amplitudes of S gradient
LgVF = (LtAmp < LsAmp) .* M0(:,:,1);
Lgm = zeros(w, h, chan);
for i = 1:chan
    Lgm(:,:,i) = LgVF(:,:);
end

% Composit Laplacian using compositing equation
Lm = nonC * Lgm + (1-nonC) * M0; % Calculates laplacian mask with influence of Guidance Vector Field Mask if nonC is true
L0 = Lm .* Ls + (1-Lm) .* Lt;

% Clean up
clear Lm

%% Build sparse matrix
% Create holder for unique indexed pixels within masked region
Iindex = zeros(w,h);
Icount = 0;
for y = 1:h
    for x = 1:w
        if M0(x,y,1) == 1
            % Increase index count, and store it at that pixel
            Icount = Icount + 1;
            Iindex(x,y) = Icount;
        end
    end
end


% Check length of matrix
MaskPixelCount = size(find(M0(:, :, 1)), 1);
if Icount ~= MaskPixelCount
    % Something went wrong
    fprintf('====\n'); 
    fprintf('Error - Poisson_dallen2ndAttempt:\n'); 
    fprintf('--\n'); 
    fprintf('index count does not match number of 1''s in mask\n'); 
    fprintf('====\n'); 
end

% Create large, sparse, linear array from indexes
A = delsq(padarray(Iindex, [1,1]));


%% Create boundary term
% Create container for new image
I0 = T0;
% Populate it with appropriate values from T0 and L0
for i = 1:chan
    b = zeros(1, MaskPixelCount); % Holder for boundary term 
    Icount = 0;
    for y = 1:h
        for x = 1:w
            
            % If the mask is on for this pixel
            if M0(x, y, 1) == 1
                
                % Increase index count
                Icount = Icount + 1;
                
                %TODO:
                %  1) Create special coefficients for situations where mask
                % goes to boundary of the image (current fix has too much
                % color bleeding into S from T around edges of image)
                %  2) nonC version of this is still not working properly.
                % Seems to be coming out negative. I might have to derive
                % the combined edges some other way
                
                % Adjust b for boundary pixels
                if y > 1
                    if M0(x, y-1, i) == 0 % Top
                        b(Icount) = b(Icount) + T0(x, y-1, i);
                    end
                else
                    b(Icount) = b(Icount) + S0(x, y, i);
                end
                
                if x > 1
                    if M0(x-1, y, i) == 0 % Left
                        b(Icount) = b(Icount) + T0(x-1, y, i);
                    end
                else
                    b(Icount) = b(Icount) + S0(x, y, i);
                end
                
                if y < h
                    if M0(x, y+1, i) == 0 % Bottom
                        b(Icount) = b(Icount) + T0(x, y+1, i);
                    end
                else
                    b(Icount) = b(Icount) + S0(x, y, i);
                end
                
                if x < w
                    if M0(x+1, y, i) == 0 % Right
                        b(Icount) = b(Icount) + T0(x+1, y, i);
                    end
                else
                    b(Icount) = b(Icount) + S0(x, y, i);
                end
                % Construct the guidance field
                IgV = L0(x, y, i);
                b(Icount) = b(Icount) + IgV;
            end
        end
    end
    
    % Calculate intensity of unknown image pixels in this channel
    Ix = A\b'; 
    % Alternate method: using biconjugate gradients
    %     Ix = bicg(A, b', [], 400); % Is it faster in any way?
    
    % Store new source values into composite image
    for y = 1:h
        for x = 1:w
            if M0(x,y,i) == 1 % Copy new intensities for pixels within the mask
                IxIndex = Iindex(x, y); % Find unique index in unique index location array
                I0(x, y, i) = Ix(IxIndex); % Assign value at that index to composite pixel
            end
        end
    end
    
    % Debugging
    %     imwrite(uint8( I0(:,:,i) ),[ iFileName '-chan' num2str(i) '.png'   ]);
end

% Clean up
clear Iindex Icount MaskPixelCount
clear A b x 

%% Save
% Debugging images
if false % These never execute
    imwrite(uint8( M0*255   ),[ iFileName '-M0.png'         ]);
    imwrite(uint8( Lt       ),[ iFileName '-Lt.png'         ]);
    imwrite(uint8( Ls       ),[ iFileName '-Ls.png'         ]);
    imwrite(uint8( LtAmp    ),[ iFileName '-LtAmp.png'      ]);
    imwrite(uint8( LsAmp    ),[ iFileName '-LsAmp.png'      ]);
    imwrite(uint8( LgVF*255 ),[ iFileName '-LgVF.png'       ]);
    imwrite(uint8( Lgm*255  ),[ iFileName '-Lgm.png'        ]);
end

% Move into debugging if statement as they begin working as expected
imwrite(uint8( L0       ),[ iFileName '-L0.png'         ]);

% Save final image
imwrite(uint8( I0        ),[ iFileName '.png'        ]);

% Report timestamp
fprintf('====\n');
fprintf('Poisson_dallen2ndAttempt :: End\n');
toc
fprintf('====\n');

%% Cleanup Variables (some of these could be cleared earlier to save memory)
clear iFileName nonC iterations dimReturnsThreshold verbose
clear w h chan 
clear S0 M0 T0 K Lt Ls LtAmp LsAmp LgVF Lgm Lm L0 I0

end

