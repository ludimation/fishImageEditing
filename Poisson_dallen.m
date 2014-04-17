function [error] = Poisson_dallen(tFileName, sFileName, mFileName, iFileName, nonC, iterations, dimReturnsThreshold, verbose)
%UNTITLED Summary of this function goes here
%   This algorithm is similar to a Jacobi method I found online
%   TODO: put full credits here


%% Handle default arguments
for i = 1 % for loop is for code collapsing only (so I don't have to look at these)
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
        iFileName = 'Poisson_dallen_test_I'; % set for unique final composit output names
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
% TODO: generalize this into a GUI that allows you to select images
t = double(imread(tFileName)); % Target image that makes up largest area of composited image
s = double(imread(sFileName)); % Source image for desired elements to be composited into target image
m = double(imread(mFileName)); % Mask image defining area ? where composited image should "look like" S0 (source image)


%% Create Image Arrays
% TODO: 
% check dimensions of all images to make sure they match
% check bit depth of S0 & T0 to make sure they match
% - throw an error if they don't

w = size(t, 1);
h = size(t, 2);
chan = size(t, 3);

T0 = t;
S0 = s;

m = (m(:,:,1) > 127); % Convert Mask to 1-bit based on first channel to avoid error
if chan > 1 % Expand back out to match channels in T for matrix math
    M0 = zeros(w, h, chan); 
    for i = 1:(chan)
        M0(:,:,i) = m(:,:);
    end
else
    M0 = m;
end


%% Create Laplacian0 from T0 & S0 Lapacians
% Calculate Laplacians for T & S images
K = [0,1,0; 1,-4,1; 0,1,0]; %Lapacian filter Kernel
Lt = imfilter(T0,K,'replicate'); %Target Laplacian
Ls = imfilter(S0,K,'replicate'); %Source Laplacian

% Create a guidance vector field 
% to keep gradients that are 
% stronger in target image
LtAmp = (Lt(:,:,1) .* Lt(:,:,2) .* Lt(:,:,3)) .^(1/3); % Greyscale amplitudes of T gradient 
LsAmp = (Ls(:,:,1) .* Ls(:,:,2) .* Ls(:,:,3)) .^(1/3); % Greyscale amplitudes of S gradient
LgVF = (LtAmp < LsAmp) .* M0(:,:,1);
Lgm = zeros(w,h,chan);
for i = 1:chan
    Lgm(:,:,i) = LgVF(:,:);
end

% Composit Laplacian using compositing equation
Lm = nonC * Lgm + (1-nonC) * M0; % Calculates laplacian mask with influence of Guidance Vector Field Mask if nonC is true
L0 = Lm .* Ls + (1-Lm) .* Lt;


%% Create initial composit image and its Laplacian

% Create I0 with pixels from S0 & T0 with alpha defined by M0
I0 = M0 .* S0 + (1-M0) .* T0;
% Calculate I0's laplacian
Li0 = imfilter(I0,K,'replicate'); %Composite Laplacian


%% Minimize the Error between the Laplacians (Li0 and L0)

Ifinal = I0;
In = Ifinal;
LiN = Li0;
Progress0 = 1E64;
for i = 1:iterations
    % TODO:
    % 1) Guidance vector field version doesn't seem to be working right
    % - image within M goes towards black rather quickly for some reason
    % - progress calculations go crazy
    % - math seems off somewhere
    % 2) This seems to be executing slowly
    % - is there a way to only execute all these operations on pixels where
    % M = 1?

    % Calculate and apply error between Laplacian filtered images
    Lerror = (L0 - LiN) / 4;
    Ifinal = M0 .* (Ifinal - Lerror) + (1-M0) .* Ifinal;
    
    % Calculate diminishing return
    Progress = abs(Ifinal - In);
    ProgressMax = max(max(max(Progress))); % might want to calculate this as an amplitude of each RGB vector instead
    % would this work? ProgressMax = max(max((Progress(:,:,1).*Progress(:,:,2).*Progress(:,:,3)).^(1/3))); % might want to calculate this as an amplitude of each RGB vector instead
    ProgressD = (Progress0 - ProgressMax)/Progress0;

    if( verbose )
        fprintf('i = %d; Progress0 = %g; ProgressMax = %g; ProgressD = %g < %g\n',i, Progress0, ProgressMax, ProgressD, dimReturnsThreshold);
    end
    
    if (ProgressD < dimReturnsThreshold )
        break;
    end
    
    % Update variables for next iteration
    Progress0 = ProgressMax;
    In = Ifinal;
    LiN = imfilter(In,K,'replicate'); %Composite Laplacian
end


%% Save
% Debugging images
if false % these never execute
    imwrite(uint8( M0*255   ),[ iFileName '-M0.png'         ]);
    imwrite(uint8( Lt       ),[ iFileName '-Lt.png'         ]);
    imwrite(uint8( Ls       ),[ iFileName '-Ls.png'         ]);
    imwrite(uint8( LtAmp    ),[ iFileName '-LtAmp.png'      ]);
    imwrite(uint8( LsAmp    ),[ iFileName '-LsAmp.png'      ]);
    imwrite(uint8( LgVF*255 ),[ iFileName '-LgVF.png'       ]);
    imwrite(uint8( Lgm*255  ),[ iFileName '-Lgm.png'        ]);
    imwrite(uint8( Progress ),[ iFileName '-Progress.png'   ]);
    imwrite(uint8( Lerror   ),[ iFileName '-Lerror.png'     ]);
end

% Move into debugging if statement as they begin working as expected
imwrite(uint8( Li0      ),[ iFileName '-Li0.png'        ]);
imwrite(uint8( L0       ),[ iFileName '-L0.png'         ]);
imwrite(uint8( I0       ),[ iFileName '-I0.png'         ]);

% Save final image
imwrite(uint8( Ifinal   ),[ iFileName '.png'        ]);

%% Cleanup Variables (some of these could be cleared earlier to save memory)
clear tFileName sFileName mFileName iFileName nonC iterations dimReturnsThreshold verbose
clear t s m w h chan 
clear S0 M0 T0 K Lt Ls LtAmp LsAmp LgVF Lgm Lm L0 I 
clear Li0 LiN Ifinal In Progress0 Lerror Progress ProgressMax ProgressD

end

