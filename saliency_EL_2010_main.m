function saliency_EL_2010_main
%
% This is a demo program of the paper L. Ma, J. Tian, and W. Yu,
% "Visual saliency detection in image using ant colony optimisation and local phase coherence," 
% Electronics Letters, Vol. 46, Jul. 2010, pp. 1066-1068.
%

clear all; close all; clc;

% Load test image
img_in = double(imread('test.jpg'));

% Perform saliency detection
mask = func_saliency_aco_phase(img_in);

% Write the output saliency map
imwrite(uint8(mask),'test_saliency_map.bmp','bmp');

%-------------------------------------------------------------------------
%------------------------------Inner Function ----------------------------
%-------------------------------------------------------------------------
function p = func_saliency_aco_phase(img)

[nrow, ncol] = size(img);
img_1D = func_2D_LexicoOrder(img);

% initialization
p = 0.001 .* ones(size(img));
p_1D = func_2D_LexicoOrder(p);

delta_p = zeros(size(img));
delta_p_1D = func_2D_LexicoOrder(delta_p);

%paramete setting
alpha = 1;
beta = 2;       
phi = 0.3;      

ant_total_num = round(sqrt(nrow*ncol));
ant_current_row = zeros(ant_total_num, 1); % record the location of ant
ant_current_col = zeros(ant_total_num, 1); % record the location of ant
rand('state', sum(clock));
temp = rand(ant_total_num, 2);
ant_current_row = round(1 + (nrow-1) * temp(:,1)); %row index
ant_current_col = round(1 + (ncol-1) * temp(:,2)); %column index
ant_current_val = img_1D((ant_current_row-1).*ncol+ant_current_col);

% build v matrix

v = zeros(size(img));
v_norm = 0;
for rr =1:nrow
    for cc=1:ncol
        %defination of clique
        temp1 = [rr-2 cc-1; rr-2 cc+1; rr-1 cc-2; rr-1 cc-1; rr-1 cc; rr-1 cc+1; rr-1 cc+2; rr cc-1];
        temp2 = [rr+2 cc+1; rr+2 cc-1; rr+1 cc+2; rr+1 cc+1; rr+1 cc; rr+1 cc-1; rr+1 cc-2; rr cc+1];
        temp0 = find(temp1(:,1)>=1 & temp1(:,1)<=nrow & temp1(:,2)>=1 & temp1(:,2)<=ncol & temp2(:,1)>=1 & temp2(:,1)<=nrow & temp2(:,2)>=1 & temp2(:,2)<=ncol);
        temp11 = temp1(temp0, :);
        temp22 = temp2(temp0, :);
        temp00 = zeros(size(temp11,1));
        for kk = 1:size(temp11,1)
            temp00(kk) = abs(img(temp11(kk,1), temp11(kk,2))-img(temp22(kk,1), temp22(kk,2)));
        end
        if size(temp11,1) == 0
            v(rr, cc) = 0;
            v_norm = v_norm + v(rr, cc);
        else
            v(rr, cc) = max(max(temp00));
        end
    end
end
% Calculate the phase
v = func_phasecong(img);
v_1D = func_2D_LexicoOrder(v);

% System setup
ant_move_step_within_iteration = 300;
total_iteration_num = 3;
search_clique_mode = 8;  

for nIteration = 1: total_iteration_num               

        ant_current_path_row = zeros(ant_total_num,ant_move_step_within_iteration);
        ant_current_path_col = zeros(ant_total_num,ant_move_step_within_iteration);
        ant_current_path_val = zeros(ant_total_num,ant_move_step_within_iteration);
        ant_current_path_row(:,1) = ant_current_row;                
        ant_current_path_col(:,1) = ant_current_col;                
        ant_current_path_val(:,1) = ant_current_val;            

        for nMoveStep = 1: ant_move_step_within_iteration-1                

            if search_clique_mode == 4
                ant_search_range_row = [ant_current_row-1, ant_current_row, ant_current_row+1, ant_current_row];
                ant_search_range_col = [ant_current_col, ant_current_col+1, ant_current_col, ant_current_col-1];                    
            elseif search_clique_mode == 8
                ant_search_range_row = [ant_current_row-1, ant_current_row-1, ant_current_row-1, ant_current_row, ant_current_row,ant_current_row+1, ant_current_row+1, ant_current_row+1];
                ant_search_range_col = [ant_current_col-1, ant_current_col, ant_current_col+1, ant_current_col-1, ant_current_col+1, ant_current_col-1, ant_current_col, ant_current_col+1];
            end

            ant_current_row_extend = padarray(ant_current_row, [0 search_clique_mode-1],'replicate','post');
            ant_current_col_extend = padarray(ant_current_col, [0 search_clique_mode-1],'replicate','post');
            ant_search_range_val = zeros(ant_total_num,search_clique_mode);

            %replace the positions our of the image's range
            temp = (ant_search_range_row>=1) & (ant_search_range_row<=nrow) & (ant_search_range_col>=1) & (ant_search_range_col<=ncol);
            ant_search_range_row = temp.*ant_search_range_row + (~temp).*ant_current_row_extend;
            ant_search_range_col = temp.*ant_search_range_col + (~temp).*ant_current_col_extend;

            ant_search_range_transit_prob_h = zeros(size(ant_search_range_val));
            ant_search_range_transit_prob_p = zeros(size(ant_search_range_val));                                

            for ii=1:search_clique_mode
                ant_search_range_transit_prob_h(:,ii) = v_1D((ant_search_range_row(:,ii)-1).*ncol+ant_search_range_col(:,ii));
                ant_search_range_transit_prob_p(:,ii) = p_1D((ant_search_range_row(:,ii)-1).*ncol+ant_search_range_col(:,ii));

            end
            temp = (ant_search_range_transit_prob_h.^alpha) .* (ant_search_range_transit_prob_p.^beta);

            temp_sum = ((sum(temp'))');
            temp_sum = padarray(temp_sum, [0 search_clique_mode-1],'replicate','post');

            ant_search_range_transit_prob = temp ./ temp_sum;

            % generate a random number to determine the next position.
            rand('state', sum(100*clock));
            temp = rand(ant_total_num,1);
            temp = padarray(temp, [0 search_clique_mode-1],'replicate','post');
            temp = cumsum(ant_search_range_transit_prob,2)>=temp;
            temp = padarray(temp, [0 1],'pre');
            temp = (diff(temp,1,2)==1);

            temp_row = (ant_search_range_row .* temp)';
            [ii, jj, vv] = find(temp_row);
            ant_next_row = vv;
            temp_col = (ant_search_range_col .* temp)';
            [ii, jj, vv] = find(temp_col);
            ant_next_col = vv;

            ant_current_path_val(:,nMoveStep+1) = img_1D((ant_next_row-1).*ncol+ant_next_col);
            ant_current_path_row(:,nMoveStep+1) = ant_next_row;
            ant_current_path_col(:,nMoveStep+1) = ant_next_col;

            ant_current_row = ant_next_row;
            ant_current_col = ant_next_col;                
            ant_current_val = img_1D((ant_current_row-1).*ncol+ant_current_col);

            rr = ant_current_row;
            cc = ant_current_col;
            delta_p_1D((rr-1).*ncol+cc,1) = 1;
            delta_p = func_LexicoOrder_2D(delta_p_1D, nrow, ncol);
        end % end of nMoveStep

        p = (1-phi).*p + delta_p.*v.*(v>mean(v(:)));          
        p_1D = func_2D_LexicoOrder(p);
        delta_p = zeros(size(img));
        delta_p_1D = func_2D_LexicoOrder(delta_p);

end % end of nIteration

% Perform normalization to be range [0,255];
p = p./(sum(sum(p)));
p = (p-min(p(:)))./(max(p(:))-min(p(:))).*255;

%-------------------------------------------------------------------------
%------------------------------Inner Function ----------------------------
%-------------------------------------------------------------------------
% Convert data from 2D format to 1D format
function result = func_2D_LexicoOrder(x)
temp = x';
result = temp(:);

%-------------------------------------------------------------------------
%------------------------------Inner Function ----------------------------
%-------------------------------------------------------------------------
% Convert data from 1D format to 2D format
function result = func_LexicoOrder_2D(x, nRow, nColumn)
result = reshape(x, nColumn, nRow)';

%-------------------------------------------------------------------------
%------------------------------Inner Function ----------------------------
%-------------------------------------------------------------------------
% Calculate local phase coherence at each pixel position
function phaseCongruency = func_phasecong(im)                    

sze = size(im);

if nargin < 2
    nscale          = 3;     % Number of wavelet scales.
end
if nargin < 3
    norient         = 6;     % Number of filter orientations.
end
if nargin < 4
    minWaveLength   = 3;     % Wavelength of smallest scale filter.
end
if nargin < 5
    mult            = 2;     % Scaling factor between successive filters.
end
if nargin < 6
    sigmaOnf        = 0.55;  % Ratio of the standard deviation of the
                             % Gaussian describing the log Gabor filter's transfer function 
			     % in the frequency domain to the filter center frequency.
end
if nargin < 7
    dThetaOnSigma   = 1.7;   % Ratio of angular interval between filter orientations
			     % and the standard deviation of the angular Gaussian
			     % function used to construct filters in the
                             % freq. plane.
end
if nargin < 8
    k               = 3.0;   % No of standard deviations of the noise energy beyond the
			     % mean at which we set the noise threshold point.
			     % standard deviation to its maximum effect
                             % on Energy.
end
if nargin < 9
    cutOff          = 0.4;   % The fractional measure of frequency spread
                             % below which phase congruency values get penalized.
end
   
g               = 10;    % Controls the sharpness of the transition in the sigmoid
                         % function used to weight phase congruency for frequency
                         % spread.
epsilon         = .0001; % Used to prevent division by zero.


thetaSigma = pi/norient/dThetaOnSigma;  % Calculate the standard deviation of the
                                        % angular Gaussian function used to
                                        % construct filters in the freq. plane.

imagefft = fft2(im);                    % Fourier transform of image
sze = size(imagefft);
rows = sze(1);
cols = sze(2);
zero = zeros(sze);

totalEnergy = zero;                     % Matrix for accumulating weighted phase 
                                        % congruency values (energy).
totalSumAn  = zero;                     % Matrix for accumulating filter response
                                        % amplitude values.
orientation = zero;                     % Matrix storing orientation with greatest
                                        % energy for each pixel.
estMeanE2n = [];

% Pre-compute some stuff to speed up filter construction

x = ones(rows,1) * (-cols/2 : (cols/2 - 1))/(cols/2);  
y = (-rows/2 : (rows/2 - 1))' * ones(1,cols)/(rows/2);
radius = sqrt(x.^2 + y.^2);       % Matrix values contain *normalised* radius from centre.
radius(round(rows/2+1),round(cols/2+1)) = 1; % Get rid of the 0 radius value in the middle 
                                             % so that taking the log of the radius will 
                                             % not cause trouble.
theta = atan2(-y,x);              % Matrix values contain polar angle.
                                  % (note -ve y is used to give +ve
                                  % anti-clockwise angles)
sintheta = sin(theta);
costheta = cos(theta);
clear x; clear y; clear theta;      % save a little memory

% The main loop...

for o = 1:norient,                   % For each orientation.
%   disp(['Processing orientation ' num2str(o)]);
  angl = (o-1)*pi/norient;           % Calculate filter angle.
  wavelength = minWaveLength;        % Initialize filter wavelength.
  sumE_ThisOrient   = zero;          % Initialize accumulator matrices.
  sumO_ThisOrient   = zero;       
  sumAn_ThisOrient  = zero;      
  Energy_ThisOrient = zero;      
  EOArray = [];          % Array of complex convolution images - one for each scale.
  ifftFilterArray = [];  % Array of inverse FFTs of filters

  ds = sintheta * cos(angl) - costheta * sin(angl); % Difference in sine.
  dc = costheta * cos(angl) + sintheta * sin(angl); % Difference in cosine.
  dtheta = abs(atan2(ds,dc));                           % Absolute angular distance.
  spread = exp((-dtheta.^2) / (2 * thetaSigma^2));      % Calculate the angular filter component.

  for s = 1:nscale,                  % For each scale.

    % Construct the filter - first calculate the radial filter component.
    fo = 1.0/wavelength;                  % Centre frequency of filter.
    rfo = fo/0.5;                         % Normalised radius from centre of frequency plane 
                                          % corresponding to fo.
    logGabor = exp((-(log(radius/rfo)).^2) / (2 * log(sigmaOnf)^2));  
    logGabor(round(rows/2+1),round(cols/2+1)) = 0; % Set the value at the center of the filter
                                                   % back to zero (undo the radius fudge).

    filter = logGabor .* spread;          % Multiply by the angular spread to get the filter.
    filter = fftshift(filter);            % Swap quadrants to move zero frequency 
                                          % to the corners.

    ifftFilt = real(ifft2(filter))*sqrt(rows*cols);  % Note rescaling to match power
    ifftFilterArray = [ifftFilterArray ifftFilt];    % record ifft2 of filter

    % Convolve image with even and odd filters returning the result in EO
    EOfft = imagefft .* filter;           % Do the convolution.
    EO = ifft2(EOfft);                    % Back transform.

    EOArray = [EOArray, EO];              % Record convolution result
    An = abs(EO);                         % Amplitude of even & odd filter response.

    sumAn_ThisOrient = sumAn_ThisOrient + An;     % Sum of amplitude responses.
    sumE_ThisOrient = sumE_ThisOrient + real(EO); % Sum of even filter convolution results.
    sumO_ThisOrient = sumO_ThisOrient + imag(EO); % Sum of odd filter convolution results.

    if s == 1                             % Record the maximum An over all scales
      maxAn = An;
    else
      maxAn = max(maxAn, An);
    end
    
    if s==1
      EM_n = sum(sum(filter.^2));           % Record mean squared filter value at smallest
    end                                     % scale. This is used for noise estimation.

    wavelength = wavelength * mult;         % Finally calculate Wavelength of next filter
  end                                       % ... and process the next scale

  % Get weighted mean filter response vector, this gives the weighted mean phase angle.

  XEnergy = sqrt(sumE_ThisOrient.^2 + sumO_ThisOrient.^2) + epsilon;   
  MeanE = sumE_ThisOrient ./ XEnergy; 
  MeanO = sumO_ThisOrient ./ XEnergy; 

  for s = 1:nscale,       
      EO = submat(EOArray,s,cols);  % Extract even and odd filter 
      E = real(EO); O = imag(EO);
      Energy_ThisOrient = Energy_ThisOrient ...
        + E.*MeanE + O.*MeanO - abs(E.*MeanO - O.*MeanE);
  end

  medianE2n = median(reshape(abs(submat(EOArray,1,cols)).^2,1,rows*cols));
  meanE2n = -medianE2n/log(0.5);
  estMeanE2n = [estMeanE2n meanE2n];

  noisePower = meanE2n/EM_n;                       % Estimate of noise power.

  % Now estimate the total energy^2 due to noise
  % Estimate for sum(An^2) + sum(Ai.*Aj.*(cphi.*cphj + sphi.*sphj))

  EstSumAn2 = zero;
  for s = 1:nscale
    EstSumAn2 = EstSumAn2+submat(ifftFilterArray,s,cols).^2;
  end

  EstSumAiAj = zero;
  for si = 1:(nscale-1)
    for sj = (si+1):nscale
      EstSumAiAj = EstSumAiAj + submat(ifftFilterArray,si,cols).*submat(ifftFilterArray,sj,cols);
    end
  end

  EstNoiseEnergy2 = 2*noisePower*sum(sum(EstSumAn2)) + 4*noisePower*sum(sum(EstSumAiAj));

  tau = sqrt(EstNoiseEnergy2/2);                     % Rayleigh parameter
  EstNoiseEnergy = tau*sqrt(pi/2);                   % Expected value of noise energy
  EstNoiseEnergySigma = sqrt( (2-pi/2)*tau^2 );

  T =  EstNoiseEnergy + k*EstNoiseEnergySigma;       % Noise threshold

  T = T/1.7;        % Empirical rescaling of the estimated noise effect to 
                    % suit the PC_2 phase congruency measure

  Energy_ThisOrient = max(Energy_ThisOrient - T, zero);  % Apply noise threshold

  width = sumAn_ThisOrient ./ (maxAn + epsilon) / nscale;    

  % Now calculate the sigmoidal weighting function for this orientation.

  weight = 1.0 ./ (1 + exp( (cutOff - width)*g)); 

  % Apply weighting

  Energy_ThisOrient =   weight.*Energy_ThisOrient;

  % Update accumulator matrix for sumAn and totalEnergy

  totalSumAn  = totalSumAn + sumAn_ThisOrient;
  totalEnergy = totalEnergy + Energy_ThisOrient;

  if(o == 1),
    maxEnergy = Energy_ThisOrient;
    featType = E + i*O;
  else
    change = Energy_ThisOrient > maxEnergy;
    orientation = (o - 1).*change + orientation.*(~change);
    featType = (E+i*O).*change + featType.*(~change);
    maxEnergy = max(maxEnergy, Energy_ThisOrient);
  end

end  % For each orientation


phaseCongruency = totalEnergy ./ (totalSumAn + epsilon);

orientation = orientation * (180 / norient);

featType = featType*i;   % Rotate feature phase angles by 90deg so that 0
                         % phase corresponds to a step edge (this is a
                         % fudge I must have something the wrong way
                         % around somewhere)

phaseCongruency = exp(phaseCongruency);
phaseCongruency = phaseCongruency./(sum(sum(phaseCongruency)));

%-------------------------------------------------------------------------
%------------------------------Inner Function ----------------------------
%-------------------------------------------------------------------------
function a = submat(big,i,cols)

a = big(:,((i-1)*cols+1):(i*cols));


