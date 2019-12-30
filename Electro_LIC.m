%
% LIC - Line Integral Convolution
% 
% Author: Davide Marcantonio
% Last update: 02 June 2016
%

%% clean the environment
close all;
% clear all;
clc;

show = true;

%% setting parameters for electric (E) field computation
epsilon_zero = 8.85e-12;  % empty space dielectric constant
k = 1 / (4 * pi * epsilon_zero);  % see electric field generated form a point source for this
delta_space = 1e-1;  % physical distance [m] from one "pixel center" to the adjacent ones

%% image parameters
image_width = 600;  % number of pixel rows
image_height = 600;  % number of pixel columns

%% compute the 2D electric field 
% of equally charged 'point' sources disposed in a 2D plane

field_lin = zeros(image_width, image_height, 2);  % initialize zero electric field 
field_modulus = zeros(image_width, image_height);  % initialize the modulus only image

% create a figure in which are plotted the charges created:
% RED => positive charge
% BLUE => negative charge
if show
    figure;  % new view
    axis([0 image_width 0 image_height]);  % needed to keep the same scale of the images
    axis ij;  % flip to see the sources as in the images (otherwise the y axis is mirrored)
    title('point charges positions (RED = "+" , BLUE = "-")');  
    hold on;  % allows the points of the charges to accumulate over the image
end

tic

n_charges = 6;  % number of charges to randomly allocate
for i = 1 : n_charges  % iter over the charges
    
    % generate a random position (x,y) for the new charge 
    % the code below forces the position to be 'in the center'
    x_charge = (rand() * image_width/2 + image_width/4); 
    y_charge = (rand() * image_height/2 + image_height/4);
    
    q = 1e-6;  % fixed electric charge modulus
    % 50% probability for the charge to become negative
    if rand() > .5
        q = -q;
    end
    
    % prints the generated charge on the plot
    if show
        if q > 0
            plot(y_charge, x_charge, 'ro');
        else
            plot(y_charge, x_charge, 'bo');
        end
    end
    
    % for each pixel of the image, compute the electric field in that point
    % relatively to the currently generated charge
    for x = 1 : image_width
        for y = 1 : image_height
            
            % compute the distance of the current pixel position with
            % respect to the charge position, scaling to make them 'physical'
            dist_x = (x - x_charge) * delta_space;
            dist_y = (y - y_charge) * delta_space;
            
            x_versor = (dist_x) / (abs(dist_x) + abs(dist_y));
            y_versor = (dist_y) / (abs(dist_x) + abs(dist_y));
            
            % compute the squared linear distance between source and pixel
            r2 = dist_x ^ 2 + dist_y ^ 2;
            
            % compute the electric field value for the current pixel
            value = (k * q) / r2;
            
            % add the field computed to the overall field in this point 
            % (superposition principle)
            field_lin(x, y, 1) = field_lin(x, y, 1) + value * x_versor; 
            field_lin(x, y, 2) = field_lin(x, y, 2) + value * y_versor; 
            
            % update the modulus image
            field_modulus(x, y) =  field_modulus(x,y) + abs(value);
        end
    end
end
disp('Finished cycle.');
toc 

%% compute the logarithmic field (for visualization purposes)
MIN = min(min(field_modulus));
field_log = log((field_modulus - MIN) + 1);  % '1' is due for the log to generate positive numbers

% normalize the logarithmic field to be in [0, 1]
field_log = field_log - min(min(field_log));
field_log = field_log / max(max(field_log));

% show the log field
if show
    figure;
    imshow(field_log);
    title('Field modulus (log)');
end

%% init images
new_image = zeros(image_width, image_height, 3);
img1 = zeros(image_width,image_height);
img2 = zeros(image_width,image_height);

%% create convolutional kernel
kernel_size = 110; 
disp(ceil(sqrt(image_width/4 * image_height/4 * .25)));

kernel1 = ones(kernel_size) / kernel_size;

x = 0 : pi / (2 * kernel_size) : pi/2;
kernel2 = cos(x);
kernel2 = kernel2 / (sum(kernel2));

%% compute LIC
random_image = rand(image_width, image_height) - 0.5;
figure;
imshow(random_image + 0.5);
title('random image');

iterations = 1;  % on cos

for j = 1 : iterations
    for x = 1 : image_width
        for y = 1 : image_height
            if rand() > 0.9999
                disp([num2str(((x * (image_height-1) + y)/(image_width * image_height))*100),' %']);
            end
            current_pixel1 = random_image(x,y) * kernel1(1);
            current_pixel2 = random_image(x,y) * kernel2(1);
            bias_x = 0;
            bias_y = 0;
            xp = x;
            yp = y;
            for w = 2 : kernel_size
                theta = atan(field_lin(xp, yp, 1) / field_lin(xp, yp, 2));

                cos_theta = cos(theta);
                sin_theta = sin(theta);
                    
                thr = 0;
                if sin_theta >  thr
                    bias_x = bias_x + abs(sin_theta);
                elseif sin_theta < - thr
                    bias_x = bias_x - abs(sin_theta);
                end

                if cos_theta > thr
                    bias_y = bias_y + abs(cos_theta);
                elseif cos_theta < - thr;
                    bias_y = bias_y - abs(cos_theta);
                end
                
                % do not go beyond the matrix
                xpc = x + round(bias_x);
                if (xpc > image_width)
                    xpc = image_width;
                end
                if (xpc <= 0)
                    xpc = 1;
                end

                ypc = y + round(bias_y);
                if (ypc > image_height)
                    ypc = image_height;
                end
                if (ypc <= 0)
                    ypc = 1;
                end


                NORM = (1 - abs(cos_theta)) * (1 - abs(sin_theta));
                current_pixel1 = current_pixel1 + random_image(xp, yp) * kernel1(w) * NORM;
                current_pixel2 = current_pixel2 + random_image(xp, yp) * kernel2(w) * NORM;


                NORM = abs(cos_theta);
                current_pixel1 = current_pixel1 + random_image(xpc, yp) * kernel1(w) * NORM;
                current_pixel2 = current_pixel2 + random_image(xpc, yp) * kernel2(w) * NORM;


                NORM = abs(sin_theta);
                current_pixel1 = current_pixel1 + random_image(xp, ypc) * kernel1(w) * NORM;
                current_pixel2 = current_pixel2 + random_image(xp, ypc) * kernel2(w) * NORM;
               
                
                NORM = abs(cos_theta) * abs(sin_theta);
                current_pixel1 = current_pixel1 + random_image(xpc, ypc) * kernel1(w) * NORM;
                current_pixel2 = current_pixel2 + random_image(xpc, ypc) * kernel2(w) * NORM;

                xp = xpc;
                yp = ypc;
                    
               
            end
            
            new_image(x, y, 1) = (current_pixel2 + .5) * 0.5 + (field_log(x,y)) * 0.5;
            new_image(x, y, 2) = (current_pixel2 + .5) * 0.5 + (field_log(x,y)) * 0.5; 
            new_image(x, y, 3) = (current_pixel2 + .5) * 0.5 + (1 - field_log(x,y)) * 0.5;

            img1(x,y) = current_pixel1;
            img2(x,y) = current_pixel2;
        end
    end
    % for the iterations
    random_image = img2;
end

if show
    figure;
    imshow(img1 + 0.5);
    title('LIC (LPF)');

    figure;
    imshow(img2 + 0.5);
    title('LIC (cos)');

    figure;
    imshow(new_image);    
end

filename = sprintf('E_field_LIC_cos.tiff');             % save images of blocks in 'blocks/' folder
imwrite(new_image, filename);

filename = sprintf('E_field_LIC_cos.csv');
csvwrite(filename, img2);

% filename = 'LIC.log';
% fileID = fopen(filename, 'a');
% 
% formatSpec = 'compute_sift:\n';
% fprintf(fileID, formatSpec);
