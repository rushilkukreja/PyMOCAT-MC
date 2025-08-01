% ------------------------------------------------------------------------------
%
%                           function angl_vec
%
%  this function finds the angle between two vectors (vectorized version).
%  the output is set to zero if either vector is zero.
%
%  author        : david vallado (modified for vectorization)
%
%  inputs          description                    range / units
%    vec1        - vector number 1 (N x 3 matrix)
%    vec2        - vector number 2 (N x 3 matrix)
%
%  outputs       :
%    theta       - angle between the two vectors  -pi to pi (N x 1 vector)
%
%  locals        :
%    temp        - temporary real variable
%
%  coupling      :
%    dot           dot product of two vectors
%    sqrt          square root function
%
%  references    :
%    vallado       2007, 193, Alg 23, Ex 4-1 (vectorized)
%
% [theta] = angl_vec ( vec1,vec2 );
% ------------------------------------------------------------------------------

function theta = angl_vec ( vec1,vec2 )

% -------------------------  implementation   -----------------
% Calculate magnitude of vectors (vectorized)
magv1 = sqrt(sum(vec1.^2, 2));
magv2 = sqrt(sum(vec2.^2, 2));

% Initialize output
theta = zeros(size(magv1));

% Find valid indices where both vectors have non-zero magnitude
valid_idx = (magv1 .* magv2) > 1.0e-12;

if any(valid_idx)
    % Calculate dot product (vectorized)
    dot_product = sum(vec1(valid_idx,:) .* vec2(valid_idx,:), 2);
    
    % Calculate cosine of angle
    temp = dot_product ./ (magv1(valid_idx) .* magv2(valid_idx));
    
    % Clamp values to [-1, 1] to avoid numerical issues
    temp(temp > 1.0) = 1.0;
    temp(temp < -1.0) = -1.0;
    
    % Calculate angle
    theta(valid_idx) = acos(temp);
end

% Set invalid angles to large number
theta(~valid_idx) = 999999.1;