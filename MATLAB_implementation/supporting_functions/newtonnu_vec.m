% ------------------------------------------------------------------------------
%
%                           function newtonnu_vec
%
%  this function finds the mean anomaly given the eccentricity and true anomaly
%  (vectorized version)
%
%  author        : adapted from Vallado for vectorization
%
%  inputs          description                    range / units
%    ecc         - eccentricity                   0.0 to inf (N x 1 vector)
%    nu          - true anomaly                   0.0  to 2pi rad (N x 1 vector)
%
%  outputs       :
%    e0          - eccentric anomaly              0.0  to 2pi rad (N x 1 vector)
%    m           - mean anomaly                   0.0  to 2pi rad (N x 1 vector)
%
%  locals        :
%    small       - tolerance for near-zero values
%    sine        - sine of eccentric anomaly
%    cose        - cosine of eccentric anomaly
%
%  coupling      :
%    none
%
%  references    :
%    vallado       2007, 85-86, alg 5
%
% [e0,m] = newtonnu_vec ( ecc, nu );
% ------------------------------------------------------------------------------

function [e0, m] = newtonnu_vec(ecc, nu)

% -------------------------  implementation   -----------------
small = 0.00000001;

% Initialize output arrays
e0 = zeros(size(ecc));
m = zeros(size(ecc));

% Handle different orbit types
for i = 1:length(ecc)
    if ecc(i) < small
        % Circular orbit
        e0(i) = nu(i);
        m(i) = nu(i);
    elseif ecc(i) < 1.0 - small
        % Elliptical orbit
        sine = (sqrt(1.0 - ecc(i)*ecc(i)) * sin(nu(i))) / (1.0 + ecc(i)*cos(nu(i)));
        cose = (ecc(i) + cos(nu(i))) / (1.0 + ecc(i)*cos(nu(i)));
        e0(i) = atan2(sine, cose);
        m(i) = e0(i) - ecc(i)*sin(e0(i));
    elseif abs(ecc(i) - 1.0) < small
        % Parabolic orbit
        e0(i) = tan(nu(i) * 0.5);
        m(i) = e0(i) + (e0(i)^3)/3.0;
    else
        % Hyperbolic orbit
        sine = (sqrt(ecc(i)*ecc(i) - 1.0) * sin(nu(i))) / (1.0 + ecc(i)*cos(nu(i)));
        e0(i) = asinh(sine);
        m(i) = ecc(i)*sinh(e0(i)) - e0(i);
    end
end