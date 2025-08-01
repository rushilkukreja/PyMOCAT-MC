% ------------------------------------------------------------------------------
%
%                           function angl
%
%  this function finds the angle between two vectors.  the output is set to
%    zero if either vector is zero.
%
%  author        : david vallado                  719-573-2600   27 may 2002
%
%  revisions
%                -
%
%  inputs          description                    range / units
%    vec1        - vector number 1
%    vec2        - vector number 2
%
%  outputs       :
%    theta       - angle between the two vectors  -pi to pi
%
%  locals        :
%    temp        - temporary real variable
%
%  coupling      :
%    dot           dot product of two vectors
%    mag           magnitude of a vector
%
%  references    :
%    vallado       2007, 193, Alg 23, Ex 4-1
%
% [theta] = angl ( vec1,vec2 );
% ------------------------------------------------------------------------------

function theta = angl ( vec1,vec2 )

% -------------------------  implementation   -----------------
magv1 = mag(vec1);
magv2 = mag(vec2);

if ( magv1*magv2 > 1.0e-12 )
    temp= dot(vec1,vec2) / (magv1*magv2);
    if ( abs( temp ) > 1.0 )
        temp= sign(temp) * 1.0;
    end
    theta= acos( temp );
else
    theta= 999999.1;
end