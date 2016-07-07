function dist = greatcircledist1(varargin)
% GREATCIRCLEDIST1 Function to compute great-circle distance in degrees
%  Calculate great circle distance between points on a sphere using the
%  Haversine Formula.  LAT1, LON1, LAT2, and LON2 are in radians.  dist is
%  arclength in degrees.
[useGeodesic, lat1, lon1, lat2, lon2, ellipsoid, ...
    units, insize, useAngularDistance] = parseDistAzInputs(varargin{:});

a = sin((lat2-lat1)/2).^2 + cos(lat1) .* cos(lat2) .* sin((lon2-lon1)/2).^2;
rng = r * 2 * atan2(sqrt(max(a,0)),sqrt(max(1 - a,0)));
dist = fromRadians('degrees',rng);

end