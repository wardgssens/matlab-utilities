% Conversion from WGS84 coordinates to WGS84 Lambert 1972
function [x, y] = GPStoLam72(lat, lon)
    % Ellipsoïd HAYFORD 1924 (= International)
    a = 6378388;
    f = 1 / 297;
    
    % Projection parameters
    phi_1 = dms2degrees(49, 50, 0.00204);
    phi_2 = dms2degrees(51, 10, 0.00204);
    phi_0 = 90;
    
    lambda_0 = dms2degrees(4, 22, 2.952);
    x_0 = 150000.013;
    y_0 = 5400088.438;
   
    % Excentricity
    e2 = 2*f - f^2;
    e = sqrt(e2);
    
    % Calculate parameters
    m_1 = cosd(phi_1) / (1 - e2 * sind(phi_1)^2)^0.5;
    m_2 = cosd(phi_2) / (1 - e2 * sind(phi_2)^2)^0.5;
    
    t_i = @(phi) tand(45 - phi/2) / ((1 - e * sind(phi)) / (1 + e * sind(phi)))^(e / 2);
    t_1 = t_i(phi_1);
    t_2 = t_i(phi_2);
    t_0 = t_i(phi_0);
    
    n = (log(m_1) - log(m_2)) / (log(t_1) - log(t_2));
    
    g = m_1 / (n * t_1^n);
    
    r_0 = a * g * t_0^n;
   
    % Transformation
    t = tand(45 - lat/2) ./ ((1 - e .* sind(lat)) ./ (1 + e .* sind(lat))).^(e/2);
    r = a .* g .* t.^n;
    theta = n .* (lon - lambda_0);
    x = x_0 + r .* sind(theta);
    y = y_0 + r_0 - r .* cosd(theta);
end

function dd = dms2degrees(d, m, s)
    dd = d + m / 60 + s / 3600;
end