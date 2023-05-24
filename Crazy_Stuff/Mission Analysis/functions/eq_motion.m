function fun=eq_motion(t,r,par,p_J2,p_solar)
% Gravitational equation of motion in cartesian coordinates
%
% Inputs:
% t         : time istant for ode evaluation (sec)
% f         : vector of functions for ode evaluation
% par       : structure with data and constants
% p_J2      : effect of J2 perturbation
% p_solar   : solar radiation pressure
%
% Outputs:
% fun       : equation of motion in cartesian coordinates
%
% Functions required: (only for solar radioation pressure)
% - n_J2000()
% - solar_coord()
% - shadow_fun_cart()
%
% Contributors:
% Gaballo Paolo
% 
% Versions:
% 2021-02-13, first version


if nargin == 4
    p_solar = [];
elseif nargin == 3
    p_solar = [];
    p_J2 = [];
end

J2 = par.const.J2;
R = par.const.R;
mu = par.const.mu;

% r_norm=sqrt(r(1)^2+r(2)^2+r(3)^2);
r_vect = [r(1),r(2),r(3)];
r_norm = norm(r_vect);

%%
% J2 perturbation:
if isempty(p_J2)
    J2_pert = [0;0;0];
else
    fac=3/2*(J2*mu*R^2)/(r_norm)^4;
    J2_pert = [fac*(r(1)/r_norm*(5*r(3)^2/r_norm^2-1));
            fac*(r(2)/r_norm*(5*r(3)^2/r_norm^2-1));
            fac*(r(3)/r_norm*(5*r(3)^2/r_norm^2-3))];
end
    
% solar radiation perturbation:
if isempty(p_solar)
    solar_pert = [0;0;0];
else
    Asc_m_ratio = par.Asc_m_ratio;
    Cr          = par.Cr;
    Psr         = par.Psr;
    
    n_Jday      = n_J2000(t,par);                           %(Jdate)
    par.n_Jday  = n_Jday;
    [L,M,lambda,eps] = solar_coord(par);
    ni = shadow_fun_cart(t,r,par,M,lambda,eps);
    
    p_sr        = ni * Psr * Cr * Asc_m_ratio/1000;              %magnite of the perturbation
    solar_pert  = -p_sr .* [cos(lambda); cos(eps)*sin(lambda); sin(eps)*sin(lambda)];
end

pert = J2_pert + solar_pert;

fun=[r(4);
    r(5);
    r(6);
    -(mu/r_norm^3)*r(1)+pert(1);
    -(mu/r_norm^3)*r(2)+pert(2);
    -(mu/r_norm^3)*r(3)+pert(3)];

end