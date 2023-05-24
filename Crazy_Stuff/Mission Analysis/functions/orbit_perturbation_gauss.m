function orbit = orbit_perturbation_gauss(a,e,i,RAAN,w,TA,n_period,plot_kep,j2_pert,ref_frame,Cr,Asc_m_ratio,JD0,Ttot)
% Evaluates the propagation and perturbation of an orbit using Gauss planetary
% equations, solved by ode113 (tollerance: 'reltol', 1.e-13, 'abstol', 1.e-14)
% J2 effect and solar radiation pressure can be considered
% For J2 two different reference frame can be chosen: 'RSW' or 'TNH'
% Solar radiation pressure is considered only in 'RSW' frame
% If solar radiation pressure is considered function need a date, expressed
% in julian days or date matrix (1,3) or (1,6)
% The output is a structure with all the data computed
% 
% Example:
% - orbit_perturbation_gauss(a,e,i): orbit propagation over 10 periods using
% random value for RAAN,w,TA. Only J2 effect is considered.
% - orbit_perturbation_gauss(a,e,i,RAAN,w,TA): orbit propagation over 10 periods.
% J2 effect is considered
% 
% Inputs:
% 1.  a             : semi-major axis                       (km)
% 2.  e             : eccentricity
% 3.  i             : inclination                           (deg)
% 4.  RAAN          : right ascension of ascending node     (deg)   (default, RAAN = 360*rand)
% 5.  w             : argument of perigee                   (deg)   (default, w = 360*rand)
% 6.  TA            : true anomaly                          (deg)   (default, TA = 360*rand)
% 7.  n_period      : number of period repetions            (default, n_period = 10)
% 8.  plot_kep      : true or false                         (default, plot_kep = false)
% 9.  j2_pert       : effect of J2 perturbation             (true or false)
% 10. ref_frame     : 'RSW' or 'TNH'                        (default, ref_frame = 'RSW')
% 11. Cr            : radiation pressure coefficient        (1<Cr<2, default Cr=0, no solar radiation)
% 12. Asc_m_ratio   : satellite's Area to mass ratio        (m^2/kg)
% 13. JD0           : initial Julian day                    (days or date)
% 14. Ttot          : Total interval of time of evaluation  (sec)
% 
% Outputs:
% orbit             : structure containing all the data evaluated
% orbit.gauss       : substructure with keplerian element propagation,
%                       position vector and velocit vector from Gauss eqs
%
% Functions required:
% - config()
% - gauss_method()
% - j2_pert()
% - solar_rad()
% - solar_rad_gauss()
% - filter_orb()
% - subplot_kep_element()
% - plot_orbit_ev()
% 
% Contributors:
% Duma Francesco
% Gaballo Paolo
% Merola Pierpaolo
% Parducci Alessandro
% 
% Versions:
% new version of orbit_perturbation(): 2021-02-13, only for Gauss eqs


if nargin < 4               %no solar radiadion
    RAAN            = 360*rand;
    w               = 360*rand;
    TA              = 360*rand;         
    n_period        = 10;
    plot_kep        = false;
    j2_pert         = true;
    ref_frame       = 'RSW';
    Cr              = 0;
    Asc_m_ratio     = 0; 
    JD0             = juliandate(datetime('today'));
    Ttot            = 0;
elseif nargin < 7           %no solar radiadion
    n_period        = 10;
    plot_kep        = false;
    j2_pert         = true;
    ref_frame       = 'RSW';
    Cr              = 0;
    Asc_m_ratio     = 0; 
    JD0             = juliandate(datetime('today'));
    Ttot            = 0;
elseif nargin < 8           %no solar radiadion
    plot_kep        = false;
    j2_pert         = true;
    ref_frame       = 'RSW';
    Cr              = 0;
    Asc_m_ratio     = 0; 
    JD0             = juliandate(datetime('today'));
    Ttot            = 0;
elseif nargin < 9           %no solar radiadion
    j2_pert         = true;
    ref_frame       = 'RSW';
    Cr              = 0;
    Asc_m_ratio     = 0; 
    JD0             = juliandate(datetime('today'));
    Ttot            = 0;
elseif nargin < 10          %no solar radiadion
    ref_frame       = 'RSW';
    Cr              = 0;
    Asc_m_ratio     = 0; 
    JD0             = juliandate(datetime('today'));
    Ttot            = 0;
elseif nargin < 11
    Cr              = 0;
    Asc_m_ratio     = 0;
    JD0             = juliandate(datetime('today'));
    Ttot            = 0;
elseif nargin < 13
    JD0             = juliandate(datetime('today'));
    Ttot            = 0;
    elseif nargin < 14
    Ttot            = 0;
end
if isempty(RAAN)
    RAAN    = 360*rand;
end
if isempty(w)
    w       = 360*rand;
end
if isempty(TA)
    TA      = 360*rand;
end
if isempty(n_period)
    n_period = 10;
end
if isempty(plot_kep)
    plot_kep = false;
end
if isempty(j2_pert)
    plot_kep = false;
end
if isempty(ref_frame)
    ref_frame = 'RSW';
end
if isempty(Cr)
    Cr          = 0;
end
if isempty(Asc_m_ratio)
    Asc_m_ratio = 0;
end
if isempty(JD0)
    JD0 = juliandate(datetime('today'));
end
if isempty(Ttot)
    Ttot = 0;
end

% Configuration constants and initial values
orbit = config(a,e,RAAN,i,w,TA,n_period,Cr,Asc_m_ratio,JD0);

% Useful constants:
deg     = orbit.const.deg;          %Degrees to radians
mu      = orbit.const.mu;           %Gravitational parameter (km^3/s^2)
       
T0      = orbit.T0;
kep0    = orbit.kep0.kep0;

%% SETTING
% nout : Number of solution points to output for plotting purposes

t0 = 0;
if Ttot ~= 0
    n_period = Ttot/orbit.T0;
    orbit.n_period = n_period;
    tf      = Ttot;
    if tf > 1e5
        nout    = Ttot/T0*360;       % at least one point every 10 deg
    else
        nout    = Ttot;
    end
else
    tf      = n_period*T0;
    nout    = n_period*2000;
end

tspan       = linspace(t0, tf, nout);
orbit.tspan = tspan;
options     = odeset('reltol', 1.e-13, 'abstol', 1.e-14, 'initialstep', T0/1000);
%% GAUSS EQUATIONS
% Use ODE113 to integrate the Gauss variational equations from t0 to tf:

if Cr == 0      %no solar radiation pressure
    if j2_pert == true
        orbit = gauss_method(tspan,options,orbit,@j2_pert,[],ref_frame);
    else
        orbit = gauss_method(tspan,options,orbit,[],[],ref_frame);
    end
else
    if j2_pert == true
        orbit = gauss_method(tspan,options,orbit,@j2_pert,@solar_rad_gauss,ref_frame);
    else
        orbit = gauss_method(tspan,options,orbit,[],@solar_rad_gauss,ref_frame);
    end
end

%% FILTER
% Filter keplerian elements (a,e,RAAN,i,w,TA)
% (unwrapping should not be required, as the come from integration with Gauss equation)
% use function 'movemean' as filter

% Sample frequency
f_sample    = 1/sum(diff(orbit.tspan)/length(orbit.tspan));
orbit.f_sample = f_sample;
% Cut_off frequency
c_off_long          = 1/(10*orbit.T0);
orbit.c_off_long    = c_off_long;
c_off_sec           = 1/(2*365*24*3600); %two years
orbit.c_off_sec     = c_off_sec;
% Filter
filt_long   = filter_orb(orbit.gauss.kep(:,1:6),f_sample,c_off_long);
orbit.filt_long = filt_long;
filt_sec    = filter_orb(orbit.gauss.kep(:,1:6),f_sample,c_off_sec);
orbit.filt_sec = filt_sec;
% Plot
subplot_kep_element(orbit.tspan,orbit.gauss,'J2+SRP',[],[],[],filt_long,filt_sec)
%% WRAP True Anomaly

orbit.gauss.kep(:,6) = mod(orbit.gauss.kep(:,6),360);

%% PLOT
% Plot keplerian elements evolution

if plot_kep == true
    subplot_kep_element(orbit.tspan,orbit.gauss,'Gauss Eqs')
end


%% Check Perigee radius

orbit.gauss.rp  = orbit.gauss.kep(:,1).*(1-orbit.gauss.kep(:,2));
rp              = orbit.gauss.rp;
orbit.rp_min    = min(rp);

if orbit.rp_min < orbit.const.RE
    error("Perigee radius is lower than Earth's radius. Check initial data (semi-major axis and eccentricity)")    
elseif orbit.rp_min < orbit.const.RE*1.015
    dif = orbit.rp_min - orbit.const.RE;
    warning("The difference between Perigee radius and Earth's radius is only %.3f km. Check initial data (semi-major axis and eccentricity)",dif)
end

%% PLOT ORBIT

h_vect=cross(orbit.gauss.r_vect,orbit.gauss.v_vect);
orbit.e_vect = cross(orbit.gauss.v_vect,h_vect)/mu - orbit.gauss.r_vect./sqrt((orbit.gauss.r_vect(:,1)).^2+(orbit.gauss.r_vect(:,2)).^2+(orbit.gauss.r_vect(:,3)).^2);

a0      = orbit.kep0.a;
e0      = orbit.kep0.e;
i0      = orbit.kep0.i * 180/pi;
RAAN0   = orbit.kep0.RAAN * 180/pi;
w0      = orbit.kep0.w * 180/pi;
TA0     = orbit.kep0.TA * 180/pi;

plot_orbit_ev(a0,e0,i0,RAAN0,w0,TA0,orbit.gauss.r_vect,orbit.gauss.kep(:,6));

% figure
% comet3(x,y,z)
end
