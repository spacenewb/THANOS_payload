function [delta,alpha] = CarToSph(r)
% Returns declination and right ascension from matrix of position vector r
%
% Inputs:
% r         : matrix of position vector (km)
% 
% Outputs:
% delta     : Declination       (rad)
% alpha     : Right Ascension   (rad)
%
% Contributors:
% Duma Francesco
% Gaballo Paolo
% Merola Pierpaolo
% Parducci Alessandro
% 
% Versions:
% 2021-02-13, first version

% cosine directors
l=zeros(size(r,1),1);
m=zeros(size(r,1),1);
n=zeros(size(r,1),1);

for i=1:size(r,1)
l(i,1)=r(i,1)/norm(r(i,:)); %direction x
m(i,1)=r(i,2)/norm(r(i,:)); %direction y
n(i,1)=r(i,3)/norm(r(i,:)); %direction z
end

% Declination delta
delta = asin(n);

% Right Ascension alpha
alpha = zeros(size(m,1),1);

for i=1:size(m,1) 
    if m(i,1) > 0
        alpha(i,1) = acos(l(i,1)/cos(delta(i,1)));
    else
        alpha(i,1) = 2*pi - acos(l(i,1)/cos(delta(i,1)));
    end
end

end


