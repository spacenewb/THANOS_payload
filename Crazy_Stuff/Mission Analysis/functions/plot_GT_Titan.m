function [orbit] = plot_GT_Titan(orbit,pix_swath,transp1,transp2,colour,plan,name)
% Plots Ground tracks over the planisphere
% 
% Inputs:
% orbit : structure with data (longitude and latitude)
% colour : colour of the Ground Track
% plan : plot planisphere or not (true or false)
% name : Name of the Ground Track
%
% Outputs:
% orbit : structure with data
%
% Contributors:
% Duma Francesco
% Gaballo Paolo
% Merola Pierpaolo
% Parducci Alessandro
% 
% Versions:
% 2021-02-13, first version

if nargin < 5
    colour = 'c';
    plan = true;
    name = false;
elseif nargin < 6
    plan = true;
    name = false;
elseif nargin < 7
    name = false;
end
if isempty(colour)
    colour = 'r';
end

longitude   = orbit.coordinates.longitude;
latitude    = orbit.coordinates.latitude;

%% Plot Planisphere
if plan == true
    h = figure('Name','Ground Track');
    set(h, 'Units', 'Normalized', 'OuterPosition', [.15 .25 .7 .7]);
    hold on
    axis equal
    grid on
    set(gca,'XTick',(-180:15:180),'XTickMode','manual');
    set(gca,'YTick',(-90:10:90),'YTickMode','manual');
    xlim([-180,180]), ylim([-90,90])

    image_file  = 'TitanTexture.jpg';
    cdata       = flip(imread(image_file));
    imagesc([-180,180],[-90, 90],cdata)
end

%% (plot a different track on every for cicle, don't lose any data)
% 
% k=1;
% if orbit.kep0.i < pi/2
%     for i=1:size(longitude,1)-1
%         if longitude(i,1)>0 && longitude(i+1,1)<0
%             if longitude(i,1)>30
%                 if longitude(i,1)<170
%                     warning('check longitude, value is %.3f, next value is%.3f',longitude(i),longitude(i+1))
%                     plot(longitude(i),latitude(i),'xy','Linewidth',3,'Markersize',12)
%                     plot(longitude(i+1),latitude(i+1),'dy','Linewidth',3,'Markersize',12)
%                 end
%                 plot(longitude(k:i,1),latitude(k:i,1),'g-');
%                 k=i+1;
%             else
%                 warning('check longitude, NO INTERRUPTION IN PLOT, value is %.3f, next value is %.3f',longitude(i),longitude(i+1))
%             end
%         elseif longitude(i,1)<-160 && longitude(i+1,1)>0
%             warning('check longitude, values goes from negative to positive, value is %.3f',longitude(i))
%             plot(longitude(k:i,1),latitude(k:i,1),'g-');
%             k=i+1;
%         end
%     end
% else
%     for i=1:size(longitude,1)-1 
%         if longitude(i,1)<0 && longitude(i+1,1)>0
%             if longitude(i,1)<-30
%                 if longitude(i)>-170
%                     warning('check longitude, value is %.3f, next value is%.3f',longitude(i),longitude(i+1))
%                     plot(longitude(i),latitude(i),'xy','Linewidth',3,'Markersize',12)
%                     plot(longitude(i+1),latitude(i+1),'dy','Linewidth',3,'Markersize',12)
%                 end
%                 plot(longitude(k:i,1),latitude(k:i,1),'g-');
%                 k=i+1;
%             else
%                 warning('check longitude, NO INTERRUPTION IN PLOT, value is %.3f, next value is %.3f',longitude(i),longitude(i+1))
%             end
%         elseif longitude(i,1)>160 && longitude(i+1,1)<0
%             warning('check longitude, values goes from positive to negative, value is %.3f',longitude(i))
%             plot(longitude(k:i,1),latitude(k:i,1),'g-');
%             k=i+1;
%         end
%     end
% end
% 
% p4=plot(longitude(k:end,1),latitude(k:end,1),'g-');
% 
%% (significantly faster, but miss a value every time longitude change sign on boundaries)

if orbit.kep0.i < pi/2
    for i=1:length(longitude)-1
        if longitude(i,1)>0 && longitude(i+1,1)<0
            if longitude(i,1)>30
                if longitude(i,1)<170
                    warning('check longitude, value is %.3f, next value is%.3f',longitude(i),longitude(i+1))
                    plot(longitude(i),latitude(i),'xy','Linewidth',3,'Markersize',12)
                    plot(longitude(i+1),latitude(i+1),'dy','Linewidth',3,'Markersize',12)
                end
                longitude(i)=NaN;
            else
                warning('check longitude, NO INTERRUPTION IN PLOT, value is %.3f, next value is %.3f',longitude(i),longitude(i+1))
            end            
        elseif longitude(i,1)<-160 && longitude(i+1,1)>0
            warning('check longitude, values goes from negative to positive, value is %.3f',longitude(i))
            longitude(i)=NaN;
        end
    end
else
    for i=1:length(longitude)-1
        if longitude(i,1)<0 && longitude(i+1,1)>0
            if longitude(i,1)<-30
                if longitude(i)>-170
                    warning('check longitude, value is %.3f, next value is%.3f',longitude(i),longitude(i+1))
                    plot(longitude(i),latitude(i),'xy','Linewidth',3,'Markersize',12)
                    plot(longitude(i+1),latitude(i+1),'dy','Linewidth',3,'Markersize',12)
                end
                longitude(i)=NaN;
            else
                warning('check longitude, NO INTERRUPTION IN PLOT, value is %.3f, next value is %.3f',longitude(i),longitude(i+1))
            end
        elseif longitude(i,1)>160 && longitude(i+1,1)<0
            warning('check longitude, values goes from positive to negative, value is %.3f',longitude(i))
            longitude(i)=NaN;
        end
    end
end

%%
orb_trx = 6;
orbit_rec = 1; % number of orbit for aquisition
trx_orb_only = orb_trx - orbit_rec;
tspan = orbit.tspan;
T0 = orbit.const.T0;

if name == false
    k = 1;
    r_t = 1;
    j = 1;
    for i=1:size(tspan,2)
        if mod(r_t,2) == 1
            if tspan(i) > (orbit_rec+(k-1)*(orbit_rec+trx_orb_only))*T0 + tspan(1) 
                color = '#04cbf3';
%                 if k == 1
%                     plot(longitude(j:i),latitude(j:i),'Color',color,'Linewidth',pix_swath);%,'DisplayName','Original Ground Track')
%                 else
                    p1 = plot(longitude(j:i),latitude(j:i),'Color',color,'Linewidth',pix_swath);%,'DisplayName','Original Ground Track');
                    p1.Color = [p1.Color transp1];
%                 end
                j = i+1;
                r_t = r_t+1;
            end
        elseif mod(r_t,2) == 0
            if tspan(i) > (orbit_rec+(k-1)*(orbit_rec+trx_orb_only)+trx_orb_only)*T0 + tspan(1)
                color = 'r';
                p2 = plot(longitude(j:i),latitude(j:i),'Color',color,'Linewidth',pix_swath);%,'DisplayName','Original Ground Track');
                p2.Color = [p2.Color transp2];
                j = i+1;
                k = k+1;
                r_t = r_t+1;
            end
        end
    end
    if strcmp(color,'#04cbf3')
        color = 'r';
    elseif color == 'r'
        color = '#04cbf3';
    end
    if strcmp(color,'#04cbf3')
        pfin1 =plot(longitude(j:end),latitude(j:end),'Color',color,'Linewidth',pix_swath);%,'DisplayName','Original Ground Track');
        pfin1.Color = [pfin1.Color transp1];
    elseif color == 'r'
        p_fin2 = plot(longitude(j:end),latitude(j:end),'Color',color,'Linewidth',pix_swath);
        p_fin2.Color = [p_fin2.Color transp2];
    end
else
    p1 = plot(longitude(1:end),latitude(1:end),'Color',colour,'Linewidth',pix_swath,'DisplayName',name);
end

xlabel('Longitude \lambda   [deg]')
ylabel('Latitude \phi   [deg]')

% m1=['o',colour];
% m2=['s',colour];
% start point (yellow circle)
p3=plot(longitude(1),latitude(1),'o','Color','#fed310','Linewidth',1.5,'Markersize',9);%,'DisplayName','start');

% end point (yellow square)
p4=plot(longitude(end),latitude(end),'s','Color','#fed310','Linewidth',1.5,'Markersize',10);%,'DisplayName','end');

legend([p1 p3 p4],{'Scanned surface','start','end'},'Location','northoutside','Orientation','horizontal')
% legend('start','end','Ground Track SAR mode','Ground Track Transmission Data')
lgd = legend('Location','north','Orientation','horizontal');
lgd.NumColumns = 4;
lgd.FontSize = 15;
% legend([p1 p2],{'Ground Track SAR mode','Ground Track Transmission Data'},'Location','northoutside','Orientation','horizontal');
% % legend('start','end','Ground Track SAR mode','Ground Track Transmission Data')
% lgd = legend('Location','north','Orientation','vertical');
% lgd.NumColumns = 1;

% orbit.plot.p1 = p1;
% orbit.plot.p2 = p2;
% orbit.plot.p3 = p3;
end

