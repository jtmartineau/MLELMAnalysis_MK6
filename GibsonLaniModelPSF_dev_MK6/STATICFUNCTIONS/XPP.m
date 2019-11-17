function ret=XPP(rho, theta, xoffset, yoffset, thetaoffset,...
    liveactuators, actuatorheights)

    liveactuators=reshape(liveactuators,...
        sqrt(numel(liveactuators)),...
        sqrt(numel(liveactuators)));

    actuatorheights=reshape(actuatorheights,...
        sqrt(numel(actuatorheights)),...
        sqrt(numel(actuatorheights)));

    x_pupil=rho.*cos(theta);
    y_pupil=rho.*sin(theta);

    liveactuators=rot90(liveactuators,2);
    actuatorheights=rot90(actuatorheights,2);

    rotation_matrix=[cos(thetaoffset),-sin(thetaoffset);
                     sin(thetaoffset),cos(thetaoffset)];
    %%%%% rotate pupil coordinates %%%%%
    %                                  % 
    pupil_coords=[reshape(x_pupil,[1,numel(x_pupil)]);...
                  reshape(y_pupil,[1,numel(y_pupil)])];
    pupil_coords_rot=zeros(size(pupil_coords));
    for ii=1:size(pupil_coords,2)
        pupil_coords_rot(:,ii)=rotation_matrix*pupil_coords(:,ii);
    end
    x_pupil_rot=reshape(pupil_coords_rot(1,:),[size(x_pupil,1),size(x_pupil,2)]);
    y_pupil_rot=reshape(pupil_coords_rot(2,:),[size(y_pupil,1),size(y_pupil,2)]);
    %                                  % 
    %%%%% rotate pupil coordinates %%%%%

    K_X_actuators=size(liveactuators,2);
    X_interactuator_distance=2/K_X_actuators;
    X_actuator_positions=linspace(-1+X_interactuator_distance/2,1-X_interactuator_distance/2,K_X_actuators);
    K_Y_actuators=size(liveactuators,1);
    Y_interactuator_distance=2/K_Y_actuators;
    Y_actuator_positions=linspace(-1+Y_interactuator_distance/2,1-Y_interactuator_distance/2,K_Y_actuators);
    [X_points,Y_points]=meshgrid(X_actuator_positions,Y_actuator_positions);

    alpha_x=100;
    alpha_y=100;
    phase_plate=zeros(size(y_pupil,1),size(x_pupil,2));
    for ii=1:size(X_points,1)
        for jj=1:size(Y_points,2)
            if liveactuators(ii,jj)==1           
                if (ii==1) && (jj==size(Y_points,1)) 
                    phase_plate=phase_plate+actuatorheights(ii,jj)*...
                        custom_sigmf(-(x_pupil_rot-xoffset),[alpha_x, X_points(ii,jj)-X_interactuator_distance/2]).*...
                        custom_sigmf(y_pupil_rot-yoffset,[alpha_y, -Y_points(ii,jj)-Y_interactuator_distance/2]);
                elseif (ii==size(X_points,2)) && (jj==1)
                    phase_plate=phase_plate+actuatorheights(ii,jj)*...
                        custom_sigmf((x_pupil_rot-xoffset),[alpha_x, -X_points(ii,jj)-X_interactuator_distance/2]).*...
                        custom_sigmf(-(y_pupil_rot-yoffset),[alpha_y, Y_points(ii,jj)-Y_interactuator_distance/2]);
                else
                    phase_plate=phase_plate+actuatorheights(ii,jj)*...
                        custom_sigmf(x_pupil_rot-xoffset,[alpha_x, -X_points(ii,jj)-X_interactuator_distance/2]).*...
                        custom_sigmf(-(x_pupil_rot-xoffset),[alpha_x, X_points(ii,jj)-X_interactuator_distance/2]).*...
                        custom_sigmf(y_pupil_rot-yoffset,[alpha_y, -Y_points(ii,jj)-Y_interactuator_distance/2]).*...
                        custom_sigmf(-(y_pupil_rot-yoffset),[alpha_y, Y_points(ii,jj)-Y_interactuator_distance/2]);     
                end
            end
        end
    end
    ret=phase_plate;
end

% this is a custom sigmf function 
function ret=custom_sigmf(x,params_vec)
    a=params_vec(1);
    c=params_vec(2);
    ret=1./(1+exp(-a*(x-c)));
end