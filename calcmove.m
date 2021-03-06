function [angle,translation] = calcmove(neato_origin,neato_orientation,cones)
    neato_origin = [neato_origin; 0];
    neato_orientation = [neato_orientation; 0];
    cones(3,:) = zeros(1,length(cones));
    [~,closest_indx] = min(vecnorm(cones-neato_origin));
    ordered_cones(:,1) = cones(:,closest_indx);
    cones = remove(ordered_cones(:,1),cones);
    ordered_cones(:,2) = next_cone(neato_origin,neato_orientation,cones);
    plot(ordered_cones(1, :), ordered_cones(2, :), 'ms');
    next_pos = generate_offset(ordered_cones);
    plot(next_pos(1), next_pos(2), 'yo');
    
    k = cross(neato_orientation,next_pos-neato_origin);
    angle3d = sign(k)*atan2d(norm(k),dot(neato_orientation,next_pos-neato_origin));
    angle = angle3d(3);
    translation = next_pos(1:2)-neato_origin(1:2);
    
    v_l = 0.1;
    v_r = 0.1;
    v = (v_l+v_r)/2;
    t_forward = norm(translation)/v;
    
    d = 0.24;
    if angle < 0
        omega_l = 0.05;
        omega_r = -0.05;
    else 
        omega_l = -0.05;
        omega_r = 0.05;
    end
    omega = (omega_r - omega_l)/d;
    t_rotate = deg2rad(angle)/omega;
    
    pub = rospublisher('/raw_vel');
    sub_bump = rossubscriber('/bump');
    msg = rosmessage(pub);
    
    disp("TURNING...")
    disp(t_rotate)
    pause

    msg.Data = [omega_l,omega_r];
    send(pub, msg);
    pause(t_rotate);
    msg.Data = [v_l,v_r];
    send(pub, msg);
    pause(t_forward);
    bumpMessage = receive(sub_bump);
    if any(bumpMessage.Data)
        msg.Data = [0.0, 0.0];
        send(pub, msg);
        pause(0.1);
    end
    msg.Data = [0,0];
    send(pub,msg);

    function point = generate_offset(cones)
        vector = diff(cones,1,2);
        norm_vector = vector./vecnorm(vector);
        midangle = 90;
        rot = [0 -1 0; 1 0 0; 0 0 0];
        offset_dir = rot*0.4*norm_vector;
        point = offset_dir+cones(:,2);
    end

    function coneless = remove(cone,dataset)
        [~,indx] = max(sum(dataset==cone));
        dataset(:,indx) = []; 
        coneless = dataset;
    end

    function next = next_cone(origin,orientation,cones)
        vectors = cones-origin;
         [~, I] = sort(vecnorm(vectors));
        closest_cones = cones(:,I);
        closest_vectors = vectors(:,I);
        for i = 1:size(closest_vectors,2)
            vector = closest_vectors(:,i);
            quiver(origin(1), origin(2), vector(1), vector(2));
              if dot(vector, orientation) > 0.1
                cone = closest_cones(:,i);
                break
            end
            
        end
        next=cone;
    end
end