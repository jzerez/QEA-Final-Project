function [angle,translation] = calcmove(neato_origin,neato_orientation,cones)
    cones(3,:) = zeros(1,length(cones));
    [~,closest_indx] = min(vecnorm(cones-neato_origin));
    ordered_cones(:,1) = cones(:,closest_indx);
    cones = remove(ordered_cones(:,1),cones);
    ordered_cones(:,2) = next_cone(neato_origin,neato_orientation,cones);
    next_pos = generate_offset(ordered_cones);
    
    k = cross(neato_orientation,next_pos-neato_origin);
    angle3d = sign(k)*atan2d(norm(k),dot(neato_orientation,next_pos-neato_origin));
    angle = angle3d(3);
    translation = next_pos(1:2)-neato_origin(1:2);
    
    v_l = 0.15;
    v_r = 0.15;
    v = (v_l+v_r)/2;
    t_forward = norm(translation)/v;
    
    d = 0.24;
    if angle < 0
        omega_l = 0.15;
        omega_r = -0.15;
    else 
        omega_l = -0.15;
        omega_r = 0.15;
    end
    omega = (omega_r-omega_l)/d;
    t_rotate = angle/omega;
    
    pub = rospublisher('/raw_vel');
    sub_bump = rossubscriber('/bump');
    msg = rosmessage(pub);

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
    offset_dir = rot*-0.25*norm_vector;
    point = offset_dir+cones(:,2);
    end

    function coneless = remove(cone,dataset)
        [~,indx] = max(sum(dataset==cone));
        dataset(:,indx) = []; 
        coneless = dataset;
    end

    function next = next_cone(origin,orientation,cones)
        vectors = cones-origin;
        distances = vecnorm(vectors);
        closest_indices = distances < 2;
        closest_cones = cones(:,closest_indices);
        closest_vectors = vectors(:,closest_indices);
        dot_length = 0;
        for i = 1:size(closest_vectors,2)
            if dot(closest_vectors(:,i),orientation) > dot_length
                dot_length = dot(closest_vectors(:,i),orientation);
                cone = closest_cones(:,i);
            end
        end
        next=cone;
    end
end