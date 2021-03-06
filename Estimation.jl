using LinearAlgebra

estimator_params = (angle_random_walk = 0.06, #in deg/sqrt(hour)
                    gyro_bias_instability = 0.8, #Bias instability in deg/hour
                    velocity_random_walk = .014, #in m/sec/sqrt(hour)
                    accel_bias_instability = 6 #in microG
)

function trajectory_uncertainty(thist,xhist,uhist,params,vehicle_params)
    N = length(thist)
    Δt = thist[2]-thist[1]
    # Thruster force Jacobian
    θ_t = 5.0*pi/180 # Thruster angle
    Bf = [cos(θ_t) 0 -sin(θ_t);
            cos(θ_t) 0 -sin(θ_t);
            cos(θ_t) 0 sin(θ_t);
            cos(θ_t) 0 sin(θ_t)]'
    m = vehicle_params[:m_dry] + vehicle_params[:m_fuel]
    Qbb = ((params[:accel_bias_instability]*(9.80665/1e6))^2)/(3600^3)
    Qgg = ((params[:gyro_bias_instability]*(pi/180))^2)/(3600^3)
    Qaa = ((params[:velocity_random_walk])^2)/3600
    Qωω = ((params[:angle_random_walk]*(pi/180))^2)/3600
    Phist = zeros(15,15,N)
    P0 = Diagonal([1000; 1000; 1000; .1*pi/180; .1*pi/180; .1*pi/180; .01; .01; .01; 1e-6; 1e-6; 1e-6; 1e-8; 1e-8; 1e-8])
    Phist[:,:,1] = P0
    x = zeros(16)
    for k = 1:(N-1)
        x = [xhist[1:10,k]; x[11:13]+sqrt(Qbb)*randn(3); x[14:16]+sqrt(Qgg)*randn(3)]
        u = [(1/m)*Bf*uhist[:,k]+sqrt(Qbb)*randn(3); xhist[11:13,k]+sqrt(Qωω)*randn(3)]
        xn, Phist[:,:,k+1] = single_step_prediction(x,u,Phist[:,:,k],Δt,params)
    end

    return Phist
end

function propagate_state(x0,P0,Δt,tfinal,params)
    N = Int(ceil(tfinal/Δt)+1)
    xhist = zeros(16,N)
    Phist = zeros(15,15,N)
    xhist[:,1] = x0
    Phist[:,:,1] = P0
    for k = 1:(N-1)
        u = 0.01*randn(6)
        xhist[:,k+1], Phist[:,:,k+1] = single_step_prediction(xhist[:,k],u,Phist[:,:,k],Δt,params)
    end

    return xhist, Phist
end

function single_step_prediction(x,u,P,Δt,params)
    #This function does a single-time-step update to the state and covariance
    #using dead-reckoning (propagating with the IMU with a zero-order hold).

    #Things in state vector
    r = x[1:3] #ECI position
    q = x[4:7]  #body -> ECI quaternion
    v = x[8:10] #ECI velocity
    b = x[11:13] #accelerometer bias
    g = x[14:16] #gyro bias

    #Things in input vector
    a = u[1:3] #accelerometer measurement (body frame)
    ω = u[4:6] #gyro measurement (body frame)

    R = qtoR(q) #rotation matrix (we'll use it a bunch)

    #Zero-order-hold update on state
    rn = r + Δt*v
    qn = qmult(q,expq(Δt*(ω-g)))
    vn = v + Δt*R*(a-b)
    xn = [rn;qn;vn;b;g]

    #Noise stuff
    Qaa = ((params[:velocity_random_walk])^2)/3600
    Qωω = ((params[:angle_random_walk]*(pi/180))^2)/3600
    Quu = Diagonal([Qaa*ones(3); Qωω*ones(3)])
    Qbb = ((params[:accel_bias_instability]*(9.80665/1e6))^2)/(3600^3)
    Qgg = ((params[:gyro_bias_instability]*(pi/180))^2)/(3600^3)
    Qxx = Diagonal([zeros(9); Qbb*ones(3); Qgg*ones(3)])

    #Linearized covariance update
    A = [I zeros(3,3) Δt*I zeros(3,6);
         zeros(3,3) Expso3(Δt*(ω-g))' zeros(3,6) -Δt*dExpso3(Δt*(ω-g));
         zeros(3,3) Δt*R*hat(b-a) I -Δt*R zeros(3,3);
         zeros(6,9) I]
    B = [zeros(3,6);
         zeros(3,3) Δt*dExpso3(Δt*(ω-g));
         Δt*R zeros(3,3);
         zeros(6,6)]
    Pn = A*P*A' + B*Quu*B' + Δt*Qxx

    return xn, Pn
end

function measurement_update(x,y,P,params)
    #This would nominally contain things like star tracker and GPS updates
end

function debug_dynamics(x,u,Δt)
    #Things in state vector
    r = x[1:3] #ECI position
    q = x[4:7]  #body -> ECI quaternion
    v = x[8:10] #ECI velocity
    b = x[11:13] #accelerometer bias
    g = x[14:16] #gyro bias

    #Things in input vector
    a = u[1:3] #accelerometer measurement (body frame)
    ω = u[4:6] #gyro measurement (body frame)

    R = qtoR(q) #rotation matrix (we'll use it a bunch)

    #Zero-order-hold update on state
    rn = r + Δt*v
    qn = qmult(q,expq(Δt*(ω-g)))
    vn = v + Δt*R*(a-b)
    xn = [rn;qn;vn;b;g]

    return xn
end

function debug_jacobians(x,u,Δt)
    #Things in state vector
    r = x[1:3] #ECI position
    q = x[4:7]  #body -> ECI quaternion
    v = x[8:10] #ECI velocity
    b = x[11:13] #accelerometer bias
    g = x[14:16] #gyro bias

    #Things in input vector
    a = u[1:3] #accelerometer measurement (body frame)
    ω = u[4:6] #gyro measurement (body frame)

    R = qtoR(q) #rotation matrix (we'll use it a bunch)

    #Zero-order-hold update on state
    rn = r + Δt*v
    qn = qmult(q,expq(Δt*(ω-g)))
    vn = v + Δt*R*(a-b)
    xn = [rn;qn;vn;b;g]

    #Linearized covariance update
    A = [I zeros(3,3) Δt*I zeros(3,6);
         zeros(3,3) Expso3(Δt*(ω-g))' zeros(3,6) -Δt*dExpso3(Δt*(ω-g));
         zeros(3,3) Δt*R*hat(b-a) I -Δt*R zeros(3,3);
         zeros(6,9) I]
    B = [zeros(3,6);
         zeros(3,3) Δt*dExpso3(Δt*(ω-g));
         Δt*R zeros(3,3);
         zeros(6,6)]

    return A, B
end

function random_estimator_state()
    r = randn(3)
    q = [1.0; (10*pi/360)*randn(3)]
    q = q/norm(q)
    v = randn(3)
    b = randn(3)
    g = randn(3)
    x = [r;q;v;b;g]

    a = randn(3)
    ω = randn(3)
    u = [a;ω]

    return x, u
end

function finite_diff_check(x,u,Δt)

    A1, B1 = debug_jacobians(x,u,Δt)

    Δx = Diagonal(5e-6*ones(15))
    Δu = Diagonal(5e-6*ones(6))
    A = zeros(15,15)
    B = zeros(15,6)

    for k = 1:15
        xp = debug_dynamics(state_add(x, Δx[:,k]), u, Δt)
        xm = debug_dynamics(state_add(x, -Δx[:,k]), u, Δt)
        A[:,k] = state_diff(xp,xm)/(2*Δx[k,k])
    end

    for k = 1:6
        xp = debug_dynamics(x, u + Δu[:,k], Δt)
        xm = debug_dynamics(x, u - Δu[:,k], Δt)
        B[:,k] = state_diff(xp,xm)/(2*Δu[k,k])
    end

    return A,B,A1,B1
end

function state_diff(x1,x2)
    q1 = x1[4:7]
    q2 = x2[4:7]
    dq = qmult(qconj(q2),q1)
    phi = logq(dq)
    return [x1[1:3]-x2[1:3]; phi; x1[8:16]-x2[8:16]]
end

function state_add(x,δ)
    q = x[4:7]
    qp = qmult(q,expq(δ[4:6]))
    return [x[1:3]+δ[1:3]; qp; x[8:16]+δ[7:15]]
end
