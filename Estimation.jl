using LinearAlgebra

estimator_params = (angle_random_walk = 0.06, #in deg/sqrt(hour)
                    gyro_bias_instability = 0.8, #Bias instability in deg/hour
                    velocity_random_walk = .014, #in m/sec/sqrt(hour)
                    accel_bias_instability = 6 #in microG
)

function imu_model(x,params)

end

function state_prediction(x,u,P,Δt,params)
    #Note that this function assumes Δt is pretty small relative
    #to the vehicle dynamics (e.g. ≈100Hz IMU update rate)

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
    rn = r + Δt*v + 0.5*Δt*R*(a-b)
    qn = qmult(q,qexp(0.5*Δt*(ω-g)))
    vn = v + Δt*R*(a-b)
    xn = [rn;qn;vn;b;g]

    #Noise stuff
    Qaa = ((params[:velocity_random_walk])^2)/3600
    Qωω = ((params[:angle_random_walk]*(pi/180))^2)/3600
    Quu = Diagonal([Qaa*ones(3); Qωω*ones(3)])
    Qbb = ((params[:accel_bias_instability])^2)/(3600^3)
    Qgg = ((params[:gyro_bias_instability]*(pi/180))^2)/(3600^3)
    Qxx = Diagonal([zeros(9); Qbb*ones(3); Qgg*ones(3)])

    #Linearized covariance update
    A = [I zeros(3,3) Δt*I zeros(3,6);
         zeros(3,3) I-Δt*hat(ω-g) zeros(3,6) Δt*I+((Δt^2)/2)*hat(ω);
         zeros(3,3) Δt*R*hat(b-a) I -Δt*R zeros(3,3);
         zeros(6,9) I]
    B = [zeros(3,6);
         zeros(3,3) Δt*I;
         Δt*R zeros(3,3);
         zeros(6,6)]
    Pn = A*P*A' + B*Quu*B' + Δt*Qxx

    return xn, Pn
end

function measurement_update(x,y,P,params)

end
