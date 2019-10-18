using Pkg
Pkg.add("PyPlot")
using PyPlot

struct Gyroscope
  angle_random_walk::Float64
  bias_stability::Float64
end

struct StarTracker
  sigma::Float64
  bias::Float64
  sampling_period::Float64
end

function bayard_attitude_by_time(t, q1, q2, l, r, p11_init, p12_init, p22_init, b)
  term1 = q2 * t.^3 / 3.0
  term2 = p22_init * t.^2
  term3 = (2 * p12_init + q1) * t
  term4 = p11_init + b^2
  return (term1 + term2 + term3) .+ term4
end

function bayard_attitude(times, q1, q2, st_sampling_period, st_sigma, st_bias)
  r = st_sampling_period * st_sigma^2
  l = sqrt(q1 + 2 * sqrt(r * q2))
  p11_init = sqrt(r) * l
  p12_init = sqrt(r * q2)
  p22_init = sqrt(q2) * l
  return bayard_attitude_by_time(times, q1, q2, l, r, p11_init, p12_init, p22_init, st_bias)
end

macro Name(arg)
  string(arg)
end

#function run()
att_meas = StarTracker(333e-6, (333e-6), 0.5)
gyros = Dict(
 "g364"        => Gyroscope(((0.09 * π/180.0)^2) / 3600.0, ((2.2 * π/180.0) / 3600.0)^2 / 3600.0), # pdca
 "honeywell"   => Gyroscope(5.8761e-12, 1.8138e-18),
 "g370"        => Gyroscope(((0.06 * π/180.0)^2) / 3600.0, ((0.8 * π/180.0) / 3600.0)^2 / 3600.0), #pdc0
 "adis16495-2" => Gyroscope(((0.1 * π/180.0)^2) / 3600.0, ((1.6 * π/180.0) / 3600.0)^2 / 3600.0),
 "sdi50x-c"    => Gyroscope(((0.02 * π/180.0)^2) / 3600.0, ((2.0 * π/180.0) / 3600.0)^2 / 3600.0))

times = 0:1.0:3600 * 100
clf()

for gyro_id ∈ keys(gyros)
  gyro = get(gyros, gyro_id, 0)
  p = bayard_attitude(times, gyro.angle_random_walk, gyro.bias_stability, att_meas.sampling_period, att_meas.sigma, att_meas.bias)
  sigma = p.^0.5 * 180/π
  pygui(true)
  plot(times / 3600.0, sigma, label=gyro_id)
end
ax = gca()
ax.set_yscale("log")
ax.set_xscale("log")
ax.set_xlabel("time (h)")
ax.set_ylabel("att err (deg) per axis 1-sigma")
ax.set_ylim((0.01, 10.0))
ax.grid(true)
ax.legend()
display(gcf())
