# Information on LQR/Linearization Implementation

The `linearize` function gives the A, B matrices of the approximated linear
model for the quadrobor, linearized around a hover at [0, 0, 1], evaluated
as a given z and u.

Run lqr_test.m to see the quadrobot go to point 1,1,1 with steady-state error

Run lqr_sim_test.m to see then quadrobor follow (but fail to catch) the uav
