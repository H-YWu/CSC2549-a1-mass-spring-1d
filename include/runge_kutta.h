//Input:
//  q - generalized coordiantes for the mass-spring system
//  qdot - generalized velocity for the mass spring system
//  dt - the time step in seconds
//  mass - the mass
//  force(q, qdot) - a function that computes the force acting on the mass as a function. This takes q and qdot as parameters.
//Output:
//  q - set q to the updated generalized coordinate using Runge-Kutta time integration
//  qdot - set qdot to the updated generalized velocity using Runge-Kutta time integration

template<typename FORCE> 
inline void runge_kutta(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double mass,  FORCE &force) {
    Eigen::VectorXd f;
    force(f, q, qdot);
    Eigen::VectorXd k1 = dt * (-f);
    Eigen::VectorXd h1 = dt * qdot;
    Eigen::VectorXd q_tmp = q + 0.5 * h1;
    force(f, q_tmp, qdot);
    Eigen::VectorXd k2 = dt * (-f);
    Eigen::VectorXd h2 = dt * (qdot + 0.5 * k1);
    q_tmp = q + 0.5 * h2;
    force(f, q_tmp, qdot);
    Eigen::VectorXd k3 = dt * (-f);
    Eigen::VectorXd h3 = dt * (qdot + 0.5 * k2);
    q_tmp = q + h3;
    force(f, q_tmp, qdot);
    Eigen::VectorXd k4 = dt * (-f);
    Eigen::VectorXd h4 = dt * (qdot + k3);

    qdot = qdot + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
    q = q + (h1 + 2.0 * h2 + 2.0 * h3 + h4) / 6.0;
}