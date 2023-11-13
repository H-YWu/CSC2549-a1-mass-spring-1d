#include <d2V_spring_particle_particle_dq2.h>

void d2V_spring_particle_particle_dq2(Eigen::MatrixXd &H, const Eigen::VectorXd &q, double stiffness) {
    H.resize(q.size(), q.size());
    H.setZero();
    for (int i = 0; i < q.size(); i ++) H(i, i) = -stiffness;
}