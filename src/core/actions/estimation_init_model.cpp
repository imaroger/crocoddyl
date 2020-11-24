#include "crocoddyl/core/utils/exception.hpp"
#include "crocoddyl/core/actions/estimation_init_model.hpp"
#include <iostream>
#include <math.h>

namespace crocoddyl {

ActionInitModelEstimation::ActionInitModelEstimation() : ActionModelAbstract(boost::make_shared<StateVector>(6), 3, 5), dt_(0.1) {
  cost_weights_ << 1., 1.2, 1.7, 0.7,1. , 1.;
  current_state_ << 0., 0., 0.;
}

ActionInitModelEstimation::~ActionInitModelEstimation() {}

void ActionInitModelEstimation::calc(const boost::shared_ptr<ActionDataAbstract>& data,
                               const Eigen::Ref<const Eigen::VectorXd>& x,
                               const Eigen::Ref<const Eigen::VectorXd>& u) {
  if (static_cast<std::size_t>(x.size()) != state_->get_nx()) {
    throw_pretty("Invalid argument: "
                 << "x has wrong dimension (it should be " + std::to_string(state_->get_nx()) + ")");
  }
  if (static_cast<std::size_t>(u.size()) != nu_) {
    throw_pretty("Invalid argument: "
                 << "u has wrong dimension (it should be " + std::to_string(nu_) + ")");
  }

  ActionDataEstimationInit* d = static_cast<ActionDataEstimationInit*>(data.get());
  const double& c = std::cos(x[2]);
  const double& s = std::sin(x[2]);
  d->xnext << x[0] + (c*x[3]-s*x[4])*dt_, 
    x[1] + (c*x[5]+s*x[4])*dt_,
    x[2] + x[5]*dt_,
    x[3] + u[0]*dt_,
    x[4] + u[1]*dt_,
    x[5] + u[2]*dt_;
  d->cost = costFunction(x[0],x[1],x[2],u[0],u[1],u[2]); 

  // std::cout << "cost : " << d->cost << std::endl;
}

double ActionInitModelEstimation::costFunction(double x, double y, double theta, 
                                        double u1, double u2, double u3){
  double cost;
  double diff_theta = current_state_[2] - theta;
  double diff_theta_final = diff_theta + 
    (( diff_theta > M_PI) ? - 2*M_PI : ((diff_theta< -M_PI) ? 2*M_PI: 0));

  cost = cost_weights_[0] + cost_weights_[1]*pow(u1,2) + cost_weights_[2]*pow(u2,2) +
    cost_weights_[3]*pow(u3,2) + cost_weights_[4]*(pow(current_state_[0]-x,2) + pow(current_state_[1]-y,2)) + 
    cost_weights_[5]*pow(diff_theta_final,2);
  return cost;
}

void ActionInitModelEstimation::calcDiff(const boost::shared_ptr<ActionDataAbstract>& data,
                                   const Eigen::Ref<const Eigen::VectorXd>& x,
                                   const Eigen::Ref<const Eigen::VectorXd>& u, 
                                   const bool& recalc) {
  if (static_cast<std::size_t>(x.size()) != state_->get_nx()) {
    throw_pretty("Invalid argument: "
                 << "x has wrong dimension (it should be " + std::to_string(state_->get_nx()) + ")");
  }
  if (static_cast<std::size_t>(u.size()) != nu_) {
    throw_pretty("Invalid argument: "
                 << "u has wrong dimension (it should be " + std::to_string(nu_) + ")");
  }

  if (recalc) {
    calc(data, x, u);
  }
  ActionDataEstimationInit* d = static_cast<ActionDataEstimationInit*>(data.get());

  // Cost derivatives
  double h = 1e-5;

  // std::cout << "--- state ---" << std::endl; 
  // std::cout << "x: " << x[0] << " ,y: " << x[1] << " ,theta: " << x[2] << " ,vf: "
  //   << x[3] << " ,w: " << x[4] << " ,vo: " << x[5] << std::endl;
  // std::cout << current_state_ << std::endl; 
  // std::cout << "--- control ---" << std::endl; 
  // std::cout << "u1: " << u[0] << "u2: " << u[1] << "u3: " << u[2] << std::endl;

  d->Lx << -2*cost_weights_[4]*(current_state_[0]-x[0]),-2*cost_weights_[4]*(current_state_[1]-x[1]),
    (costFunction(x[0],x[1],x[2]+h,u[0],u[1],u[2]) - costFunction(x[0],x[1],x[2],u[0],u[1],u[2]))/h,
    0.,0.,0.;

  // std::cout << "--- Lx---" << std::endl; 
  // std::cout << d->Lx << std::endl;

  d->Lu << 2*cost_weights_[1]*u[0], 2*cost_weights_[2]*u[1], 2*cost_weights_[3]*u[2];

  // std::cout << "--- Lu---" << std::endl; 
  // std::cout << d->Lu << std::endl;   

  d->Lxx.diagonal() << 2*cost_weights_[4],2*cost_weights_[4],
    (costFunction(x[0],x[1],x[2]+h,u[0],u[1],u[2]) - 2*costFunction(x[0],x[1],x[2],u[0],u[1],u[2]) 
    + costFunction(x[0],x[1],x[2]-h,u[0],u[1],u[2]))/pow(h,2), // dl/dtheta²;
    0.,0.,0.;

  // std::cout << "--- Lxx---" << std::endl; 
  // std::cout << d->Lxx << std::endl;     

  d->Luu.diagonal() << 2*cost_weights_[1], 2*cost_weights_[2], 2*cost_weights_[3];

  // std::cout << "--- Luu---" << std::endl; 
  // std::cout << d->Luu << std::endl;     

  // Dynamic derivatives
  const double& c = std::cos(x[2]);
  const double& s = std::sin(x[2]);

  d->Fx << 1., 0., (-s*x[3]-c*x[4]) * dt_, c*dt_, -s*dt_, 0., 
    0., 1., (c*x[3]-s*x[4]) * dt_, s*dt_, c*dt_, 0., 
    0., 0., 1., 0.,0.,  dt_,
    0., 0., 0., 1., 0., 0.,
    0., 0., 0., 0., 1., 0.,
    0., 0., 0., 0., 0., 1.;

  // std::cout << "--- Fx---" << std::endl; 
  // std::cout << d->Fx << std::endl;  

  d->Fu.block(3,0,3,3).diagonal() << dt_, dt_, dt_;

  // std::cout << "--- Fu---" << std::endl; 
  // std::cout << d->Fu << std::endl;  

}

boost::shared_ptr<ActionDataAbstract> ActionInitModelEstimation::createData() {
  return boost::make_shared<ActionDataEstimationInit>(this);
}

const Eigen::VectorXd& ActionInitModelEstimation::get_cost_weights() const { return cost_weights_; }

void ActionInitModelEstimation::set_cost_weights(const Eigen::VectorXd& weights) { cost_weights_ = weights; }

const Eigen::Vector3d& ActionInitModelEstimation::get_current_state() const { return current_state_; } 

void ActionInitModelEstimation::set_current_state(const Eigen::Vector3d& state){current_state_ = state; }
}// namespace crocoddyl
