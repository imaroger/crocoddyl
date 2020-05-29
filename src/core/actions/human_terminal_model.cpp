#include "crocoddyl/core/utils/exception.hpp"
#include "crocoddyl/core/actions/human_terminal_model.hpp"
#include <iostream>
#include <math.h>

namespace crocoddyl {

ActionTerminalModelHuman::ActionTerminalModelHuman() : ActionModelAbstract(boost::make_shared<StateVector>(6), 3, 5), dt_(0.1) {
  cost_weights_ << 1. , 1.;
  final_state_ << 0., 0., 0.;
  alpha_ = 1.0;
}

ActionTerminalModelHuman::~ActionTerminalModelHuman() {}

void ActionTerminalModelHuman::calc(const boost::shared_ptr<ActionDataAbstract>& data,
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

  ActionDataHumanTerminal* d = static_cast<ActionDataHumanTerminal*>(data.get());
  const double& c = std::cos(x[2]);
  const double& s = std::sin(x[2]);
  d->xnext << x[0] + alpha_*(c*x[3]-s*x[4])*dt_, 
    x[1] + alpha_*(c*x[5]+s*x[4])*dt_,
    x[2] + x[5]*dt_,
    x[3] + u[0]*dt_,
    x[4] + u[1]*dt_,
    x[5] + u[2]*dt_;
  d->cost = costFunction(x[0],x[1],x[2]); 

  // std::cout << "cost : " << d->cost << std::endl;
}

double ActionTerminalModelHuman::costFunction(double x, double y, double theta){
  double cost;
  double diff_theta = final_state_[2] - theta;
  double diff_theta_final = diff_theta + 
    (( diff_theta > M_PI) ? - 2*M_PI : ((diff_theta< -M_PI) ? 2*M_PI: 0));

  cost = cost_weights_[0]*(pow(final_state_[0]-x,2) + pow(final_state_[1]-y,2)) + 
    cost_weights_[1]*pow(diff_theta_final,2);
  return cost;
}

void ActionTerminalModelHuman::calcDiff(const boost::shared_ptr<ActionDataAbstract>& data,
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
  ActionDataHumanTerminal* d = static_cast<ActionDataHumanTerminal*>(data.get());

  // Cost derivatives
  double h = 1e-5;

  // std::cout << "--- state ---" << std::endl; 
  // std::cout << "x: " << x[0] << " ,y: " << x[1] << " ,theta: " << x[2] << " ,vf: "
  //   << x[3] << " ,w: " << x[4] << " ,vo: " << x[5] << std::endl;
  // std::cout << final_state_ << std::endl; 
  // std::cout << "--- control ---" << std::endl; 
  // std::cout << "u1: " << u[0] << "u2: " << u[1] << "u3: " << u[2] << std::endl;

  d->Lx << (costFunction(x[0]+h,x[1],x[2]) - costFunction(x[0],x[1],x[2]))/h,
    (costFunction(x[0],x[1]+h,x[2]) - costFunction(x[0],x[1],x[2]))/h,
    (costFunction(x[0],x[1],x[2]+h) - costFunction(x[0],x[1],x[2]))/h,
    0, 0, 0;

  // std::cout << "--- Lx---" << std::endl; 
  // std::cout << d->Lx << std::endl;

  d->Lu << 0.,0.,0.;

  // std::cout << "--- Lu---" << std::endl; 
  // std::cout << d->Lu << std::endl;   

  d->Lxx.block(0,0,3,3) << (costFunction(x[0]+h,x[1],x[2]) - 2*costFunction(x[0],x[1],x[2]) 
    + costFunction(x[0]-h,x[1],x[2]))/pow(h,2), // dl/dx²
    (costFunction(x[0]+h,x[1]+h,x[2]) - costFunction(x[0]+h,x[1],x[2]) 
    - costFunction(x[0],x[1]+h,x[2]) + costFunction(x[0],x[1],x[2]))/pow(h,2), // dl/dxdy
    (costFunction(x[0]+h,x[1],x[2]+h) - costFunction(x[0]+h,x[1],x[2]) 
    - costFunction(x[0],x[1],x[2]+h) + costFunction(x[0],x[1],x[2]))/pow(h,2), // dl/dxdtheta
    (costFunction(x[0]+h,x[1]+h,x[2]) - costFunction(x[0]+h,x[1],x[2]) 
    - costFunction(x[0],x[1]+h,x[2]) + costFunction(x[0],x[1],x[2]))/pow(h,2), // dl/dydx 
    (costFunction(x[0],x[1]+h,x[2]) - 2*costFunction(x[0],x[1],x[2]) 
    + costFunction(x[0],x[1]-h,x[2]))/pow(h,2), // dl/dy² 
    (costFunction(x[0],x[1]+h,x[2]+h) - costFunction(x[0],x[1],x[2]+h) 
    - costFunction(x[0],x[1]+h,x[2]) + costFunction(x[0],x[1],x[2]))/pow(h,2), // dl/dydtheta 
    (costFunction(x[0]+h,x[1],x[2]+h) - costFunction(x[0]+h,x[1],x[2]) 
    - costFunction(x[0],x[1],x[2]+h) + costFunction(x[0],x[1],x[2]))/pow(h,2), // dl/dthetadx
    (costFunction(x[0],x[1]+h,x[2]+h) - costFunction(x[0],x[1],x[2]+h) 
    - costFunction(x[0],x[1]+h,x[2]) + costFunction(x[0],x[1],x[2]))/pow(h,2), // dl/dthetady
    (costFunction(x[0],x[1],x[2]+h) - 2*costFunction(x[0],x[1],x[2]) 
    + costFunction(x[0],x[1],x[2]-h))/pow(h,2); // dl/dtheta²;

  // std::cout << "--- Lxx---" << std::endl; 
  // std::cout << d->Lxx << std::endl;     

  d->Luu.diagonal() << 0.,0.,0.;

  // std::cout << "--- Luu---" << std::endl; 
  // std::cout << d->Luu << std::endl;     

  // Dynamic derivatives
  const double& c = std::cos(x[2]);
  const double& s = std::sin(x[2]);

  d->Fx << 1., 0., alpha_*(-s*x[3]-c*x[4]) * dt_, alpha_*c*dt_, alpha_*-s*dt_, 0., 
    0., 1., alpha_*(c*x[3]-s*x[4]) * dt_, alpha_*s*dt_, alpha_*c*dt_, 0., 
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

boost::shared_ptr<ActionDataAbstract> ActionTerminalModelHuman::createData() {
  return boost::make_shared<ActionDataHumanTerminal>(this);
}

const Eigen::VectorXd& ActionTerminalModelHuman::get_cost_weights() const { return cost_weights_; }

void ActionTerminalModelHuman::set_cost_weights(const Eigen::VectorXd& weights) { cost_weights_ = weights; }

const Eigen::Vector3d& ActionTerminalModelHuman::get_final_state() const { return final_state_; } 

void ActionTerminalModelHuman::set_final_state(const Eigen::Vector3d& statef){final_state_ = statef; }

const double ActionTerminalModelHuman::get_alpha() const { return alpha_; }

void ActionTerminalModelHuman::set_alpha(const double slowing_param){alpha_ = slowing_param; }
}  // namespace crocoddyl
