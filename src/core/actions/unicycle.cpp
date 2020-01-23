///////////////////////////////////////////////////////////////////////////////
// BSD 3-Clause License
//
// Copyright (C) 2018-2020, LAAS-CNRS, University of Edinburgh
// Copyright note valid unless otherwise stated in individual files.
// All rights reserved.
///////////////////////////////////////////////////////////////////////////////

#include "crocoddyl/core/utils/exception.hpp"
#include "crocoddyl/core/actions/unicycle.hpp"

namespace crocoddyl {

ActionModelUnicycle::ActionModelUnicycle() : ActionModelAbstract(boost::make_shared<StateVector>(3), 2, 5), dt_(0.1) {
  cost_weights_ << 10., 1.;
}

ActionModelUnicycle::~ActionModelUnicycle() {}

void ActionModelUnicycle::calc(const boost::shared_ptr<ActionDataAbstract>& data,
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

  ActionDataUnicycle* d = static_cast<ActionDataUnicycle*>(data.get());
  const double& c = std::cos(x[2]);
  const double& s = std::sin(x[2]);
  d->xnext << x[0] + c * u[0] * dt_, x[1] + s * u[0] * dt_, x[2] + u[1] * dt_;
  d->r.head<3>() = cost_weights_[0] * x;
  d->r.tail<2>() = cost_weights_[1] * u;
  d->cost = 0.5 * d->r.transpose() * d->r;
}

void ActionModelUnicycle::calcDiff(const boost::shared_ptr<ActionDataAbstract>& data,
                                   const Eigen::Ref<const Eigen::VectorXd>& x,
                                   const Eigen::Ref<const Eigen::VectorXd>& u, const bool& recalc) {
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
  ActionDataUnicycle* d = static_cast<ActionDataUnicycle*>(data.get());

  // Cost derivatives
  const double& w_x = cost_weights_[0] * cost_weights_[0];
  const double& w_u = cost_weights_[1] * cost_weights_[1];
  d->Lx = x.cwiseProduct(Eigen::VectorXd::Constant(state_->get_nx(), w_x));
  d->Lu = u.cwiseProduct(Eigen::VectorXd::Constant(nu_, w_u));
  d->Lxx.diagonal() << w_x, w_x, w_x;
  d->Luu.diagonal() << w_u, w_u;

  // Dynamic derivatives
  const double& c = std::cos(x[2]);
  const double& s = std::sin(x[2]);
  d->Fx << 1., 0., -s * u[0] * dt_, 0., 1., c * u[0] * dt_, 0., 0., 1.;
  d->Fu << c * dt_, 0., s * dt_, 0., 0., dt_;
}

// ActionModelHuman::ActionModelHuman() : ActionModelAbstract(boost::make_shared<StateVector>(3), 2, 5), dt_(0.1) {
//   cost_weights_ << 10., 1.;
// }

// ActionModelHuman::~ActionModelHuman() {}

// void ActionModelHuman::calc(const boost::shared_ptr<ActionDataAbstract>& data,
//                               const Eigen::Ref<const Eigen::VectorXd>& x,
//                               const Eigen::Ref<const Eigen::VectorXd>& u){
//   if (static_cast<std::size_t>(x.size()) != state_->get_nx()) {
//     throw_pretty("Invalid argument: "
//                  << "x has wrong dimension (it should be " + std::to_string(state_->get_nx()) + ")");
//   }
//   if (static_cast<std::size_t>(u.size()) != nu_) {
//     throw_pretty("Invalid argument: "
//                  << "u has wrong dimension (it should be " + std::to_string(nu_) + ")");
//   }

//   ActionDataHuman* d = static_cast<ActionDataHuman*>(data.get());
//   const double& c = std::cos(x[2]);
//   const double& s = std::sin(x[2]);
//   d->xnext << x[0] + c * u[0] * dt_, x[1] + s * u[0] * dt_, x[2] + u[1] * dt_;
//   d->cost = costFunction(x[0],x[1],x[2],u[0],u[1]);
// }

// double ActionModelHuman::costFunction(double x, double y, double theta, 
//                                         double u1, double u2){
//   double cost;
//   cost = (pow(cost_weights_[0]*x,2) + pow(cost_weights_[0]*y,2) + 
//     pow(cost_weights_[0]*theta,2) + pow(cost_weights_[1]*u1,2) + 
//     pow(cost_weights_[1]*u2,2))/2;
//   return cost;
// }

// void ActionModelHuman::calcDiff(const boost::shared_ptr<ActionDataAbstract>& data,
//                                    const Eigen::Ref<const Eigen::VectorXd>& x,
//                                    const Eigen::Ref<const Eigen::VectorXd>& u, const bool& recalc) {
//   if (static_cast<std::size_t>(x.size()) != state_->get_nx()) {
//     throw_pretty("Invalid argument: "
//                  << "x has wrong dimension (it should be " + std::to_string(state_->get_nx()) + ")");
//   }
//   if (static_cast<std::size_t>(u.size()) != nu_) {
//     throw_pretty("Invalid argument: "
//                  << "u has wrong dimension (it should be " + std::to_string(nu_) + ")");
//   }

//   if (recalc) {
//     calc(data, x, u);
//   }
//   ActionDataHuman* d = static_cast<ActionDataHuman*>(data.get());

//   // Cost derivatives
//   d->Lx << (costFunction(x[0]+1e-6,x[1],x[2],u[0],u[1]) - costFunction(x[0],x[1],x[2],u[0],u[1]))/1e-6,
//     (costFunction(x[0],x[1]+1e-6,x[2],u[0],u[1]) - costFunction(x[0],x[1],x[2],u[0],u[1]))/1e-6,
//     (costFunction(x[0],x[1],x[2]+1e-6,u[0],u[1]) - costFunction(x[0],x[1],x[2],u[0],u[1]))/1e-6;

//   std::cout << d->Lx << "| \n" << x.cwiseProduct(Eigen::VectorXd::Constant(state_->get_nx(), cost_weights_[0]*cost_weights_[0])) << std::endl;  
//   std::cout << "------" << std::endl;

//   d->Lu << (costFunction(x[0],x[1],x[2],u[0]+1e-6,u[1]) - costFunction(x[0],x[1],x[2],u[0],u[1]))/1e-6,
//     (costFunction(x[0],x[1],x[2],u[0],u[1]+1e-6) - costFunction(x[0],x[1],x[2],u[0],u[1]))/1e-6;
  
//   std::cout << d->Lu << "| \n" << u.cwiseProduct(Eigen::VectorXd::Constant(nu_, cost_weights_[1]*cost_weights_[1])) << std::endl;  
//   std::cout << "------" << std::endl;

//   d->Lxx << (costFunction(x[0]+1e-6,x[1],x[2],u[0],u[1]) - 2*costFunction(x[0],x[1],x[2],u[0],u[1]) 
//     + costFunction(x[0]-1e-6,x[1],x[2],u[0],u[1]))/pow(1e-6,2), // dl/dxdx
//     (costFunction(x[0]+1e-6,x[1]+1e-6,x[2],u[0],u[1]) - costFunction(x[0]+1e-6,x[1],x[2],u[0],u[1]) 
//     - costFunction(x[0],x[1]+1e-6,x[2],u[0],u[1]) + costFunction(x[0],x[1],x[2],u[0],u[1]))/pow(1e-6,2), // dl/dxdy
//     (costFunction(x[0]+1e-6,x[1],x[2]+1e-6,u[0],u[1]) - costFunction(x[0]+1e-6,x[1],x[2],u[0],u[1]) 
//     - costFunction(x[0],x[1],x[2]+1e-6,u[0],u[1]) + costFunction(x[0],x[1],x[2],u[0],u[1]))/pow(1e-6,2), // dl/dxdz
//     (costFunction(x[0]+1e-6,x[1]+1e-6,x[2],u[0],u[1]) - costFunction(x[0]+1e-6,x[1],x[2],u[0],u[1]) 
//     - costFunction(x[0],x[1]+1e-6,x[2],u[0],u[1]) + costFunction(x[0],x[1],x[2],u[0],u[1]))/pow(1e-6,2), // dl/dydx 
//     (costFunction(x[0],x[1]+1e-6,x[2],u[0],u[1]) - 2*costFunction(x[0],x[1],x[2],u[0],u[1]) 
//     + costFunction(x[0],x[1]-1e-6,x[2],u[0],u[1]))/pow(1e-6,2), // dl/dydy 
//     (costFunction(x[0],x[1]+1e-6,x[2]+1e-6,u[0],u[1]) - costFunction(x[0],x[1],x[2]+1e-6,u[0],u[1]) 
//     - costFunction(x[0],x[1]+1e-6,x[2],u[0],u[1]) + costFunction(x[0],x[1],x[2],u[0],u[1]))/pow(1e-6,2), // dl/dydz 
//     (costFunction(x[0]+1e-6,x[1],x[2]+1e-6,u[0],u[1]) - costFunction(x[0]+1e-6,x[1],x[2],u[0],u[1]) 
//     - costFunction(x[0],x[1],x[2]+1e-6,u[0],u[1]) + costFunction(x[0],x[1],x[2],u[0],u[1]))/pow(1e-6,2), // dl/dzdx
//     (costFunction(x[0],x[1]+1e-6,x[2]+1e-6,u[0],u[1]) - costFunction(x[0],x[1],x[2]+1e-6,u[0],u[1]) 
//     - costFunction(x[0],x[1]+1e-6,x[2],u[0],u[1]) + costFunction(x[0],x[1],x[2],u[0],u[1]))/pow(1e-6,2), // dl/dzdy
//     (costFunction(x[0],x[1],x[2]+1e-6,u[0],u[1]) - 2*costFunction(x[0],x[1],x[2],u[0],u[1]) 
//     + costFunction(x[0],x[1],x[2]-1e-6,u[0],u[1]))/pow(1e-6,2); // dl/dzdz

//   std::cout << d->Lxx << std::endl;
//   std::cout << "------" << std::endl;    

//   d->Luu << (costFunction(x[0],x[1],x[2],u[0]+1e-6,u[1]) - 2*costFunction(x[0],x[1],x[2],u[0],u[1]) 
//     + costFunction(x[0],x[1],x[2],u[0]-1e-6,u[1]))/pow(1e-6,2), // dl/du1du1
//     (costFunction(x[0],x[1],x[2],u[0]+1e-6,u[1]+1e-6) - costFunction(x[0],x[1],x[2],u[0]+1e-6,u[1]) 
//     - costFunction(x[0],x[1],x[2],u[0],u[1]+1e-6) + costFunction(x[0],x[1],x[2],u[0],u[1]))/pow(1e-6,2), // dl/du1du2
//     (costFunction(x[0],x[1],x[2],u[0]+1e-6,u[1]+1e-6) - costFunction(x[0],x[1],x[2],u[0]+1e-6,u[1]) 
//     - costFunction(x[0],x[1],x[2],u[0],u[1]+1e-6) + costFunction(x[0],x[1],x[2],u[0],u[1]))/pow(1e-6,2), // dl/du2du1 
//     (costFunction(x[0],x[1],x[2],u[0],u[1]+1e-6) - 2*costFunction(x[0],x[1],x[2],u[0],u[1]) 
//     + costFunction(x[0],x[1],x[2],u[0],u[1]-1e-6))/pow(1e-6,2); // dl/du2du2 
//   std::cout << d->Luu << std::endl;
//   std::cout << "------" << std::endl;  

//   // Dynamic derivatives
//   const double& c = std::cos(x[2]);
//   const double& s = std::sin(x[2]);
//   d->Fx << 1., 0., -s * u[0] * dt_, 0., 1., c * u[0] * dt_, 0., 0., 1.;
//   d->Fu << c * dt_, 0., s * dt_, 0., 0., dt_;
// }

boost::shared_ptr<ActionDataAbstract> ActionModelUnicycle::createData() {
  return boost::make_shared<ActionDataUnicycle>(this);
}

const Eigen::Vector2d& ActionModelUnicycle::get_cost_weights() const { return cost_weights_; }

void ActionModelUnicycle::set_cost_weights(const Eigen::Vector2d& weights) { cost_weights_ = weights; }

}  // namespace crocoddyl
