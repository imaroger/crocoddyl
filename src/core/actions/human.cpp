#include "crocoddyl/core/utils/exception.hpp"
#include "crocoddyl/core/actions/human.hpp"
#include <iostream>
#include <math.h>

namespace crocoddyl {

ActionModelHuman::ActionModelHuman() : ActionModelAbstract(boost::make_shared<StateVector>(6), 3, 5), dt_(0.1) {
  cost_weights_ << 1., 1.2, 1.7, 0.7, 5.2;
}

ActionModelHuman::~ActionModelHuman() {}

void ActionModelHuman::calc(const boost::shared_ptr<ActionDataAbstract>& data,
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

  std::cout << state_->get_nx() << " | " << nu_<< std::endl;

  ActionDataHuman* d = static_cast<ActionDataHuman*>(data.get());
  const double& c = std::cos(x[2]);
  const double& s = std::sin(x[2]);
  d->xnext << x[0] + (c*x[3]-s*u[5])*dt_, 
    x[1] + (c*x[5]+s*x[3])*dt_,
    x[2] + x[4]*dt_,
    x[3] + u[0]*dt_,
    x[4] + u[1]*dt_,
    x[5] + u[2]*dt_;
  d->r << cost_weights_[0],
  cost_weights_[1]*pow(u[0],2),
  cost_weights_[2]*pow(u[1],2),
  cost_weights_[3]*pow(u[2],2),
  cost_weights_[4]*pow(atan((final_state_[1]-x[1])/(final_state_[0]-x[0]))-x[2],2);
  d->cost = 0.5 * d->r.transpose() * d->r; //Here J is sum not vector
}

void ActionModelHuman::calcDiff(const boost::shared_ptr<ActionDataAbstract>& data,
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
  ActionDataHuman* d = static_cast<ActionDataHuman*>(data.get());

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

boost::shared_ptr<ActionDataAbstract> ActionModelHuman::createData() {
  return boost::make_shared<ActionDataHuman>(this);
}

const Eigen::VectorXd& ActionModelHuman::get_cost_weights() const { return cost_weights_; }

void ActionModelHuman::set_cost_weights(const Eigen::VectorXd& weights) { cost_weights_ = weights; }

const Eigen::Vector2d& ActionModelHuman::get_cost_weights() const { return final_state_; } 

void ActionModelHuman::set_final_state(const Eigen::Vector2d& statef){final_state_ = statef; }

}  // namespace crocoddyl
