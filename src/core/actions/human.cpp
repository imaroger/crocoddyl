#include "crocoddyl/core/utils/exception.hpp"
#include "crocoddyl/core/actions/human.hpp"
#include <iostream>
#include <math.h>

namespace crocoddyl {


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
//   // d->r.head<3>() = cost_weights_[0] * x;
//   // d->r.tail<2>() = cost_weights_[1] * u;
//   //d->cost = 0.5 * d->r.transpose() * d->r;
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

ActionModelHuman::ActionModelHuman() : ActionModelAbstract(boost::make_shared<StateVector>(6), 3, 5), dt_(0.01) {
  cost_weights_ << 1., 1.2, 1.7, 0.7, 5.2, 100;
  final_state_ << 0.1, 0.1, 0.1;
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

  ActionDataHuman* d = static_cast<ActionDataHuman*>(data.get());
  const double& c = std::cos(x[2]);
  const double& s = std::sin(x[2]);
  d->xnext << x[0] + (c*x[3]-s*x[5])*dt_, 
    x[1] + (c*x[5]+s*x[3])*dt_,
    x[2] + x[4]*dt_,
    x[3] + u[0]*dt_,
    x[4] + u[1]*dt_,
    x[5] + u[2]*dt_;
  // d->r << sqrt(cost_weights_[0]),
  //   sqrt(cost_weights_[1])*u[0],
  //   sqrt(cost_weights_[2])*u[1],
  //   sqrt(cost_weights_[3])*u[2],
  //   sqrt(cost_weights_[4])*atan((final_state_[1]-x[1])/(final_state_[0]-x[0]))-x[2];
  d->cost = costFunction(x[0],x[1],x[2],u[0],u[1],u[2]); 

  std::cout << "cost : " << d->cost << std::endl;
}

double ActionModelHuman::costFunction(double x, double y, double theta, 
                                        double u1, double u2, double u3){
  double cost;
  cost = 
    cost_weights_[0] + cost_weights_[1]*pow(u1,2) + cost_weights_[2]*pow(u2,2) +
    cost_weights_[3]*pow(u3,2) + cost_weights_[4]*
    pow(atan2(final_state_[1]-y,final_state_[0]-x)-theta,2) + 
    cost_weights_[5]*(pow(final_state_[0]-x,2) + pow(final_state_[1]-y,2) + pow(final_state_[2]-theta,2));
  return cost;
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

  //const double& psi = atan((final_state_[1]-x[1])/(final_state_[0]-x[0]))-x[2];
  // Eigen::MatrixXd Lx2;
  // Eigen::MatrixXd Lu2;  
  // Lx2 << (final_state_[1]-x[1])*cost_weights_[4]/(pow(final_state_[1]-x[1],2)+
  //   (pow(final_state_[0]-x[0],2)))*2*psi,
  //   (final_state_[0]-x[0])*cost_weights_[4]/(pow(final_state_[1]-x[1],2)+
  //   (pow(final_state_[0]-x[0],2)))*2*psi,
  //   -2*cost_weights_[4]*psi, 0, 0, 0;
  //Lu2 << 2*cost_weights_[1]*u[1],2*cost_weights_[2]*u[2],2*cost_weights_[3]*u[3];
  
  double h = 1e-6;

  std::cout << "--- state ---" << std::endl; 
  std::cout << "x: " << x[0] << " ,y: " << x[1] << " ,theta: " << x[2] << " ,vf: "
    << x[3] << " ,w: " << x[4] << " ,vo: " << x[5] << std::endl;
  std::cout << final_state_ << std::endl; 
  std::cout << "--- control ---" << std::endl; 
  std::cout << "u1: " << u[0] << "u2: " << u[1] << "u3: " << u[2] << std::endl;

  d->Lx << (costFunction(x[0]+h,x[1],x[2],u[0],u[1],u[2]) - costFunction(x[0],x[1],x[2],u[0],u[1],u[2]))/h,
    (costFunction(x[0],x[1]+h,x[2],u[0],u[1],u[2]) - costFunction(x[0],x[1],x[2],u[0],u[1],u[2]))/h,
    (costFunction(x[0],x[1],x[2]+h,u[0],u[1],u[2]) - costFunction(x[0],x[1],x[2],u[0],u[1],u[2]))/h,
    0, 0, 0;

  std::cout << "--- Lx---" << std::endl; 
  std::cout << d->Lx << std::endl;

  // std::cout << "--- Lx derivate ---" << std::endl; 
  // std::cout << (final_state_[1]-x[1])*cost_weights_[4]/(pow(final_state_[1]-x[1],2)+
  //   (pow(final_state_[0]-x[0],2)))*2*psi << std::endl; 

  d->Lu << (costFunction(x[0],x[1],x[2],u[0]+h,u[1],u[2]) - costFunction(x[0],x[1],x[2],u[0],u[1],u[2]))/h,
    (costFunction(x[0],x[1],x[2],u[0],u[1]+h,u[2]) - costFunction(x[0],x[1],x[2],u[0],u[1],u[2]))/h,
    (costFunction(x[0],x[1],x[2],u[0],u[1],u[2]+h) - costFunction(x[0],x[1],x[2],u[0],u[1],u[2]))/h;
  
  // d->Lu << 2*cost_weights_[1], 2*cost_weights_[2], 2*cost_weights_[3];

  std::cout << "--- Lu---" << std::endl; 
  std::cout << d->Lu << std::endl;   

  // std::cout << "--- Lu derivate ---" << std::endl; 
  // std::cout << Lu2 << std::endl; 

  d->Lxx.block(0,0,3,3) << (costFunction(x[0]+h,x[1],x[2],u[0],u[1],u[2]) - 2*costFunction(x[0],x[1],x[2],u[0],u[1],u[2]) 
    + costFunction(x[0]-h,x[1],x[2],u[0],u[1],u[2]))/pow(h,2), // dl/dx²
    (costFunction(x[0]+h,x[1]+h,x[2],u[0],u[1],u[2]) - costFunction(x[0]+h,x[1],x[2],u[0],u[1],u[2]) 
    - costFunction(x[0],x[1]+h,x[2],u[0],u[1],u[2]) + costFunction(x[0],x[1],x[2],u[0],u[1],u[2]))/pow(h,2), // dl/dxdy
    (costFunction(x[0]+h,x[1],x[2]+h,u[0],u[1],u[2]) - costFunction(x[0]+h,x[1],x[2],u[0],u[1],u[2]) 
    - costFunction(x[0],x[1],x[2]+h,u[0],u[1],u[2]) + costFunction(x[0],x[1],x[2],u[0],u[1],u[2]))/pow(h,2), // dl/dxdtheta
    (costFunction(x[0]+h,x[1]+h,x[2],u[0],u[1],u[2]) - costFunction(x[0]+h,x[1],x[2],u[0],u[1],u[2]) 
    - costFunction(x[0],x[1]+h,x[2],u[0],u[1],u[2]) + costFunction(x[0],x[1],x[2],u[0],u[1],u[2]))/pow(h,2), // dl/dydx 
    (costFunction(x[0],x[1]+h,x[2],u[0],u[1],u[2]) - 2*costFunction(x[0],x[1],x[2],u[0],u[1],u[2]) 
    + costFunction(x[0],x[1]-h,x[2],u[0],u[1],u[2]))/pow(h,2), // dl/dy² 
    (costFunction(x[0],x[1]+h,x[2]+h,u[0],u[1],u[2]) - costFunction(x[0],x[1],x[2]+h,u[0],u[1],u[2]) 
    - costFunction(x[0],x[1]+h,x[2],u[0],u[1],u[2]) + costFunction(x[0],x[1],x[2],u[0],u[1],u[2]))/pow(h,2), // dl/dydtheta 
    (costFunction(x[0]+h,x[1],x[2]+h,u[0],u[1],u[2]) - costFunction(x[0]+h,x[1],x[2],u[0],u[1],u[2]) 
    - costFunction(x[0],x[1],x[2]+h,u[0],u[1],u[2]) + costFunction(x[0],x[1],x[2],u[0],u[1],u[2]))/pow(h,2), // dl/dthetadx
    (costFunction(x[0],x[1]+h,x[2]+h,u[0],u[1],u[2]) - costFunction(x[0],x[1],x[2]+h,u[0],u[1],u[2]) 
    - costFunction(x[0],x[1]+h,x[2],u[0],u[1],u[2]) + costFunction(x[0],x[1],x[2],u[0],u[1],u[2]))/pow(h,2), // dl/dthetady
    (costFunction(x[0],x[1],x[2]+h,u[0],u[1],u[2]) - 2*costFunction(x[0],x[1],x[2],u[0],u[1],u[2]) 
    + costFunction(x[0],x[1],x[2]-h,u[0],u[1],u[2]))/pow(h,2); // dl/dtheta²;

  std::cout << "--- Lxx---" << std::endl; 
  std::cout << d->Lxx << std::endl;     

  d->Luu << (costFunction(x[0],x[1],x[2],u[0]+h,u[1],u[2]) - 2*costFunction(x[0],x[1],x[2],u[0],u[1],u[2]) 
    + costFunction(x[0],x[1],x[2],u[0]-h,u[1],u[2]))/pow(h,2), // dl/du1²
    (costFunction(x[0],x[1],x[2],u[0]+h,u[1]+h,u[2]) - costFunction(x[0],x[1],x[2],u[0]+h,u[1],u[2]) 
    - costFunction(x[0],x[1],x[2],u[0],u[1]+h,u[2]) + costFunction(x[0],x[1],x[2],u[0],u[1],u[2]))/pow(h,2), // dl/du1du2
    (costFunction(x[0],x[1],x[2],u[0]+h,u[1],u[2]+h) - costFunction(x[0],x[1],x[2],u[0]+h,u[1],u[2]) 
    - costFunction(x[0],x[1],x[2],u[0],u[1],u[2]+h) + costFunction(x[0],x[1],x[2],u[0],u[1],u[2]))/pow(h,2), // dl/du1du3    
    (costFunction(x[0],x[1],x[2],u[0]+h,u[1]+h,u[2]) - costFunction(x[0],x[1],x[2],u[0]+h,u[1],u[2]) 
    - costFunction(x[0],x[1],x[2],u[0],u[1]+h,u[2]) + costFunction(x[0],x[1],x[2],u[0],u[1],u[2]))/pow(h,2), // dl/du2du1 
    (costFunction(x[0],x[1],x[2],u[0],u[1]+h,u[2]) - 2*costFunction(x[0],x[1],x[2],u[0],u[1],u[2]) 
    + costFunction(x[0],x[1],x[2],u[0],u[1]-h,u[2]))/pow(h,2), // dl/du2² 
    (costFunction(x[0],x[1],x[2],u[0],u[1]+h,u[2]+h) - costFunction(x[0],x[1],x[2],u[0],u[1],u[2]+h) 
    - costFunction(x[0],x[1],x[2],u[0],u[1]+h,u[2]) + costFunction(x[0],x[1],x[2],u[0],u[1],u[2]))/pow(h,2), // dl/du2du3     
    (costFunction(x[0],x[1],x[2],u[0]+h,u[1],u[2]+h) - costFunction(x[0],x[1],x[2],u[0]+h,u[1],u[2]) 
    - costFunction(x[0],x[1],x[2],u[0],u[1],u[2]+h) + costFunction(x[0],x[1],x[2],u[0],u[1],u[2]))/pow(h,2), // dl/du3du1 
    (costFunction(x[0],x[1],x[2],u[0],u[1]+h,u[2]+h) - costFunction(x[0],x[1],x[2],u[0],u[1],u[2]+h) 
    - costFunction(x[0],x[1],x[2],u[0],u[1]+h,u[2]) + costFunction(x[0],x[1],x[2],u[0],u[1],u[2]))/pow(h,2), // dl/du3du2 
    (costFunction(x[0],x[1],x[2],u[0],u[1],u[2]+h) - 2*costFunction(x[0],x[1],x[2],u[0],u[1],u[2]) 
    + costFunction(x[0],x[1],x[2],u[0],u[1],u[2]-h))/pow(h,2), // dl/du3² 

  // d->Luu.diagonal() << 2*cost_weights_[1], 2*cost_weights_[2], 2*cost_weights_[3];

  std::cout << "--- Luu---" << std::endl; 
  std::cout << d->Luu << std::endl;     

  // Dynamic derivatives
  const double& c = std::cos(x[2]);
  const double& s = std::sin(x[2]);

  d->Fx << 1., 0., (-s*x[3]-c*x[5]) * dt_, c*dt_, 0., -s*dt_,
    0., 1., (c*x[3]-s*x[5]) * dt_, s*dt_, 0., c*dt_,
    0., 0., 1., 0., dt_, 0.,
    0., 0., 0., 1., 0., 0.,
    0., 0., 0., 0., 1., 0.,
    0., 0., 0., 0., 0., 1.;

  std::cout << "--- Fx---" << std::endl; 
  std::cout << d->Fx << std::endl;  

  d->Fu.block(3,0,3,3).diagonal() << dt_, dt_, dt_;

  std::cout << "--- Fu---" << std::endl; 
  std::cout << d->Fu << std::endl;  

}

boost::shared_ptr<ActionDataAbstract> ActionModelHuman::createData() {
  return boost::make_shared<ActionDataHuman>(this);
}

const Eigen::VectorXd& ActionModelHuman::get_cost_weights() const { return cost_weights_; }

void ActionModelHuman::set_cost_weights(const Eigen::VectorXd& weights) { cost_weights_ = weights; }

const Eigen::Vector3d& ActionModelHuman::get_final_state() const { return final_state_; } 

void ActionModelHuman::set_final_state(const Eigen::Vector3d& statef){final_state_ = statef; }

}  // namespace crocoddyl
