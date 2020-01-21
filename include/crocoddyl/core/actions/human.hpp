///////////////////////////////////////////////////////////////////////////////
// BSD 3-Clause License
//
// Copyright (C) 2018-2020, LAAS-CNRS, University of Edinburgh
// Copyright note valid unless otherwise stated in individual files.
// All rights reserved.
///////////////////////////////////////////////////////////////////////////////

#ifndef CROCODDYL_CORE_ACTIONS_Human_HPP_
#define CROCODDYL_CORE_ACTIONS_Human_HPP_

#include "crocoddyl/core/action-base.hpp"
#include "crocoddyl/core/states/euclidean.hpp"
#include <stdexcept>

namespace crocoddyl {

class ActionModelHuman : public ActionModelAbstract {
 public: 
  ActionModelHuman();
  ~ActionModelHuman();

  void calc(const boost::shared_ptr<ActionDataAbstract>& data, const Eigen::Ref<const Eigen::VectorXd>& x,
            const Eigen::Ref<const Eigen::VectorXd>& u);
  void calcDiff(const boost::shared_ptr<ActionDataAbstract>& data, const Eigen::Ref<const Eigen::VectorXd>& x,
                const Eigen::Ref<const Eigen::VectorXd>& u, const bool& recalc = true);
  boost::shared_ptr<ActionDataAbstract> createData();

  double costFunction(double x, double y, double theta, 
                        double u1, double u2, double u3);

  const Eigen::VectorXd& get_cost_weights() const;
  void set_cost_weights(const Eigen::VectorXd& weights);

  const Eigen::Vector3d& get_final_state() const;
  void set_final_state(const Eigen::Vector3d& statef);

 private:
  Eigen::VectorXd cost_weights_ = Eigen::VectorXd(6);
  Eigen::Vector3d final_state_;
  double dt_;
};

struct ActionDataHuman : public ActionDataAbstract {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  template <typename Model>
  explicit ActionDataHuman(Model* const model) : ActionDataAbstract(model) {}
};

}  // namespace crocoddyl

#endif  // CROCODDYL_CORE_ACTIONS_Human_HPP_
