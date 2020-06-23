///////////////////////////////////////////////////////////////////////////////
// BSD 3-Clause License
//
// Copyright (C) 2018-2019, LAAS-CNRS
// Copyright note valid unless otherwise stated in individual files.
// All rights reserved.
///////////////////////////////////////////////////////////////////////////////

#include "crocoddyl/core/utils/exception.hpp"
#include "crocoddyl/core/solver-base.hpp"

namespace crocoddyl {

SolverAbstract::SolverAbstract(boost::shared_ptr<ShootingProblem> problem)
    : problem_(problem),
      is_feasible_(false),
      cost_(0.),
      stop_(0.),
      xreg_(NAN),
      ureg_(NAN),
      steplength_(1.),
      dV_(0.),
      dVexp_(0.),
      th_acceptstep_(0.1),
      th_stop_(1e-4),
      iter_(0) {
  // Allocate common data
  const std::size_t& T = problem_->get_T();
  xs_.resize(T + 1);
  us_.resize(T);
  models_.resize(T + 1);
  datas_.resize(T + 1);
  for (std::size_t t = 0; t < T; ++t) {
    const boost::shared_ptr<ActionModelAbstract>& model = problem_->get_runningModels()[t];
    const boost::shared_ptr<ActionDataAbstract>& data = problem_->get_runningDatas()[t];
    const std::size_t& nu = model->get_nu();

    xs_[t] = model->get_state()->zero();
    us_[t] = Eigen::VectorXd::Zero(nu);
    models_[t] = model;
    datas_[t] = data;
  }
  xs_.back() = problem_->get_terminalModel()->get_state()->zero();
  models_.back() = problem_->get_terminalModel();
  datas_.back() = problem_->get_terminalData();
}

SolverAbstract::~SolverAbstract() {}

void SolverAbstract::setCandidate(const std::vector<Eigen::VectorXd>& xs_warm,
                                  const std::vector<Eigen::VectorXd>& us_warm, const bool& is_feasible) {
  const std::size_t& T = problem_->get_T();

  if (xs_warm.size() == 0) {
    for (std::size_t t = 0; t < T; ++t) {
      xs_[t] = problem_->get_runningModels()[t]->get_state()->zero();
    }
    xs_.back() = problem_->get_terminalModel()->get_state()->zero();
  } else {
    assert_pretty(xs_warm.size() == T + 1,
                  "Warm start state has wrong dimension, got " << xs_warm.size() << " expecting " << (T + 1));
    std::copy(xs_warm.begin(), xs_warm.end(), xs_.begin());
  }

  if (us_warm.size() == 0) {
    for (std::size_t t = 0; t < T; ++t) {
      const std::size_t& nu = problem_->get_runningModels()[t]->get_nu();
      us_[t] = Eigen::VectorXd::Zero(nu);
    }
  } else {
    assert_pretty(us_warm.size() == T,
                  "Warm start control has wrong dimension, got " << us_warm.size() << " expecting " << T);
    std::copy(us_warm.begin(), us_warm.end(), us_.begin());
  }
  is_feasible_ = is_feasible;
}

void SolverAbstract::setCallbacks(const std::vector<boost::shared_ptr<CallbackAbstract> >& callbacks) {
  callbacks_ = callbacks;
}

const std::vector<boost::shared_ptr<CallbackAbstract> >& SolverAbstract::getCallbacks() const { return callbacks_; }

const boost::shared_ptr<ShootingProblem>& SolverAbstract::get_problem() const { return problem_; }

const std::vector<boost::shared_ptr<ActionModelAbstract> >& SolverAbstract::get_models() const { return models_; }

const std::vector<boost::shared_ptr<ActionDataAbstract> >& SolverAbstract::get_datas() const { return datas_; }

const std::vector<Eigen::VectorXd>& SolverAbstract::get_xs() const { return xs_; }

const std::vector<Eigen::VectorXd>& SolverAbstract::get_us() const { return us_; }

const bool& SolverAbstract::get_isFeasible() const { return is_feasible_; }

const std::size_t& SolverAbstract::get_iter() const { return iter_; }

const double& SolverAbstract::get_cost() const { return cost_; }

const double& SolverAbstract::get_stop() const { return stop_; }

const Eigen::Vector2d& SolverAbstract::get_d() const { return d_; }

const double& SolverAbstract::get_xreg() const { return xreg_; }

const double& SolverAbstract::get_ureg() const { return ureg_; }

const double& SolverAbstract::get_stepLength() const { return steplength_; }

const double& SolverAbstract::get_dV() const { return dV_; }

const double& SolverAbstract::get_dVexp() const { return dVexp_; }

bool raiseIfNaN(const double& value) {
  if (std::isnan(value) || std::isinf(value) || value >= 1e30) {
    return true;
  } else {
    return false;
  }
}

}  // namespace crocoddyl
