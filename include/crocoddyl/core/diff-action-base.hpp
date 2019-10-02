///////////////////////////////////////////////////////////////////////////////
// BSD 3-Clause License
//
// Copyright (C) 2018-2019, LAAS-CNRS
// Copyright note valid unless otherwise stated in individual files.
// All rights reserved.
///////////////////////////////////////////////////////////////////////////////

#ifndef CROCODDYL_CORE_DIFF_ACTION_BASE_HPP_
#define CROCODDYL_CORE_DIFF_ACTION_BASE_HPP_

#include "crocoddyl/core/state-base.hpp"
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

namespace crocoddyl {

struct DifferentialActionDataAbstract;  // forward declaration

/**
 * @brief This class DifferentialActionModelAbstract represents a first-order
 * ODE, i.e.
 * \f[
 * \mathbf{\dot{v}} = \mathbf{f}(\mathbf{q}, \mathbf{v}, \boldsymbol{\tau})
 * \f]
 * where \f$ xout = \mathbf{\dot{v}} \f$ and represents the  acceleration of the
 * system. Note that Jacobians Fx and Fu in the
 * DifferentialActionDataAbstract are in \f$ \mathbb{R}^{nv\times ndx} \f$ and
 * \f$ \mathbb{R}^{nv\times nu} \f$, respectively.
 *
 * Then we use the acceleration to integrate the system, and as consequence we
 * obtain:
 * \f[
 * \mathbf{\dot{x}} = (\mathbf{v}, \mathbf{\dot{v}}) = \mathbf{f}(\mathbf{x},\mathbf{u})
 * \f]
 * where this \f$ f \f$ function is different to the other one.
 * So \f$ xout \f$ is interpreted here as \f$ vdout \f$ or \f$ aout \f$.
 */
class DifferentialActionModelAbstract {
 public:
  DifferentialActionModelAbstract(StateAbstract& state, unsigned int const& nu, unsigned int const& nr = 0);
  virtual ~DifferentialActionModelAbstract();

  virtual void calc(const boost::shared_ptr<DifferentialActionDataAbstract>& data,
                    const Eigen::Ref<const Eigen::VectorXd>& x, const Eigen::Ref<const Eigen::VectorXd>& u) = 0;
  virtual void calcDiff(const boost::shared_ptr<DifferentialActionDataAbstract>& data,
                        const Eigen::Ref<const Eigen::VectorXd>& x, const Eigen::Ref<const Eigen::VectorXd>& u,
                        const bool& recalc = true) = 0;
  virtual boost::shared_ptr<DifferentialActionDataAbstract> createData();

  void calc(const boost::shared_ptr<DifferentialActionDataAbstract>& data, const Eigen::Ref<const Eigen::VectorXd>& x);
  void calcDiff(const boost::shared_ptr<DifferentialActionDataAbstract>& data,
                const Eigen::Ref<const Eigen::VectorXd>& x);

  unsigned int const& get_nu() const;
  unsigned int const& get_nr() const;
  StateAbstract& get_state() const;

 protected:
  unsigned int nu_;
  unsigned int nr_;
  StateAbstract& state_;
  Eigen::VectorXd unone_;

#ifdef PYTHON_BINDINGS

 public:
  void calc_wrap(const boost::shared_ptr<DifferentialActionDataAbstract>& data, const Eigen::VectorXd& x,
                 const Eigen::VectorXd& u = Eigen::VectorXd()) {
    if (u.size() == 0) {
      calc(data, x);
    } else {
      calc(data, x, u);
    }
  }

  void calcDiff_wrap(const boost::shared_ptr<DifferentialActionDataAbstract>& data, const Eigen::VectorXd& x,
                     const Eigen::VectorXd& u, const bool& recalc) {
    calcDiff(data, x, u, recalc);
  }
  void calcDiff_wrap(const boost::shared_ptr<DifferentialActionDataAbstract>& data, const Eigen::VectorXd& x,
                     const Eigen::VectorXd& u) {
    calcDiff(data, x, u, true);
  }
  void calcDiff_wrap(const boost::shared_ptr<DifferentialActionDataAbstract>& data, const Eigen::VectorXd& x) {
    calcDiff(data, x, unone_, true);
  }
  void calcDiff_wrap(const boost::shared_ptr<DifferentialActionDataAbstract>& data, const Eigen::VectorXd& x,
                     const bool& recalc) {
    calcDiff(data, x, unone_, recalc);
  }

#endif
};

struct DifferentialActionDataAbstract {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  template <typename Model>
  explicit DifferentialActionDataAbstract(Model* const model)
      : cost(0.),
        xout(model->get_state().get_nv()),
        Fx(model->get_state().get_nv(), model->get_state().get_ndx()),
        Fu(model->get_state().get_nv(), model->get_nu()),
        r(model->get_nr()),
        Lx(model->get_state().get_ndx()),
        Lu(model->get_nu()),
        Lxx(model->get_state().get_ndx(), model->get_state().get_ndx()),
        Lxu(model->get_state().get_ndx(), model->get_nu()),
        Luu(model->get_nu(), model->get_nu()) {
    xout.setZero();
    r.setZero();
    Fx.setZero();
    Fu.setZero();
    Lx.setZero();
    Lu.setZero();
    Lxx.setZero();
    Lxu.setZero();
    Luu.setZero();
  }

  double cost;
  Eigen::VectorXd xout;
  Eigen::MatrixXd Fx;
  Eigen::MatrixXd Fu;
  Eigen::VectorXd r;
  Eigen::VectorXd Lx;
  Eigen::VectorXd Lu;
  Eigen::MatrixXd Lxx;
  Eigen::MatrixXd Lxu;
  Eigen::MatrixXd Luu;
};

}  // namespace crocoddyl

#endif  // CROCODDYL_CORE_DIFF_ACTION_BASE_HPP_