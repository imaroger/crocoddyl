///////////////////////////////////////////////////////////////////////////////
// BSD 3-Clause License
//
// Copyright (C) 2018-2019, LAAS-CNRS
// Copyright note valid unless otherwise stated in individual files.
// All rights reserved.
///////////////////////////////////////////////////////////////////////////////

#ifndef BINDINGS_PYTHON_CROCODDYL_CORE_ACTIONS_Human_Terminal_HPP_
#define BINDINGS_PYTHON_CROCODDYL_CORE_ACTIONS_Human_Terminal_HPP_

#include "crocoddyl/core/actions/human_terminal_model.hpp"

namespace crocoddyl {
namespace python {

namespace bp = boost::python;

void exposeActionHumanTerminal() {
  bp::class_<ActionTerminalModelHuman, bp::bases<ActionModelAbstract> >(
      "ActionTerminalModelHuman",
      "Human action model.\n\n"
      "The transition model of an Human system is described as\n"
      "    xnext = [v*cos(theta); v*sin(theta); w],\n"
      "where the position is defined by (x, y, theta) and the control input\n"
      "by (v,w). Note that the state is defined only with the position. On the\n"
      "other hand, we define the quadratic cost functions for the state and\n"
      "control.",
      bp::init<>(bp::args("self"), "Initialize the Human action model."))
      .def("calc", &ActionTerminalModelHuman::calc_wrap,
           ActionModel_calc_wraps(bp::args("self", "data", "x", "u"),
                                  "Compute the next state and cost value.\n\n"
                                  "It describes the time-discrete evolution of the Human system.\n"
                                  "Additionally it computes the cost value associated to this discrete\n"
                                  "state and control pair.\n"
                                  ":param data: action data\n"
                                  ":param x: time-discrete state vector\n"
                                  ":param u: time-discrete control input"))
      .def<void (ActionTerminalModelHuman::*)(const boost::shared_ptr<ActionDataAbstract>&, const Eigen::VectorXd&,
                                         const Eigen::VectorXd&, const bool&)>(
          "calcDiff", &ActionTerminalModelHuman::calcDiff_wrap, bp::args("self", "data", "x", "u", "recalc"),
          "Compute the derivatives of the Human dynamics and cost functions.\n\n"
          "It computes the partial derivatives of the Human system and the\n"
          "cost function. If recalc == True, it first updates the state evolution\n"
          "and cost value. This function builds a quadratic approximation of the\n"
          "action model (i.e. dynamical system and cost function).\n"
          ":param data: action data\n"
          ":param x: time-discrete state vector\n"
          ":param u: time-discrete control input\n"
          ":param recalc: If true, it updates the state evolution and the cost value (default True).")
      .def<void (ActionTerminalModelHuman::*)(const boost::shared_ptr<ActionDataAbstract>&, const Eigen::VectorXd&,
                                         const Eigen::VectorXd&)>("calcDiff", &ActionTerminalModelHuman::calcDiff_wrap,
                                                                  bp::args("self", "data", "x", "u"))
      .def<void (ActionTerminalModelHuman::*)(const boost::shared_ptr<ActionDataAbstract>&, const Eigen::VectorXd&)>(
          "calcDiff", &ActionTerminalModelHuman::calcDiff_wrap, bp::args("self", "data", "x"))
      .def<void (ActionTerminalModelHuman::*)(const boost::shared_ptr<ActionDataAbstract>&, const Eigen::VectorXd&,
                                         const bool&)>("calcDiff", &ActionTerminalModelHuman::calcDiff_wrap,
                                                       bp::args("self", "data", "x", "recalc"))
      .def("createData", &ActionTerminalModelHuman::createData, bp::args("self"), "Create the Human action data.")
      .add_property(
          "costWeights",
          bp::make_function(&ActionTerminalModelHuman::get_cost_weights, bp::return_value_policy<bp::return_by_value>()),
          bp::make_function(&ActionTerminalModelHuman::set_cost_weights), "cost weights")
      .add_property(
          "finalState",
          bp::make_function(&ActionTerminalModelHuman::get_final_state, bp::return_value_policy<bp::return_by_value>()),
          bp::make_function(&ActionTerminalModelHuman::set_final_state), "final state")
      .add_property(
          "alpha",
          bp::make_function(&ActionTerminalModelHuman::get_alpha, bp::return_value_policy<bp::return_by_value>()),
          bp::make_function(&ActionTerminalModelHuman::set_alpha), "slowing parameter alpha");
  bp::register_ptr_to_python<boost::shared_ptr<ActionDataHumanTerminal> >();

  bp::class_<ActionDataHumanTerminal, bp::bases<ActionDataAbstract> >(
      "ActionDataHumanTerminal",
      "Action data for the Human system.\n\n"
      "The Human data, apart of common one, contains the cost residuals used\n"
      "for the computation of calc and calcDiff.",
      bp::init<ActionTerminalModelHuman*>(bp::args("self", "model"),
                                     "Create Human data.\n\n"
                                     ":param model: Human action model"));
}

}  // namespace python
}  // namespace crocoddyl

#endif  // BINDINGS_PYTHON_CROCODDYL_CORE_ACTIONS_Human_HPP_
