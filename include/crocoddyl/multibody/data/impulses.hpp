///////////////////////////////////////////////////////////////////////////////
// BSD 3-Clause License
//
// Copyright (C) 2018-2020, University of Edinburgh
// Copyright note valid unless otherwise stated in individual files.
// All rights reserved.
///////////////////////////////////////////////////////////////////////////////

#ifndef CROCODDYL_CORE_DATA_IMPULSES_HPP_
#define CROCODDYL_CORE_DATA_IMPULSES_HPP_

#include <boost/shared_ptr.hpp>
#include "crocoddyl/multibody/data/multibody.hpp"
#include "crocoddyl/multibody/impulses/multiple-impulses.hpp"

namespace crocoddyl {

struct DataCollectorImpulse : virtual DataCollectorAbstract {
  DataCollectorImpulse(boost::shared_ptr<ImpulseDataMultiple> impulses)
      : DataCollectorAbstract(), impulses(impulses) {}
  virtual ~DataCollectorImpulse() {}

  boost::shared_ptr<ImpulseDataMultiple> impulses;
};

struct DataCollectorMultibodyInImpulse : DataCollectorMultibody, DataCollectorImpulse {
  DataCollectorMultibodyInImpulse(pinocchio::Data* const pinocchio, boost::shared_ptr<ImpulseDataMultiple> impulses)
      : DataCollectorMultibody(pinocchio), DataCollectorImpulse(impulses) {}
  virtual ~DataCollectorMultibodyInImpulse() {}
};

}  // namespace crocoddyl

#endif  // CROCODDYL_CORE_DATA_MULTIBODY_IN_IMPULSE_HPP_