import os
import sys

import crocoddyl
import numpy as np
import example_robot_data
import pinocchio

WITHDISPLAY = 'display' in sys.argv or 'CROCODDYL_DISPLAY' in os.environ

crocoddyl.switchToNumpyMatrix()

# Load robot
robot = example_robot_data.loadTalos()
rmodel = robot.model

# Create data structures
rdata = rmodel.createData()
state = crocoddyl.StateMultibody(rmodel)
actuation = crocoddyl.ActuationModelFloatingBase(state)

# Set integration time
DT = 5e-2
T = 60
target = np.array([0.4, 0, 1.2])

# Initialize reference state, target and reference CoM
rightFoot = 'right_sole_link'
leftFoot = 'left_sole_link'
endEffector = 'gripper_left_joint'
endEffectorId = rmodel.getFrameId(endEffector)
rightFootId = rmodel.getFrameId(rightFoot)
leftFootId = rmodel.getFrameId(leftFoot)
q0 = rmodel.referenceConfigurations["half_sitting"]
rmodel.defaultState = np.concatenate([q0, np.zeros((rmodel.nv, 1))])
pinocchio.forwardKinematics(rmodel, rdata, q0)
pinocchio.updateFramePlacements(rmodel, rdata)
rfPos0 = rdata.oMf[rightFootId].translation
lfPos0 = rdata.oMf[leftFootId].translation
refGripper = rdata.oMf[rmodel.getFrameId("gripper_left_joint")].translation
comRef = (rfPos0 + lfPos0) / 2
comRef[2] = np.asscalar(pinocchio.centerOfMass(rmodel, rdata, q0)[2])

# Initialize Gepetto viewer
if WITHDISPLAY:
    display = crocoddyl.GepettoDisplay(robot, frameNames=[rightFoot, leftFoot])
    display.robot.viewer.gui.addSphere('world/point', .05, [1., 0., 0., 1.])  # radius = .1, RGBA=1001
    display.robot.viewer.gui.applyConfiguration('world/point', target.tolist() + [0., 0., 0., 1.])  # xyz+quaternion

# Add contact to the model
contactModel = crocoddyl.ContactModelMultiple(state, actuation.nu)
framePlacementLeft = crocoddyl.FramePlacement(leftFootId, pinocchio.SE3.Identity())
supportContactModelLeft = crocoddyl.ContactModel6D(state, framePlacementLeft, actuation.nu, np.matrix([0, 0]).T)
contactModel.addContact(leftFoot + "_contact", supportContactModelLeft)
framePlacementRight = crocoddyl.FramePlacement(rightFootId, pinocchio.SE3.Identity())
supportContactModelRight = crocoddyl.ContactModel6D(state, framePlacementRight, actuation.nu, np.matrix([0, 0]).T)
contactModel.addContact(rightFoot + "_contact", supportContactModelRight)
contactData = contactModel.createData(rdata)

# Cost for self-collision
maxfloat = sys.float_info.max
xlb = np.vstack([
    -maxfloat * np.matrix(np.ones((6, 1))),  # dimension of the SE(3) manifold
    rmodel.lowerPositionLimit[7:],
    -maxfloat * np.matrix(np.ones((state.nv, 1)))
])
xub = np.vstack([
    maxfloat * np.matrix(np.ones((6, 1))),  # dimension of the SE(3) manifold
    rmodel.upperPositionLimit[7:],
    maxfloat * np.matrix(np.ones((state.nv, 1)))
])
bounds = crocoddyl.ActivationBounds(xlb, xub, 1.)
limitCost = crocoddyl.CostModelState(state, crocoddyl.ActivationModelQuadraticBarrier(bounds), rmodel.defaultState,
                                     actuation.nu)

# Cost for state and control
stateWeights = np.array([0] * 3 + [10.] * 3 + [0.01] * (state.nv - 6) + [10] * state.nv)
stateWeightsTerm = np.array([0] * 3 + [10.] * 3 + [0.01] * (state.nv - 6) + [100] * state.nv)
xRegCost = crocoddyl.CostModelState(state, crocoddyl.ActivationModelWeightedQuad(np.matrix(stateWeights**2).T),
                                    rmodel.defaultState, actuation.nu)
uRegCost = crocoddyl.CostModelControl(state, actuation.nu)
xRegTermCost = crocoddyl.CostModelState(state, crocoddyl.ActivationModelWeightedQuad(np.matrix(stateWeightsTerm**2).T),
                                        rmodel.defaultState, actuation.nu)

# Cost for target reaching
goaltrackingWeights = np.array([1] * 3 + [0.0001] * 3)
framePoseEff = pinocchio.SE3.Identity()
framePoseEff.translation = np.matrix(target).T
Pref = crocoddyl.FramePlacement(endEffectorId, framePoseEff)
goalTrackingCost = crocoddyl.CostModelFramePlacement(
    state, crocoddyl.ActivationModelWeightedQuad(np.matrix(goaltrackingWeights**2).T), Pref, actuation.nu)

# Cost for CoM reference
comTrack = crocoddyl.CostModelCoMPosition(state, comRef, actuation.nu)

# Create cost model per each action model
runningCostModel = crocoddyl.CostModelSum(state, actuation.nu)
terminalCostModel = crocoddyl.CostModelSum(state, actuation.nu)

# Then let's added the running and terminal cost functions
runningCostModel.addCost("gripperPose", goalTrackingCost, 1e2)
runningCostModel.addCost("stateReg", xRegCost, 1e-3)
runningCostModel.addCost("ctrlReg", uRegCost, 1e-4)
runningCostModel.addCost("limitCost", limitCost, 1e3)

terminalCostModel.addCost("gripperPose", goalTrackingCost, 1e2)
terminalCostModel.addCost("stateReg", xRegTermCost, 1e-3)
terminalCostModel.addCost("limitCost", limitCost, 1e3)

# Create the action model
dmodelRunning = crocoddyl.DifferentialActionModelContactFwdDynamics(state, actuation, contactModel, runningCostModel)
dmodelTerminal = crocoddyl.DifferentialActionModelContactFwdDynamics(state, actuation, contactModel, terminalCostModel)
runningModel = crocoddyl.IntegratedActionModelEuler(dmodelRunning, DT)
terminalModel = crocoddyl.IntegratedActionModelEuler(dmodelTerminal, 0)

# Problem definition
x0 = np.concatenate([q0, pinocchio.utils.zero(state.nv)])
problem = crocoddyl.ShootingProblem(x0, [runningModel] * T, terminalModel)

# Creating the DDP solver for this OC problem, defining a logger
ddp = crocoddyl.SolverFDDP(problem)
if WITHDISPLAY:
    ddp.setCallbacks([
        crocoddyl.CallbackVerbose(),
        crocoddyl.CallbackDisplay(crocoddyl.GepettoDisplay(robot, 4, 4, frameNames=[rightFoot, leftFoot]))
    ])
else:
    ddp.setCallbacks([crocoddyl.CallbackVerbose()])

# Solving it with the DDP algorithm
xs = [rmodel.defaultState] * len(ddp.models())
us = [m.quasiStatic(d, rmodel.defaultState) for m, d in list(zip(ddp.models(), ddp.datas()))[:-1]]
ddp.solve(xs, us, 500, False, 0.1)
ddp.calc()

# Visualizing the solution in gepetto-viewer
if WITHDISPLAY:
    display.displayFromSolver(ddp)

# Get final state and end effector position
xT = ddp.xs[-1]
pinocchio.forwardKinematics(rmodel, rdata, xT[:state.nq])
pinocchio.updateFramePlacements(rmodel, rdata)
com = pinocchio.centerOfMass(rmodel, rdata, xT[:state.nq])
finalPosEff = np.array(rdata.oMf[rmodel.getFrameId("gripper_left_joint")].translation.T.flat)

print('Finally reached = ', finalPosEff)
print('Distance between hand and target = ', np.linalg.norm(finalPosEff - target))
print('Distance to default state = ', np.linalg.norm(rmodel.defaultState - np.array(xT.flat)))
print('XY distance to CoM reference = ', np.linalg.norm(com[:2] - comRef[:2]))
