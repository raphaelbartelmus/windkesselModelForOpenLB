/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  by Raphael Bartelmus
 *
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
 */

#include "olb3D.h"
#include "olb3D.hh"
#include "wkfct3d.h"

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

using T            = FLOATING_POINT_TYPE;
using DESCRIPTOR   = D3Q27<>;
using BulkDynamics = BGKdynamics<T, DESCRIPTOR>;

// PARAMETERS TO CHANGE ----------------------------------------

// --- Choose only one flow profile! ---
// #define pedrizettiFlowProfile
#define wang2023FlowProfileVOLUME
// #define wang2023FlowProfileVELOCITY

// GRID/GEOMETRY
#define resolution (int)93     // resolution: number of voxels per charPhysL
#define charPhysLength (T)0.03 // reference length of simulation geometry in m
#define noOfCuboids (int)32    // number of cuboids

// VELOCITY
#define maxFlowVelocity (T)0.05    // maximum flow velocity in m/s
#define latticeVelocity (T)0.01732 // lattice velocity

// OUTPUT/SAVING
#define startUpTime (T)4        // time for smooth start up in s
#define noHeartBeatPeriods (T)2 // defines simulation length
#define heartBeatPeriod (T)1    // 1/HR (for pedrizetti flow)
#define maxPhysT (T) startUpTime + noHeartBeatPeriods* heartBeatPeriod
#define deltaVTK (T)0.04        // time between vtk outputs in phys Simulation s
#define deltaVTKstartup (T)0.04 // time between vtk outputs during startup

// 3-Element-Windkessel-Model parameters
#define windkesselFactor                                                       \
  (T)5e-03 // prefactor for calculated pressure values for numerical stability
#define latticeDeltaT                                                          \
  (int)1 // deltaT for 3WK model in lattice time steps, usually 1
#define order                                                                  \
  (int)3 // order of backward finite difference for 1st and 2nd derivative
// -------------------------------------------------------------

// (Maybe) improves speed
static T   flowRateVector[5]; // [0]= summed FlowRates, [1]=Area, [2:]=Vector
static T   pressureVector[5]; // [0]= summed Pressures, [1]=Area, [2:]=Vector
static int input[1];          // irrelevant

// ------------------------ Global Variables -----------------------
#define vectorLength 6 // length of the velocity vector
T pressureOutletDesc[2] = {0.0, 0.0};
T velocityOutletDesc[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
T pressureOutletBCT[2]  = {0.0, 0.0};
T velocityOutletBCT[6]  = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
T pressureOutletLCCA[2] = {0.0, 0.0};
T velocityOutletLCCA[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
T pressureOutletLSA[2]  = {0.0, 0.0};
T velocityOutletLSA[6]  = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

#pragma region windkesselParameters
// Descending aorta
const T RP_AO    = 2.01e07;
const T RC_AO    = 26.6e07;
const T C_AO     = 57e-10;
const T ALPHA_AO = (RC_AO + RP_AO) / (RP_AO * C_AO);
const T BETA_AO  = RC_AO;
const T GAMMA_AO = (1 / (RP_AO * C_AO));
// BCT
const T RP_BCT    = 4.3e07;
const T RC_BCT    = 115e07;
const T C_BCT     = 15e-10;
const T ALPHA_BCT = (RC_BCT + RP_BCT) / (RP_BCT * C_BCT);
const T BETA_BCT  = RC_BCT;
const T GAMMA_BCT = (1 / (RP_BCT * C_BCT));
// LCCA
const T RP_LCCA    = 5e07;
const T RC_LCCA    = 240e07;
const T C_LCCA     = 8e-10;
const T ALPHA_LCCA = (RC_LCCA + RP_LCCA) / (RP_LCCA * C_LCCA);
const T BETA_LCCA  = RC_LCCA;
const T GAMMA_LCCA = (1 / (RP_LCCA * C_LCCA));
// LSA
const T RP_LSA    = 6e07;
const T RC_LSA    = 210e07;
const T C_LSA     = 9e-10;
const T ALPHA_LSA = (RC_LSA + RP_LSA) / (RP_LSA * C_LSA);
const T BETA_LSA  = RC_LSA;
const T GAMMA_LSA = (1 / (RP_LSA * C_LSA));
#pragma endregion

// Stores data from .stl file in geometry in form of material numbers
void prepareGeometry(UnitConverter<T, DESCRIPTOR> const& converter,
                     IndicatorF3D<T>& indicator, STLreader<T>& stlReader,
                     SuperGeometry<T, 3>& superGeometry)
{

  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename(0, 2, indicator);
  superGeometry.rename(2, 1, stlReader);

  superGeometry.clean();

  // IndicatorCircle3D<T>        (3x center, 3x normal, radius)

  // Set material number for inflow
  IndicatorCircle3D<T>   inflow(0., 0., 0., 0., 0., 1., 0.017);
  IndicatorCylinder3D<T> layerInflow(
      inflow, 2. * converter.getConversionFactorLength());
  superGeometry.rename(2, 3, 1, layerInflow);

  // Set material number for 1st coronary artery
  IndicatorCircle3D<T>   coronary1(0.0343083, 0.0259034, 0., 0., 0., 1.,
                                   0.0017647);
  IndicatorCylinder3D<T> coronary1layer(
      coronary1, 2. * converter.getConversionFactorLength());
  superGeometry.rename(2, 4, 1, coronary1layer);

  // // Set material number for 2nd coronary artery
  IndicatorCircle3D<T>   coronary2(-0.00529951, -0.04390055, 0., 0., 0., 1.,
                                   0.0017647);
  IndicatorCylinder3D<T> coronary2layer(
      coronary2, 2. * converter.getConversionFactorLength());
  superGeometry.rename(2, 5, 1, coronary2layer);

  // // Set material number for descending outflow
  IndicatorCircle3D<T>   outflowDescending(-0.004515225, 0.07217795, 0., 0., 0.,
                                           1., 0.0109504);
  IndicatorCylinder3D<T> layerOutflowDescending(
      outflowDescending, 2. * converter.getConversionFactorLength());
  superGeometry.rename(2, 6, 1, layerOutflowDescending);

  // // Set material number for BCT
  IndicatorCircle3D<T>   outflowBCT(0.0517593, 0.024923, 0.200954, 0., 0., 1.,
                                    0.009411768);
  IndicatorCylinder3D<T> layerOutflowBCT(
      outflowBCT, 2. * converter.getConversionFactorLength());
  superGeometry.rename(2, 7, 1, layerOutflowBCT);

  // // Set material number for LCCA
  IndicatorCircle3D<T>   outflowLCCA(0.0797742, 0.0311975, 0.200954, 0., 0., 1.,
                                     0.004705884);
  IndicatorCylinder3D<T> layerOutflowLCCA(
      outflowLCCA, 2. * converter.getConversionFactorLength());
  superGeometry.rename(2, 8, 1, layerOutflowLCCA);

  // // Set material number for LSA
  IndicatorCircle3D<T> outflowLSA(0.08215115, 0.04629555, 0.200954, 0., 0., 1.,
                                  0.00490195);
  IndicatorCylinder3D<T> layerOutflowLSA(
      outflowLSA, 2. * converter.getConversionFactorLength());
  superGeometry.rename(2, 9, 1, layerOutflowLSA);

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  // Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean(3);
  superGeometry.checkForErrors();

  superGeometry.print();
  clout << "Prepare Geometry ... OK" << std::endl;
}

// Set up the geometry of the simulation
void prepareLattice(SuperLattice<T, DESCRIPTOR>&        lattice,
                    UnitConverter<T, DESCRIPTOR> const& converter,
                    STLreader<T>& stlReader, SuperGeometry<T, 3>& superGeometry)
{

  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // material=1 --> bulk dynamics
  lattice.defineDynamics<BulkDynamics>(superGeometry, 1);

  // material=2 --> Bouzidi boundary (wall)
  setBouzidiBoundary<T, DESCRIPTOR>(lattice, superGeometry, 2, stlReader);

  // Inlet: material=3 --> bulk dynamics + velocity (inflow)
  lattice.defineDynamics<BulkDynamics>(superGeometry, 3);
  setInterpolatedVelocityBoundary<T, DESCRIPTOR>(lattice, omega, superGeometry,
                                                 3);

  // CorAs: material=4,5 --> bulk dynamics + velocity (outflow)
  lattice.defineDynamics<BulkDynamics>(superGeometry, 4);
  setInterpolatedVelocityBoundary<T, DESCRIPTOR>(lattice, omega, superGeometry,
                                                 4);
  lattice.defineDynamics<BulkDynamics>(superGeometry, 5);
  setInterpolatedVelocityBoundary<T, DESCRIPTOR>(lattice, omega, superGeometry,
                                                 5);

  // material=6,7,8,9 --> bulk dynamics + pressure (outflow)
  lattice.defineDynamics<BulkDynamics>(
      superGeometry.getMaterialIndicator({6, 7, 8, 9}));
  setInterpolatedPressureBoundary<T, DESCRIPTOR>(
      lattice, omega, superGeometry.getMaterialIndicator({6, 7, 8, 9}));

  // Initial conditions
  AnalyticalConst3D<T, T> rhoF(1);
  std::vector<T>          velocity(3, T());
  AnalyticalConst3D<T, T> uF(velocity);

  // Initialize all values of distribution functions to their local equilibrium
  lattice.defineRhoU(superGeometry.getMaterialIndicator({1, 3}), rhoF, uF);
  lattice.iniEquilibrium(superGeometry.getMaterialIndicator({1, 3}), rhoF, uF);

  lattice.setParameter<descriptors::OMEGA>(omega);
  lattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setBoundaryValues(SuperLattice<T, DESCRIPTOR>&        sLattice,
                       UnitConverter<T, DESCRIPTOR> const& converter, int iT,
                       SuperGeometry<T, 3>& superGeometry)
{
  //  INITIAL SMOOTH START-UP OF FLOWS
  static int iTstartUp = converter.getLatticeTime(startUpTime);
  static int iTperiod  = converter.getLatticeTime(heartBeatPeriod);
  T          startUpFactor =
      (int(iT) > int(iTstartUp)) ? 1.0 : linearScale(iT, iTstartUp);

// SET VELOCITY BOUNDARY CONDITIONS
#if defined(pedrizettiFlowProfile)
  T poiseuilleMaxU =
      pedrizettiFlow(iT, iTperiod) * startUpFactor * maxFlowVelocity * 1.5;
#elif defined(wang2023FlowProfileVOLUME)
  T poiseuilleMaxU = wang2023Flow(iT, iTperiod) * startUpFactor *
                     maxFlowVelocity * 1.5 * 4.355;
#elif defined(wang2023FlowProfileVELOCITY)
  T poiseuilleMaxU =
      wang2023Flow(iT, iTperiod) * startUpFactor * maxFlowVelocity * 1.5 / 0.35;
#else
#error "No flow profile defined!"
#endif
  CirclePoiseuille3D<T> velocityIn(
      superGeometry, 3, converter.getLatticeVelocity(poiseuilleMaxU), 0);
  sLattice.defineU(superGeometry, 3, velocityIn);

// SCALING needed to adjust the calculated flow rate to the actual flow rate
// set by maxFlowVelocity and the defined FlowProfile
#if defined(pedrizettiFlowProfile)
  T uCorA =
      corAflow(iT, iTperiod) * startUpFactor * (maxFlowVelocity / 0.2) * 1.5;
#elif defined(wang2023FlowProfileVOLUME)
  T uCorA =
      corAflow(iT, iTperiod) * startUpFactor * (maxFlowVelocity / 0.2) * 1.5;
#elif defined(wang2023FlowProfileVELOCITY)
  T uCorA = corAflow(iT, iTperiod) * startUpFactor * (maxFlowVelocity / 0.2) *
            ((0.1 / 0.35) / 0.4355) * 1.5;
#endif
  CirclePoiseuille3D<T> velocityCorA1(superGeometry, 4,
                                      converter.getLatticeVelocity(uCorA), 0);
  CirclePoiseuille3D<T> velocityCorA2(superGeometry, 5,
                                      converter.getLatticeVelocity(uCorA), 0);
  sLattice.defineU(superGeometry, 4, velocityCorA1);
  sLattice.defineU(superGeometry, 5, velocityCorA2);

  sLattice.setProcessingContext<
      Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
      ProcessingContext::Simulation);

  //  SET PRESSURE BOUNDARY CONDITIONS
  // 3WK-Model for AoD, BCT, LCCA, LSA
  AnalyticalConst3D<T, T> pressureOutDesc(
      converter.getLatticeDensityFromPhysPressure(pressureOutletDesc[0] *
                                                  windkesselFactor));
  sLattice.defineRho(superGeometry, 6, pressureOutDesc);
  AnalyticalConst3D<T, T> pressureOutBCT(
      converter.getLatticeDensityFromPhysPressure(pressureOutletBCT[0] *
                                                  windkesselFactor));
  sLattice.defineRho(superGeometry, 7, pressureOutBCT);
  AnalyticalConst3D<T, T> pressureOutLCCA(
      converter.getLatticeDensityFromPhysPressure(pressureOutletLCCA[0] *
                                                  windkesselFactor));
  sLattice.defineRho(superGeometry, 8, pressureOutLCCA);
  AnalyticalConst3D<T, T> pressureOutLSA(
      converter.getLatticeDensityFromPhysPressure(pressureOutletLSA[0] *
                                                  windkesselFactor));
  sLattice.defineRho(superGeometry, 9, pressureOutLSA);

  sLattice.setProcessingContext<Array<momenta::FixedDensity::RHO>>(
      ProcessingContext::Simulation);
}

void getResults(SuperLattice<T, DESCRIPTOR>&  sLattice,
                UnitConverter<T, DESCRIPTOR>& converter, int iT,
                SuperGeometry<T, 3>& superGeometry, util::Timer<T>& timer)
{

  OstreamManager clout(std::cout, "getResults");

  sLattice.setProcessingContext(ProcessingContext::Evaluation);
  sLattice.scheduleBackgroundOutputVTK([&, iT](auto task) {
    SuperVTMwriter3D<T>        vtmWriter("windkessel3d");
    SuperLatticePhysVelocity3D velocity(sLattice, converter);
    SuperLatticePhysPressure3D pressure(sLattice, converter);
    vtmWriter.addFunctor(velocity);
    vtmWriter.addFunctor(pressure);
    task(vtmWriter, iT);
  });

  // Timer console output
  timer.update(iT);
  timer.print(iT);
}

void windkesselUpdater(SuperPlaneIntegralF3D<T>& velocityIntegralDesc,
               SuperPlaneIntegralF3D<T>& pressureIntegralDesc,
               SuperPlaneIntegralF3D<T>& velocityIntegralBCT,
               SuperPlaneIntegralF3D<T>& pressureIntegralBCT,
               SuperPlaneIntegralF3D<T>& velocityIntegralLCCA,
               SuperPlaneIntegralF3D<T>& pressureIntegralLCCA,
               SuperPlaneIntegralF3D<T>& velocityIntegralLSA,
               SuperPlaneIntegralF3D<T>& pressureIntegralLSA, T deltaT)
{
  // AoD ----------------------------------------------------------------------
  velocityIntegralDesc(flowRateVector, input);
  pressureIntegralDesc(pressureVector, input);
  // shift the velocity vector by one position
  for (int i = vectorLength - 1; i > 0; --i) {
    velocityOutletDesc[i] = velocityOutletDesc[i - 1];
  }
  velocityOutletDesc[0] = flowRateVector[0];
  pressureOutletDesc[1] = pressureVector[0] / pressureVector[1];
  // calculate 1st and 2nd derivative
  T deltaFlow = backwardFiniteDifference(1, order, deltaT, velocityOutletDesc,
                                         vectorLength);
  T deltaDeltaFlow = backwardFiniteDifference(2, order, deltaT,
                                              velocityOutletDesc, vectorLength);
  // calculate pressure BC for next timestep with 3WK model
  pressureOutletDesc[0] =
      rungeKutta(pressureOutletDesc[1], velocityOutletDesc[0], deltaFlow,
                 deltaDeltaFlow, deltaT, ALPHA_AO, BETA_AO, GAMMA_AO);

  // BCT  -------------------------------------------------------------------
  velocityIntegralBCT(flowRateVector, input);
  pressureIntegralBCT(pressureVector, input);
  // shift the velocity vector by one position
  for (int i = vectorLength - 1; i > 0; --i) {
    velocityOutletBCT[i] = velocityOutletBCT[i - 1];
  }
  velocityOutletBCT[0] = flowRateVector[0];
  pressureOutletBCT[1] = pressureVector[0] / pressureVector[1];
  // calculate 1st and 2nd derivative
  deltaFlow      = backwardFiniteDifference(1, order, deltaT, velocityOutletBCT,
                                            vectorLength);
  deltaDeltaFlow = backwardFiniteDifference(2, order, deltaT, velocityOutletBCT,
                                            vectorLength);
  // calculate pressure BC for next timestep with 3WK model
  pressureOutletBCT[0] =
      rungeKutta(pressureOutletBCT[1], velocityOutletBCT[0], deltaFlow,
                 deltaDeltaFlow, deltaT, ALPHA_BCT, BETA_BCT, GAMMA_BCT);

  // LCCA -------------------------------------------------------------------
  velocityIntegralLCCA(flowRateVector, input);
  pressureIntegralLCCA(pressureVector, input);
  // shift the velocity vector by one position
  for (int i = vectorLength - 1; i > 0; --i) {
    velocityOutletLCCA[i] = velocityOutletLCCA[i - 1];
  }
  velocityOutletLCCA[0] = flowRateVector[0];
  pressureOutletLCCA[1] = pressureVector[0] / pressureVector[1];
  // calculate 1st and 2nd derivative
  deltaFlow = backwardFiniteDifference(1, order, deltaT, velocityOutletLCCA,
                                       vectorLength);
  deltaDeltaFlow = backwardFiniteDifference(2, order, deltaT,
                                            velocityOutletLCCA, vectorLength);
  // calculate pressure BC for next timestep with 3WK model
  pressureOutletLCCA[0] =
      rungeKutta(pressureOutletLCCA[1], velocityOutletLCCA[0], deltaFlow,
                 deltaDeltaFlow, deltaT, ALPHA_LCCA, BETA_LCCA, GAMMA_LCCA);

  // LSA -------------------------------------------------------------------
  velocityIntegralLSA(flowRateVector, input);
  pressureIntegralLSA(pressureVector, input);
  // shift the velocity vector by one position
  for (int i = vectorLength - 1; i > 0; --i) {
    velocityOutletLSA[i] = velocityOutletLSA[i - 1];
  }
  velocityOutletLSA[0] = flowRateVector[0];
  pressureOutletLSA[1] = pressureVector[0] / pressureVector[1];
  // calculate 1st and 2nd derivative
  deltaFlow      = backwardFiniteDifference(1, order, deltaT, velocityOutletLSA,
                                            vectorLength);
  deltaDeltaFlow = backwardFiniteDifference(2, order, deltaT, velocityOutletLSA,
                                            vectorLength);
  // calculate pressure BC for next timestep with 3WK model
  pressureOutletLSA[0] =
      rungeKutta(pressureOutletLSA[1], velocityOutletLSA[0], deltaFlow,
                 deltaDeltaFlow, deltaT, ALPHA_LSA, BETA_LSA, GAMMA_LSA);
}

int main(int argc, char* argv[])
{
  // === 1st Step: Initialization ===
  olbInit(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");
  OstreamManager clout(std::cout, "main");
  // display messages from every single mpi process
  clout.setMultiOutput(false);
  UnitConverterFromResolutionAndLatticeVelocity<T, DESCRIPTOR> converter(
      int {resolution},   // resolution: number of voxels per charPhysL
      (T)latticeVelocity, // latticeVelocity: lattice velocity
      (T)charPhysLength,  // charPhysLength: reference length of simulation
                          // geometry
      (T)maxFlowVelocity, // charPhysVelocity: maximal/highest expected
                          // velocity during simulation in m / s
      (T)3.3e-06,         // physical kinematic viscosity in m^2 /s
      (T)1060.0           // physical density in kg / m^3
  );
  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("windkessel3d");

  // print 3WK parameters in log file
  std::cout << std::endl << std::endl;
  std::cout << "Aorta: " << ALPHA_AO << ", " << BETA_AO << ", " << GAMMA_AO
            << std::endl;
  std::cout << "BCT: " << ALPHA_BCT << ", " << BETA_BCT << ", " << GAMMA_BCT
            << std::endl;
  std::cout << "LCCA : " << ALPHA_LCCA << ", " << BETA_LCCA << ", "
            << GAMMA_LCCA << std::endl;
  std::cout << "LSA  : " << ALPHA_LSA << ", " << BETA_LSA << ", " << GAMMA_LSA
            << std::endl;
  std::cout << std::endl << std::endl;

  T   latticeL      = converter.getPhysDeltaX();
  T   deltaT        = converter.getPhysTime(latticeDeltaT);
  int updateTimevtk = converter.getLatticeTime(deltaVTK);

  // === 2nd Step: Prepare Geometry ===
  // Instantiation of the STLreader class
  // file name, voxel size in meter, stl unit in meter
  STLreader<T>        stlReader("AortaShapeFinal.stl",
                                converter.getConversionFactorLength(), 0.001, 0, true);
  IndicatorLayer3D<T> extendedDomain(stlReader,
                                     converter.getConversionFactorLength());
  CuboidGeometry3D<T> cuboidGeometry(extendedDomain,
                                     converter.getConversionFactorLength(),
                                     noOfCuboids, "volume");
  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);
  // Instantiation of a superGeometry
  SuperGeometry<T, 3> superGeometry(cuboidGeometry, loadBalancer);
  prepareGeometry(converter, extendedDomain, stlReader, superGeometry);

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T, DESCRIPTOR> sLattice(superGeometry);
  prepareLattice(sLattice, converter, stlReader, superGeometry);

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer(converter.getLatticeTime(maxPhysT),
                       superGeometry.getStatistics().getNvoxel());
  timer.start();

  // Writes the geometry, cuboid no. and rank no. as vti file for visualization
#pragma region VTKWriterInitialization
  SuperVTMwriter3D<T>                   vtmWriter("windkessel3d");
  SuperLatticeGeometry3D<T, DESCRIPTOR> geometry(sLattice, superGeometry);
  SuperLatticeCuboid3D<T, DESCRIPTOR>   cuboid(sLattice);
  SuperLatticeRank3D<T, DESCRIPTOR>     rank(sLattice);
  vtmWriter.write(geometry);
  vtmWriter.write(cuboid);
  vtmWriter.write(rank);
  vtmWriter.createMasterFile();
#pragma endregion

// definition of velocity and pressure integrals for inlet and outlets
#pragma region measurementIntegralDefinition
  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity(sLattice, converter);
  SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure(sLattice, converter);
  AnalyticalFfromSuperF3D<T>                intpolatePressure(pressure, true);
  AnalyticalFfromSuperF3D<T>                intpolateVelocity(velocity, true);

  // Inlet: Currently not used
  IndicatorCircle3D<T> inflowPressure(
      0., 0., 0. + 10 * converter.getPhysDeltaX(), 0., 0., 1., 0.017);
  SuperPlaneIntegralF3D<T> velocityPlaneInflow(velocity, superGeometry,
                                               inflowPressure, {1, 3});
  SuperPlaneIntegralF3D<T> pressurePlaneInflow(pressure, superGeometry,
                                               inflowPressure, {1, 3});

  // Coronary Arteries: Currently not used
  IndicatorCircle3D<T>     coronary1(0.0343083, 0.0259034,
                                     0. + 10 * converter.getPhysDeltaX(), 0., 0.,
                                     -1., 0.0017647);
  SuperPlaneIntegralF3D<T> velocityPlaneCoronary1(velocity, superGeometry,
                                                  coronary1, {1, 4});
  SuperPlaneIntegralF3D<T> pressurePlaneCoronary1(pressure, superGeometry,
                                                  coronary1, {1, 4});

  IndicatorCircle3D<T>     coronary2(-0.00529951, -0.04390055,
                                     0. + 10 * converter.getPhysDeltaX(), 0., 0.,
                                     -1., 0.0017647);
  SuperPlaneIntegralF3D<T> velocityPlaneCoronary2(velocity, superGeometry,
                                                  coronary2, {1, 5});
  SuperPlaneIntegralF3D<T> pressurePlaneCoronary2(pressure, superGeometry,
                                                  coronary2, {1, 5});

  // Descending outflow
  IndicatorCircle3D<T>     outflowDescending(-0.004515225, 0.07217795,
                                             0. + 10 * converter.getPhysDeltaX(),
                                             0., 0., -1., 0.0109504);
  SuperPlaneIntegralF3D<T> velocityPlaneDescending(velocity, superGeometry,
                                                   outflowDescending, {1, 6});
  SuperPlaneIntegralF3D<T> pressurePlaneDescending(pressure, superGeometry,
                                                   outflowDescending, {1, 6});

  // Ascending outflows
  IndicatorCircle3D<T>     outflowBCT(0.0517593, 0.024923,
                                      0.200954 - 10 * converter.getPhysDeltaX(), 0.,
                                      0., 1., 0.009411768);
  SuperPlaneIntegralF3D<T> velocityPlaneBCT(velocity, superGeometry, outflowBCT,
                                            {1, 7});
  SuperPlaneIntegralF3D<T> pressurePlaneBCT(pressure, superGeometry, outflowBCT,
                                            {1, 7});
  IndicatorCircle3D<T>     outflowLCCA(0.0797742, 0.0311975,
                                       0.200954 - 10 * converter.getPhysDeltaX(),
                                       0., 0., 1., 0.004705884);
  SuperPlaneIntegralF3D<T> velocityPlaneLCCA(velocity, superGeometry,
                                             outflowLCCA, {1, 8});
  SuperPlaneIntegralF3D<T> pressurePlaneLCCA(pressure, superGeometry,
                                             outflowLCCA, {1, 8});
  IndicatorCircle3D<T>     outflowLSA(0.08215115, 0.04629555,
                                      0.200954 - 10 * converter.getPhysDeltaX(), 0.,
                                      0., 1., 0.00490195);
  SuperPlaneIntegralF3D<T> velocityPlaneLSA(velocity, superGeometry, outflowLSA,
                                            {1, 9});
  SuperPlaneIntegralF3D<T> pressurePlaneLSA(pressure, superGeometry, outflowLSA,
                                            {1, 9});
#pragma endregion

  for (std::size_t iT = 0; iT <= converter.getLatticeTime(startUpTime); iT++) {
    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues(sLattice, converter, iT, superGeometry);

    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    // === 7th Step: Output of the results ===
    if (iT % converter.getLatticeTime(deltaVTKstartup) == 0) {
      getResults(sLattice, converter, iT, superGeometry, timer);
    }

    // === 8th Step: update 3WK model
    windkesselUpdater(velocityPlaneDescending, pressurePlaneDescending,
              velocityPlaneBCT, pressurePlaneBCT, velocityPlaneLCCA,
              pressurePlaneLCCA, velocityPlaneLSA, pressurePlaneLSA, deltaT);

    if (std::isnan(velocityOutletDesc[0])) {
      clout << "##### ERROR at TimeStep: " << iT << " #####" << std::endl;
      break;
    }
  }

  for (std::size_t iT = converter.getLatticeTime(startUpTime);
       iT <= converter.getLatticeTime(maxPhysT); iT++) {
    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues(sLattice, converter, iT, superGeometry);

    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    // === 7th Step: Output of the results ===
    if (iT % updateTimevtk == 0) {
      getResults(sLattice, converter, iT, superGeometry, timer);
    }

    // === 8th Step: update 3WK model
    windkesselUpdater(velocityPlaneDescending, pressurePlaneDescending,
              velocityPlaneBCT, pressurePlaneBCT, velocityPlaneLCCA,
              pressurePlaneLCCA, velocityPlaneLSA, pressurePlaneLSA, deltaT);

    if (std::isnan(velocityOutletDesc[0])) {
      clout << "##### ERROR at TimeStep: " << iT << " #####" << std::endl;
      break;
    }
  }

  timer.stop();
  timer.printSummary();
}
