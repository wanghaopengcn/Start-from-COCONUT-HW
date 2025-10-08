// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NewtonMethod/NewtonMethodMHD.hh"


#include "NewtonMethod/UpdateSolMHD.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Common/CFLog.hh"
#include "Framework/State.hh"
#include "Common/BadValueException.hh"

#include "Framework/SubSystemStatus.hh"
#include "Framework/SpaceMethod.hh"
#include "Framework/SpaceMethodData.hh"
#include "Framework/ConvectiveVarSet.hh"

#include "Framework/ConvergenceMethod.hh"
#include "Framework/ConvergenceMethodData.hh"
#include "Framework/ConvergenceStatus.hh"
#include "NewtonMethod/NewtonIteratorData.hh"
#include "NewtonMethod/NewtonIterator.hh"

#include "Framework/PathAppender.hh"
#include "Framework/GeometricEntity.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<UpdateSolMHD, NewtonIteratorData, NewtonMethodMHDModule> 
updateSolMHDProvider("UpdateSolMHD");

//////////////////////////////////////////////////////////////////////////////

void UpdateSolMHD::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption<vector<CFreal>,Config::DynamicOption<> >("Relaxation","Relaxation factor");
  options.addConfigOption< bool >("Validate","Check that each update creates variables with physical meaning");
  options.addConfigOption< CFreal >("pressureCorrectionValue","the correction value for the pressure in negative pressure occurring states, should be a very small positive number.");
}

//////////////////////////////////////////////////////////////////////////////

UpdateSolMHD::UpdateSolMHD(const std::string& name) : 
  NewtonIteratorCom(name),
  socket_states("states"),
  socket_rhs("rhs"),
  socket_updateCoeff("updateCoeff"),
  _model(),
  _pressureCorrectionVal()
{
  addConfigOptionsTo(this);
   
  m_alpha = vector<CFreal>();
  setParameter("Relaxation",&m_alpha);

  m_validate = false;
  setParameter("Validate",&m_validate);

  _pressureCorrectionVal = MathTools::MathConsts::CFrealEps();
  setParameter("pressureCorrectionValue",&_pressureCorrectionVal);
}

//////////////////////////////////////////////////////////////////////////////

void UpdateSolMHD::setup()
{
  NewtonIteratorCom::setup();

  _model = PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm().d_castTo<MHDTerm>();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
 
  if (m_alpha.size() == 0) {
    m_alpha.resize(nbEqs);
    for (CFuint i = 0; i < nbEqs; ++i) {
      m_alpha[i] = 1.;
    }
  }

  if (m_alpha.size() == 1) {
    const CFreal value = m_alpha[0];
    m_alpha.resize(nbEqs);
    for (CFuint i = 0; i < nbEqs; ++i) {
      m_alpha[i] = value;
    }
  }

  if (m_alpha.size() != nbEqs) {
    throw BadValueException (FromHere(),"UpdateSolMHD::setup() : m_alpha.size() != nbEqs");
  }
}

//////////////////////////////////////////////////////////////////////////////

void UpdateSolMHD::execute()
{
  CFAUTOTRACE;
  
  DataHandle < Framework::State*, Framework::GLOBAL > states  = socket_states.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  
  // rhs is the temporary placeholder for the dU
  DataHandle<CFreal>& dU = rhs;

  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint states_size = states.size();

  CFuint nbStatesWithNegPressure = 0;
  CFuint nbStatesWithNegrho = 0;

  SafePtr<FilterState> filterState = getMethodData().getFilterState();
  SafePtr<FilterRHS> filterRHS     = getMethodData().getFilterRHS();
  
  CFreal T_Sun = 1.5e6;   // K
  CFreal mu_cor = 1.27;   // Mean molecular weight
  CFreal mH = 1.67e-27;   // Mass hydrogen
  CFreal Bref = 2.2e-4;   // T
  CFreal mu0 = 1.2566e-6;
  CFreal rhoref = 1.67e-13;
  CFreal kB = 1.38e-23;
  CFreal pref = 0.03851;
  CFreal vref = 4.8e5; 
  CFreal smooth_factor = 0.0;
  
  for (CFuint iState = 0; iState < states_size; ++iState)
  {
    State& cur_state = *states[iState];
    // do the update only if the state is parallel updatable
    if (cur_state.isParUpdatable())
    {
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
      {
	filterRHS->filter(iEq, dU(iState, iEq, nbEqs));
	cur_state[iEq] += m_alpha[iEq] * dU(iState, iEq, nbEqs);
      }

      // checking the pressure value in the state and correction in case of negative pressure 
      // (especially necessary for solar wind/planet magnetosphere interaction)  

/*
  smooth_factor = 0.0;
  CFreal dist = (minPlasmaBeta - plasmaBeta)/2.0e-6; //if current beta is 0.019, then dist is positive, smooth factor is close to 1
  smooth_factor = 0.5 + 0.5 * std::tanh(dist * 3.141592654);  //if the current beta is  0.3, then dist is negative and so smooth factor is 0
  cur_state[7] = (smooth_factor * pmag * minPlasmaBeta + pth * (1.0-smooth_factor))/ pref;

  PressureBoundary_dimless = cur_state[7]; 
  pth = PressureBoundary_dimless * pref;
  plasmaBeta = pth/pmag;
  
  if ((plasmaBeta < 0.02) || (pth < 1.e-9)) {
		cur_state[7] -= m_alpha[7] * dU(iState, 7, nbEqs);
		nbStatesWithNegPressure += 1;
      } 
*/
/*	  
      //if (cur_state[7] < 0.0) {
	  if (cur_state[7] < 1.d-9) {  
        const CFreal pNegative = cur_state[7];
	// the pressure value is corrected if it is initially negative
        //cur_state[7] = _pressureCorrectionVal;
//		cur_state[7] -= m_alpha[7] * dU(iState, 7, nbEqs);

//        nbStatesWithNegPressure += 1;

        // position of the negative pressure occurring state is also important

        //RealVector stateCoord(PhysicalModelStack::getActive()->getDim());
        stateCoord = cur_state.getCoordinates();
        if (PhysicalModelStack::getActive()->getDim() == 2)
          cout << "Pressure was " << pNegative << " in " << iState << ". state with coordinates ("
             << stateCoord[0] << "," << stateCoord[1] << ") and is corrected to " << _pressureCorrectionVal << "." << endl;
        else
          cout << "Pressure was " << pNegative << " in " << iState << ". state with coordinates ("
             << stateCoord[0] << "," << stateCoord[1] << "," << stateCoord[2] << ") and is corrected to " << _pressureCorrectionVal << "." << endl;
      }
*/	 
	  
  CFreal x = cur_state.getCoordinates()[XX];
  CFreal y = cur_state.getCoordinates()[YY];
  CFreal z = cur_state.getCoordinates()[ZZ];
  CFreal r2 = x*x + y*y + z*z;
  CFreal r = std::sqrt(r2);
  
  if (r < 1.15){
	  CFreal Bmag = std::sqrt(cur_state[4] * cur_state[4] + cur_state[5] * cur_state[5] + cur_state[6] * cur_state[6]);
	  CFreal minPlasmaBeta = 0.02;
	  //CFreal minPlasmaBeta = 0.0001;
	  CFreal maxPlasmaBeta = 1.0;
	  CFreal pmag = (Bmag*Bref)*(Bmag*Bref) / 2. / mu0;

	  CFreal PressureBoundary_dimless = cur_state[7];
	  CFreal pth = PressureBoundary_dimless * pref;
	  CFreal plasmaBeta = pth / pmag;

	  //CFreal Vr = x / r*cur_state[1] + y / r*cur_state[2] + z / r*cur_state[3];
	  //CFreal Vmagmax = -10.0e3 / vref;
	  //if (Vr < Vmagmax) {
	//	  cur_state[1] -= m_alpha[1] * dU(iState, 1, nbEqs);
	//	  cur_state[2] -= m_alpha[2] * dU(iState, 2, nbEqs);
	//	  cur_state[3] -= m_alpha[3] * dU(iState, 3, nbEqs);
	 // }
	  //if (std::abs(dU(iState, 0, nbEqs)) > std::abs(cur_state[0])){
	  //}

	  CFreal densityBoundary_dimless = cur_state[0];
	  CFreal maxvA = 3.0e6;
	  CFreal vA = Bmag*Bref / pow((densityBoundary_dimless * rhoref * mu0), 0.5);

	  smooth_factor = 0.0;
	  CFreal dist = (vA - maxvA) / 2.0e3;
	  smooth_factor = 0.5 + 0.5 * std::tanh(dist * 3.141592654);

	  CFreal rho_temp = (smooth_factor * (Bmag*Bref)*(Bmag*Bref) / (maxvA*maxvA) / mu0 + densityBoundary_dimless * rhoref * (1.0 - smooth_factor)) / rhoref;
	  if (rho_temp < _pressureCorrectionVal) {
		  cur_state[0] -= m_alpha[0] * dU(iState, 0, nbEqs);
		  rho_temp = cur_state[0];
		  nbStatesWithNegrho += 1;
	  }
	  cur_state[0] = rho_temp;


	  CFreal p_temp = cur_state[7];
	  //(pth < 1.e-9)
	  if ((plasmaBeta < 0.0001) || (pth < 1.e-5)){
	  //if (p_temp < 1.e-3){
		 cur_state[7] -= m_alpha[7] * dU(iState, 7, nbEqs);
		 p_temp = cur_state[7];
		 nbStatesWithNegrho += 1;
	  }

	  smooth_factor = 0.0;
	  pth = p_temp*pref;
	  plasmaBeta = pth / pmag;
	  dist = (minPlasmaBeta - plasmaBeta) / 2.0e-6; //if current beta is 0.019, then dist is positive, smooth factor is close to 1
	  smooth_factor = 0.5 + 0.5 * std::tanh(dist * 3.141592654);  //if the current beta is  0.3, then dist is negative and so smooth factor is 0
	  p_temp = (smooth_factor * pmag * minPlasmaBeta + pth * (1.0 - smooth_factor)) / pref;
	  cur_state[7] = p_temp;
  }

      // apply a polymorphic filter to the state
      filterState->filter(cur_state);
    }
    else {
      // reset to 0 the RHS for ghost states in order to avoid 
      // inconsistencies in the parallel L2 norm computation
      cf_assert(!cur_state.isParUpdatable());
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
      {
	dU(iState, iEq, nbEqs) = 0.;
      }
    }
  }

  if (nbStatesWithNegPressure > 0)
    cout << "There were " << nbStatesWithNegPressure << " states with negative pressure." << endl;
  if (nbStatesWithNegrho > 0)
    cout << "There were " << nbStatesWithNegrho << " states with negative rho." << endl;
  
  ///This vector will temporally contain the ID of any vector considered not valid.
  std::vector<CFuint> badStatesIDs;

  ///Name of the file where to store invalid states in.
  boost::filesystem::path filepath ( "unphysical_states.plt");

  filepath = PathAppender::getInstance().appendParallel( filepath );

  // loop for unphysicalness check:
  if ( m_validate )
  {
    Common::SafePtr<SpaceMethod> theSpaceMethod = getMethodData().getCollaborator<SpaceMethod>();
    Common::SafePtr<SpaceMethodData> theSpaceMethodData = theSpaceMethod->getSpaceMethodData();
    Common::SafePtr<Framework::ConvectiveVarSet> theVarSet = theSpaceMethodData->getUpdateVar();

    ///Reset the vector:
    badStatesIDs.resize(0);

    for (CFuint iState = 0; iState < states_size; ++iState)
    {
      State& cur_state = *states[iState];
      if ( !( theVarSet->isValid(cur_state) ) )
      {
        badStatesIDs.push_back( iState );
      }

    } // for all states
  } // if validate

  if ( badStatesIDs.size() > 0 )
  {
    ///Gets the global iteration
    CFuint cur_global_iter = SubSystemStatusStack::getActive()->getNbIter();
    ///Gets the Newton method procedure iteration
    CFuint cur_iter = getMethodData().getConvergenceStatus().iter;

    ///opens a file to put the states which have unphysical variables
    ofstream fout ( filepath.string().c_str(), ios::app );

    if ( !fout.is_open() )
    {
      cout << endl << "File [" << filepath.string().c_str() << "]  is not open.\n";
    }
    ///
    const std::vector<std::string>& theVarNames=getMethodData().getCollaborator<SpaceMethod>()->getSpaceMethodData()->getUpdateVar()->getVarNames();
    CFuint dim = PhysicalModelStack::getActive()->getDim();
    ///

    ///Zone header:
    fout << "VARIABLES = \"X\"\n ";
    fout << "\"Y\"" << endl;
    if (dim == 3 ){
        fout << "\"Z\"" << endl;
    }
    for (CFuint i = 0; i < theVarNames.size(); i++){
      fout << "\"" << theVarNames[i] << "\"" << endl;
    }
    ///


    CFuint badStates_size = badStatesIDs.size();

    ///Line  
    fout << "ZONE T = \"Global iteration "<< cur_global_iter << ", Newton iteration " << cur_iter << "\"" << endl;
    ///Line
    fout << "I="<< badStates_size << ", J=1, K=1, ZONETYPE=Ordered"<< endl;
    ///Line
    fout << "DATAPACKING=POINT" << endl;
    ///Line
    fout << "DT=(";
    for ( CFuint i = 0 ; i < dim + theVarNames.size(); i++){
      fout << "DOUBLE ";
    }

    fout << ")" << endl;
    ///

    for (CFuint index = 0; index < badStates_size; ++index)
    {
     State& invalid_state  = *states[ badStatesIDs[index] ];
     fout << invalid_state.getCoordinates() << " " << invalid_state << std::endl;
    } // for all invalid states

    fout.close();

  } // writing badStates

  // reset to 0 the update coefficient
  updateCoeff = 0.0;
}

//////////////////////////////////////////////////////////////////////

void UpdateSolMHD::configure ( Config::ConfigArgs& args )
{
  NewtonIteratorCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > UpdateSolMHD::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;
  
  result.push_back(&socket_states);
  result.push_back(&socket_rhs);
  result.push_back(&socket_updateCoeff);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD
