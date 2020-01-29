// Copyright 2012-2017 CNRS-UM LIRMM, CNRS-AIST JRL
//
// This file is part of RBDyn.
//
// RBDyn is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// RBDyn is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with RBDyn.  If not, see <http://www.gnu.org/licenses/>.

// associated header
#include "RBDyn/Hessian.h"

// includes
// std
#include <algorithm>
#include <stdexcept>

// RBDyn
#include "RBDyn/MultiBodyConfig.h"

#include <chrono>
#include <ctime>

#include <iostream>

namespace rbd
{

Hessian::Hessian()
{}

Hessian::Hessian(const MultiBody& mb, const std::string& bodyName,
		const Eigen::Vector3d& point):
  jointsPath_(),
  point_(point),
  jac_(),
	hes_()
{
  int index = mb.sBodyIndexByName(bodyName);

	int dof = 0;
	while(index != -1)
	{
		jointsPath_.insert(jointsPath_.begin(), index);
		dof += mb.joint(index).dof();

		index = mb.parent(index);
	}

	jac_.resize(6, dof);
	hes_.resize(6);
	for (size_t dim = 0; dim<6; dim++)
		hes_[dim] = Eigen::MatrixXd::Zero(dof,dof);
}

/// private implementation of the Jacobian computation
/// We use the Transform template allow Eigen3 to
/// remove the Matrix3d from the computation.
template<typename Transform>
static inline const Eigen::MatrixXd&
jacobian_(const MultiBody& mb, const MultiBodyConfig& mbc,
	const Transform& Trans_0_p, const std::vector<int>& jointsPath,
	Eigen::MatrixXd& jac)
{
	const std::vector<Joint>& joints = mb.joints();
	int curJ = 0;

	sva::PTransformd X_0_p(Trans_0_p);

	for(std::size_t index = 0; index < jointsPath.size(); ++index)
	{
		int i = jointsPath[index];

		sva::PTransformd X_i_N = X_0_p*mbc.bodyPosW[i].inv();

		for(int dof = 0; dof < joints[i].dof(); ++dof)
		{
			jac.col(curJ + dof).noalias() =
				(X_i_N*(sva::MotionVecd(mbc.motionSubspace[i].col(dof)))).vector();
		}

		curJ += joints[i].dof();
	}

	return jac;
}

/// private implementation of the Hessian computation
/// We use the Transform template allow Eigen3 to
/// remove the Matrix3d from the computation.
static inline const std::vector<Eigen::MatrixXd>&
hessian_(const MultiBody& mb, const MultiBodyConfig& mbc,
	const std::vector<int>& jointsPath,
	const Eigen::MatrixXd& jac, std::vector<Eigen::MatrixXd> hes_)
{
  /* Hessian computation according to Erleben 2017 */
  // we go through all dof's on the kinematic chain in two loops i (parent) and j (child)
  // we take the corresponding joint axis in the base frame 0_u_j = J.col(j).head(3)
  // the geometric jacobian is composed as J = [0_u; 0_u x deltap] = [J_rot; J_lin]
  // the axes in the world frame are thus 0_u = J.head(3)
  // we calculate the cross product 0_u_j x J.col(i).head(3) and 0_u_j x J.col(i).tail(3) which are the hessian entries 1 to 6
  // FIXME: handling prismatic joints more elegantly like in the jacobian computation

  // FIXME: quaternions are not handled so far
  auto start = std::chrono::high_resolution_clock::now();

  bool print = false;
  if (print)
  std::cout<<"hessian"<<std::endl;
  if (print)
  std::cout<<"jac:\n"<<jac<<std::endl;
	for (size_t dim = 0; dim<6; dim++)
		hes_[dim] = Eigen::MatrixXd::Zero(jac.cols(),jac.cols()); // since we want H and not B = JTJ + sum(H)

  Eigen::Vector3d jacCol;
  Eigen::Vector3d uxJ;
  Eigen::Vector3d axis;
  Eigen::Vector3d jacCol_w;
  Eigen::Vector3d uxJ_w;
  
	const std::vector<Joint>& joints = mb.joints();

	int curParent = 0;
	for(std::size_t indexParent = 0; indexParent < jointsPath.size(); indexParent++) // outer joints, rowwise
	{
    if (print)
    {
      // std::cout<<"indexParent: "<<indexParent<<std::endl;
      std::cout<<"=============curParent: "<<curParent<<std::endl;
    }
		int i = jointsPath[indexParent];
    if (print)
      std::cout<<"joints[i].type(): "<<joints[i].type()<<" revolute: "<<Joint::Rev<<" prismatic: "<<Joint::Prism<<std::endl;

    for (int dof_i = 0; dof_i < joints[i].dof(); ++dof_i) // colwise
    {
      // handle prismatic joints
      if (joints[i].type() == Joint::Prism || (joints[i].type() == Joint::Free && dof_i > 2))
      {
        for (int dim=0; dim<6; dim++)
        {
          hes_[dim].coeffRef(curParent+dof_i,curParent+dof_i) = 1;
        }
        continue;
      }
      // revolute joints
      else
      {
        jacCol_w = jac.col(curParent+dof_i).head(3);
        jacCol = jac.col(curParent+dof_i).tail(3);

		    int curChild = 0;
		    for(std::size_t indexChild = 0; indexChild <= indexParent; indexChild++) // inner joints, columnwise
		    {
          if (print)
          {
            // std::cout<<"indexChild: "<<indexChild<<std::endl;
            std::cout<<"-------------curChild: "<<curChild<<std::endl;
          }

		      int j = jointsPath[indexChild];

		        for(int dof_j = 0; dof_j < joints[j].dof(); ++dof_j) // rowwise
		        {

              if (joints[j].type() == Joint::Prism || (joints[j].type() == Joint::Free && dof_j > 2))
                continue;
              {
                if (joints[j].type() == Joint::Free) // FIXME: why is this required while it is not necessary for the Jacobian calculation
                {
                  switch (dof_j) // according to xyz-Euler convention of RBDYN
                  {
                    case 0:{
                      axis = mbc.bodyPosW[j].rotation().transpose() * mbc.motionSubspace[j].col(dof_j).head(3);
                      break;
                    }
                    case 1:{
                      axis =
                        mbc.bodyPosW[j].rotation().transpose() * sva::RotX(mbc.q[j][1]) * mbc.motionSubspace[j].col(dof_j).head(3);
                      break;
                    }
                    case 2:{
                      axis =
                        mbc.bodyPosW[j].rotation().transpose() * sva::RotX(mbc.q[j][1]) * sva::RotY(mbc.q[j][2]) * mbc.motionSubspace[j].col(dof_j).head(3);
                      break;
                    }
                  }
                }
                else
                {
                  axis = jac.col(curChild+dof_j).head(3);
                }

                if (print)
                {
                  std::cout<<"untransformed axis: "<<mbc.motionSubspace[j].col(dof_j).head(3).transpose()<<std::endl;
                  std::cout<<"transformed axis: "<<axis.transpose()<<std::endl;
                }

                uxJ_w = axis.cross(jacCol_w);
                uxJ = axis.cross(jacCol); // + mbc.motionSubspace[j].col(dof_j).tail(3).cwiseProduct(mbc.motionSubspace[i].col(dof_i).tail(3));
                // std::cout<<"axis "<<mbc.motionSubspace[j].col(dof_j).tail(3).transpose()<<" axis: "<<mbc.motionSubspace[i].col(dof_i).tail(3).transpose()<<": "<<(mbc.motionSubspace[j].col(dof_j).tail(3).cwiseProduct(mbc.motionSubspace[i].col(dof_i).tail(3))).transpose()<<std::endl;

                if (print){
                  std::cout<<"jacCol: "<<jacCol.transpose()<<std::endl;
                  std::cout<<"uxJ: "<<uxJ.transpose()<<std::endl;
                  std::cout<<"jacCol_w: "<<jacCol_w.transpose()<<std::endl;
                  std::cout<<"uxJ_w: "<<uxJ_w.transpose()<<std::endl;
                }

                if (print)
                  std::cout<<">>>>H[ "<<curParent + dof_i<<", "<<curChild + dof_j<<" ]"<<std::endl;
		            for (int dim = 0; dim < 3; dim++)
		            {
	 	            	hes_[dim](curParent + dof_i, curChild + dof_j) += uxJ_w[dim]; // -uxJ if mulitplied directly by error; + because we multiply by w ~ -e
                  if (curParent + dof_i != curChild + dof_j)
	 	            	  hes_[dim](curChild + dof_j, curParent + dof_i) = hes_[dim](curParent + dof_i, curChild + dof_j);
                }
                for (int dim = 3; dim < 6; dim++)
		            {
	 	            	hes_[dim](curParent + dof_i, curChild + dof_j) += uxJ[dim-3]; // -uxJ if mulitplied directly by error; + because we multiply by w ~ -e
                  if (curParent + dof_i != curChild + dof_j)
	 	            	  hes_[dim](curChild + dof_j, curParent + dof_i) = hes_[dim](curParent + dof_i, curChild + dof_j);
                }
              }
            }
		      curChild += joints[j].dof();
        }
		  }
    }
		curParent += joints[i].dof();
	}

  if (print){
    std::cout<<"hes in RBDyn:\n";
    std::cout<<"H[0] = \n"<<hes_[0]<<std::endl;
    std::cout<<"H[1] = \n"<<hes_[1]<<std::endl;
    std::cout<<"H[2] = \n"<<hes_[2]<<std::endl;
    std::cout<<"H[3] = \n"<<hes_[3]<<std::endl;
    std::cout<<"H[4] = \n"<<hes_[4]<<std::endl;
    std::cout<<"H[5] = \n"<<hes_[5]<<std::endl;
  }
  // std::cout<<"hes sum:\n"<<hes_[0]+hes_[1]+hes_[2]+hes_[3]+hes_[4]+hes_[5]<<std::endl;

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> difference = end - start;
  // std::cout<<"duration task hessian: "<<difference.count()<<"[s]"<<std::endl;

	return hes_;
}

const std::vector<Eigen::MatrixXd>&
Hessian::hessian(const MultiBody& mb, const MultiBodyConfig& mbc)
{
	int N = jointsPath_.back();

	// the transformation must be read {}^0E_p {}^pT_N {}^NX_0
  Eigen::Vector3d T_0_Np((point_ * mbc.bodyPosW[N]).translation());
  jac_ = jacobian_(mb, mbc, T_0_Np, jointsPath_, jac_); /* defined in Jacobian.cpp */
	hes_ = hessian_(mb, mbc, jointsPath_, jac_, hes_);
  return hes_;
}

const std::vector<Eigen::MatrixXd>&
Hessian::hessian_preCompJ(const MultiBody& mb, const MultiBodyConfig& mbc, const Eigen::MatrixXd& jac)
{
	int N = jointsPath_.back();

	// the transformation must be read {}^0E_p {}^pT_N {}^NX_0
  Eigen::Vector3d T_0_Np((point_ * mbc.bodyPosW[N]).translation());
	hes_ = hessian_(mb, mbc, jointsPath_, jac, hes_);
  return hes_;
}

const std::vector<Eigen::MatrixXd>&
Hessian::hessian_complete(const MultiBody& mb, const MultiBodyConfig& mbc)
{
	int N = jointsPath_.back();

	// the transformation must be read {}^0E_p {}^pT_N {}^NX_0
  Eigen::Vector3d T_0_Np((point_ * mbc.bodyPosW[N]).translation());
  jac_ = jacobian_(mb, mbc, T_0_Np, jointsPath_, jac_); /* defined in Jacobian.cpp */
	hes_ = hessian_(mb, mbc, jointsPath_, jac_, hes_);
  return hes_;
}

const std::vector<Eigen::MatrixXd>&
Hessian::hessian_complete_preCompJ(const MultiBody& mb, const MultiBodyConfig& mbc, const Eigen::MatrixXd& jac)
{
	int N = jointsPath_.back();

	// the transformation must be read {}^0E_p {}^pT_N {}^NX_0
  Eigen::Vector3d T_0_Np((point_ * mbc.bodyPosW[N]).translation());
	hes_ = hessian_(mb, mbc, jointsPath_, jac, hes_);
  return hes_;
}

void Hessian::fullHessian(const MultiBody& mb,
	const Eigen::Ref<const Eigen::MatrixXd>& hes,
	Eigen::MatrixXd& res) const
{
  res.block(0, 0, mb.nrDof(), mb.nrDof()).setZero();
  int hesPos_i = 0;
  for(std::size_t index = 0; index < jointsPath_.size(); ++index)
  {
    int i = jointsPath_[index];
    int dof_i = mb.joint(i).dof();
    int hesPos_j = 0;
    for(std::size_t jndex = 0; jndex < jointsPath_.size(); ++jndex)
    {
      int j = jointsPath_[jndex];
      int dof_j = mb.joint(j).dof();
    	res.block(mb.jointPosInDof(i), mb.jointPosInDof(j), dof_i, dof_j) =
              hes.block(hesPos_i, hesPos_j, dof_i, dof_j);
    	hesPos_j += dof_j;
    }
    hesPos_i += dof_i;
  }
}

} // namespace rbd
