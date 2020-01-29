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

#pragma once

// includes
// RBDyn
#include "MultiBody.h"

#include <rbdyn/config.hh>

namespace rbd
{
struct MultiBodyConfig;

/**
	* Algorithm to compute the jacobian of a specified body.
	*/
class RBDYN_DLLAPI Hessian
{
public:
	Hessian();

	/**
		* Create a hessian from the root body to the specified body.
		* @param mb Multibody where bodyId is in.
		* @param bodyName Specified body.
		* @param point Point in the body exprimed in body coordinate.
		* @throw std::out_of_range If bodyId don't exist.
		*/
	Hessian(const MultiBody& mb, const std::string& bodyName,
		const Eigen::Vector3d& point=Eigen::Vector3d::Zero());

	  /**
		* Compute the hessian in world frame.
		* @param mb MultiBody used has model.
		* @param mbc Use bodyPosW and motionSubspace.
		* @return Returns H of B = JTJ + sum_i e_i H_i with i=1,2,3 of mb with mbc configuration.
		* @return which corresponds to the second order derivatives of f
		*/
	const std::vector<Eigen::MatrixXd>&
		hessian(const MultiBody& mb, const MultiBodyConfig& mbc);
  /* with precomputed jacobian */
	const std::vector<Eigen::MatrixXd>&
		hessian_preCompJ(const MultiBody& mb, const MultiBodyConfig& mbc,
			 		  const Eigen::MatrixXd& jac);

	/**
		* Compute the hessian in world frame.
		* @param mb MultiBody used has model.
		* @param mbc Use bodyPosW and motionSubspace.
		* @return Returns B of B = JTJ + sum_i e_i H_i with i=1,2,3 of mb with mbc configuration.
		*/
	const std::vector<Eigen::MatrixXd>&
		hessian_complete(const MultiBody& mb, const MultiBodyConfig& mbc);
  /* with precomputed jacobian */
	const std::vector<Eigen::MatrixXd>&
		hessian_complete_preCompJ(const MultiBody& mb, const MultiBodyConfig& mbc,
			 		  const Eigen::MatrixXd& jac);

	  /**
		* Project the hessian in the full robot parameters vector.
		* @param mb MuliBody used has model.
		* @param hes Hessian to project.
		* @param res Projected Hessian (must be allocated).
		*/
	void fullHessian(const MultiBody& mb,
		const Eigen::Ref<const Eigen::MatrixXd>& hes,
		Eigen::MatrixXd& res) const;

	/// @return The number of degree of freedom in the joint path
	int dof() const
	{
		return static_cast<int>(hes_[0].cols());
	}

	/// @param point Static translation in the body exprimed in body coordinate.
	void point(const Eigen::Vector3d& point)
	{
		point_ = sva::PTransformd(point);
	}

private:
	std::vector<int> jointsPath_;
	sva::PTransformd point_;

	Eigen::MatrixXd jac_;
	std::vector<Eigen::MatrixXd> hes_;
};

} // namespace rbd
