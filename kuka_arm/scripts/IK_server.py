#!/usr/bin/env python

# Copyright (C) 2017 Udacity Inc.
#
# This file is part of Robotic Arm: Pick and Place project for Udacity
# Robotics nano-degree program
#
# All Rights Reserved.

# Author: Harsh Pandya

# import modules
import rospy
import tf
from kuka_arm.srv import *
from trajectory_msgs.msg import JointTrajectory, JointTrajectoryPoint
from geometry_msgs.msg import Pose
from mpmath import *
from sympy import *


def handle_calculate_IK(req):
    rospy.loginfo("Received %s eef-poses from the plan" % len(req.poses))
    if len(req.poses) < 1:
        print "No valid poses received"
        return -1
    else:

        ### Your FK code here
        # Create symbols
	alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6 = symbols('alpha0:7')
	a0, a1, a2, a3, a4, a5, a6 = symbols('a0:7')
	d1, d2, d3, d4, d5, d6, d7 = symbols('d1:8')
	q1, q2, q3, q4, q5, q6, q7 = symbols('q1:8') # theta variables
	
	# Create Modified DH parameters
	dh_table = {alpha0:      0, a0:      0, d1:  0.75, q1:        q1, 
		alpha1: -pi/2., a1:   0.35, d2:     0, q2: -pi/2.+q2, 
                alpha2:      0, a2:   1.25, d3:     0, q3:        q3, 
                alpha3: -pi/2., a3: -0.054, d4:   1.5, q4:        q4, 
                alpha4:  pi/2., a4:      0, d5:     0, q5:        q5, 
                alpha5: -pi/2., a5:      0, d6:     0, q6:        q6, 
                alpha6:      0, a6:      0, d7: 0.303, q7:         0} 

	# Define Modified DH Transformation matrix
	def tf_matrix(alpha, a, d, q):
		tf = Matrix([[            cos(q),           -sin(q),           0,             a],
			[ sin(q)*cos(alpha), cos(q)*cos(alpha), -sin(alpha), -sin(alpha)*d],
			[ sin(q)*sin(alpha), cos(q)*sin(alpha),  cos(alpha),  cos(alpha)*d],
			[                 0,                 0,           0,             1]])
		return tf


	# Create individual transformation matrices
	T0_1 = tf_matrix(alpha0, a0, d1, q1).subs(dh_table)
	T1_2 = tf_matrix(alpha1, a1, d2, q2).subs(dh_table)
	T2_3 = tf_matrix(alpha2, a2, d3, q3).subs(dh_table)
	T3_4 = tf_matrix(alpha3, a3, d4, q4).subs(dh_table)
	T4_5 = tf_matrix(alpha4, a4, d5, q5).subs(dh_table)
	T5_6 = tf_matrix(alpha5, a5, d6, q6).subs(dh_table)
	T6_E = tf_matrix(alpha6, a6, d7, q7).subs(dh_table)
	# print("T0_1=", T0_1)
	# print("T1_2=", T1_2)
	# print("T2_1=", T2_3)
	# print("T3_4=", T3_4)
	# print("T4_5=", T4_5)
	# print("T5_6=", T5_6)
	# print("T6_E=", T6_E)

	# Generate homogeneous transform matrix from base_link to gripper_link
	T0_E = T0_1 * T1_2 * T2_3 * T3_4 * T4_5 * T5_6 * T6_E
	# print("T0_E=", T0_E)

	# Extract rotation matrices from the transformation matrices
	r, p, y = symbols('r p y')

	R_x = Matrix([[       1,       0,       0],
			[       0,  cos(r), -sin(r)],
                        [       0,  sin(r),  cos(r)]])

	R_y = Matrix([[  cos(p),       0,  sin(p)],
                        [       0,       1,       0],
                        [ -sin(p),       0,  cos(p)]])

	R_z = Matrix([[  cos(y), -sin(y),       0],
                        [  sin(y),  cos(y),       0],
                        [       0,       0,       1]])


	R_correct = R_z.evalf(subs={'y': pi}) * R_y.evalf(subs={'p': -pi/2})

	Rrpy = R_z * R_y * R_x * R_correct


	R0_3 = T0_1[0:3,0:3] * T1_2[0:3,0:3] * T2_3[0:3,0:3]


        # Initialize service response
        joint_trajectory_list = []
        for x in xrange(0, len(req.poses)):
		# IK code starts here
		joint_trajectory_point = JointTrajectoryPoint()

		# Extract end-effector position and orientation from request
		# px,py,pz = end-effector position
		# roll, pitch, yaw = end-effector orientation
		px = req.poses[x].position.x
		py = req.poses[x].position.y
		pz = req.poses[x].position.z

		(roll, pitch, yaw) = tf.transformations.euler_from_quaternion(
		[req.poses[x].orientation.x, req.poses[x].orientation.y,
		    req.poses[x].orientation.z, req.poses[x].orientation.w])

		### Your IK code here
		R_EE = Rrpy.subs({'r': roll, 'p': pitch, 'y': yaw})
		P_WC = Matrix([[px], [py], [pz]]) - (0.303) * R_EE[:,2]

		# Calculate joint angles using Geometric IK method
		theta1 = atan2(P_WC[1], P_WC[0])

		A = sqrt(1.5**2 + 0.054**2)
		B = sqrt((sqrt(P_WC[0]**2 + P_WC[1]**2) - 0.35)**2 + (P_WC[2] - 0.75)**2)
		C = 1.25

		a = acos(( C**2 + B**2 - A**2 ) / (2*A*B))
		b = acos(( A**2 + C**2 - B**2 ) / (2*A*C))
		c = acos(( A**2 + B**2 - C**2 ) / (2*B*C))

		theta2 = pi/2 - a - atan2(P_WC[2]-0.75, sqrt(P_WC[0]**2+P_WC[1]**2)-0.35)
		theta3 = pi/2 - b - atan2(0.054, 1.5) 

		# print(R0_3.inv("LU") * R_EE)

		R0_3 = R0_3.evalf(subs={q1: theta1, q2: theta2, q3:theta3})
		R3_6 = R0_3.inv("LU") * R_EE

		# print(R3_6) # debug end effector rotation matrix

		theta5 = atan2(sqrt(R3_6[0,2]**2 + R3_6[2,2]**2),R3_6[1,2])
		
		if (sin(theta5) > 0):
		 	theta4 = atan2(-R3_6[2,2], R3_6[0,2])
		 	theta6 = atan2(R3_6[1,1], -R3_6[1,0])
		else:
		 	theta4 = atan2(R3_6[2,2], -R3_6[0,2])
		 	theta6 = atan2(-R3_6[1,1], R3_6[1,0])

		# Populate response for the IK request
		# In the next line replace theta1,theta2...,theta6 by your joint angle variables
		joint_trajectory_point.positions = [theta1, theta2, theta3, theta4, theta5, theta6]
		joint_trajectory_list.append(joint_trajectory_point)

        rospy.loginfo("length of Joint Trajectory List: %s" % len(joint_trajectory_list))
        return CalculateIKResponse(joint_trajectory_list)


def IK_server():
    # initialize node and declare calculate_ik service
    rospy.init_node('IK_server')
    s = rospy.Service('calculate_ik', CalculateIK, handle_calculate_IK)
    print "Ready to receive an IK request"
    rospy.spin()

if __name__ == "__main__":
    IK_server()
