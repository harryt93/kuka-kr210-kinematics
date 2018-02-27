#!/usr/bin/env python

# Copyright (C) 2017 Electric Movement Inc.
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

        # Define D-H parameters
        #######################

        d1, d2, d3, d4, d5, d6, d7 = symbols('d1:8')   
        q1, q2, q3, q4, q5, q6, q7 = symbols('q1:8') 
        q1, q2, q3, q4, q5, q6, q7 = symbols('q1:8')
        a0, a1, a2, a3, a4, a5, a6 = symbols('a0:7')
        alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6 = symbols('alpha0:7')

        # Distance and Angle (in meters and radians)
        d0_1 = 0.75
        d3_4 = 1.5
        d4_7 = 0.303

        a1_2 = 0.35
        a2_3 = 1.25
        a3_4 = -0.054

        alpha1_2 = -pi/2
        alpha3_4 = -pi/2
        alpha4_5 =  pi/2
        alpha5_6 = -pi/2


        # Construct DH Parameter Table
        ##############################

        DH_Table = {    alpha0:      0, a0:      0,   d1:  0.75,   q1:         q1,
                        alpha1:  -pi/2, a1:   0.35,   d2:     0,   q2: -pi/2 + q2,
                        alpha2:      0, a2:   1.25,   d3:     0,   q3:         q3,
                        alpha3:  -pi/2, a3: -0.054,   d4:   1.5,   q4:         q4,
                        alpha4:   pi/2, a4:      0,   d5:     0,   q5:         q5,
                        alpha5:  -pi/2, a5:      0,   d6:     0,   q6:         q6,
                        alpha6:      0, a6:      0,   d7: 0.303,   q7:          0,}


        # Define a function for quickly generating homogeneous transforms
        #################################################################

        def generate_homogeneous_transform(A, a, d, q):

                m = Matrix([[          cos(q),       -sin(q),       0,         a],
                          [ sin(q)*cos(A), cos(q)*cos(A), -sin(A), -sin(A)*d],
                          [ sin(q)*sin(A), cos(q)*sin(A),  cos(A),  cos(A)*d],
                          [             0,             0,       0,         1]])
                return m


        # Construct individual transformation matrices
        ###########################################
    
        T0_1 = generate_homogeneous_transform(alpha0,a0,d1,q1).subs(DH_Table)
        T1_2 = generate_homogeneous_transform(alpha1,a1,d2,q2).subs(DH_Table)
        T2_3 = generate_homogeneous_transform(alpha2,a2,d3,q3).subs(DH_Table)
        T3_4 = generate_homogeneous_transform(alpha3,a3,d4,q4).subs(DH_Table)
        T4_5 = generate_homogeneous_transform(alpha4,a4,d5,q5).subs(DH_Table)
        T5_6 = generate_homogeneous_transform(alpha5,a5,d6,q6).subs(DH_Table)
        T6_7 = generate_homogeneous_transform(alpha6,a6,d7,q7).subs(DH_Table)


        # Construct the full transformation from base link to end effector
        ##################################################################

        T0_2 = T0_1 * T1_2
        T0_3 = T0_2 * T2_3
        T0_4 = T0_3 * T3_4
        T0_5 = T0_4 * T4_5
        T0_6 = T0_5 * T5_6
        T0_7 = T0_6 * T6_7


        # Extract rotation matrices
        ###########################

        R0_1 = T0_1[0:3,0:3]
        R0_2 = T0_2[0:3,0:3]
        R0_3 = T0_3[0:3,0:3]
        R0_4 = T0_4[0:3,0:3]
        R0_5 = T0_5[0:3,0:3]
        R0_6 = T0_6[0:3,0:3]
        R0_7 = T0_7[0:3,0:3]

        # Initialize service response
        joint_trajectory_list = []
        for x in xrange(0, len(req.poses)):
            joint_trajectory_point = JointTrajectoryPoint()


            # Extract end effector positions
            #####################################

            px = req.poses[x].position.x
            py = req.poses[x].position.y
            pz = req.poses[x].position.z


            # Extract R, P, and Y
            #####################################

            (r, p, y) = tf.transformations.euler_from_quaternion(
                [req.poses[x].orientation.x, req.poses[x].orientation.y,
                req.poses[x].orientation.z, req.poses[x].orientation.w])

            # R_rpy = roll, pitch, yaw
            R_rpy =  Matrix([[cos(y)*cos(p), cos(y)*sin(p)*sin(r)-sin(y)*cos(r), cos(y)*sin(p)*cos(r)+sin(y)*sin(r)],
                            [sin(y)*cos(p), sin(y)*sin(p)*sin(r)+cos(y)*cos(r), sin(y)*sin(p)*cos(r)-cos(y)*sin(r)],
                            [-sin(p),         cos(p)*sin(r),                             cos(p)*cos(r)]])

            # Compute Correction Matrix
            #####################################

            R_corr = Matrix([[0,0,1],
                            [0,-1,0],
                             [1,0,0]])
            
            # Calculate the origin matrix by the rpy and correlation transpose
            Rrpy = R_rpy*(R_corr.T)

            # Compute EE Matrix
            #####################################

            EE = Matrix([[px],
                         [py],
                         [pz]])

            # Compute Wrist Center
            #####################################

            wc = EE - d4_7 * Rrpy[:,2]

            # Theta 1
            #####################################
            theta1 = atan2(wc[1],wc[0])


            # Calculate a, b, c 
            #####################################

            a = a2_3
            c = sqrt(pow((sqrt(wc[0] * wc[0]+wc[1] * wc[1]) - 0.35),2)+pow((wc[2] - 0.75), 2))
            b = sqrt(d3_4**2+a3_4**2)

            # Take the inverse cosine to get the angle (phi)
            phi1 = acos((c * c + a * a -b * b) / (2 * c * a))
            phi2 = acos((b * b + a * a -c * c) / (2 * b * a))
            phi3 = acos((c * c + b * b -a * a) / (2 * c * b))

            # Theta 2
            #####################################
            
            theta2 = pi/2 - phi1 - atan2((wc[2]-0.75),(sqrt(wc[0] * wc[0]+wc[1] * wc[1]) - 0.35))
            

            # Theta 3
            #####################################
            
            theta3 = pi/2 - phi2 - 0.036

            # Calculate rotational matrices
            R0_3_solved = R0_3.evalf(subs={q1: theta1, q2: theta2, q3: theta3})
            R3_6 = R0_3_solved.T * Rrpy

            # Use the previous information to acquire the remaining angles (theta)
            theta4 = atan2(R3_6[2,2], -R3_6[0,2])
            theta5 = atan2(sqrt(R3_6[0,2]*R3_6[0,2] + R3_6[2,2]*R3_6[2,2]),R3_6[1,2])
            theta6 = atan2(-R3_6[1,1],R3_6[1,0])

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
