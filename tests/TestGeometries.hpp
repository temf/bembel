// This file is part of Bembel, the higher order C++ boundary element library.
//
// Copyright (C) 2022 see <http://www.bembel.eu>
//
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
//
#ifndef TESTS_TESTGEOMETRIES_HPP_
#define TESTS_TESTGEOMETRIES_HPP_

#include "tests/Test.hpp"

namespace Test {

struct TestGeometryWriter {
  static void writeEdgeCase1() {
    std::ofstream test_file;
    test_file.open("test_EdgeCase1.dat");
    test_file
        << "# nurbs mesh v.2.1\n"
           "# test_EdgeCase1.dat\n"
           "# Generated by BEMBEL, see www.bembel.eu\n"
           "#\n"
           "2 3 5 0 0 \n"
           "PATCH 0 \n"
           "1 1 \n"
           "2 2 \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "0.000000000000000   1.000000000000000   0.000000000000000   "
           "1.000000000000000   \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "0.000000000000000   0.000000000000000   0.000000000000000   "
           "0.000000000000000   \n"
           "1.000000000000000   1.000000000000000   1.000000000000000   "
           "1.000000000000000  \n"
           "PATCH 1 \n"
           "1 1 \n"
           "2 2 \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "1.000000000000000   2.000000000000000   1.000000000000000   "
           "2.000000000000000   \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "0.000000000000000   0.000000000000000   0.000000000000000   "
           "0.000000000000000   \n"
           "1.000000000000000   1.000000000000000   1.000000000000000   "
           "1.000000000000000  \n"
           "PATCH 2 \n"
           "1 1 \n"
           "2 2 \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "-1.000000000000000   0.000000000000000   -1.000000000000000   "
           "0.000000000000000   \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "0.000000000000000   0.000000000000000   0.000000000000000   "
           "0.000000000000000   \n"
           "1.000000000000000   1.000000000000000   1.000000000000000   "
           "1.000000000000000  \n"
           "PATCH 3 \n"
           "1 1 \n"
           "2 2 \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "0.000000000000000   1.000000000000000   0.000000000000000   "
           "1.000000000000000   \n"
           "-1.000000000000000   -1.000000000000000   0.000000000000000   "
           "0.000000000000000   \n"
           "0.000000000000000   0.000000000000000   0.000000000000000   "
           "0.000000000000000   \n"
           "1.000000000000000   1.000000000000000   1.000000000000000   "
           "1.000000000000000  \n"
           "PATCH 4 \n"
           "1 1 \n"
           "2 2 \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "0.000000000000000   1.000000000000000   0.000000000000000   "
           "1.000000000000000   \n"
           "1.000000000000000   1.000000000000000   2.000000000000000   "
           "2.000000000000000   \n"
           "0.000000000000000   0.000000000000000   0.000000000000000   "
           "0.000000000000000   \n"
           "1.000000000000000   1.000000000000000   1.000000000000000   "
           "1.000000000000000  \n";
    test_file.close();
  }
  static void writeEdgeCase2() {
    std::ofstream test_file;
    test_file.open("test_EdgeCase2.dat");
    test_file
        << "# nurbs mesh v.2.1\n"
           "# test_EdgeCase2.dat\n"
           "# Generated by BEMBEL, see www.bembel.eu\n"
           "#\n"
           "2 3 5 0 0 \n"
           "PATCH 0 \n"
           "1 1 \n"
           "2 2 \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "1.000000000000000   0.000000000000000   1.000000000000000   "
           "0.000000000000000   \n"
           "1.000000000000000   1.000000000000000   0.000000000000000   "
           "0.000000000000000   \n"
           "0.000000000000000   0.000000000000000   0.000000000000000   "
           "0.000000000000000   \n"
           "1.000000000000000   1.000000000000000   1.000000000000000   "
           "1.000000000000000  \n"
           "PATCH 1 \n"
           "1 1 \n"
           "2 2 \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "1.000000000000000   2.000000000000000   1.000000000000000   "
           "2.000000000000000   \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "0.000000000000000   0.000000000000000   0.000000000000000   "
           "0.000000000000000   \n"
           "1.000000000000000   1.000000000000000   1.000000000000000   "
           "1.000000000000000  \n"
           "PATCH 2 \n"
           "1 1 \n"
           "2 2 \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "-1.000000000000000   0.000000000000000   -1.000000000000000   "
           "0.000000000000000   \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "0.000000000000000   0.000000000000000   0.000000000000000   "
           "0.000000000000000   \n"
           "1.000000000000000   1.000000000000000   1.000000000000000   "
           "1.000000000000000  \n"
           "PATCH 3 \n"
           "1 1 \n"
           "2 2 \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "0.000000000000000   1.000000000000000   0.000000000000000   "
           "1.000000000000000   \n"
           "-1.000000000000000   -1.000000000000000   0.000000000000000   "
           "0.000000000000000   \n"
           "0.000000000000000   0.000000000000000   0.000000000000000   "
           "0.000000000000000   \n"
           "1.000000000000000   1.000000000000000   1.000000000000000   "
           "1.000000000000000  \n"
           "PATCH 4 \n"
           "1 1 \n"
           "2 2 \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "0.000000000000000   1.000000000000000   0.000000000000000   "
           "1.000000000000000   \n"
           "1.000000000000000   1.000000000000000   2.000000000000000   "
           "2.000000000000000   \n"
           "0.000000000000000   0.000000000000000   0.000000000000000   "
           "0.000000000000000   \n"
           "1.000000000000000   1.000000000000000   1.000000000000000   "
           "1.000000000000000  \n";
    test_file.close();
  }
  static void writeEdgeCase3() {
    std::ofstream test_file;
    test_file.open("test_EdgeCase3.dat");
    test_file
        << "# nurbs mesh v.2.1\n"
           "# test_EdgeCase3.dat\n"
           "# Generated by BEMBEL, see www.bembel.eu\n"
           "#\n"
           "2 3 5 0 0 \n"
           "PATCH 0 \n"
           "1 1 \n"
           "2 2 \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "1.000000000000000   1.000000000000000   0.000000000000000   "
           "0.000000000000000   \n"
           "0.000000000000000   1.000000000000000   0.000000000000000   "
           "1.000000000000000   \n"
           "0.000000000000000   0.000000000000000   0.000000000000000   "
           "0.000000000000000   \n"
           "1.000000000000000   1.000000000000000   1.000000000000000   "
           "1.000000000000000  \n"
           "PATCH 1 \n"
           "1 1 \n"
           "2 2 \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "1.000000000000000   2.000000000000000   1.000000000000000   "
           "2.000000000000000   \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "0.000000000000000   0.000000000000000   0.000000000000000   "
           "0.000000000000000   \n"
           "1.000000000000000   1.000000000000000   1.000000000000000   "
           "1.000000000000000  \n"
           "PATCH 2 \n"
           "1 1 \n"
           "2 2 \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "-1.000000000000000   0.000000000000000   -1.000000000000000   "
           "0.000000000000000   \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "0.000000000000000   0.000000000000000   0.000000000000000   "
           "0.000000000000000   \n"
           "1.000000000000000   1.000000000000000   1.000000000000000   "
           "1.000000000000000  \n"
           "PATCH 3 \n"
           "1 1 \n"
           "2 2 \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "0.000000000000000   1.000000000000000   0.000000000000000   "
           "1.000000000000000   \n"
           "-1.000000000000000   -1.000000000000000   0.000000000000000   "
           "0.000000000000000   \n"
           "0.000000000000000   0.000000000000000   0.000000000000000   "
           "0.000000000000000   \n"
           "1.000000000000000   1.000000000000000   1.000000000000000   "
           "1.000000000000000  \n"
           "PATCH 4 \n"
           "1 1 \n"
           "2 2 \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "0.000000000000000   1.000000000000000   0.000000000000000   "
           "1.000000000000000   \n"
           "1.000000000000000   1.000000000000000   2.000000000000000   "
           "2.000000000000000   \n"
           "0.000000000000000   0.000000000000000   0.000000000000000   "
           "0.000000000000000   \n"
           "1.000000000000000   1.000000000000000   1.000000000000000   "
           "1.000000000000000  \n";
    test_file.close();
  }
  static void writeEdgeCase4() {
    std::ofstream test_file;
    test_file.open("test_EdgeCase4.dat");
    test_file
        << "# nurbs mesh v.2.1\n"
           "# test_EdgeCase4.dat\n"
           "# Generated by BEMBEL, see www.bembel.eu\n"
           "#\n"
           "2 3 5 0 0 \n"
           "PATCH 0 \n"
           "1 1 \n"
           "2 2 \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "1.000000000000000   0.000000000000000   1.000000000000000   "
           "0.000000000000000   \n"
           "0.000000000000000   0.000000000000000   0.000000000000000   "
           "0.000000000000000   \n"
           "1.000000000000000   1.000000000000000   1.000000000000000   "
           "1.000000000000000  \n"
           "PATCH 1 \n"
           "1 1 \n"
           "2 2 \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "1.000000000000000   2.000000000000000   1.000000000000000   "
           "2.000000000000000   \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "0.000000000000000   0.000000000000000   0.000000000000000   "
           "0.000000000000000   \n"
           "1.000000000000000   1.000000000000000   1.000000000000000   "
           "1.000000000000000  \n"
           "PATCH 2 \n"
           "1 1 \n"
           "2 2 \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "-1.000000000000000   0.000000000000000   -1.000000000000000   "
           "0.000000000000000   \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "0.000000000000000   0.000000000000000   0.000000000000000   "
           "0.000000000000000   \n"
           "1.000000000000000   1.000000000000000   1.000000000000000   "
           "1.000000000000000  \n"
           "PATCH 3 \n"
           "1 1 \n"
           "2 2 \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "0.000000000000000   1.000000000000000   0.000000000000000   "
           "1.000000000000000   \n"
           "-1.000000000000000   -1.000000000000000   0.000000000000000   "
           "0.000000000000000   \n"
           "0.000000000000000   0.000000000000000   0.000000000000000   "
           "0.000000000000000   \n"
           "1.000000000000000   1.000000000000000   1.000000000000000   "
           "1.000000000000000  \n"
           "PATCH 4 \n"
           "1 1 \n"
           "2 2 \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "0.000000000000000   1.000000000000000   0.000000000000000   "
           "1.000000000000000   \n"
           "1.000000000000000   1.000000000000000   2.000000000000000   "
           "2.000000000000000   \n"
           "0.000000000000000   0.000000000000000   0.000000000000000   "
           "0.000000000000000   \n"
           "1.000000000000000   1.000000000000000   1.000000000000000   "
           "1.000000000000000  \n";
    test_file.close();
  }
  static void writeSpherePanel() {
    std::ofstream test_file;
    test_file.open("test_SpherePanel.dat");
    test_file
        << "# nurbs mesh v.2.1\n"
           "# test_Sphere.dat\n"
           "# Generated by BEMBEL, see www.bembel.eu\n"
           "#\n"
           "2 3 1 0 0 \n"
           "PATCH 0 \n"
           "4 4 \n"
           "5 5 \n"
           "0.000000000000000   0.000000000000000   0.000000000000000   "
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   1.000000000000000   1.000000000000000   "
           "1.000000000000000  \n"
           "0.000000000000000   0.000000000000000   0.000000000000000   "
           "0.000000000000000   0.000000000000000   1.000000000000000   "
           "1.000000000000000   1.000000000000000   1.000000000000000   "
           "1.000000000000000   \n"
           "0.577350269189626   0.278838767912603   0.000000000000000   "
           "-0.278838767912603   -0.577350269189626   0.632392158505876   "
           "0.315090742770461   0.000000000000000   -0.315090742770461   "
           "-0.632392158505876   0.647791890991355   0.328648516366383   "
           "0.000000000000000   -0.328648516366383   -0.647791890991355   "
           "0.632392158505876   0.315090742770461   0.000000000000000   "
           "-0.315090742770461   -0.632392158505876   0.577350269189626   "
           "0.278838767912603   0.000000000000000   -0.278838767912603   "
           "-0.577350269189626   \n"
           "0.577350269189626   0.632392158505876   0.647791890991355   "
           "0.632392158505876   0.577350269189626   0.632392158505876   "
           "0.762259526419164   0.804938188574224   0.762259526419164   "
           "0.632392158505876   0.647791890991355   0.804938188574224   "
           "0.859116756396542   0.804938188574224   0.647791890991355   "
           "0.632392158505876   0.762259526419164   0.804938188574224   "
           "0.762259526419164   0.632392158505876   0.577350269189626   "
           "0.632392158505876   0.647791890991355   0.632392158505876   "
           "0.577350269189626   \n"
           "-0.577350269189626   -0.632392158505876   -0.647791890991355   "
           "-0.632392158505876   -0.577350269189626   -0.278838767912603   "
           "-0.315090742770461   -0.328648516366383   -0.315090742770461   "
           "-0.278838767912603   -0.000000000000000   -0.000000000000000   "
           "-0.000000000000000   -0.000000000000000   -0.000000000000000   "
           "0.278838767912603   0.315090742770461   0.328648516366383   "
           "0.315090742770461   0.278838767912603   0.577350269189626   "
           "0.632392158505876   0.647791890991355   0.632392158505876   "
           "0.577350269189626  \n "
           "1.000000000000000   0.891211203608397   0.859116756396542   "
           "0.891211203608397   1.000000000000000   0.891211203608397   "
           "0.762259526419164   0.718665173540050   0.762259526419164   "
           "0.891211203608397   0.859116756396542   0.718665173540050   "
           "0.671272431591931   0.718665173540050   0.859116756396542   "
           "0.891211203608397   0.762259526419164   0.718665173540050   "
           "0.762259526419164   0.891211203608397   1.000000000000000   "
           "0.891211203608397   0.859116756396542   0.891211203608397   "
           "1.000000000000000   ";
    test_file.close();
  }

  static void writeScaledSpherePanel() {
    std::ofstream test_file;
    test_file.open("test_ScaledSpherePanel.dat");
    test_file << "# nurbs mesh v.2.1\n"
                 "# test_ScaledSphere.dat\n"
                 "# Generated by BEMBEL, see www.bembel.eu\n"
                 "#\n"
                 "2 3 1 0 0 \n"
                 "PATCH 0 \n"
                 "4 4 \n"
                 "5 5 \n"
                 "0.000000000000000   0.000000000000000   0.000000000000000   "
                 "0.000000000000000   0.000000000000000   1.000000000000000   "
                 "1.000000000000000   1.000000000000000   1.000000000000000   "
                 "1.000000000000000  \n"
                 "0.000000000000000   0.000000000000000   0.000000000000000   "
                 "0.000000000000000   0.000000000000000   1.000000000000000   "
                 "1.000000000000000   1.000000000000000   1.000000000000000   "
                 "1.000000000000000   \n"
                 "5.77350269189626   2.78838767912603   0.00000000000000   "
                 "-2.78838767912603   -5.77350269189626   6.32392158505876   "
                 "3.15090742770461   0.00000000000000   -3.15090742770461   "
                 "-6.32392158505876   6.47791890991355   3.28648516366383   "
                 "0.00000000000000   -3.28648516366383   -6.47791890991355   "
                 "6.32392158505876   3.15090742770461   0.00000000000000   "
                 "-3.15090742770461   -6.32392158505876   5.77350269189626   "
                 "2.78838767912603   0.00000000000000   -2.78838767912603   "
                 "-5.77350269189626   \n"
                 "5.77350269189626   6.32392158505876   6.47791890991355   "
                 "6.32392158505876   5.77350269189626   6.32392158505876   "
                 "7.62259526419164   8.04938188574224   7.62259526419164   "
                 "6.32392158505876   6.47791890991355   8.04938188574224   "
                 "8.59116756396542   8.04938188574224   6.47791890991355   "
                 "6.32392158505876   7.62259526419164   8.04938188574224   "
                 "7.62259526419164   6.32392158505876   5.77350269189626   "
                 "6.32392158505876   6.47791890991355   6.32392158505876   "
                 "5.77350269189626   \n"
                 "-5.77350269189626   -6.32392158505876   -6.47791890991355   "
                 "-6.32392158505876   -5.77350269189626   -2.78838767912603   "
                 "-3.15090742770461   -3.28648516366383   -3.15090742770461   "
                 "-2.78838767912603   -0.00000000000000   -0.00000000000000   "
                 "-0.00000000000000   -0.00000000000000   -0.00000000000000   "
                 "2.78838767912603   3.15090742770461   3.28648516366383   "
                 "3.15090742770461   2.78838767912603   5.77350269189626   "
                 "6.32392158505876   6.47791890991355   6.32392158505876   "
                 "5.77350269189626  \n "
                 "10.00000000000000   8.91211203608397   8.59116756396542   "
                 "8.91211203608397   10.00000000000000   8.91211203608397   "
                 "7.62259526419164   7.18665173540050   7.62259526419164   "
                 "8.91211203608397   8.59116756396542   7.18665173540050   "
                 "6.71272431591931   7.18665173540050   8.59116756396542   "
                 "8.91211203608397   7.62259526419164   7.18665173540050   "
                 "7.62259526419164   8.91211203608397   10.00000000000000   "
                 "8.91211203608397   8.59116756396542   8.91211203608397   "
                 "10.00000000000000   ";

    test_file.close();
  }
  static void writeScreen() {
    std::ofstream test_file;
    test_file.open("test_Screen.dat");
    test_file << "# nurbs mesh v.2.1\n"
                 "# test_Screen.dat\n"
                 "# Generated by BEMBEL, see www.bembel.eu\n"
                 "#\n"
                 "2 3 1 0 0 \n"
                 "PATCH 0 \n"
                 "1 1 \n"
                 "2 2 \n"
                 "0.000000000000000   0.000000000000000   1.000000000000000   "
                 "1.000000000000000 \n"
                 "0.000000000000000   0.000000000000000   1.000000000000000   "
                 "1.000000000000000 \n"
                 "0.000000000000000   1.000000000000000   0.000000000000000   "
                 "1.000000000000000 \n"
                 "0.000000000000000   0.000000000000000   1.000000000000000   "
                 "1.000000000000000 \n"
                 "0.000000000000000   0.000000000000000   0.000000000000000   "
                 "0.000000000000000 \n"
                 "1.000000000000000   1.000000000000000   1.000000000000000   "
                 "1.000000000000000 \n";
    test_file.close();
  }
};

}  // namespace Test
#endif  // TESTS_TESTGEOMETRIES_HPP_
