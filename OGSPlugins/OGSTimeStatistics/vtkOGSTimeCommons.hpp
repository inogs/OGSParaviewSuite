/*=========================================================================

  Program:   OGSTimeCommons
  Module:    vtkOGSTimeCommons.hpp

  Copyright (c) 2020 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkOGSTimeCommons_h
#define vtkOGSTimeCommons_h

void BuildTimeList(Time::TimeList &TL, vtkInformation *Info);
void RecoverMasterFileName(std::string &fname, vtkDataSet *input);
void strsplit(std::string str, std::string splitBy, std::vector<std::string> &tokens);

#endif