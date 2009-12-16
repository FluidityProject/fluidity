/*  Copyright (C) 2006 Imperial College London and others.
    
    Please see the AUTHORS file in the main source directory for a full list
    of copyright holders.

    Prof. C Pain
    Applied Modelling and Computation Group
    Department of Earth Science and Engineering
    Imperial College London

    C.Pain@Imperial.ac.uk
    
    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation,
    version 2.1 of the License.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
    USA
*/

#include "FLComms_IO.h"
#include "Unittest_Tools_Cpp.h"

using namespace Fluidity;
using namespace std;

extern "C"{
#define test_flcomms_io_fc F77_FUNC(test_flcomms_io, TEST_FLCOMMS_IO)
  void test_flcomms_io_fc(){
    // Clear the existing output file
    halo_manager.Reset();
    halo_manager.SetNProcs(1);
        
    report_test("[Write halos]", WriteHalos(string("data/test_flcomms_io_out"), false), false, "Failed to write halo");
    halo_manager.Reset();
    report_test("[Read halos]", ReadHalos(string("data/test_flcomms_io_out")), false, "Failed to read halo");
    
    ostringstream buffer;
    
    set<int> tagsOut = halo_manager.GetTags();
    buffer << tagsOut.size();
    report_test("[Correct number of tags read]", tagsOut.size() != 0, false, buffer.str() + " tags found");
    buffer.str("");
  
    // Add some data to the halo
    halo_manager.Reset();
    halo_manager.SetNProcs(1);
    
    // Test with ordered list and out of order list, as halo node order matters    
    vector<vector<int> > sendIn, recvIn;
    
    sendIn.push_back(vector<int>());
    sendIn[0].push_back(2);
    sendIn[0].push_back(3);
    sendIn[0].push_back(5);
    sendIn[0].push_back(7);
    sendIn[0].push_back(11);
    sendIn[0].push_back(13);
    sendIn[0].push_back(17);
    sendIn[0].push_back(19);
    sendIn[0].push_back(23);
    sendIn[0].push_back(29);
    
    recvIn.push_back(vector<int>());
    recvIn[0].push_back(2);
    recvIn[0].push_back(1);
    recvIn[0].push_back(4);
    recvIn[0].push_back(3);
    recvIn[0].push_back(6);
    recvIn[0].push_back(5);
    recvIn[0].push_back(8);
    recvIn[0].push_back(7);
    recvIn[0].push_back(10);
    recvIn[0].push_back(9);
    
    set<int> tagsIn;
    tagsIn.insert(42);
    int npnodesIn = 666;
    report_test("[Register halo]", halo_manager.RegisterHalo(*tagsIn.begin(), npnodesIn, sendIn, recvIn), false, "Failed to register halo");
    
    report_test("[Write halos]", WriteHalos(string("data/test_flcomms_io_out"), false), false, "Failed to write halo");
    halo_manager.Reset();
    report_test("[Read halos]", ReadHalos(string("data/test_flcomms_io_out")), false, "Failed to read halo");
    
    tagsOut = halo_manager.GetTags();
    buffer << tagsOut.size();
    report_test("[Correct number of tags read]", tagsOut.size() != tagsIn.size(), false, buffer.str() + " tags found");
    buffer.str("");

    bool fail;
    if(tagsOut.size() != tagsIn.size()){
      fail = true;
    }else{
      fail = false;
      for(set<int>::const_iterator iter = tagsOut.begin();iter != tagsOut.end();iter++){
        if(tagsIn.count(*iter) == 0){
          fail = true;
          break;
        }
      }
    }
    report_test("[Correct tags read]", fail, false, "Incorrect tags");
    
    int npnodesOut = halo_manager.GetNOwnedNodes(*tagsIn.begin());
    buffer << npnodesOut;
    report_test("[Correct number of private nodes read]", npnodesOut != npnodesIn, false, buffer.str() + " found");
    buffer.str("");
    
    vector<vector<int> > sendOut, recvOut;
    report_test("[Export halo]", halo_manager.ExportHalo(*tagsIn.begin(), sendOut, recvOut), false, "Failed to export halo");

    buffer << sendOut.size();
    report_test("[Correct number of processes with sends read]", sendOut.size() != sendIn.size(), false, buffer.str() + " found");
    buffer.str("");
    buffer << recvOut.size();
    report_test("[Correct number of processes with receives read]", recvOut.size() != recvIn.size(), false, buffer.str() + " found");
    buffer.str("");

    buffer << sendOut[0].size();
    report_test("[Correct number of sends read]", sendOut[0].size() != sendIn[0].size(), false, buffer.str() + " found");
    buffer.str("");
    buffer << sendOut[0].size();
    report_test("[Correct number of receives read]", recvOut[0].size() != recvIn[0].size(), false, buffer.str() + " found");
    buffer.str("");

    if(sendOut[0].size() != sendIn[0].size()){
      fail = true;
    }else{
      fail = false;
      for(size_t i = 0;i < sendOut[0].size();i++){
        if(sendOut[0][i] != sendIn[0][i]){
          fail = true;
          break;
        }
      }
    }
    report_test("[Correct sends read]", fail, false, "Incorrect sends read");
    if(recvOut[0].size() != recvIn[0].size()){
      fail = true;
    }else{
      fail = false;
      for(size_t i = 0;i < recvOut[0].size();i++){
        if(recvOut[0][i] != recvIn[0][i]){
          fail = true;
          break;
        }
      }
    }
    report_test("[Correct receives read]", fail, false, "Incorrect receives read");
    
    return;
  }
}
