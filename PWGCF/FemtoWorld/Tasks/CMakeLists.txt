# Copyright 2019-2020 CERN and copyright holders of ALICE O2.
# See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
# All rights not expressly granted are reserved.
#
# This software is distributed under the terms of the GNU General Public
# License v3 (GPL Version 3), copied verbatim in the file "COPYING".
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization
# or submit itself to any jurisdiction.

#o2physics_add_dpl_workflow(femto-world-debug-track
#          SOURCES femtoWorldDebugTrack.cxx
#          PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
#          COMPONENT_NAME Analysis)

#o2physics_add_dpl_workflow(femto-world-debug-v0
#          SOURCES femtoWorldDebugV0.cxx
#          PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
#          COMPONENT_NAME Analysis)

#o2physics_add_dpl_workflow(femto-world-hash
#          SOURCES femtoWorldHashTask.cxx
#          PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
#          COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(femto-world-pair-track-track
          SOURCES femtoWorldPairTaskTrackTrack.cxx
          PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
          COMPONENT_NAME Analysis)

#o2physics_add_dpl_workflow(femto-world-pair-track-v0
#          SOURCES femtoWorldPairTaskTrackV0.cxx
#          PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
#          COMPONENT_NAME Analysis)

