//$8    16.01.2002    Matthias D.  Changed structure for MEG-Sensor
//$7    28.02.2001  Matthias D.  Added method to remove Electrodes
//$6    21.01.2001  Matthias D.  Removed Capital letters from file names
//$5    16.11.2000  Matthias D.  Added Constructor containing dataHandler and subsequent functionality
//$4    23.10.2000  Matthias D.  Added std:: to declaration of string and vector variables
//$3    12.10.2000  Matthias D.  Added methods to adress member variables for General Description and EEG Electrodes, set member variables for General description and Electordes as protected
//$2    20.09.2000  Matthias D.  changed Electrode Positions to double
//$1    31.08.2000  Matthias D.  created

#ifndef __IP_SIM_SENSORCONFIGURATION_C_H__
#define __IP_SIM_SENSORCONFIGURATION_C_H__

#include "CON_Vector_t.h"
#include "CON_Matrix_t.h"
#include "CON_Block_t.h"

#include "AnalysisDef.h" // added by Steinstraeter
#include "IP_SIG_AbstractDatahandler_c.h"

//
///////////////////////////////////////////////////////////////////////////////
//
//    NeuroFEM license:
//    =================
//    Copyright (c) 2007 by
//    Dr.Carsten Wolters, Dr.Alfred Anwander, Dr.Matthias Duempelmann,
//    Dr.Thomas Knoesche.
//
//    Permission is hereby granted, free of charge, to any person obtaining
//    a copy of the NeuroFEM software and associated documentation files
//    (the "Software"), to deal in the Software without restrictions, including
//    without limitations the rights to use, copy, modify, merge, publish,
//    distribute, sublicense, and/or sell copies of the Software, and to permit
//    persons to whom the Software is furnished to do so, subject to the following
//    conditions:
//
//    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
//    THE SOFTWARE. THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE
//    SOFTWARE IS WITH YOU.
//
//    The above copyright notice, permission and warranty notices shall be
//    included in all copies of the Software. A copyright notice of the
//    contributing authors and the above permission and warranty notices
//    shall be included in all source code files of the Software.
//
///////////////////////////////////////////////////////////////////////////////
// IP_SIM_SensorConfiguration_c.h: interface for the IP_SIM_SensorConfiguration_c class.
//
//////////////////////////////////////////////////////////////////////

enum sensortype_e {sensortype_unknown, sensortype_EEG, sensortype_MEG, sensortype_EEGMEG};


// Identification

struct SensorConfiguraton_Identification
    {
    std::string            SubjectName;               // Subject (Name of patient or volunteer)
    std::string            Study_Date_Time;           // Date and time of EEG/MEG Study
    };

// EEG

class ANALYSIS_EXPORT EEG_Electrode_Description_c
{
public:
      EEG_Electrode_Description_c();
      ~EEG_Electrode_Description_c();

      void   setElectrodeLabel(std::string inLabel);    // Set Channel Label (Fp1, Fp2,...)
      std::string getElectrodeLabel();                  // Get Channel Label (Fp1, Fp2,...)

      bool   setElectrodePosition(CON_Vector_t<double> inPosition);   // Set Electrode Position x,y,z in mm
      CON_Vector_t<double>& getElectrodePosition();                    // Get Electrode Position x,y,z in mm

protected:
std::string             m_Label;                     // Channel Label (Fp1, Fp2,...)
CON_Vector_t<double>   m_ElectrodePosition;         // Electrode Positions x,y,z in mm
};


class ANALYSIS_EXPORT anEEGSensorConfiguration_c
{
public:
    anEEGSensorConfiguration_c();
    anEEGSensorConfiguration_c(IP_SIG_AbstractDatahandler_c& inDataHandler);  // Constructor using data handler to set sensor positions
    virtual ~anEEGSensorConfiguration_c();

    // DataHandler

    bool                     setDataHandler(IP_SIG_AbstractDatahandler_c& inDatahandler);
    bool                     setPositionsbyDataHandler();

    // General
    bool setSensorConfiguraton_Identification(std::string inSubject,std::string inStudy_Date_Time);   //      Identification of SensorConfiguration
    std::string getSubject();
    std::string getStudy_Date_Time();

    bool setSensorType(sensortype_e inSensorType);              //   Set Type of sensor EEG, MEG or EEG+MEG
    sensortype_e getSensorType();                                //   Get Type of sensor EEG, MEG or EEG+MEG

    int  getNumberEEGElectrodes();                                                 // Get the number of EEG-Electrodes
    bool addEEGElectrode(std::string inLabel,CON_Vector_t<double> inPosition);       // Add EEG Electrode
    bool removeEEGElectrode(int Index);                                               // Remove EEG Electrode at Index Position

    void   changeElectrodeLabel(int Index,std::string inLabel);                    // Set Channel Label (Fp1, Fp2,...)
    std::string getElectrodeLabel(int Index);                                      // Get Channel Label (Fp1, Fp2,...)

    CON_Vector_t<double> getEEGElectrodePosition(int Index);                         // Get EEG Electrode Position
    bool changeEEGElectrodePosition(int Index, CON_Vector_t<double> inPosition);     // Change the Position of an EEG-Electrode

    void   setEEGReferenceLabel(std::string inLabel);                                   // Set Label of Reference Elecrode
    std::string getEEGReferenceLabel();                                                   // Get Label of Reference Elecrode

protected:
IP_SIG_AbstractDatahandler_c*                     m_pDataHandler;        //  Data Handler, can be used to get access to Sensorpositions

int                                             m_nEEGElectrodes;            //   Number of EEG Electrodes
std::vector<EEG_Electrode_Description_c>     m_EEGElectrodes;            //   Electrode labels and positions
std::string                                  m_RefLabel;                //   Label of reference electrode
sensortype_e                                 m_SensorType;                //   Type of sensor EEG, MEG or EEG+MEG
SensorConfiguraton_Identification             m_SenSor_Ident;            //   Identification of SensorConfiguration
};

// MEG
struct anMEGSensor_s
{
   // data
   CON_Matrix_t<double> positions;             // 3D positions of sensors
   CON_Matrix_t<double> directions;            // 3D directions of sensors
   CON_Vector_t<double> coil_positions;        // positions of coils
   CON_Vector_t<double> coil_areas;            // area of coils
   CON_Vector_t<int>    coil_windings;         // windings of coils
   CON_Matrix_t<double> integration_points;    // 2D positions of integration points
   CON_Vector_t<double> integration_weights;   // weights of integrations points

   // check consistency of sensor
   bool checkConsistency ()
   {
      if (positions.length  () != directions.length ()) return false;
      if (positions.height  () != 3)                    return false;
      if (directions.height () != 3)                    return false;

      if (coil_positions.length () != coil_areas.length ()) return false;
      if (coil_areas.length () != coil_windings.length ())  return false;

      if (integration_points.height () != integration_weights.length ()) return false;
      if (integration_points.length () != 2)                             return false;

      return true;
   };

   // Transfer local 2D integration coordinates to global 3D coordinates
   CON_Block_t<double> computeIntegrationPoints(int Chan);
   // Transform result from Weber to Tesla
   double WeberToTesla();
};

class ANALYSIS_EXPORT anMEGSensorConfiguration_c
{
 public:
 anMEGSensorConfiguration_c();
 bool setMatrices(CON_Matrix_t<double>&  ChannelCoilMatrix,  CON_Matrix_t<double>& correctionMatrix, CON_Matrix_t<double>& CoilPositions, CON_Matrix_t<double>& CoilNormal, CON_Vector_t<double>& CoilAreas, CON_Vector_t<int>& CoilWindings, CON_Matrix_t<int>& Polygons,CON_Vector_t<int>& numberintegrationPoints ,CON_Matrix_t<double>& integrationPoints, CON_Matrix_t<double>& integrationWeights, CON_Vector_t<double>& CircleRadius, CON_Matrix_t<double>& CircleNormals,CON_Vector_t<double> weber2tesla
 ,CON_Matrix_t<double>& PolygonsPoints,CON_Vector_t<int>& PolygonCoilMatch);
 //virtual ~anMEGSensorConfiguration_c();

// POINT Integration
            bool  computeIntegrationPointsNonWeighted(int coil_number, CON_Matrix_t<double>& Positions, CON_Matrix_t<double>& Normals,  CON_Vector_t<double>& weights) ;
            //  give back all needed information to make point integration
            bool computeIntegrationPointsWeighted(int channel_number, CON_Matrix_t<double>& Positions, CON_Matrix_t<double>& Normals,  CON_Vector_t<double>& weights);
           // like above, integration weight is weighted by estimated noise correction value

 // LINE Integration
           bool computeIntegrationPolygonNonWeighted(int channel_number,  CON_Block_t<double>& PolygonsAtChannel,  CON_Matrix_t<double>& weights);
          //Block PolygonsAtChannel is constructed like this:  [coil][PolgonPiece][Polygon-Coordinates]
          //Matrix weights= [coil][weight]

           bool computeIntegrationPolygonWeighted(int channel_number,  CON_Block_t<double>& PolygonsAtChannel,
                      CON_Block_t<double>& NormalAtChannel, CON_Matrix_t<double>& weights);
          // same structure like above, weights are weighted by estimated noise correction value

           bool getCoilCirclePickup(int coil_number, CON_Vector_t<double>& Position);
           bool getCoilCircleNormal(int coil_number, CON_Vector_t<double>& Normal);
           double getCoilCircleRadius(int coil_number);
     CON_Vector_t<double> getPositionOfCoil(int number);
     CON_Matrix_t<double> getPositionMatrix();
     CON_Matrix_t<double> getCorrectionMatrix();
         int getNumberCoilsPerChannel(int channel,CON_Vector_t<int>& numbers);
     int getNumberIntegrationPointsCoil(int Coilnumber);
     double getIntegrationPointWeights(int Coilnumber,int IntegrationPoint);
     double getCoilWindings(int Coilnumber);
     double getCoilArea(int Coilnumber);
     int getNumberOfChannels();
     int getNumberOfCoils();
         //int getNumberOfSensors();
         bool checkConsistency();
         bool isActivePointIntegration();
         bool isActiveCircleIntegration();
         bool isActivePolygonIntegration();
     double WeberToTesla(long coil);

protected:

struct anCoil_s
{
  CON_Vector_t<double> PositionOfCoil;  //next 3 arrays contain coil-information
  CON_Vector_t<double> Normal;

 double coil_area;            // area of coil
 int coil_winding;         // winding of coil

 bool point_integration;
 CON_Matrix_t<double> integration_points;    // 2D positions of integration points
 CON_Vector_t<double> integration_weights;   // weights of integration points

 bool polygon_integration;  //switch to polygon integration
 CON_Matrix_t<int> polygons; //polygon integration
 CON_Matrix_t<double> polygons_points;

  //Switch to circle integration mode?
  bool circle_integration;  // for circle integration
  double circle_radius;

   // check consistency of coil
   bool checkConsistency ()
   {
      if  ((point_integration && polygon_integration)  || (point_integration && circle_integration ) || (polygon_integration && circle_integration)) return false;

      if (polygons.height () != integration_weights.length ()) return false;

      if (integration_points.height () != integration_weights.length ()) return false;

     // if ((xDirectionsOfCoil.length () != 3) ||  (NormalOfCoil.length()!=3)  ||  (PositionOfCoil.length()!=3) ) return false;

      return true;
   };
};

CON_Vector_t<double> m_weber2tesla;
anCoil_s *m_Coils; //contain all coils
CON_Matrix_t<double> m_Channel; //sensor = combination of coils
CON_Matrix_t<double> m_ChannelSensor; // for fem, CTF

CON_Matrix_t<double> m_weightedChannel;
CON_Matrix_t<double> m_correctionMatrix; // noise correction
bool m_PointIntegration;
bool m_CircleIntegration;
bool m_PolygonIntegration;
int m_NumberofChannels;
};

#endif // __anSensorConfiguration_c_H__
