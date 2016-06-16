//$13    10.04.2008  Dannhauer, remake of the class
//$12    01.07.2003  Anwander A.  Remoced bug in removeEEGElectrode() for Linux
//$11    19.03.2002  Matthias D.  Removed bugs in removeEEGElectrode() for Linux
//$10    19.03.2002  Matthias D.  Adapted removeEEGElectrode() for Linux
//$9    22.01.2002  Matthias D.  Removed bugs from anMEGSensor_s::computeIntegrationPoints() and anMEGSensor_s::WeberToTesla()
//$8    16.01.2002    Matthias D.  Changed structure for MEG-Sensor
//$7    28.02.2001  Matthias D.  Added Method to remove Electrodes
//$6    16.11.2000  Matthias D.  Added Constructor containing dataHandler and subsequent functionality
//$5    23.10.2000  Matthias D.  Added std:: to declaration of string and vector variables
//$4    12.10.2000  Matthias D.  Added methods to adress member variables for General Description and EEG Electrodes
//$3    22.09.2000  Matthias D.  Changed from \ to / for directory names
//$2    20.09.2000  Matthias D.  changed Electrode Positions to double
//$1    31.08.2000  Matthias D.  created

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
// IP_SIM_SensorConfiguration_c.cpp: implementation of the IP_SIM_SensorConfiguration_c class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAnalysisPch.h"

using namespace std;

#include "LIN_Auxil.h"

#include "IP_SIM_SensorConfiguration_c.h"

#include <math.h>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

anEEGSensorConfiguration_c::anEEGSensorConfiguration_c():
m_pDataHandler(NULL),
m_nEEGElectrodes(0)
{
   m_SenSor_Ident.SubjectName=" ";
   m_SenSor_Ident.Study_Date_Time = "00.00.0000";
   m_SensorType = sensortype_unknown;
}



anEEGSensorConfiguration_c::anEEGSensorConfiguration_c(IP_SIG_AbstractDatahandler_c& inDataHandler):
m_pDataHandler(&inDataHandler),
m_nEEGElectrodes(0)
{
 m_SenSor_Ident.SubjectName=" ";
 m_SenSor_Ident.Study_Date_Time = "00.00.0000";
 m_SensorType = sensortype_unknown;

 setPositionsbyDataHandler(); // Set Sensor Positions
}


anEEGSensorConfiguration_c::~anEEGSensorConfiguration_c()
{
 m_pDataHandler = NULL;
}

// Set Datahandler

bool anEEGSensorConfiguration_c::setDataHandler(IP_SIG_AbstractDatahandler_c& inDatahandler)
{
  m_pDataHandler= &inDatahandler;
  return true;
}

bool anEEGSensorConfiguration_c::setPositionsbyDataHandler()
{
  CON_Matrix_t<float>  SensorPositions, Directions;
  CON_Vector_t<bool>   EnabledSensor;
  std::string        Label=" ";
  CON_Vector_t<double> Pos(3);
  int i,k;

  m_pDataHandler->GetSignalSources(SensorPositions,Directions,EnabledSensor);

  for (i=0; i < SensorPositions.height(); i++)
    {
      if(EnabledSensor[i])
        {
         for(k=0;k<3;k++) Pos[k] = SensorPositions[i][k];
         addEEGElectrode(Label,Pos);
        }
    }
  return true;
}

// Members for General Properties

bool anEEGSensorConfiguration_c::setSensorConfiguraton_Identification(std::string inSubject,std::string inStudy_Date_Time)   //      Identification of SensorConfiguration
{
m_SenSor_Ident.SubjectName = inSubject;
m_SenSor_Ident.Study_Date_Time = inStudy_Date_Time;
return true;
}

std::string anEEGSensorConfiguration_c::getSubject()
{
return m_SenSor_Ident.SubjectName;
}

std::string anEEGSensorConfiguration_c::getStudy_Date_Time()
{
return m_SenSor_Ident.Study_Date_Time;
}

bool anEEGSensorConfiguration_c::setSensorType(sensortype_e inSensorType)            //   Set Type of sensor EEG, MEG or EEG+MEG
{
m_SensorType = inSensorType;
return true;
}

sensortype_e anEEGSensorConfiguration_c::getSensorType()                                //   Get Type of sensor EEG, MEG or EEG+MEG
{
return m_SensorType;
}

// Members for EEG Electrodes

int anEEGSensorConfiguration_c::getNumberEEGElectrodes()
{
return m_EEGElectrodes.size();
}

bool anEEGSensorConfiguration_c::addEEGElectrode(std::string inLabel,CON_Vector_t<double> inPosition)            // Add EEG Electrode
{
if (inPosition.length() != 3) return false;
m_EEGElectrodes.resize(m_nEEGElectrodes+1);
m_EEGElectrodes[m_nEEGElectrodes].setElectrodeLabel(inLabel);
m_EEGElectrodes[m_nEEGElectrodes].setElectrodePosition(inPosition);
m_nEEGElectrodes++;
return true;
}


bool anEEGSensorConfiguration_c::removeEEGElectrode(int Index)                                               // Remove EEG Electrode at Index Position
{
if((Index < 0)|| Index > (getNumberEEGElectrodes()-1)) return false;

#ifdef WIN32
m_EEGElectrodes.erase(&m_EEGElectrodes[Index]);
#else
int i;
for (i=Index;i<(int)(m_EEGElectrodes.size())-1;i++)
    m_EEGElectrodes[i] = m_EEGElectrodes[i+1];

m_nEEGElectrodes = m_EEGElectrodes.size()-1;
m_EEGElectrodes.resize(m_nEEGElectrodes);
#endif

return true;
}

void anEEGSensorConfiguration_c::changeElectrodeLabel(int Index, std::string inLabel)                         // Set Channel Label (Fp1, Fp2,...)
{
m_EEGElectrodes[Index].setElectrodeLabel(inLabel);
}

std::string anEEGSensorConfiguration_c::getElectrodeLabel(int Index)                                          // Get Channel Label (Fp1, Fp2,...)
{
if ((Index < 0)||(Index>m_nEEGElectrodes-1))
    {
     return "";
    }
return m_EEGElectrodes[Index].getElectrodeLabel();
}

CON_Vector_t<double> anEEGSensorConfiguration_c::getEEGElectrodePosition(int Index)                        // Get EEG Electrode Position
{
if ((Index < 0)||(Index>m_nEEGElectrodes-1))
    {CON_Vector_t<double> NULLPOS(3,0);
     return NULLPOS;
    }
return m_EEGElectrodes[Index].getElectrodePosition();
}

bool anEEGSensorConfiguration_c::changeEEGElectrodePosition(int Index,CON_Vector_t<double> inPosition)      // Change the Position of an EEG-Electrode
{
if ((Index < 0)||(Index>m_nEEGElectrodes-1)) return false;
m_EEGElectrodes[Index].setElectrodePosition(inPosition);
return true;
}

void anEEGSensorConfiguration_c::setEEGReferenceLabel(std::string inLabel)                                     // Set Label of Reference Elecrode
{
m_RefLabel = inLabel;
}

std::string anEEGSensorConfiguration_c::getEEGReferenceLabel()                                                  // Get Label of Reference Elecrode
{
return m_RefLabel;
}

// EEG electrode description

EEG_Electrode_Description_c::EEG_Electrode_Description_c()
{
}

EEG_Electrode_Description_c::~EEG_Electrode_Description_c()
{

}


void EEG_Electrode_Description_c::setElectrodeLabel(std::string inLabel)
{
 m_Label = inLabel;
}

std::string EEG_Electrode_Description_c::getElectrodeLabel()
{
 return m_Label;
}

bool EEG_Electrode_Description_c::setElectrodePosition(CON_Vector_t<double> inPosition)
{
if (inPosition.length() != 3) return false;
m_ElectrodePosition = inPosition;
return false;
}

CON_Vector_t<double>& EEG_Electrode_Description_c::getElectrodePosition()
{
    return m_ElectrodePosition;
}

//////////////////////////////////////////////////////////////////////////
// MEG
// Transfer local 2D integration coordinates to global 3D coordinates

CON_Block_t<double> anMEGSensor_s::computeIntegrationPoints(int Chan) // Transfer local 2D integration coordinates to global 3D coordinates
{
    CON_Block_t<double>    IntegrationPoints3D;    // Block containing return values
    double                xl,yl,zl;                // local coordinates of the current integration point
    CON_Vector_t<double>    Dir(3);
    CON_Vector_t<double>    v1(3);
    CON_Vector_t<double>    v2(3);
    double                P;
    double                lendir;
    int                    c,p,k;
    CON_Matrix_t<double>    AMatrix;
    CON_Vector_t<float>    xMatrix(4);
    CON_Matrix_t<double>    BMatrix;


      for(k=0;k<3;k++) Dir[k] = directions[k][Chan];
    lendir = Dir.norm();
    if ( lendir < 1e-30 )
        {
            return IntegrationPoints3D;
        }

    Dir[0] = Dir[0] / lendir;
    Dir[1] = Dir[1] / lendir;
    Dir[2] = Dir[2] / lendir;

    //  Form right-hand coordinate-system [v1,v2,v3], where..
    //  v1 and v2 make up the plane rectangular to dir
     if (((Dir[0] != 0) || (Dir[2] != 0)))
        {
         //  v1 |- to Dir[]
         v1[0] = Dir[2];v1[1] = 0.0;v1[2] = -Dir[0];
        }
        else   //  Dir[] equals (0,1,0)
        {
         v1[0] = Dir[1];v1[1] = 0.0;v1[2] = 0.0;
        }

     CON_Vector_t<double> Direction(3);
     CON_Vector_t<double> XDirection(3);

     for (k=0;k<3;k++)
        {
         Direction[k]= directions[0][Chan];
         XDirection[k]= directions[0][Chan];
        }

/*     if(Scalar(Direction,XDirection)!=0)
        //If directions and RedXDirections are not perpendicular, they should be made perpendicular
        {
         AMatrix.allocate(4,3);
         BMatrix.allocate(4,1);


         AMatrix[0][0]= directions[0][Chan];
         AMatrix[0][1]= directions[1][Chan];
         AMatrix[0][2]= directions[2][Chan];

            AMatrix[1][0]= 0;
            AMatrix[1][1]= directions[2][Chan];
            AMatrix[1][2]=-directions[1][Chan];

            AMatrix[2][0]=-directions[2][Chan];
            AMatrix[2][1]= 0;
            AMatrix[2][2]= directions[0][Chan];

            AMatrix[3][0]= directions[1][Chan];
            AMatrix[3][1]=-directions[0][Chan];
            AMatrix[3][2]= 0;

         BMatrix[0][0]=0;
            BMatrix[1][0]=(XDirection[1]*directions[2][Chan])-(XDirection[2]*directions[1][Chan]);
            BMatrix[2][0]=(XDirection[2]*directions[0][Chan])-(XDirection[0]*directions[2][Chan]);
            BMatrix[3][0]=(XDirection[0]*directions[1][Chan])-(XDirection[1]*directions[0][Chan]);

         CON_Matrix_t<double> Res;

         PseudoInverse(4,3,1, AMatrix,BMatrix,Res);

          for (int i=0;i<3;i++) xMatrix[i] = Res[i][0];

         for(i=0;i<3;i++) v1[i]=xMatrix[i];
        }
     else
        {
         for(int i=0;i<3;i++) v1[i]=XDirection[i];
        } */


    // allocate returning Block
    int nNumberOfSensorPositions = positions.length();
    int nNumberOfCoils = coil_positions.length();
    IntegrationPoints3D.allocate(nNumberOfCoils,7,3);

    P= v1.norm();

    v1[0] /= P;
    v1[1] /= P;
    v1[2] /= P;

    Crossprod(v2,Dir,v1);   //  v2 = Dir x v1

    for (c=0;c<nNumberOfCoils;c++)
        for(p=0;p<7;p++)
            {
             xl = integration_points[p][0];
             yl = integration_points[p][1];
             zl = coil_positions[c];
             for (k=0;k<3;k++)
                 IntegrationPoints3D[c][p][k] = positions[k][Chan] + xl*v1[k]+yl*v2[k]+zl*Dir[k];
            }

    return IntegrationPoints3D;
}

double anMEGSensor_s::WeberToTesla()
{
     int i;
    double N= 0;

    int nNumberOfCoils = coil_positions.length();

    // a) find coil with lowest z amongst those with positive sense
    float MinZ = 1e30f;
    for (i=0;i<nNumberOfCoils;i++)
    if (coil_windings[i]<0)
    if (coil_positions[i] < MinZ)
      MinZ = coil_positions[i];
    // b) sum up numbers of windings of all coils with negative sense and z smaller than a)
    for (i=0;i<nNumberOfCoils;i++)
    if (coil_windings[i]>0 && coil_positions[i] < MinZ)
      N += coil_windings[i]*coil_areas[i];
    // c) if N=0 negative senses are summed up
    if (!N)
    for (i=0;i<nNumberOfCoils;i++)
        N += fabs(coil_windings[i]*coil_areas[i]);
    return N;
}

anMEGSensorConfiguration_c::anMEGSensorConfiguration_c()
{

}

CON_Vector_t<double> anMEGSensorConfiguration_c::getPositionOfCoil(int number)
{
 CON_Vector_t<double> tmp;
 if (number < getNumberOfChannels()) return m_Coils[number].PositionOfCoil;

 return tmp;
}

CON_Matrix_t<double> anMEGSensorConfiguration_c::getPositionMatrix()
{
 CON_Matrix_t<double> tmp;
 tmp.allocate(getNumberOfCoils(),3);

 for (int i=0;i<getNumberOfCoils();i++)
   tmp[i]=getPositionOfCoil(i);

 return tmp;
}

/*
int anMEGSensorConfiguration_c::getNumberOfSensors()
{
 return m_Channel.height();
}
*/

CON_Matrix_t<double> anMEGSensorConfiguration_c::getCorrectionMatrix()
{

 return m_correctionMatrix;
}

bool anMEGSensorConfiguration_c::setMatrices(CON_Matrix_t<double>&  ChannelCoilMatrix,  CON_Matrix_t<double>&
correctionMatrix, CON_Matrix_t<double>& CoilPositions, CON_Matrix_t<double>& CoilNormal, CON_Vector_t<double>& CoilAreas, CON_Vector_t<int>& CoilWindings, CON_Matrix_t<int>& Polygons,CON_Vector_t<int>& numberintegrationPoints ,CON_Matrix_t<double>& integrationPoints, CON_Matrix_t<double>& integrationWeights,
 CON_Vector_t<double>& CircleRadius, CON_Matrix_t<double>& CircleNormals,CON_Vector_t<double> weber2tesla, CON_Matrix_t<double>& PolygonsPoints,CON_Vector_t<int>& PolygonCoilMatch)
{
 m_PointIntegration=false;
 m_CircleIntegration=false;
 m_PolygonIntegration=false;

 int count=0; //consistent integration scheme
 if (Polygons.length()!=0)  {m_PolygonIntegration=true; count++;}
 if ((integrationPoints.length()!=0) && (integrationPoints.height()!=0)) {m_PointIntegration=true;count++;}
 if ((CircleRadius.length()!=0)  && (CircleNormals.length()!=0) && (CircleNormals.height()!=0)) {m_CircleIntegration=true; count++;}
 //if (count > 1) {printf("ERROR: inconsistent integration scheme \n"); return false;}
 count=1;
 if (m_PointIntegration && weber2tesla.length()>0) m_weber2tesla=weber2tesla;
// if (m_PolygonIntegration && ( ) && () )

 m_Channel=ChannelCoilMatrix;

 m_Coils = (anCoil_s *) malloc(sizeof(anCoil_s)  * getNumberOfCoils());

 CON_Vector_t<double> Dir,v1,v2;
 m_NumberofChannels=ChannelCoilMatrix.height();

 if(correctionMatrix.length()==m_NumberofChannels && correctionMatrix.height()==m_NumberofChannels) m_correctionMatrix=correctionMatrix; else
  {
    m_correctionMatrix.allocate(m_NumberofChannels,m_NumberofChannels);
    for(int i=0;i<m_NumberofChannels;i++)
    {
      m_correctionMatrix[i][i]=1;
    }
  }

 int tmp=getNumberOfCoils();
 int t=0;
 for(int i=0;i<tmp;i++)
 {

   m_Coils[i].PositionOfCoil.allocate(CoilPositions.length());
   m_Coils[i].Normal.allocate(CoilPositions.length());

   m_Coils[i].PositionOfCoil[0]=CoilPositions[i][0] ; m_Coils[i].PositionOfCoil[1]=CoilPositions[i][1]; m_Coils[i].PositionOfCoil[2]=CoilPositions[i][2];
   m_Coils[i].Normal[0]=CoilNormal[i][0] ; m_Coils[i].Normal[1]=CoilNormal[i][1]; m_Coils[i].Normal[2]=CoilNormal[i][2];

   if(CoilAreas.length()>=i+1)
     m_Coils[i].coil_area=CoilAreas[i];
   else  m_Coils[i].coil_area=CoilAreas[0];

   if(CoilWindings.length()>=i+1) m_Coils[i].coil_winding=CoilWindings[i];
   else m_Coils[i].coil_winding=CoilWindings[0];

   //point integration
   if (m_PointIntegration)
   {
     int number_integration_points=numberintegrationPoints[i];
     m_Coils[i].point_integration=true;
     m_Coils[i].integration_points.allocate(number_integration_points,3);
     m_Coils[i].integration_weights.allocate(number_integration_points);
     for (int k=0;k<number_integration_points;k++)
     {
      double x=integrationPoints[t+k][0];
      m_Coils[i].integration_points[k][0]=x;
      double y=integrationPoints[t+k][1];
      m_Coils[i].integration_points[k][1]=y;
      double z=integrationPoints[t+k][2];
      m_Coils[i].integration_points[k][2]=z;
      m_Coils[i].integration_weights[k]=integrationWeights[i][k];
     }
     t=t+number_integration_points;
   }


   if (m_PolygonIntegration)
   {
     int tmpp=0;
     for (int k=0;k<PolygonCoilMatch.length();k++)
       if(PolygonCoilMatch[k]==i) tmpp++;

     m_Coils[i].polygons.allocate(tmpp,2);
     m_Coils[i].polygons_points.allocate(tmpp,3);

     tmpp=0;
     for (int k=0;k<PolygonCoilMatch.length();k++)
       if(PolygonCoilMatch[k]==i)
          {
        m_Coils[i].polygons[tmpp][0]=Polygons[k][0];
        m_Coils[i].polygons[tmpp][1]=Polygons[k][1];
        m_Coils[i].polygons_points[tmpp][0]=PolygonsPoints[k][0];
        m_Coils[i].polygons_points[tmpp][1]=PolygonsPoints[k][1];
        m_Coils[i].polygons_points[tmpp][2]=PolygonsPoints[k][2];
        tmpp++;
      }

   }
   /*
   else
   if (m_CircleIntegration)
   {


   }*/

 }

 return true;
}

int anMEGSensorConfiguration_c::getNumberOfChannels()
{
  return m_NumberofChannels;
}

int anMEGSensorConfiguration_c::getNumberOfCoils() //number of total coils
{
  return m_Channel.length();
}


int anMEGSensorConfiguration_c::getNumberCoilsPerChannel(int channel,CON_Vector_t<int>& numbers)
{
int count=0;

 for (int i=0;i<getNumberOfCoils();i++)
  if(m_Channel[channel][i]!=0) {count++;}

numbers.allocate(count);
 count=0;
 for (int i=0;i<getNumberOfCoils();i++)
  if(m_Channel[channel][i]!=0) {numbers[count]=i;count++;}

  return count;
}


bool anMEGSensorConfiguration_c::computeIntegrationPointsNonWeighted(int coil_number, CON_Matrix_t<double>& Positions, CON_Matrix_t<double>& Normals,  CON_Vector_t<double>& IntegrationWeights)
{

  CON_Vector_t<int> coilnumbers;

  int maxintpoints=0;

  for(int i=0;i<m_Coils[coil_number].integration_points.height();i++)
   if( m_Coils[coil_number].integration_points.height()>maxintpoints) maxintpoints=m_Coils[coil_number].integration_points.height();

  Positions.allocate(maxintpoints,3);Normals.allocate(1,3);


   for(int k=0;k<m_Coils[coil_number].integration_points.height();k++)
    {
     Positions[k][0]=m_Coils[coil_number].integration_points[k][0];
     Positions[k][1]=m_Coils[coil_number].integration_points[k][1];
     Positions[k][2]=m_Coils[coil_number].integration_points[k][2];
    }

      Normals[0][0]=m_Coils[coil_number].Normal[0]; //direction [2][0..2] is the normal
      Normals[0][1]=m_Coils[coil_number].Normal[1];
      Normals[0][2]=m_Coils[coil_number].Normal[2];
         IntegrationWeights=m_Coils[coil_number].integration_weights;


  return true;
}

bool anMEGSensorConfiguration_c::computeIntegrationPointsWeighted(int channel_number, CON_Matrix_t<double>& Positions, CON_Matrix_t<double>& Normals,  CON_Vector_t<double>& weights)
{
/*
  bool result=computeIntegrationPointsNonWeighted(channel_number,Positions,Normals,weights);
  if (!result) return false;

  CON_Vector_t<int> coilnumbers;
  int coils=getNumberCoilsPerChannel(channel_number,coilnumbers);
  if(coilnumbers.getMaxi()>getNumberOfCoils()) {printf("Error in Sensordefinition \n"); return false;}

  for(int i=0;i<weights.length();i++)
   weights[i] = m_weightedChannel[channel_number][coilnumbers[i]];
*/
return true;
}

double anMEGSensorConfiguration_c::WeberToTesla(long coil)
{
 return m_weber2tesla[coil];
}

int anMEGSensorConfiguration_c::getNumberIntegrationPointsCoil(int Coilnumber)
{
  return m_Coils[Coilnumber].integration_points.height();
}

double anMEGSensorConfiguration_c::getIntegrationPointWeights(int Coilnumber,int IntegrationPoint)
{
  if(Coilnumber>getNumberOfCoils()) {printf("Error: coil does not exist! \n"); return 0;}
  if(IntegrationPoint>m_Coils[Coilnumber].integration_points.height()) {printf("Error: integration point does not exist! \n"); return false;}
  if(IntegrationPoint>m_Coils[Coilnumber].integration_weights.length()) {printf("Error: integration weight does not exist! \n"); return false;}

  return m_Coils[Coilnumber].integration_weights[IntegrationPoint];

}

double anMEGSensorConfiguration_c::getCoilWindings(int Coilnumber)
{
  if(Coilnumber>getNumberOfCoils()) {printf("Error: coil does not exist! \n"); return 0;}
  return m_Coils[Coilnumber].coil_winding;
}


double anMEGSensorConfiguration_c::getCoilArea(int Coilnumber)
{
 if(Coilnumber>getNumberOfCoils()) {printf("Error: coil does not exist! \n"); return 0;}
 return m_Coils[Coilnumber].coil_area;
}

bool anMEGSensorConfiguration_c::computeIntegrationPolygonNonWeighted(int channel_number,  CON_Block_t<double>& PolygonsAtChannel,  CON_Matrix_t<double>& weights)
{
return true;
}

bool computeIntegrationPolygonWeighted(int channel_number,  CON_Block_t<double>& PolygonsAtChannel,CON_Block_t<double>& NormalAtChannel, CON_Matrix_t<double>& weights)
{
return true;
}

bool getCoilCirclePickup(int coil_number, CON_Vector_t<double>& Position)
{
return true;
}

bool getCoilCircleNormal(int coil_number, CON_Vector_t<double>& Normal)
{
return true;
}

double getCoilCircleRadius(int coil_number)
{
return 1.0;
}


bool checkConsistency()
{
return true;
}

bool isActivePointIntegration()
{
return true;
}

bool isActiveCircleIntegration()
{
return true;
}


bool isActivePolygonIntegration()
{
return true;
}
