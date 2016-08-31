//$15   04.06.2003  Anwander A.  Extention of the simulator to 2 or 3 source directions
//$14    11.03.2002  Matthias D.  Removed spelline error DoSetParamterValues();
//$13    07.02.2002  Matthias D.  Made ANSI compliant
//$12     06.02.2002    Matthias D.  Corrected spelling error
//$11    14.01.2002    Matthias D.  Changed IP_SIM_SensorConfiguration_c to anEEGSensorConfiguration_c
//$10    03.08.2001  Matthias D.  Added DoSetParamterValues()
//$9    31.07.2001  Matthias D.  m_SimParameters are filled with conductivity values. !!! Currently m_SimParameters are not used for calculations !!
//$8    04.04.2001  Matthias D.  Corrected electrode positions: Projection on outer sphere
//$7    04.04.2001  Matthias D.  Removed bug: Bias is removed from output results
//$6    14.03.2001    Matthias D.  Removed bug: Adjustement of source positions due to center of sphere
//$5    14.11.2000  Matthias D.  Changed design for the various types of source directions
//$4    12.10.2000  Matthias D.  member functions are used now for the access of EEG electrode variables
//$3    22.09.2000  Matthias D.  Changed from \ to / for directory names
//$2    20.09.2000  Matthias D.  copied method from fweegforwardsimulator to IP_SIM_SimulatorEEGSpheres_c
//$1    28.08.2000  Matthias D.  created

// IP_SIM_SimulatorEEGSpheres_c.cpp: implementation of the IP_SIM_SimulatorEEGSpheres_c class.
//
//////////////////////////////////////////////////////////////////////


#include <stdio.h>
#include "StdAnalysisPch.h"

using namespace std;

#include "IP_SIM_SimulatorEEGSpheres_c.h"


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


IP_SIM_SimulatorEEGSpheres_c::IP_SIM_SimulatorEEGSpheres_c(anEEGSensorConfiguration_c& inSensorconfiguration,CON_Vector_t<double>& inRadii,CON_Vector_t<double>& inCenter,
                            CON_Vector_t<double>& inConductivities):
IP_SIM_AbstractSimulatorEEGMEG_c(),
m_SensorConfiguration(inSensorconfiguration),
m_Radii(inRadii),
m_Center(inCenter),
m_Conductivities(inConductivities)
{
  for (int i= 0; i < m_Conductivities.length(); i++)
    {
      std::string name = "Conductivity Sphere: ";
      char buffer[20];
      sprintf (buffer, "%10i", i); //itoa(i,buffer,10);
      name += buffer;

      std::string unit = "S/m";
      IP_SIM_ParameterSimulator_c Conductivity(name, unit, 0.0, 10e35, m_Conductivities[i], 0.9*m_Conductivities[i], m_Conductivities[i]);

      m_SimParameters.push_back(Conductivity);
    }
}

IP_SIM_SimulatorEEGSpheres_c::~IP_SIM_SimulatorEEGSpheres_c()
{

}


bool IP_SIM_SimulatorEEGSpheres_c::isValidSimulator()
{
return true;
}

// Set and Get Radii and Conductivities

bool IP_SIM_SimulatorEEGSpheres_c::setRadii(CON_Vector_t<double>& inRadii)
{
if ((inRadii.length()<1)||(inRadii.length()<4)) return false;
m_Radii = inRadii;
return true;
}

bool IP_SIM_SimulatorEEGSpheres_c::getRadii(CON_Vector_t<double>& outRadii)
{
if (m_Radii.length() <1) return false;
outRadii = m_Radii;
return true;
}

int  IP_SIM_SimulatorEEGSpheres_c::getNumberofSpheres()
{
return m_Radii.length();
}

bool IP_SIM_SimulatorEEGSpheres_c::setCenter(CON_Vector_t<double>& inCenter)
{
if (inCenter.length()!=3) return false;
m_Center = inCenter;
return true;
}

bool IP_SIM_SimulatorEEGSpheres_c::getCenter(CON_Vector_t<double>& outCenter)
{
if (m_Center.length()!=3) return false;
outCenter = m_Center;
return true;
}

bool IP_SIM_SimulatorEEGSpheres_c::setConductivities(CON_Vector_t<double>& inConductivities)
{
if ((inConductivities.length()<1)||(inConductivities.length()>4)) return false;
m_Conductivities = inConductivities;
return true;
}

bool IP_SIM_SimulatorEEGSpheres_c::getConductivities(CON_Vector_t<double>& outConductivities)
{
if (m_Conductivities.length()<1) return false;
outConductivities = m_Conductivities;
return true;
}



bool IP_SIM_SimulatorEEGSpheres_c::DoSetParameterValues()
{

    for (int ParameterNumber = 0; ParameterNumber < getNumberofParameters(); ParameterNumber++)
        {
         std::string ParameterName;
         ParameterNgetName(  ParameterNumber, ParameterName);

        // Conducitvity values
        if ( ParameterName == "Conductivity Sphere: 0")
            {
              double value;
            ParameterNgetValue(  ParameterNumber, value);
            m_Conductivities[0] = value;
            }

        if ( ParameterName == "Conductivity Sphere: 1")
            {
              double value;
            ParameterNgetValue(  ParameterNumber, value);
            m_Conductivities[1] = value;
            }

        if ( ParameterName == "Conductivity Sphere: 2")
            {
              double value;
            ParameterNgetValue(  ParameterNumber, value);
            m_Conductivities[2] = value;
            }

        if ( ParameterName == "Conductivity Sphere: 3")
            {
              double value;
            ParameterNgetValue(  ParameterNumber, value);
            m_Conductivities[3] = value;
            }

        }

 return true;

}
/////////////Compute Gain Matrix/Estimated Potentials at Sensor Positions
// This function should return a gain matrix, each row referring to a channel and each column
// to a dipole. If no directions are given, a gain matrix is computed, containing 3 columns
// for each dipole, representing the standard directions (x,y,z). If only one direction is
// passed, this one is used for all dipoles.

void IP_SIM_SimulatorEEGSpheres_c::computeGainMatrix(int NumberDip, const CON_Matrix_t<double>& inPos, const CON_Matrix_t<double>& inDir, CON_Matrix_t<double>& outSimData)
{
    CON_Vector_t<double>    rpp(3,0);                 // Polar coordinates of dipole location
    double                r0,r0aux;                 // r of dipole (polar coordinates)
    CON_Vector_t<double>    q(3,0);                     // rotated dipole moments
    double                cost,sint,cosp,sinp;
    double                epsilon[4] = {0,0,0,0};     // Will be filled with conductivity values
    CON_Vector_t<double>    elc_pos(3,0);             // Electrode position
    CON_Vector_t<double>    auxpos(3,0);             // Corrected electrode position
    CON_Vector_t<double>    r_k(3,0);
    double                r[4];                    // radii ot the spheres
    CON_Vector_t<double>    pdip(3,0);
    int i,j,k;                                     // Counter variables


    int numberSegments = 0;

    switch (inDir.height())
    {
        case 0:
            numberSegments = 3;
            break;
        case 3:
            numberSegments = 1;
            break;
        case 6:
            numberSegments = 2;
            break;
        case 9:
            numberSegments = 3;
            break;
    }

    // Allocate Result Matrix
    outSimData.allocate(m_SensorConfiguration.getNumberEEGElectrodes(),numberSegments * NumberDip);

    // Get number of spheres and conductivities
    for (j=0;j<m_Conductivities.length();j++)
        epsilon[j] = m_Conductivities[j];
    // fill remaining shells with conductivity of innermost radius
    for (;j>0 && j<4;j++)
        epsilon[j] = epsilon[j-1];

    // Compute for each Dipole ...
    for (int s = 0; s < numberSegments;s ++)
        for (i = 0; i < NumberDip; i ++)
        {
            CON_Vector_t<double>   Pos(3,0.0);
            CON_Vector_t<double>   Dir(3,0.0);
            for(j=0;j<3;j++)
            {
                Pos[j]=inPos[j][i]-m_Center[j];
                // if no directions are given, use standard ones
                if (!inDir.exist() )
                    Dir[j] = (j==s);
                // if fewer directions are given than dipoles, first direction is taken for all dipoles
                if (inDir.length() < NumberDip && inDir.length() > 0)
                    Dir[j] = inDir[j+3*s][0];
                // if there is a direction for every dipole, just use them
                // the directions for the different segments are stored in lines 1-3 / 4-6 / 7-9 of inDir
                if (inDir.length() == NumberDip)
                    Dir[j]=inDir[j+3*s][i];
            }
            Cart2Sphere(Pos, rpp); // rpp:  r,theta phi of the dipole location
            r0 = rpp[0];
            cost = cos(rpp[1]); sint = sin(rpp[1]);
            cosp = cos(rpp[2]); sinp = sin(rpp[2]);

            // rotation of the dipole moment along theta and phi of the dipole location
            q[0] =  Dir[0] * cosp * cost + Dir[1] * sinp * cost - Dir[2] * sint;
            q[1] = -Dir[0] * sinp        + Dir[1] * cosp;
            q[2] =  Dir[0] * cosp * sint + Dir[1] * sinp * sint + Dir[2] * cost;

            // Compute for each Electrode ...
            for (k=0;k < m_SensorConfiguration.getNumberEEGElectrodes();k++)
            {
                elc_pos = m_SensorConfiguration.getEEGElectrodePosition(k);
                auxpos = elc_pos - m_Center;
                // Correct auxpos, that auxpos is on surface of outer sphere
                double auxposradius=sqrt(auxpos[0]*auxpos[0]+auxpos[1]*auxpos[1]+auxpos[2]*auxpos[2]);
                for (j=0;j<3;j++) auxpos[j] = auxpos[j]*m_Radii[0]/auxposradius;

                // rotation of the field point along theta and phi of the dipole location
                r_k[0] =  auxpos[0] * cosp * cost + auxpos[1] * sinp * cost - auxpos[2] * sint;
                r_k[1] = -auxpos[0] * sinp        + auxpos[1] * cosp;
                r_k[2] =  auxpos[0] * cosp * sint + auxpos[1] * sinp * sint + auxpos[2] * cost;

                Cart2Sphere(r_k,rpp); // rpp: r,theta,phi of the electrode position

                for (j=0;j<m_Radii.length();j++)
                    r[j] = m_Radii[j];
                // fill remaining shells with innermost radius
                for (;j>0 && j<4;j++)
                    r[j] = r[j-1];

                r0aux = r0;
                if ( r0aux < 1E-10 )
                    r0aux = 1E-10;
                // look at the constraint for the dipole pos
                if ( r0aux/r[3] > 0.9999 )
                    r0aux = 0.9999*r[3];

                computeEEGFourLayerSphere(rpp[1],rpp[2],r0aux,epsilon,epsilon,r,pdip);
                // Arguments: rpp[1] = theta of electrode, rpp2(2) = phi of electrode position,
                // r0aux = corrected r of dipole position, epsilon = conductivity values,
                // r = radii of spheres, pdip = result

                outSimData[k][s*NumberDip+i] = Scalar(pdip,q);
            }    // loop over electrodes

        // remove bias
        double dNorm = 0;
        for(k=0;k<outSimData.height();k++) dNorm += outSimData[k][s*NumberDip+i];
        dNorm =dNorm/outSimData.height();
        for(k=0;k<outSimData.height();k++)
                outSimData[k][s*NumberDip+i] -= dNorm;
        } // loop over dipoles
}


void IP_SIM_SimulatorEEGSpheres_c::computeGainMatrix(const CON_Vector_t<double>& inParameters,CON_Matrix_t<double>& outSimData)
{

}


#define EPSILON 1E-6
#define MAXTERM 3000
void IP_SIM_SimulatorEEGSpheres_c::computeEEGFourLayerSphere(double cta, double fi, double r0,double epslon[4], double eta[4], double r[4], CON_Vector_t<double>& pdip)
 {
    int         i, j, k, n,nterm;
    double        dom, cdom, dx, cosfi, sinfi,a;
    double      d11, d12, d13, d14, d21, d22, d23, d24, e1, e2, e3, e4, err1, err2;
    double        sumc[2];
    double      v[4];
    double                    c1[MAXTERM+1], c2[MAXTERM+1];
    double                    c[4][2][2];
    double                    ep = EPSILON;

    n = 0;
    sumc[0] = 0;
    sumc[1] = 0;

    for (;;)
    {
        n++;
        for (j = 0; j < 4; j++)
            v[j] = 0.5 * (-1 + sqrt(1 + 4 * (n + 1) * n * eta[j] / epslon[j]));
        for (j = 1; j < 4; j++)
        {
            dom = -(2 * v[j] + 1) * epslon[j];
            c[j][0][0] = (-(epslon[j - 1] * v[j - 1] + epslon[j] * (v[j] + 1))) / dom;
            c[j][0][1] = (epslon[j - 1] * (v[j - 1] + 1) - epslon[j] * (v[j] + 1)) / dom;
            c[j][1][0] = (epslon[j - 1] * v[j - 1] - epslon[j] * v[j]) / dom;
            c[j][1][1] = (-(epslon[j - 1] * (v[j - 1] + 1) + epslon[j] * v[j])) / dom;
        }

        d11 = c[3][1][0] * c[2][0][0] * c[1][0][0];
        d12 = c[3][1][0] * c[2][0][1] * c[1][1][0];
        d13 = c[3][1][1] * c[2][1][0] * c[1][0][0];
        d14 = c[3][1][1] * c[2][1][1] * c[1][1][0];

        d21 = c[3][1][0] * c[2][0][0] * c[1][0][1];
        d22 = c[3][1][0] * c[2][0][1] * c[1][1][1];
        d23 = c[3][1][1] * c[2][1][0] * c[1][0][1];
        d24 = c[3][1][1] * c[2][1][1] * c[1][1][1];

        e1 = pow(r[1] / r[0], v[0]);
        e2 = pow(r[2] / r[1], v[1]);
        e3 = pow(r[3] / r[2], v[2]);
        e4 = pow(r0 / r[3], v[3]);

        cdom = d11 * e1 * e1 * e2 * e2 * e3 * e3 * r[3] / r[0];
        cdom += d12 * e1 * e1 * e3 * e3 * r[3] / r[2] * r[1] / r[0];
        cdom += d13 * e2 * e2 * e1 * e1 * r[2] / r[0];
        cdom += d14 * e1 * e1 * r[1] / r[0];
        cdom *= (v[0] + 1) / v[0];
        cdom += d21 * e3 * e3 * e2 * e2 * r[3] / r[1];
        cdom += d22 * e3 * e3 * r[3] / r[2];
        cdom += d23 * e2 * e2 * r[2] / r[1];
        cdom += d24;

        c1[n] = (2 * n + 1) * (2 * v[0] + 1) * v[3] / ((2 * v[3] + 1) * v[0] * epslon[3])
        * e1 * e2 * e3 * e4 / (cdom * r0 * r[0]);
        c2[n] = (2 * n + 1) * (2 * v[0] + 1) / ((2 * v[3] + 1) * v[0] * epslon[3])
        * e1 * e2 * e3 * e4 / (cdom * r[0]);

        sumc[0] += c1[n];
        sumc[1] += c2[n];
        err1 = fabs(c1[n] / sumc[0]);
        err2 = fabs(c2[n] / sumc[1]);

        if (err1 < ep && err2 < ep || n >= MAXTERM )
        {
            nterm = n;
            if ( n >= MAXTERM )
            {
            }
            break;
        }
    }
    dx = cos(cta);
    cosfi = cos(fi);
    sinfi = sin(fi);
    for (k = 0; k < 3; k++)
        pdip[k] = 0;
    for (n = 1; n <= nterm; n++)
    {
        a= -c2[n] * plgndr(n, 1, dx) / r0;
        pdip[0] += a * cosfi;
        pdip[1] += a * sinfi;
        pdip[2] += c1[n] * plgndr(n, 0, dx);
    }
    for (k = 0; k < 3; k++)
        pdip[k] *= 0.25 / M_PI;
}

double IP_SIM_SimulatorEEGSpheres_c:: plgndr(int l, int m, double x)   // Compute Result of Legendre Series
{
    int        i, ll;
    double    pmm, pmmp1, somx2, fact, pll;

    if (m < 0 || m > l || fabs(x) > 1)
        return 0;
    pmm = 1;
    if (m > 0)
    {
        somx2 = sqrt((1 - x) * (1 + x));
        fact = 1;
        for (i = 0; i < m; i++)
        {
            pmm = -pmm * fact * somx2;
            fact += 2;
        }
    }
    if (l == m)
        return pmm;
    pmmp1 = x * (2 * m + 1) * pmm;
    if (l == m + 1)
        return pmmp1;
    for (ll = m + 2; ll <= l; ll++)
    {
        pll = (x * (2 * ll - 1) * pmmp1 - (ll + m - 1) * pmm) / (ll - m);
        pmm = pmmp1;
        pmmp1 = pll;
    }
    return pll;
 }
