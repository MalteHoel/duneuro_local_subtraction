//$18    06.02.2002    Matthias D. Removed #include <AnalysisDef.h>
//$17    09.01.2002    Nico        made forward decl. of class anAnalyzerAveragerEventParameters_c
//$16    29.10.2001    Matthias D. Set second parameter in GetData() to long
//$14    16.05.2001    Joachim        Added GetData() and GetDataForEventList() overrides using a matrix of float
//$13    23.04.2001    Frank N.    Changed parameters of GetDataForEventList()
//$12    12.01.2001    Matthias D. Added parameters to GetDataforEventList() to deal with uncomplete epochs
//$11    08.01.2001    Matthias D. Added GetDataforEventList();
//$10   21.12.2000    Matthias D. Added GetDataForSpikeEvents();
//$9    13.10.2000    Frank N.    Removed argument bSignalSourcesWithPosition to GetData()
//$8    29.08.2000    Frank N.    Added argument bSignalSourcesWithPosition to GetData()
//$7    15.08.2000    Frank N.    Made GetData() return a bool
//$6    04.04.2000    Frank N.    Added GetSampleCount() as pure virtual
//$5    30.11.1999    Frank Z.    changed declaration of getsignalsources
//$4    10.11.1999    Frank Z.    inserted virtual keyword
//$3    03.11.1999    Frank Z.    inserted GetSamplingFrequency()
//$2    27.08.1999    Frank Z.    Took out includes
//$1    20.08.1999    Frank Z.    Created

#ifndef __IP_SIG_ABSTRACTDATAHANDLER_C_H__
#define __IP_SIG_ABSTRACTDATAHANDLER_C_H__

#include "CON_Vector_t.h"
#include "CON_Matrix_t.h"
#include "CON_Block_t.h"

class anAnalyzerAveragerEventParameters_c;


////////////////////// --- abstract base classes for EEG/MEG analysis --- ///////////////////////////////////

// class IP_SIG_AbstractDatahandler_c serves as a base class for all methods used in detection methods
// IP_SIG_AbstractDatahandler_c defines the interfaces needed to implement a spike or feature detection method.

class ANALYSIS_EXPORT IP_SIG_AbstractDatahandler_c
{
protected:
    IP_SIG_AbstractDatahandler_c();
public:
    virtual ~IP_SIG_AbstractDatahandler_c();

    virtual bool isValidDataHandler() const = 0;

    virtual float GetSamplingFrequency() const = 0;
    virtual long GetSampleCount() const = 0;

    virtual bool GetData(long lSample, long lSamples, CON_Matrix_t<double>& inData) const = 0;
    virtual bool GetData(long lSample, long lSamples, CON_Matrix_t<float>& inData) const = 0;

    virtual void GetSignalSources(CON_Matrix_t<float>& pos, CON_Matrix_t<float>& ori,CON_Vector_t<bool>&) const = 0;

    // Methods to average data
    virtual bool GetDataForEventList(const std::vector<anAnalyzerAveragerEventParameters_c>& eventlist, CON_Vector_t<int>& ZeroPos, CON_Vector_t<int>& NumSamples, CON_Vector_t<bool>& bEpochComplete, CON_Block_t<float>& data) const = 0;
    virtual bool GetDataForEventList(const std::vector<anAnalyzerAveragerEventParameters_c>& eventlist, CON_Vector_t<int>& ZeroPos, CON_Vector_t<int>& NumSamples, CON_Vector_t<bool>& bEpochComplete, CON_Block_t<double>& data) const = 0;
};

#endif //__anAbstractDataHandler_c_H__
