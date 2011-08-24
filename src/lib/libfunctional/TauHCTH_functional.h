#ifndef TauHCTH_functional_h
#define TauHCTH_functional_h
/**********************************************************
* TauHCTH_functional.h: declarations for TauHCTH_functional for KS-DFT
* Robert Parrish, robparrish@gmail.com
* Autogenerated by MATLAB Script on 25-May-2011
*
***********************************************************/
#include "functional.h"

namespace psi { namespace functional {

class TauHCTH_Functional : public Functional {
public:
    TauHCTH_Functional(int npoints, int deriv);
    virtual ~TauHCTH_Functional();
    virtual void computeRKSFunctional(boost::shared_ptr<RKSFunctions> prop);
    virtual void computeUKSFunctional(boost::shared_ptr<UKSFunctions> prop);
};
}}
#endif

