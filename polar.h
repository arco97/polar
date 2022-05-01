/*****************************************************************//**
 * @file   polar.h
 * @brief  Polar is a template class to handl complex data in polar format <br>
 * The Magnitue mag[] and the Angle ang&deg;[rad]
 *
 * @author Andre Hoehne
 * @date   Maerz 2022
***********************************************************************/

#ifndef POLAR_H
#define POLAR_H
#include <cmath>
#include <complex>
#include <sstream>
#include <iostream>
#include <ostream>
#include <math.h>
#include <vector>
namespace usr{
  template<typename _Tp> class polar;
  /**
   * .
   */
  template<typename _Tp = double>
  class polar{
  public:
      typedef enum {Rad,Grad180}phase;
    /**
     * .
     */
    polar(_Tp mag = _Tp(),_Tp ang = _Tp())
      :_m_mag(mag),_m_ang(ang){}
    /**
     * .
     */
    polar(std::complex<_Tp> comp)
      :_m_mag(std::abs(comp)),_m_ang(std::arg(comp)){}
    /**
    @param mag set the magnitude r[ ].
    @return void.
    */
    void mag(_Tp mag){_m_mag=mag;}
    /**
     * @param ang set the angel &deg;[rad].
     */
    void ang(_Tp ang){_m_ang=ang;}
    /**
     * @return get the Magnitued r[ ].
     */
    const _Tp mag(){return _m_mag;}
    /**
     * @return get the Angle &deg;[rad].
     */
    const _Tp ang(){return _m_ang;}
    /**
     * @return get the imaginary part of the Magnitude and the Angle $$ z \in \mathbb{C} \newline \Im(z) = mag*\sin(ang)$$
     */
    const _Tp imag(){return _m_mag*sin(_m_ang);}
    /**
     * @return get the real part of the Magnitude and the Angle $$ z \in \mathbb{C} \newline \Re(z) = mag*\cos(ang)$$
     */
    const _Tp real(){return _m_mag*cos(_m_ang);}
    /**
    * @param factor is the factor befor the logarytmical
     * @return get the Magnitude logarytmical  $$ mag = factor*\log_{10}(mag) dB$$
     */
    const _Tp decibel(_Tp factor=20){return factor*log10(_m_mag);}
    /**
    * 
    */
    static const _Tp mag(const std::complex<_Tp> comp) { return std::abs(comp);}
    /**
    * 
    */
    static const _Tp ang(const std::complex<_Tp> comp, const phase pha) {
        if (pha == phase::Rad)return std::arg(comp);
        else if (pha == phase::Grad180) return std::arg(comp) * (180 / M_PI);
    }
    /**
    * 
    */
    static const _Tp mag(const polar<_Tp> pol) { return pol.mag(); }
    /**
    * 
    */
    static const _Tp ang(const polar<_Tp> pol,const phase pha) {
        if (pha == phase::Rad)return pol.ang();
        else if (pha == phase::Grad180) return pol.ang() * (180 / M_PI);
    }
    /**
    * 
   */
    std::vector<polar<_Tp>> unwrapped(const std::vector<polar<_Tp>> vpol,bool underZero = 0,phase pha = phase::Rad) {
        double offset =0;
        double angle = 0;
        std::vector<polar<_Tp>> result = vpol;
        switch (pha) {
        case phase::Rad:
            angle = M_PI;
            break;
        case phase::Grad180:
            angle = 180;
            break;
        }
           
        for (size_t i = 1; i < result.size(); i++) {
            double difference = result.at(i).ang() - result.at(i - 1).ang();
            if (difference > angle) offset -= 2 * angle;
            else if (difference < -angle) offset += 2 * angle;
            result.at(i).ang() += offset;
        }
        if(underZero&&result.at(0)>0)
            for(size_t i = 1; i < result.size(); i++)result.at(i).ang() -= 360;
    return std::vector<polar<_Tp>>(result);
    }

    std::ostream& operator <<(const usr::polar<_Tp> pol) {
        std::ostream *os;
        *os<<'('<<pol.mag()<<';'<<pol.ang()<<')';
        return *os;
    }

    private:
    _Tp _m_mag,_m_ang;
  };
}

template<typename _Tp>
std::ostream& operator <<(std::ostream &os ,const usr::polar<_Tp> &pol) {
    os << '(' << pol.mag() << ';' << pol.ang() << ')';
    return os;
}
///  Insertion operator for complex values.
template<typename _Tp, typename _CharT, class _Traits>
std::basic_ostream<_CharT, _Traits>&
operator<<(std::basic_ostream<_CharT, _Traits>& __os, const usr::polar<_Tp>& __x)
{
    std::basic_ostringstream<_CharT, _Traits> __s;
    __s.flags(__os.flags());
    __s.imbue(__os.getloc());
    __s.precision(__os.precision());
    __s << '(' << __x.mag() << ';' << __x.ang() << ')';
    return __os << __s.str();
}
#endif // POLAR_H
