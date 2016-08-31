// ### declaration and inline methods of class UTI_ByteOrder_c ###
// Initial author: Olaf Steinstraeter, 2007

#ifndef UTI_BYTEORDER_C_H
#define UTI_BYTEORDER_C_H


// ### Included files ###{{{
#include <string>
//###}}}


// ### Class UTI_ByteOrder_c ###//{{{
// Converts basic data types from external representations
// to native representations, e.g.
//   '16 bit integer, MSB format to 'long'.
class UTI_ByteOrder_c
{
  protected:
    // protected data fields//{{{
    // =====================

    // native byte order: LSB or MSB ?
    //   1 : LSB
    //   0 : MSB
    //  -1 : Not set yet.
    //       Initialized to -1. Will be set by the constructor of the first instance of this class.
    static int byteOrderIsLSB;
    //}}}


    // protected methods//{{{
    // =================

    // determine native data format, especially the byte order
    void getNativeFormat();
    //}}}

  public:
    // public methods//{{{
    // ==============

    // default constructor
    //   The following is tested by the first constructor of this class:
    //     sizeof(char)  == 1 (this may be the standard),
    //     sizeof(short) == 2 (not standard but used in the most implementations and expected by a lot of people, int16 == short),
    //     sizeof(float) == 4 (a hint that the implementation uses single-precission floation point numbers according
    //                        to IEEE Standard for Binary Floating-Point Arithmetic (IEEE 754), float == float32),
    //     sizeof(double) == 8 (a hint that the implementation uses double-precission floation point numbers according
    //                         to IEEE Standard for Binary Floating-Point Arithmetic (IEEE 754), double = float64).
    //   If any of the tests failed the program is aborted.
    inline UTI_ByteOrder_c();


    // information: byte order
    //   native_format_is_lsb():
    //     inline instance method:
    //        native_format_is_lsb == true:
    //           The Least Significant Byte of a multibyte C++ data type (like int) is stored at the lowest memory address.
    //           Other name: Little-Endian.
    //        native_format_is_lsb == false:
    //           The Most Significant Byte of a multibyte C++ data type (like int) is stored at the lowest memory address.
    //           Other names: Big-Endian, Network Byte Order.
    //   get_host_byte_order():
    //     inline class method:
    //       Returns 'l' if host byte order is LSB first, returns 'm' if host byte order is MSB first,
    //       and '?' if the byte order could not be determined (e.g. sizeof(short) != 2).
    //       Can be used even if the constructor aborts the program.
    inline bool native_format_is_lsb() const;
    static inline char get_host_byte_order();


    // swap byte order for unstructured data: exchange two bytes
    static inline void swapBytes(char *b1, char *b2);

    // swap byte order for unstructured data: 16 bit value
    static inline void swap16(char *ptrToNumber);

    // swap byte order for unstructured data: 32 bit value
    static inline void swap32(char *ptrToNumber);

    // swap byte order for unstructured data: 64 bit value
    static inline void swap64(char *ptrToNumber);


    // in-line conversion of single values: MSB <-> native byte order
    //   value in native byte order  ->  value in MSB byte order
    //   and
    //   value in MSB byte order     ->  value in native byte order
    inline void msb(unsigned short &value) const;
    inline void msb(  signed short &value) const;
    inline void msb(         float &value) const;
    inline void msb(        double &value) const;

    // in-line conversion of single values: LSB <-> native byte order
    //   value in native byte order  ->  value in LSB byte order
    //   and
    //   value in LSB byte order     ->  value in native byte order
    inline void lsb(unsigned short &value) const;
    inline void lsb(  signed short &value) const;
    inline void lsb(         float &value) const;
    inline void lsb(        double &value) const;


    // in-line conversion of arrays: MSB <-> native byte order
    //   arrays of n numbers (not n bytes!)
    inline void msb(unsigned short *array, unsigned long n) const;
    inline void msb(  signed short *array, unsigned long n) const;
    inline void msb(         float *array, unsigned long n) const;
    inline void msb(        double *array, unsigned long n) const;

    // in-line conversion of arrays: LSB <-> native byte order
    //   arrays of n numbers (not n bytes!)
    inline void lsb(unsigned short *array, unsigned long n) const;
    inline void lsb(  signed short *array, unsigned long n) const;
    inline void lsb(         float *array, unsigned long n) const;
    inline void lsb(        double *array, unsigned long n) const;


    // string conversion: hexadecimal
    static std::string toHex(const unsigned char *byteArray, int numberOfBytes);
    template<class T> static inline std::string toHex(const T &value);

    // string conversion: binary
    static std::string toBin(const unsigned char *byteArray, int numberOfBytes);
    template<class T> static inline std::string toBin(const T &value);
    //}}}
};
//}}}


// [inline] UTI_ByteOrder_c::UTI_ByteOrder_c(): default constructor//{{{
inline UTI_ByteOrder_c::UTI_ByteOrder_c()
{
  if (byteOrderIsLSB == -1)
    getNativeFormat();
}
//}}}

// [inline] UTI_ByteOrder_c::native_format_is_lsb(): information: byte order//{{{
inline bool UTI_ByteOrder_c::native_format_is_lsb() const
{
  if (byteOrderIsLSB)
    return true;
  else
    return false;
}
//}}}

// [inline class method] UTI_ByteOrder_c::swapBytes() swap byte order for unstructured data: exchange two bytes//{{{
inline void UTI_ByteOrder_c::swapBytes(char *b1, char *b2)
{
  static char h;
  h = *b1;
  *b1 = *b2;
  *b2 = h;
}
//}}}

// [inline class method] UTI_ByteOrder_c::swap16() swap byte order for unstructured data: 16 bit value//{{{
inline void UTI_ByteOrder_c::swap16(char *ptrToNumber)
{
  swapBytes(ptrToNumber, ptrToNumber+1);
}
//}}}

// [inline class method] UTI_ByteOrder_c::swap32(): swap byte order for unstructured data: 32 bit value//{{{
inline void UTI_ByteOrder_c::swap32(char *ptrToNumber)
{
  swapBytes(ptrToNumber,   ptrToNumber+3);
  swapBytes(ptrToNumber+1, ptrToNumber+2);
}
//}}}

// [inline class method] UTI_ByteOrder_c::swap64(): swap byte order for unstructured data: 64 bit value//{{{
inline void UTI_ByteOrder_c::swap64(char *ptrToNumber)
{
  swapBytes(ptrToNumber,   ptrToNumber+7);
  swapBytes(ptrToNumber+1, ptrToNumber+6);
  swapBytes(ptrToNumber+2, ptrToNumber+5);
  swapBytes(ptrToNumber+3, ptrToNumber+4);
}
//}}}

// [inline] UTI_ByteOrder_c::msb(unsigned short&): in-line conversion of single values: MSB <-> native byte order//{{{
inline void UTI_ByteOrder_c::msb(unsigned short &value) const
{
  if (byteOrderIsLSB) swap16(reinterpret_cast<char*>(&value));
}
//}}}

// [inline] UTI_ByteOrder_c::msb(signed short&): in-line conversion of single values: MSB <-> native byte order//{{{
inline void UTI_ByteOrder_c::msb(signed short &value) const
{
  if (byteOrderIsLSB) swap16(reinterpret_cast<char*>(&value));
}
//}}}

// [inline] UTI_ByteOrder_c::msb(float&): in-line conversion of single values: MSB <-> native byte order//{{{
inline void UTI_ByteOrder_c::msb(float &value) const
{
  if (byteOrderIsLSB) swap32(reinterpret_cast<char*>(&value));
}
//}}}

// [inline] UTI_ByteOrder_c::msb(double&): in-line conversion of single values: MSB <-> native byte order//{{{
inline void UTI_ByteOrder_c::msb(double &value) const
{
  if (byteOrderIsLSB) swap64(reinterpret_cast<char*>(&value));
}
//}}}

// [inline] UTI_ByteOrder_c::lsb(unsigned short&): in-line conversion of single values: LSB <-> native byte order//{{{
inline void UTI_ByteOrder_c::lsb(unsigned short &value) const
{
  if (!byteOrderIsLSB) swap16(reinterpret_cast<char*>(&value));
}
//}}}

// [inline] UTI_ByteOrder_c::lsb(signed short&): in-line conversion of single values: LSB <-> native byte order//{{{
inline void UTI_ByteOrder_c::lsb(signed short &value) const
{
  if (!byteOrderIsLSB) swap16(reinterpret_cast<char*>(&value));
}
//}}}

// [inline] UTI_ByteOrder_c::lsb(float&): in-line conversion of single values: LSB <-> native byte order//{{{
inline void UTI_ByteOrder_c::lsb(float &value) const
{
  if (!byteOrderIsLSB) swap32(reinterpret_cast<char*>(&value));
}
//}}}

// [inline] UTI_ByteOrder_c::lsb(double&): in-line conversion of single values: LSB <-> native byte order//{{{
inline void UTI_ByteOrder_c::lsb(double &value) const
{
  if (!byteOrderIsLSB) swap64(reinterpret_cast<char*>(&value));
}
//}}}

// [inline] UTI_ByteOrder_c::msb(unsigned short*, unsigned long): in-line conversion of arrays: MSB <-> native byte order//{{{
inline void UTI_ByteOrder_c::msb(unsigned short *array, unsigned long n) const
{
  if (byteOrderIsLSB)
    {
      unsigned short *arrayEnd = array + n;
      for (; array != arrayEnd; array++)
        swap16(reinterpret_cast<char*>(array));
    }
}
//}}}

// [inline] UTI_ByteOrder_c::msb(signed short*, unsigned long): in-line conversion of arrays: MSB <-> native byte order//{{{
inline void UTI_ByteOrder_c::msb(signed short *array, unsigned long n) const
{
  if (byteOrderIsLSB)
    {
      signed short *arrayEnd = array + n;
      for (; array != arrayEnd; array++)
        swap16(reinterpret_cast<char*>(array));
    }
}
//}}}

// [inline] UTI_ByteOrder_c::msb(float*, unsigned long): in-line conversion of arrays: MSB <-> native byte order//{{{
inline void UTI_ByteOrder_c::msb(float *array, unsigned long n) const
{
  if (byteOrderIsLSB)
    {
      float *arrayEnd = array + n;
      for (; array != arrayEnd; array++)
        swap32(reinterpret_cast<char*>(array));
    }
}
//}}}

// [inline] UTI_ByteOrder_c::msb(double*, unsigned long): in-line conversion of arrays: MSB <-> native byte order//{{{
inline void UTI_ByteOrder_c::msb(double *array, unsigned long n) const
{
  if (byteOrderIsLSB)
    {
      double *arrayEnd = array + n;
      for (; array != arrayEnd; array++)
        swap64(reinterpret_cast<char*>(array));
    }
}
//}}}

// [inline] UTI_ByteOrder_c::lsb(unsigned short*, unsigned long): in-line conversion of arrays: LSB <-> native byte order//{{{
inline void UTI_ByteOrder_c::lsb(unsigned short *array, unsigned long n) const
{
  if (!byteOrderIsLSB)
    {
      unsigned short *arrayEnd = array + n;
      for (; array != arrayEnd; array++)
        swap16(reinterpret_cast<char*>(array));
    }
}
//}}}

// [inline] UTI_ByteOrder_c::lsb(signed short*, unsigned long): in-line conversion of arrays: LSB <-> native byte order//{{{
inline void UTI_ByteOrder_c::lsb(signed short *array, unsigned long n) const
{
  if (!byteOrderIsLSB)
    {
      signed short *arrayEnd = array + n;
      for (; array != arrayEnd; array++)
        swap16(reinterpret_cast<char*>(array));
    }
}
//}}}

// [inline] UTI_ByteOrder_c::lsb(float*, unsigned long): in-line conversion of arrays: LSB <-> native byte order//{{{
inline void UTI_ByteOrder_c::lsb(float *array, unsigned long n) const
{
  if (!byteOrderIsLSB)
    {
      float *arrayEnd = array + n;
      for (; array != arrayEnd; array++)
        swap32(reinterpret_cast<char*>(array));
    }
}
//}}}

// [inline] UTI_ByteOrder_c::lsb(double*, unsigned long): in-line conversion of arrays: LSB <-> native byte order//{{{
inline void UTI_ByteOrder_c::lsb(double *array, unsigned long n) const
{
  if (!byteOrderIsLSB)
    {
      double *arrayEnd = array + n;
      for (; array != arrayEnd; array++)
        swap64(reinterpret_cast<char*>(array));
    }
}
//}}}

// [inline] UTI_ByteOrder_c::toHex(): string conversion: hexadecimal//{{{
template<class T>
inline std::string UTI_ByteOrder_c::toHex(const T &value)
{
  return toHex(reinterpret_cast<const unsigned char*>(&value), sizeof(value));
}
//}}}

// [inline] UTI_ByteOrder_c::toBin(): string conversion: binary//{{{
template<class T>
inline std::string UTI_ByteOrder_c::toBin(const T &value)
{
  return toBin(reinterpret_cast<const unsigned char*>(&value), sizeof(value));
}
//}}}

// [inline class method] UTI_ByteOrder_c::get_host_byte_order(): byte order of the host//{{{
inline char UTI_ByteOrder_c::get_host_byte_order()
{
  if (UTI_ByteOrder_c::byteOrderIsLSB == -1)
    {
      short s = 1;
      if (sizeof(s) != 2)
        return '?';

      char *c = reinterpret_cast<char *>(&s);
      if (c[0] == 1)
        {
          UTI_ByteOrder_c::byteOrderIsLSB = 1;
          return 'l';
        }
      else if (c[1] == 1)
        {
          UTI_ByteOrder_c::byteOrderIsLSB = 1;
          return 'm';
        }
      else
        return '?';
    }
  else
    {
      if (UTI_ByteOrder_c::byteOrderIsLSB)
        return 'l';
      else
        return 'm';
    }
}
///}}}
#endif
