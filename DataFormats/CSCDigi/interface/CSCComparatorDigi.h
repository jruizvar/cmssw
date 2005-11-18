#ifndef CSCComparatorDigi_CSCComparatorDigi_h
#define CSCComparatorDigi_CSCComparatorDigi_h

/** \class CSCComparatorDigi
 *
 * Digi for CSC Comparators.
 *  
 *  $Date: 2005/11/17 13:04:41 $
 *  $Revision: 1.1 $
 *
 * \author M. Schmitt, Northwestern
 *
 */

class CSCComparatorDigi{

public:

  //
  // The definition of the CSC Comparator Digi.
  // This should never be used directly, only by calling data().
  // Made public to be able to generate lcgdict, SA, 27/4/05
  //
  struct theComparatorDigi {
    unsigned int strip;
    int comparator;
  };

  // Construct from the strip number and the ADC readings.
  explicit CSCComparatorDigi (int strip, int comparator);

  // Copy constructor
  CSCComparatorDigi (const CSCComparatorDigi& digi);

  // Default construction.
  CSCComparatorDigi ();

  // Assignment operator
  CSCComparatorDigi& operator=(const CSCComparatorDigi& digi);

  // Digis are equal if they are on the same strip and have same Comparator data
  bool operator==(const CSCComparatorDigi& digi) const;

  // Get the strip number
  int getStrip() const;

  // Get Comparator readings
  int getComparator() const;

  // Set the strip number
  void setStrip(int strip);

  // Set Comparator data
  void setComparator (int comparator);

  // Print content of digi
  void print() const;

  // Print the binary representation of the digi
  void dump() const;

private:
  friend class testCSCComparatorDigis;

  // Set data words
  void set(int strip, int comparator);

  // access
  theComparatorDigi* data(); 
  const theComparatorDigi* data() const;
  theComparatorDigi aComparatorDigi;
};

#include<iostream>
// needed by COBRA
//inline std::ostream & operator<<(std::ostream & o, const CSCComparatorDigi& digi) {
//  return o << " " << digi.strip()
//	   << " " << digi.ADCCounts();
//}
#endif


